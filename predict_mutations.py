#!/usr/bin/env python3
"""In-silico saturation mutagenesis with batched prediction on SLURM.

For every position in a region, write one mutant genome per alternative base,
submit a prediction job for each, wait for the batch, then delete the FASTAs
before starting the next batch. Batching bounds peak disk usage (every mutant
holds a full copy of the chromosome) while still letting many jobs queue at once.

    python ism_pipeline.py chr9:45260000-45261900 --experiment TMPRSS13
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, NamedTuple, Sequence

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BASES = ("A", "C", "G", "T")
LOG = logging.getLogger("ism")
REGION_RE = re.compile(r"^(?P<chrom>.+):(?P<start>[\d,_]+)(?:-(?P<end>[\d,_]+))?$")

# Prediction settings that were constants in the original JSON writer. Edit here
# if the pipeline ever needs them changed; they are not worth a CLI flag each.
# Timing for the SLURM wait loop. Deliberately not CLI flags: these are tuning
# details, not per-run choices.
POLL_SECONDS = 30       # gap between squeue checks
TIMEOUT_SECONDS = 24 * 3600   # give up on a batch after this long
GRACE_SECONDS = 120     # extra wait for output to land after a job leaves the queue

JSON_CONSTANTS = {
    "get_peaks": True,
    "merged": False,
    "is_explained": False,
    "save_tensors": False,
    "save_sequences": False,
}


class UserError(Exception):
    """A problem the caller can fix; reported as a message, not a traceback."""


class Mutant(NamedTuple):
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str

    @property
    def name(self) -> str:
        return f"{self.chrom}_{self.pos}_{self.ref}-{self.alt}"


@dataclass
class Config:
    """Everything the run needs. Built from the CLI, easy to fake in tests."""

    genome: Path
    output_dir: Path
    json_dir: Path
    model_py: Path
    pipeline_script: Path
    hf_model: str
    cell_type: int
    tag: str
    pad: int
    batch_size: int
    threads: int
    poll_seconds: int = POLL_SECONDS
    timeout_seconds: int = TIMEOUT_SECONDS
    grace_seconds: int = GRACE_SECONDS
    index: bool = True
    keep_fasta: bool = False
    dry_run: bool = False

    def dir_for(self, mutant: Mutant) -> Path:
        return self.output_dir / mutant.name

    def predictions_dir_for(self, mutant: Mutant) -> Path:
        return self.dir_for(mutant) / "Predictions"

    def marker_for(self, mutant: Mutant) -> Path:
        """The bigwig that means this mutant is done. Derived, so it cannot
        drift from the output_folder handed to the pipeline in the JSON."""
        return (
            self.predictions_dir_for(mutant)
            / "predictions"
            / str(self.cell_type)
            / "prediction_ref.bw"
        )

    def json_for(self, mutant: Mutant) -> Path:
        return self.json_dir / f"{mutant.name}.json"


# --------------------------------------------------------------------------- #
# Region and mutants
# --------------------------------------------------------------------------- #


def parse_region(text: str) -> tuple[str, int, int]:
    """'chr9:45260000-45261900' -> ('chr9', 45260000, 45261900), 1-based inclusive."""
    match = REGION_RE.match(text.strip())
    if match is None:
        raise UserError(
            f"Could not parse region {text!r}. Expected 'chrom:start-end' or "
            "'chrom:pos', 1-based inclusive, e.g. 'chr9:45260000-45261900'."
        )
    clean = lambda s: int(s.replace(",", "").replace("_", ""))
    start = clean(match["start"])
    end = clean(match["end"]) if match["end"] else start
    if start < 1:
        raise UserError(f"Region start must be >= 1 (1-based), got {start}.")
    if end < start:
        raise UserError(f"Region end ({end}) is before its start ({start}).")
    return match["chrom"], start, end


def parse_cell_type(text: str) -> int:
    """A single cell type id, e.g. '3'."""
    try:
        value = int(text.strip())
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Cell type must be a single integer id, got {text!r}."
        ) from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"Cell type id must be >= 0, got {value}.")
    return value


def load_chromosome(fasta: Path, chrom: str) -> SeqRecord:
    """Return the record named `chrom`, reading the FASTA only once."""
    if not fasta.is_file():
        raise UserError(f"Genome FASTA not found: {fasta}")
    seen: list[str] = []
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id == chrom:
            return record
        seen.append(record.id)
    raise UserError(
        f"Sequence {chrom!r} not found in {fasta}. "
        f"Available: {', '.join(seen[:10]) if seen else 'none - is this a FASTA?'}"
    )


def plan_mutants(record: SeqRecord, start: int, end: int) -> list[Mutant]:
    """Every substitution in the region, skipping non-ACGT reference bases."""
    seq = str(record.seq)
    if end > len(seq):
        raise UserError(
            f"Region ends at {end} but {record.id} is only {len(seq):,} bp long."
        )
    out = []
    for pos in range(start, end + 1):
        ref = seq[pos - 1].upper()
        if ref in BASES:
            out += [Mutant(record.id, pos, ref, a) for a in BASES if a != ref]
    return out


def check_windows(mutants: Sequence[Mutant], pad: int, chrom_length: int) -> None:
    """Warn if the model's input window runs off the chromosome."""
    bad = [m for m in mutants if m.pos - pad < 1 or m.pos + pad > chrom_length]
    if bad:
        LOG.warning(
            "%d position(s) sit within %d bp of a chromosome end, so their "
            "prediction window (e.g. %d-%d) falls outside %s. The pipeline may "
            "pad or fail on these.",
            len(bad), pad, bad[0].pos - pad, bad[0].pos + pad, bad[0].chrom,
        )


# --------------------------------------------------------------------------- #
# Preparing one mutant
# --------------------------------------------------------------------------- #


def write_json(path: Path, mutant: Mutant, config: Config) -> None:
    """Write the prediction config. Paths follow output_dir, not a hardcoded tree."""
    folder = config.dir_for(mutant)
    payload = {
        "chromosome": mutant.chrom,
        # The pipeline reads this as an array; one id still goes in a list.
        "cell_types": [config.cell_type],
        "start": mutant.pos - config.pad,
        "end": mutant.pos + config.pad,
        **JSON_CONSTANTS,
        "output_folder": str(config.predictions_dir_for(mutant)),
        "input_folder": str(folder),
        "tag": config.tag,
        "HFmodel": config.hf_model,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def prepare_mutant(mutant: Mutant, template: str, record: SeqRecord, config: Config) -> None:
    """Write the mutant FASTA, index it, copy the model, write the JSON."""
    folder = config.dir_for(mutant)
    folder.mkdir(parents=True, exist_ok=True)
    fasta = folder / "genome.fa"

    index = mutant.pos - 1
    if template[index].upper() != mutant.ref:
        raise ValueError(f"reference mismatch at {mutant.chrom}:{mutant.pos}")
    mutated = SeqRecord(
        Seq(template[:index] + mutant.alt + template[index + 1 :]),
        id=record.id,
        description=record.description,
    )
    # Write then rename, so an interrupted run leaves no half-written FASTA that
    # a later run would mistake for finished work.
    tmp = folder / "genome.fa.tmp"
    SeqIO.write(mutated, tmp, "fasta")
    os.replace(tmp, fasta)

    if config.index and not (folder / "genome.fa.fai").is_file():
        subprocess.run(
            ["samtools", "faidx", str(fasta)],
            check=True, capture_output=True, text=True,
        )
    shutil.copy2(config.model_py, folder / config.model_py.name)
    write_json(config.json_for(mutant), mutant, config)


def cleanup_mutant(mutant: Mutant, config: Config) -> None:
    """Drop the big intermediates, keeping the predictions."""
    if config.keep_fasta:
        return
    folder = config.dir_for(mutant)
    for name in ("genome.fa", "genome.fa.fai", "genome.fa.tmp", config.model_py.name):
        (folder / name).unlink(missing_ok=True)
    shutil.rmtree(folder / "__pycache__", ignore_errors=True)


# --------------------------------------------------------------------------- #
# SLURM
# --------------------------------------------------------------------------- #


def submit(mutant: Mutant, config: Config) -> str:
    """sbatch one prediction, returning its job id."""
    result = subprocess.run(
        ["sbatch", "--parsable", str(config.pipeline_script), str(config.json_for(mutant))],
        check=True, capture_output=True, text=True,
    )
    return result.stdout.strip().split(";")[0]


def still_running(job_ids: set[str]) -> set[str]:
    """Which of these jobs squeue still knows about."""
    if not job_ids:
        return set()
    result = subprocess.run(
        ["squeue", "-h", "-o", "%i", "-j", ",".join(sorted(job_ids))],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        # squeue errors when every id has been purged, which means they finished.
        LOG.debug("squeue returned %d: %s", result.returncode, result.stderr.strip())
        return set()
    listed = {line.strip().split("_")[0] for line in result.stdout.split() if line.strip()}
    return job_ids & listed


def wait_for_batch(jobs: Sequence[tuple[str, Mutant]], config: Config) -> list[Mutant]:
    """Block until every job leaves the queue. Returns mutants with no output."""
    pending = {job_id for job_id, _ in jobs}
    deadline = time.monotonic() + config.timeout_seconds
    while pending:
        if time.monotonic() > deadline:
            LOG.error(
                "Timed out after %ds with %d job(s) still queued: %s",
                config.timeout_seconds, len(pending), ", ".join(sorted(pending)),
            )
            break
        time.sleep(config.poll_seconds)
        pending = still_running(pending)
        LOG.debug("%d/%d jobs still in the queue", len(pending), len(jobs))

    # A job can leave the queue slightly before its output lands on a shared
    # filesystem, so give the markers a moment to appear before calling failure.
    missing = [m for _, m in jobs if not config.marker_for(m).is_file()]
    grace = time.monotonic() + config.grace_seconds
    while missing and time.monotonic() < grace:
        time.sleep(min(config.poll_seconds, 10))
        missing = [m for m in missing if not config.marker_for(m).is_file()]
    return missing


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #


def batched(items: Sequence[Mutant], size: int) -> Iterator[list[Mutant]]:
    for i in range(0, len(items), size):
        yield list(items[i : i + size])


def run(mutants: Sequence[Mutant], record: SeqRecord, config: Config) -> int:
    """Prepare, submit, wait and clean up, one batch at a time. Returns failures."""
    template = str(record.seq)
    todo = [m for m in mutants if not config.marker_for(m).is_file()]
    done = len(mutants) - len(todo)
    if done:
        LOG.info("%d mutant(s) already predicted, skipping.", done)
    if not todo:
        return 0

    failures = 0
    batches = list(batched(todo, config.batch_size))
    for number, batch in enumerate(batches, start=1):
        LOG.info("Batch %d/%d: preparing %d mutant(s)", number, len(batches), len(batch))
        if config.dry_run:
            for mutant in batch:
                LOG.info("  would submit %s", config.json_for(mutant))
            continue

        def prepare(mutant: Mutant) -> Mutant | None:
            try:
                prepare_mutant(mutant, template, record, config)
                return mutant
            except Exception as exc:  # don't lose the whole batch to one mutant
                LOG.error("Could not prepare %s: %s", mutant.name, exc)
                return None

        with ThreadPoolExecutor(max_workers=config.threads) as pool:
            ready = [m for m in pool.map(prepare, batch) if m is not None]
        failures += len(batch) - len(ready)

        jobs: list[tuple[str, Mutant]] = []
        for mutant in ready:
            try:
                jobs.append((submit(mutant, config), mutant))
            except subprocess.CalledProcessError as exc:
                LOG.error("sbatch failed for %s: %s", mutant.name, exc.stderr.strip())
                failures += 1
        LOG.info("Batch %d: submitted %d job(s), waiting", number, len(jobs))

        missing = wait_for_batch(jobs, config)
        for mutant in missing:
            LOG.error("No prediction for %s (expected %s)", mutant.name, config.marker_for(mutant))
        failures += len(missing)

        succeeded = [m for m in ready if m not in missing]
        for mutant in succeeded:
            cleanup_mutant(mutant, config)
        LOG.info("Batch %d done: %d ok, %d missing", number, len(succeeded), len(missing))

    return failures


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("region", help="1-based inclusive, e.g. chr9:45260000-45261900")

    names = parser.add_argument_group("names and locations")
    names.add_argument("-e", "--experiment", default="Locus",
                       help="Experiment name, used to build the default paths "
                            "(default: %(default)s)")
    names.add_argument("--base-dir", type=Path, default=Path("08_Pipeline_predictions"),
                       help="Root of the prediction tree (default: %(default)s)")
    names.add_argument("--genome", type=Path,
                       help="Reference FASTA (default: <base-dir>/genome.fa)")
    names.add_argument("--output-dir", type=Path,
                       help="Mutant genomes and predictions "
                            "(default: <base-dir>/<experiment>)")
    names.add_argument("--json-dir", type=Path,
                       help="Prediction configs (default: 07_JSONs/<experiment>)")
    names.add_argument("--model-py", type=Path,
                       help="model.py copied into each mutant folder "
                            "(default: <base-dir>/model.py)")
    names.add_argument("--pipeline-script", type=Path, default=Path("i_Run_Pipeline.sh"),
                       help="Script passed to sbatch (default: %(default)s)")

    model = parser.add_argument_group("prediction")
    model.add_argument("-m", "--model", dest="hf_model", default="GFIO/Human_Ery",
                       help="HFmodel field in the JSON (default: %(default)s)")
    model.add_argument("-c", "--cell-type", type=parse_cell_type, required=True,
                       help="Cell type id to predict, e.g. 3")
    model.add_argument("--tag", default="pred", help="tag field (default: %(default)s)")
    model.add_argument("--pad", type=int, default=98304,
                       help="Half-width of the model input window "
                            "(default: %(default)s)")

    sched = parser.add_argument_group("scheduling")
    sched.add_argument("-b", "--batch-size", type=int, default=30,
                       help="Mutants in flight at once. Each holds a full copy of "
                            "the chromosome on disk until its batch finishes "
                            "(default: %(default)s)")
    sched.add_argument("-t", "--threads", type=int, default=4,
                       help="Parallel FASTA writers (default: %(default)s)")

    other = parser.add_argument_group("other")
    other.add_argument("--no-index", dest="index", action="store_false",
                       help="Skip samtools faidx")
    other.add_argument("--keep-fasta", action="store_true",
                       help="Do not delete mutant genomes after prediction")
    other.add_argument("-n", "--dry-run", action="store_true",
                       help="Report what would be submitted, then stop")
    other.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    return parser


def config_from_args(args: argparse.Namespace) -> Config:
    output_dir = args.output_dir or args.base_dir / args.experiment
    return Config(
        genome=args.genome or args.base_dir / "genome.fa",
        output_dir=output_dir,
        json_dir=args.json_dir or Path("07_JSONs") / args.experiment,
        model_py=args.model_py or args.base_dir / "model.py",
        pipeline_script=args.pipeline_script,
        hf_model=args.hf_model,
        cell_type=args.cell_type,
        tag=args.tag,
        pad=args.pad,
        batch_size=args.batch_size,
        threads=args.threads,
        index=args.index,
        keep_fasta=args.keep_fasta,
        dry_run=args.dry_run,
    )


def check_environment(config: Config) -> None:
    for value, flag in ((config.batch_size, "--batch-size"), (config.threads, "--threads")):
        if value < 1:
            raise UserError(f"{flag} must be >= 1.")
    if not config.model_py.is_file():
        raise UserError(f"model.py not found: {config.model_py}")
    if config.dry_run:
        return
    if not config.pipeline_script.is_file():
        raise UserError(f"Pipeline script not found: {config.pipeline_script}")
    for tool in (["sbatch", "squeue"] + (["samtools"] if config.index else [])):
        if shutil.which(tool) is None:
            raise UserError(
                f"{tool} is not on PATH. Run this from a submit node"
                + (", or pass --no-index." if tool == "samtools" else ".")
            )


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S",
    )
    try:
        chrom, start, end = parse_region(args.region)
        config = config_from_args(args)
        check_environment(config)
        record = load_chromosome(config.genome, chrom)
        mutants = plan_mutants(record, start, end)
    except UserError as exc:
        LOG.error("%s", exc)
        return 2

    if not mutants:
        LOG.warning("No mutable positions in %s:%d-%d.", chrom, start, end)
        return 0

    check_windows(mutants, config.pad, len(record.seq))
    LOG.info(
        "%d mutants over %d position(s), %d per batch, model %s, cell type %d",
        len(mutants), end - start + 1, config.batch_size, config.hf_model,
        config.cell_type,
    )
    try:
        failures = run(mutants, record, config)
    except KeyboardInterrupt:
        LOG.error("Interrupted. Submitted jobs keep running; scancel them if needed.")
        return 130
    if failures:
        LOG.error("%d mutant(s) have no prediction. Re-run to retry just those.", failures)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
