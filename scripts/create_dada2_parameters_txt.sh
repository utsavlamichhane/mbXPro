#!/usr/bin/env bash
set -euo pipefail

# Standalone DADA2 parameter recommender.
#
# This file does not depend on any other local script. It reads a QIIME 2
# SampleData[PairedEndSequencesWithQuality] .qza, chooses DADA2 paired-end
# trim/truncation parameters, and writes a concise .txt file.
#
# Primer behavior:
#   - If --forward-primer is provided, its base count is used for
#     --p-trim-left-f.
#   - If --reverse-primer is provided, its base count is used for
#     --p-trim-left-r.
#   - If either primer is not provided, that side defaults to 20 bases.
#
# Examples:
#   ./create_dada2_parameters_txt.sh
#   ./create_dada2_parameters_txt.sh Paired_End_artifact.qza
#   ./create_dada2_parameters_txt.sh Paired_End_artifact.qza dada2_parameters.txt
#   ./create_dada2_parameters_txt.sh -i Paired_End_artifact.qza -o params.txt \
#     --forward-primer GTGYCAGCMGCCGCGGTAA \
#     --reverse-primer GGACTACNVGGGTWTCTAAT
#
# Useful optional settings:
#   --assume-primer-length 20
#   --min-overlap 12
#   --reads-per-sample 1000
#   --amplicon-length 253
#   --min-pair-pass-rate 0.80
#   --min-forward-pass-rate 0.90
#   --min-reverse-pass-rate 0.84
#   --length-slack 1

usage() {
  sed -n '2,38p' "$0" >&2
}

DEMUX_QZA=""
PARAM_TXT="dada2_parameters.txt"
FORWARD_PRIMER="${FORWARD_PRIMER:-}"
REVERSE_PRIMER="${REVERSE_PRIMER:-}"
ASSUME_PRIMER_LENGTH="${ASSUME_PRIMER_LENGTH:-20}"
MIN_OVERLAP="${MIN_OVERLAP:-12}"
READS_PER_SAMPLE="${READS_PER_SAMPLE:-1000}"
INSERT_SAMPLE_PAIRS="${INSERT_SAMPLE_PAIRS:-2500}"
AMPLICON_LENGTH="${AMPLICON_LENGTH:-}"
MAX_EE_F="${MAX_EE_F:-2.0}"
MAX_EE_R="${MAX_EE_R:-2.0}"
MIN_PAIR_PASS_RATE="${MIN_PAIR_PASS_RATE:-0.80}"
MIN_FORWARD_PASS_RATE="${MIN_FORWARD_PASS_RATE:-0.90}"
MIN_REVERSE_PASS_RATE="${MIN_REVERSE_PASS_RATE:-0.84}"
LENGTH_SLACK="${LENGTH_SLACK:-1}"
TOP_CANDIDATES="${TOP_CANDIDATES:-20}"

positionals=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input|--demux-qza)
      DEMUX_QZA="${2:?Missing value for $1}"
      shift 2
      ;;
    -o|--output)
      PARAM_TXT="${2:?Missing value for $1}"
      shift 2
      ;;
    --forward-primer|-f)
      FORWARD_PRIMER="${2:?Missing value for $1}"
      shift 2
      ;;
    --reverse-primer|-r)
      REVERSE_PRIMER="${2:?Missing value for $1}"
      shift 2
      ;;
    --assume-primer-length)
      ASSUME_PRIMER_LENGTH="${2:?Missing value for $1}"
      shift 2
      ;;
    --min-overlap)
      MIN_OVERLAP="${2:?Missing value for $1}"
      shift 2
      ;;
    --reads-per-sample)
      READS_PER_SAMPLE="${2:?Missing value for $1}"
      shift 2
      ;;
    --insert-sample-pairs)
      INSERT_SAMPLE_PAIRS="${2:?Missing value for $1}"
      shift 2
      ;;
    --amplicon-length)
      AMPLICON_LENGTH="${2:?Missing value for $1}"
      shift 2
      ;;
    --max-ee-f)
      MAX_EE_F="${2:?Missing value for $1}"
      shift 2
      ;;
    --max-ee-r)
      MAX_EE_R="${2:?Missing value for $1}"
      shift 2
      ;;
    --min-pair-pass-rate)
      MIN_PAIR_PASS_RATE="${2:?Missing value for $1}"
      shift 2
      ;;
    --min-forward-pass-rate)
      MIN_FORWARD_PASS_RATE="${2:?Missing value for $1}"
      shift 2
      ;;
    --min-reverse-pass-rate)
      MIN_REVERSE_PASS_RATE="${2:?Missing value for $1}"
      shift 2
      ;;
    --length-slack)
      LENGTH_SLACK="${2:?Missing value for $1}"
      shift 2
      ;;
    --top-candidates)
      TOP_CANDIDATES="${2:?Missing value for $1}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        positionals+=("$1")
        shift
      done
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      usage
      exit 1
      ;;
    *)
      positionals+=("$1")
      shift
      ;;
  esac
done

if [[ "${#positionals[@]}" -ge 1 && -z "${DEMUX_QZA}" ]]; then
  DEMUX_QZA="${positionals[0]}"
fi

if [[ "${#positionals[@]}" -ge 2 && "${PARAM_TXT}" == "dada2_parameters.txt" ]]; then
  PARAM_TXT="${positionals[1]}"
fi

export DEMUX_QZA
export PARAM_TXT
export FORWARD_PRIMER
export REVERSE_PRIMER
export ASSUME_PRIMER_LENGTH
export MIN_OVERLAP
export READS_PER_SAMPLE
export INSERT_SAMPLE_PAIRS
export AMPLICON_LENGTH
export MAX_EE_F
export MAX_EE_R
export MIN_PAIR_PASS_RATE
export MIN_FORWARD_PASS_RATE
export MIN_REVERSE_PASS_RATE
export LENGTH_SLACK
export TOP_CANDIDATES

python3 - <<'PY'
from __future__ import annotations

import csv
import gzip
import math
import os
import statistics
import sys
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}
RC_TABLE = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


@dataclass(frozen=True)
class ManifestEntry:
    sample_id: str
    filename: str
    direction: str
    member: str


@dataclass
class ReadPair:
    forward_seq: str
    reverse_seq: str
    forward_qual: bytes
    reverse_qual: bytes


@dataclass
class InsertEstimate:
    used: int
    aligned: int
    median: int
    p95: int
    p99: int
    requirement: int
    identity_median: float | None


@dataclass
class TrimInfo:
    forward: int
    reverse: int
    forward_source: str
    reverse_source: str
    forward_primer: str
    reverse_primer: str
    forward_primer_matches: tuple[int, int] | None
    reverse_primer_matches: tuple[int, int] | None


@dataclass
class Candidate:
    trunc_f: int
    trunc_r: int
    trim_f: int
    trim_r: int
    overlap: int
    pair_pass_rate: float
    forward_pass_rate: float
    reverse_pass_rate: float
    effective_f: int
    effective_r: int


def die(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def env_int(name: str, default: int | None = None) -> int | None:
    value = os.environ.get(name, "")
    if value == "":
        return default
    try:
        return int(value)
    except ValueError:
        die(f"{name} must be an integer, got {value!r}")


def env_float(name: str, default: float) -> float:
    value = os.environ.get(name, "")
    if value == "":
        return default
    try:
        return float(value)
    except ValueError:
        die(f"{name} must be numeric, got {value!r}")


def normalize_primer(primer: str) -> str:
    primer = "".join(ch for ch in primer.upper().replace("U", "T") if not ch.isspace() and ch != "-")
    if primer and any(ch not in IUPAC for ch in primer):
        bad = sorted({ch for ch in primer if ch not in IUPAC})
        die(f"Primer contains unsupported character(s): {''.join(bad)}")
    return primer


def percentile(values: list[int] | list[float], pct: float) -> int:
    if not values:
        die("Cannot calculate percentile on an empty list.")
    ordered = sorted(values)
    idx = math.ceil((pct / 100.0) * len(ordered)) - 1
    idx = min(max(idx, 0), len(ordered) - 1)
    return int(round(ordered[idx]))


class QzaFastqSource:
    def __init__(self, qza_path: Path):
        self.qza_path = qza_path
        self.zip_file = zipfile.ZipFile(qza_path)
        self.entries = self._read_manifest()

    def close(self) -> None:
        self.zip_file.close()

    def _read_text_member(self, suffix: str) -> tuple[str, str] | None:
        for name in self.zip_file.namelist():
            if name.endswith(suffix):
                with self.zip_file.open(name) as handle:
                    return name, handle.read().decode("utf-8")
        return None

    def metadata_text(self) -> str:
        result = self._read_text_member("/metadata.yaml")
        return result[1] if result else ""

    def _read_manifest(self) -> list[ManifestEntry]:
        result = self._read_text_member("/data/MANIFEST")
        if result is None:
            raise ValueError("No data/MANIFEST found.")
        manifest_member, text = result
        data_prefix = manifest_member.rsplit("/", 1)[0]
        rows = list(csv.DictReader(text.splitlines()))
        entries: list[ManifestEntry] = []
        for row in rows:
            sample_id = row.get("sample-id") or row.get("sample id") or row.get("sample_id")
            filename = row.get("filename")
            direction = row.get("direction")
            if not sample_id or not filename or not direction:
                raise ValueError(f"Unexpected MANIFEST row: {row}")
            entries.append(ManifestEntry(sample_id, filename, direction, f"{data_prefix}/{filename}"))
        return entries

    def open_fastq(self, entry: ManifestEntry) -> Iterator[tuple[str, str]]:
        with self.zip_file.open(entry.member) as raw:
            with gzip.open(raw, "rt") as fastq:
                while True:
                    header = fastq.readline()
                    if not header:
                        break
                    seq = fastq.readline().strip().upper()
                    fastq.readline()
                    qual = fastq.readline().strip()
                    if not qual:
                        break
                    yield seq, qual


def find_default_qza() -> Path:
    qza_env = os.environ.get("DEMUX_QZA", "").strip()
    if qza_env:
        return Path(qza_env)

    paired: list[Path] = []
    for path in sorted(Path.cwd().glob("*.qza")):
        try:
            source = QzaFastqSource(path)
            metadata = source.metadata_text()
            source.close()
        except Exception:
            continue
        if "SampleData[PairedEndSequencesWithQuality]" in metadata:
            paired.append(path)

    if not paired:
        die("No paired-end demultiplexed .qza found. Pass one with -i.")
    if len(paired) > 1:
        die("Multiple paired-end demultiplexed .qza files found. Pass one with -i.")
    return paired[0]


def sample_pairs(source: QzaFastqSource, reads_per_sample: int) -> tuple[list[ReadPair], dict[str, int]]:
    by_sample: dict[str, dict[str, ManifestEntry]] = defaultdict(dict)
    for entry in source.entries:
        by_sample[entry.sample_id][entry.direction] = entry

    pairs: list[ReadPair] = []
    sample_counts: dict[str, int] = {}
    for sample_id in sorted(by_sample):
        directions = by_sample[sample_id]
        if "forward" not in directions or "reverse" not in directions:
            continue
        count = 0
        f_iter = source.open_fastq(directions["forward"])
        r_iter = source.open_fastq(directions["reverse"])
        for (f_seq, f_qual), (r_seq, r_qual) in zip(f_iter, r_iter):
            pairs.append(
                ReadPair(
                    f_seq,
                    r_seq,
                    bytes(ord(ch) - 33 for ch in f_qual),
                    bytes(ord(ch) - 33 for ch in r_qual),
                )
            )
            count += 1
            if count >= reads_per_sample:
                break
        sample_counts[sample_id] = count

    if not pairs:
        die("No paired FASTQ reads could be sampled from the .qza.")
    return pairs, sample_counts


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper().replace("U", "T")


def best_pair_alignment(forward: str, reverse: str, kmer: int = 13, min_overlap: int = 40):
    rc_reverse = reverse_complement(reverse)
    if len(forward) < min_overlap or len(rc_reverse) < min_overlap:
        return None

    positions: dict[str, list[int]] = defaultdict(list)
    for i in range(0, len(forward) - kmer + 1):
        word = forward[i : i + kmer]
        if "N" not in word:
            positions[word].append(i)

    shift_counts: Counter[int] = Counter()
    for j in range(0, len(rc_reverse) - kmer + 1):
        word = rc_reverse[j : j + kmer]
        if "N" in word:
            continue
        for i in positions.get(word, []):
            shift_counts[i - j] += 1

    best = None
    for shift, _ in shift_counts.most_common(12):
        start = max(0, shift)
        end = min(len(forward), shift + len(rc_reverse))
        overlap = end - start
        if overlap < min_overlap:
            continue
        matches = 0
        comparable = 0
        for i in range(start, end):
            f_base = forward[i]
            r_base = rc_reverse[i - shift]
            if f_base == "N" or r_base == "N":
                continue
            comparable += 1
            if f_base == r_base:
                matches += 1
        if comparable < min_overlap:
            continue
        identity = matches / comparable
        score = identity * comparable
        insert_len = max(len(forward), shift + len(rc_reverse)) - min(0, shift)
        candidate = (score, insert_len, identity, overlap)
        if best is None or candidate > best:
            best = candidate

    if best is None:
        return None
    _, insert_len, identity, overlap = best
    return insert_len, identity, overlap


def estimate_insert_length(pairs: list[ReadPair], max_pairs: int, amplicon_length: int | None) -> InsertEstimate:
    if amplicon_length is not None:
        return InsertEstimate(0, 0, amplicon_length, amplicon_length, amplicon_length, amplicon_length, None)

    inserts: list[int] = []
    identities: list[float] = []
    used = min(max_pairs, len(pairs))
    for pair in pairs[:used]:
        alignment = best_pair_alignment(pair.forward_seq, pair.reverse_seq)
        if alignment is None:
            continue
        insert_len, identity, _overlap = alignment
        if identity >= 0.90:
            inserts.append(insert_len)
            identities.append(identity)

    if not inserts:
        die("Could not infer insert length. Pass --amplicon-length if raw reads do not overlap.")

    median = int(round(statistics.median(inserts)))
    if len(inserts) >= 8:
        q1, q3 = statistics.quantiles(inserts, n=4)[0], statistics.quantiles(inserts, n=4)[2]
        iqr = q3 - q1
        clean = [x for x in inserts if (q1 - 1.5 * iqr) <= x <= (q3 + 1.5 * iqr)]
    else:
        clean = inserts[:]
    if len(clean) < max(20, len(inserts) // 2):
        clean = [x for x in inserts if abs(x - median) <= 50] or inserts[:]

    return InsertEstimate(
        used=used,
        aligned=len(inserts),
        median=int(round(statistics.median(clean))),
        p95=percentile(clean, 95),
        p99=percentile(clean, 99),
        requirement=percentile(clean, 99),
        identity_median=round(float(statistics.median(identities)), 4) if identities else None,
    )


def primer_match_counts(pairs: list[ReadPair], direction: str, primer: str) -> tuple[int, int]:
    if not primer:
        return (0, 0)
    matched = 0
    total = 0
    for pair in pairs:
        seq = pair.forward_seq if direction == "forward" else pair.reverse_seq
        if len(seq) < len(primer):
            continue
        total += 1
        mismatches = 0
        comparable = 0
        for observed, expected in zip(seq[: len(primer)], primer):
            if observed == "N":
                continue
            comparable += 1
            if observed not in IUPAC[expected]:
                mismatches += 1
        if comparable and mismatches / comparable <= 0.15:
            matched += 1
    return matched, total


def decide_trim(pairs: list[ReadPair], default_len: int, forward_primer: str, reverse_primer: str) -> TrimInfo:
    if default_len < 0:
        die("--assume-primer-length must be >= 0")

    f_primer = normalize_primer(forward_primer)
    r_primer = normalize_primer(reverse_primer)

    if f_primer:
        trim_f = len(f_primer)
        f_source = f"forward primer length ({len(f_primer)} bases)"
        f_matches = primer_match_counts(pairs, "forward", f_primer)
    else:
        trim_f = default_len
        f_source = f"default unknown-primer assumption ({default_len} bases)"
        f_matches = None

    if r_primer:
        trim_r = len(r_primer)
        r_source = f"reverse primer length ({len(r_primer)} bases)"
        r_matches = primer_match_counts(pairs, "reverse", r_primer)
    else:
        trim_r = default_len
        r_source = f"default unknown-primer assumption ({default_len} bases)"
        r_matches = None

    return TrimInfo(trim_f, trim_r, f_source, r_source, f_primer, r_primer, f_matches, r_matches)


def max_trunc_from_read(seq: str, qual: bytes, trim_left: int, max_ee: float) -> int:
    ee = 0.0
    max_trunc = trim_left
    for pos in range(trim_left, len(qual)):
        if seq[pos] == "N":
            break
        ee += 10 ** (-(qual[pos] / 10))
        trunc_len = pos + 1
        if ee <= max_ee:
            max_trunc = trunc_len
    return max_trunc


def build_suffix_matrix(
    pairs: list[ReadPair],
    trim_f: int,
    trim_r: int,
    max_ee_f: float,
    max_ee_r: float,
    max_f_len: int,
    max_r_len: int,
):
    matrix = [[0 for _ in range(max_r_len + 1)] for _ in range(max_f_len + 1)]
    f_pass_counts = [0 for _ in range(max_f_len + 1)]
    r_pass_counts = [0 for _ in range(max_r_len + 1)]

    for pair in pairs:
        f_max = min(max_trunc_from_read(pair.forward_seq, pair.forward_qual, trim_f, max_ee_f), max_f_len)
        r_max = min(max_trunc_from_read(pair.reverse_seq, pair.reverse_qual, trim_r, max_ee_r), max_r_len)
        f_pass_counts[f_max] += 1
        r_pass_counts[r_max] += 1
        matrix[f_max][r_max] += 1

    suffix = [row[:] for row in matrix]
    for i in range(max_f_len - 1, -1, -1):
        row = suffix[i]
        below = suffix[i + 1]
        for j in range(max_r_len + 1):
            row[j] += below[j]
    for i in range(max_f_len + 1):
        row = suffix[i]
        for j in range(max_r_len - 1, -1, -1):
            row[j] += row[j + 1]
    for i in range(max_f_len - 1, -1, -1):
        f_pass_counts[i] += f_pass_counts[i + 1]
    for j in range(max_r_len - 1, -1, -1):
        r_pass_counts[j] += r_pass_counts[j + 1]

    return suffix, f_pass_counts, r_pass_counts


def rank_candidates(
    pairs: list[ReadPair],
    trim: TrimInfo,
    insert_requirement: int,
    min_overlap: int,
    max_ee_f: float,
    max_ee_r: float,
    min_pair_pass_rate: float,
    min_forward_pass_rate: float,
    min_reverse_pass_rate: float,
    length_slack: int,
) -> list[Candidate]:
    max_f_len = max(len(pair.forward_qual) for pair in pairs)
    max_r_len = max(len(pair.reverse_qual) for pair in pairs)
    suffix, f_pass_counts, r_pass_counts = build_suffix_matrix(
        pairs, trim.forward, trim.reverse, max_ee_f, max_ee_r, max_f_len, max_r_len
    )

    n = len(pairs)
    candidates: list[Candidate] = []
    for trunc_f in range(trim.forward + 50, max_f_len + 1):
        for trunc_r in range(trim.reverse + 50, max_r_len + 1):
            overlap = trunc_f + trunc_r - insert_requirement
            if overlap < min_overlap:
                continue
            pair_pass = suffix[trunc_f][trunc_r] / n
            if pair_pass == 0:
                continue
            candidates.append(
                Candidate(
                    trunc_f=trunc_f,
                    trunc_r=trunc_r,
                    trim_f=trim.forward,
                    trim_r=trim.reverse,
                    overlap=overlap,
                    pair_pass_rate=pair_pass,
                    forward_pass_rate=f_pass_counts[trunc_f] / n,
                    reverse_pass_rate=r_pass_counts[trunc_r] / n,
                    effective_f=trunc_f - trim.forward,
                    effective_r=trunc_r - trim.reverse,
                )
            )

    if not candidates:
        die("No candidate truncation lengths satisfy the constraints.")

    eligible = [
        c
        for c in candidates
        if c.pair_pass_rate >= min_pair_pass_rate
        and c.forward_pass_rate >= min_forward_pass_rate
        and c.reverse_pass_rate >= min_reverse_pass_rate
    ]
    if not eligible:
        eligible = [c for c in candidates if c.pair_pass_rate >= min_pair_pass_rate]
    if not eligible:
        eligible = candidates[:]

    max_effective_len = max(c.effective_f + c.effective_r for c in eligible)
    near_longest = [c for c in eligible if c.effective_f + c.effective_r >= max_effective_len - length_slack]
    near_longest.sort(
        key=lambda c: (
            c.pair_pass_rate,
            -abs(c.effective_f - c.effective_r),
            c.overlap,
            c.effective_f + c.effective_r,
        ),
        reverse=True,
    )
    selected = near_longest[0]

    eligible_ids = {id(c) for c in eligible}
    remaining = [c for c in candidates if c is not selected]
    remaining.sort(
        key=lambda c: (
            id(c) in eligible_ids,
            c.effective_f + c.effective_r,
            c.pair_pass_rate,
            -abs(c.effective_f - c.effective_r),
        ),
        reverse=True,
    )
    return [selected] + remaining


def format_match_counts(counts: tuple[int, int] | None) -> str:
    if counts is None:
        return "not checked; primer not provided"
    matched, total = counts
    pct = (matched / total * 100.0) if total else 0.0
    return f"{matched}/{total} sampled reads ({pct:.1f}%)"


def write_output(
    output_path: Path,
    qza_path: Path,
    sample_counts: dict[str, int],
    insert: InsertEstimate,
    trim: TrimInfo,
    candidates: list[Candidate],
    top_n: int,
) -> None:
    best = candidates[0]
    output_path.parent.mkdir(parents=True, exist_ok=True) if output_path.parent != Path(".") else None

    with output_path.open("w") as handle:
        handle.write("Recommended DADA2 parameters\n\n")
        handle.write(f"Input artifact: {qza_path}\n")
        handle.write(f"Samples used: {len(sample_counts)}\n")
        handle.write(f"Sampled read pairs: {sum(sample_counts.values())}\n\n")

        handle.write(f"--p-trunc-len-f {best.trunc_f}\n")
        handle.write(f"--p-trunc-len-r {best.trunc_r}\n")
        handle.write(f"--p-trim-left-f {best.trim_f}\n")
        handle.write(f"--p-trim-left-r {best.trim_r}\n\n")

        handle.write("Primer handling\n")
        handle.write(f"Forward primer provided: {trim.forward_primer or 'no'}\n")
        handle.write(f"Reverse primer provided: {trim.reverse_primer or 'no'}\n")
        handle.write(f"Forward trim-left source: {trim.forward_source}\n")
        handle.write(f"Reverse trim-left source: {trim.reverse_source}\n")
        handle.write(f"Forward primer match count: {format_match_counts(trim.forward_primer_matches)}\n")
        handle.write(f"Reverse primer match count: {format_match_counts(trim.reverse_primer_matches)}\n\n")

        handle.write("Quality/overlap estimates\n")
        handle.write(f"Insert aligned pairs: {insert.aligned}/{insert.used}\n")
        handle.write(f"Insert median/p95/p99: {insert.median}/{insert.p95}/{insert.p99}\n")
        handle.write(f"Insert requirement used: {insert.requirement}\n")
        handle.write(f"Median alignment identity: {insert.identity_median}\n")
        handle.write(f"Estimated overlap buffer: {best.overlap} bp\n")
        handle.write(f"Estimated paired-read pass rate: {best.pair_pass_rate:.4f}\n")
        handle.write(f"Estimated forward-read pass rate: {best.forward_pass_rate:.4f}\n")
        handle.write(f"Estimated reverse-read pass rate: {best.reverse_pass_rate:.4f}\n")
        handle.write(f"Effective retained forward length: {best.effective_f} bp\n")
        handle.write(f"Effective retained reverse length: {best.effective_r} bp\n\n")

        handle.write("QIIME 2 command\n")
        handle.write("qiime dada2 denoise-paired \\\n")
        handle.write(f"  --i-demultiplexed-seqs {qza_path.name} \\\n")
        handle.write(f"  --p-trunc-len-f {best.trunc_f} \\\n")
        handle.write(f"  --p-trunc-len-r {best.trunc_r} \\\n")
        handle.write(f"  --p-trim-left-f {best.trim_f} \\\n")
        handle.write(f"  --p-trim-left-r {best.trim_r} \\\n")
        handle.write("  --o-table feature_table.qza \\\n")
        handle.write("  --o-representative-sequences representative_sequences.qza \\\n")
        handle.write("  --o-denoising-stats dada2_stats.qza\n\n")

        handle.write("Top candidate table\n")
        handle.write(
            "rank\ttrunc_len_f\ttrunc_len_r\ttrim_left_f\ttrim_left_r\t"
            "estimated_overlap\testimated_pair_pass_rate\testimated_forward_pass_rate\t"
            "estimated_reverse_pass_rate\teffective_forward_length\teffective_reverse_length\n"
        )
        for rank, candidate in enumerate(candidates[:top_n], start=1):
            handle.write(
                f"{rank}\t{candidate.trunc_f}\t{candidate.trunc_r}\t{candidate.trim_f}\t"
                f"{candidate.trim_r}\t{candidate.overlap}\t{candidate.pair_pass_rate:.6f}\t"
                f"{candidate.forward_pass_rate:.6f}\t{candidate.reverse_pass_rate:.6f}\t"
                f"{candidate.effective_f}\t{candidate.effective_r}\n"
            )


def main() -> int:
    qza_path = find_default_qza()
    if not qza_path.exists():
        die(f"Input artifact does not exist: {qza_path}")

    reads_per_sample = env_int("READS_PER_SAMPLE", 1000)
    insert_sample_pairs = env_int("INSERT_SAMPLE_PAIRS", 2500)
    assume_primer_length = env_int("ASSUME_PRIMER_LENGTH", 20)
    min_overlap = env_int("MIN_OVERLAP", 12)
    top_candidates = env_int("TOP_CANDIDATES", 20)
    amplicon_length = env_int("AMPLICON_LENGTH", None)
    max_ee_f = env_float("MAX_EE_F", 2.0)
    max_ee_r = env_float("MAX_EE_R", 2.0)
    min_pair_pass_rate = env_float("MIN_PAIR_PASS_RATE", 0.80)
    min_forward_pass_rate = env_float("MIN_FORWARD_PASS_RATE", 0.90)
    min_reverse_pass_rate = env_float("MIN_REVERSE_PASS_RATE", 0.84)
    length_slack = env_int("LENGTH_SLACK", 1)

    source = QzaFastqSource(qza_path)
    try:
        metadata = source.metadata_text()
        if "SampleData[PairedEndSequencesWithQuality]" not in metadata:
            die(f"{qza_path} is not SampleData[PairedEndSequencesWithQuality].")
        pairs, sample_counts = sample_pairs(source, reads_per_sample)
    finally:
        source.close()

    trim = decide_trim(
        pairs,
        assume_primer_length if assume_primer_length is not None else 20,
        os.environ.get("FORWARD_PRIMER", ""),
        os.environ.get("REVERSE_PRIMER", ""),
    )
    insert = estimate_insert_length(pairs, insert_sample_pairs or 2500, amplicon_length)
    candidates = rank_candidates(
        pairs,
        trim,
        insert.requirement,
        min_overlap or 12,
        max_ee_f,
        max_ee_r,
        min_pair_pass_rate,
        min_forward_pass_rate,
        min_reverse_pass_rate,
        length_slack or 1,
    )

    output_path = Path(os.environ.get("PARAM_TXT", "dada2_parameters.txt"))
    write_output(output_path, qza_path, sample_counts, insert, trim, candidates, top_candidates or 20)

    best = candidates[0]
    print("Recommended DADA2 parameters")
    print(f"  --p-trunc-len-f {best.trunc_f}")
    print(f"  --p-trunc-len-r {best.trunc_r}")
    print(f"  --p-trim-left-f {best.trim_f}")
    print(f"  --p-trim-left-r {best.trim_r}")
    print(f"Wrote {output_path}")
    return 0


raise SystemExit(main())
PY
