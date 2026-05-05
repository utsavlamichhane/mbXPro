#!/usr/bin/env bash
# =============================================================================
#  mbx_primer_identifier.sh  v2.0
#  Identify which 16S/18S/ITS primer pair was used in a set of FASTQ files,
#  OR detect that primers have been already trimmed and infer the V-region.
#
#  Compatible with bash 3.2+ (macOS default shell).
#  Requires: python3 (standard with QIIME2 conda env), gzip.
#
#  THREE-TIER DETECTION ALGORITHM
#  ──────────────────────────────
#    TIER 1  Direct primer detection (Cutadapt-style):
#              IUPAC-aware sliding window, all primers × 3 orientations,
#              best by combined match rate.  N (in primer or read) treated as
#              wild per IUPAC convention.
#    TIER 2  V-region motif anchor detection:
#              When TIER 1 finds nothing, scan the 5'-end of reads for
#              conserved 16S motifs (e.g. V4 starts with TACG.AGG right after
#              515F).  If matched, report DETECTION_STATUS=TRIMMED with the
#              inferred region.
#    TIER 3  Failure with rich diagnostics:
#              Read-length distribution, position-wise base composition,
#              top-5 candidate primers below threshold, top-3 most frequent
#              5'-prefixes, and concrete next-step suggestions.
#
#  Output: <out_dir>/0_primer_handling/mbx_primer_info.txt   (key=value pairs)
# =============================================================================

set -euo pipefail

# ── Helpers ───────────────────────────────────────────────────────────────────
err() {
  echo "" >&2
  echo "╔══════════════════════════════════════════════════════════════╗" >&2
  echo "║  ERROR                                                       ║" >&2
  echo "╚══════════════════════════════════════════════════════════════╝" >&2
  echo "[ERROR] $*" >&2
  echo "" >&2
  exit 1
}
warn() { echo "[WARN]  $*" >&2; }
info() { echo "[INFO]  $*"; }
ok()   { echo "[OK]    $*"; }
sep()  { echo "────────────────────────────────────────────────────────────────"; }

_abspath() {
  if [[ -d "$1" ]]; then cd "$1" && pwd
  elif [[ -f "$1" ]]; then echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  else return 1; fi
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_primer_identifier.sh v2.0 — Identify 16S primers from FASTQ files,
                                 OR detect that primers were trimmed off.

USAGE:
  mbx_primer_identifier.sh <fastq_dir> [OPTIONS]

DESCRIPTION:
  Samples reads from .fastq[.gz] files and runs a three-tier detection:
    1. Direct primer detection (IUPAC sliding window, all orientations)
    2. V-region motif anchors (detects "primers were already trimmed")
    3. Rich failure diagnostics with concrete suggestions

  Search order: <fastq_dir> itself, its parent, and one level of subdirectories
  (matches create_manifest.sh).

  Output:  <parent>/mbX_pro_outputs_<timestamp>/0_primer_handling/mbx_primer_info.txt
           (or under MBX_OUT_DIR if exported by the orchestrator)

BUILT-IN PRIMER DATABASE (35+ pairs, compiled from published protocols):
  V1-V2:  27F/338R,  8F/341R
  V1-V3:  27F/519R,  27F/534R,  8F/534R,  68F/518R
  V1-V9:  27F/1492R
  V3:     338F/537R
  V3-V4:  341F/785R, 341F/805R, PRK341F/PRK806R, Bakt_341F/Bakt_805R
  V3-V5:  341F/926R
  V4:     515F/806R (Caporaso),  515F-Parada/806R-Apprill (EMP)
  V4-V5:  515F-Parada/926R
  V5-V7:  799F/1193R, 799F/1391R   (chloroplast-skipping, plant-friendly)
  V6-V7:  967F/1391R
  V6-V8:  926F/1392R, B969F/BA1406R, 968F/1401R
  Plus Earth Microbiome Project, Nadkarni, and Klindworth-2013 variants.

OPTIONS:
  --forward-primer <SEQ> Skip detection — record this 5'->3' forward primer.
  --reverse-primer <SEQ> Skip detection — record this 5'->3' reverse primer.
                         (Both flags together produce a primer_info.txt with
                         DETECTION_STATUS=USER_SUPPLIED.)
  --samples <N>          Reads to sample per orientation             (default: 10000)
  --mismatches <N>       Max IUPAC mismatches per primer per read    (default: 3)
  --min-rate <0-1>       Minimum match-rate to call a hit            (default: 0.05)
  --offset <N>           Max 5'-offset to slide primer over          (default: 25)
  --extra-primers FILE   Additional primers (TSV: name<TAB>seq<TAB>region<TAB>F|R)
  --report-best          Report best candidate even below --min-rate (tentative)
  --debug-csv            Dump full primer×orientation rate matrix to CSV
  --strict               Fail on any non-ASCII/non-IUPAC char in DB
  --no-search-parent     Skip searching the parent directory
  --no-search-child      Skip searching child directories
  --dry-run              Show what would run without doing it
  -h, --help             Show this help message and exit

EXAMPLES:
  # Standard run
  mbx_primer_identifier.sh /path/to/FASTQ

  # Lower thresholds for hard-to-detect or partial primers
  mbx_primer_identifier.sh /path/to/FASTQ --mismatches 4 --min-rate 0.02 --offset 30

  # Add a custom primer not in the built-in DB
  mbx_primer_identifier.sh /path/to/FASTQ --extra-primers my_primers.tsv

  # Always emit a tentative best-candidate even when no primer passes threshold
  mbx_primer_identifier.sh /path/to/FASTQ --report-best

OUTPUT FILE FORMAT (mbx_primer_info.txt):
  READ_TYPE=paired
  DETECTION_STATUS=DETECTED          # DETECTED | TRIMMED | UNKNOWN
  CONFIDENCE_LEVEL=HIGH              # HIGH (≥80%) | MEDIUM (30-80%) | LOW (5-30%) | NONE
  INFERRED_REGION=V4                 # populated when STATUS=TRIMMED
  DETECTION_NOTE=Primer detected with high confidence.
  FORWARD_PRIMER_NAME=515F (Parada)
  FORWARD_PRIMER_SEQUENCE=GTGYCAGCMGCCGCGGTAA
  ...

  Backwards-compatible: every key the previous version emitted is still here.
  Downstream consumers (mbx_dada2_parameter_finder.sh, mbx_classifier_arranger.sh)
  read FORWARD_PRIMER_SEQUENCE / REVERSE_PRIMER_SEQUENCE — those keys are unchanged.

COMMON ERRORS:
  "python3 not found"          → activate your QIIME2 conda env
  "No .fastq[.gz] files found" → check the directory path
  "DETECTION_STATUS=TRIMMED"   → primers were already cut by the sequencing
                                 facility; pipeline can still proceed.
                                 Pass --no-trim or set trim-left=0 downstream.
  "DETECTION_STATUS=UNKNOWN"   → run again with --mismatches 4 --offset 30,
                                 or supply --forward-primer / --reverse-primer
                                 directly to mbx_dada2_parameter_finder.sh.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
FASTQ_DIR=""
N_SAMPLES=10000
MAX_MISMATCHES=3
MAX_OFFSET=25
MIN_MATCH_RATE="0.05"
SEARCH_PARENT=true
SEARCH_CHILD=true
DRY_RUN=false
EXTRA_PRIMERS=""
REPORT_BEST=false
DEBUG_CSV=false
STRICT=false
USER_FWD_PRIMER=""
USER_REV_PRIMER=""

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)             usage ;;
    --dry-run)             DRY_RUN=true; shift ;;
    --no-search-parent)    SEARCH_PARENT=false; shift ;;
    --no-search-child)     SEARCH_CHILD=false; shift ;;
    --report-best)         REPORT_BEST=true; shift ;;
    --debug-csv)           DEBUG_CSV=true; shift ;;
    --strict)              STRICT=true; shift ;;
    --samples)             N_SAMPLES="${2:?--samples needs a number}"; shift 2 ;;
    --mismatches)          MAX_MISMATCHES="${2:?--mismatches needs a number}"; shift 2 ;;
    --offset)              MAX_OFFSET="${2:?--offset needs a number}"; shift 2 ;;
    --min-rate)            MIN_MATCH_RATE="${2:?--min-rate needs a 0-1 fraction}"; shift 2 ;;
    --extra-primers)       EXTRA_PRIMERS="${2:?--extra-primers needs a TSV path}"; shift 2 ;;
    --forward-primer)      USER_FWD_PRIMER="${2:?--forward-primer needs a sequence}"; shift 2 ;;
    --reverse-primer)      USER_REV_PRIMER="${2:?--reverse-primer needs a sequence}"; shift 2 ;;
    -*)                    err "Unknown option: '$1'.  Run with --help." ;;
    *)
      if [[ -z "$FASTQ_DIR" ]]; then FASTQ_DIR="$1"
      else err "Unexpected extra argument: '$1'"; fi
      shift ;;
  esac
done

[[ -z "$FASTQ_DIR" ]] && err "No FASTQ directory provided.  Run with --help."
[[ -d "$FASTQ_DIR" ]] || err "Directory does not exist: '${FASTQ_DIR}'"

command -v python3 >/dev/null 2>&1 || err "python3 not found.
  → Activate your QIIME2 conda environment first:
    conda activate qiime2-amplicon-2025.4"

if [[ -n "$EXTRA_PRIMERS" ]]; then
  [[ -f "$EXTRA_PRIMERS" ]] || err "--extra-primers file not found: '$EXTRA_PRIMERS'"
  EXTRA_PRIMERS="$(_abspath "$EXTRA_PRIMERS")"
fi

FASTQ_DIR="$(_abspath "$FASTQ_DIR")"

# ── Output directory ──────────────────────────────────────────────────────────
# Standalone   : create timestamped sibling dir of FASTQ.
# Orchestrator : reuse MBX_OUT_DIR (single shared output dir for the whole run).
if [[ -n "${MBX_OUT_DIR:-}" ]]; then
  OUT_ROOT="$MBX_OUT_DIR"
  TIMESTAMP="$(basename "$OUT_ROOT" | sed -E 's/^mbX_pro_outputs_//')"
  [[ -z "$TIMESTAMP" ]] && TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
else
  TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
  OUT_ROOT="$(dirname "$FASTQ_DIR")/mbX_pro_outputs_${TIMESTAMP}"
fi
PRIMER_DIR="${OUT_ROOT}/0_primer_handling"
OUTPUT_TXT="${PRIMER_DIR}/mbx_primer_info.txt"
DEBUG_OUT="${PRIMER_DIR}/mbx_primer_debug_matrix.csv"

mkdir -p "$PRIMER_DIR" || err "Could not create: $PRIMER_DIR  —  permissions?"

sep
info "FASTQ directory      : $FASTQ_DIR"
info "Output directory     : $PRIMER_DIR"
if [[ -n "$USER_FWD_PRIMER" || -n "$USER_REV_PRIMER" ]]; then
  info "Mode                 : USER_SUPPLIED (skipping read-scan detection)"
  [[ -n "$USER_FWD_PRIMER" ]] && info "Forward primer       : $USER_FWD_PRIMER"
  [[ -n "$USER_REV_PRIMER" ]] && info "Reverse primer       : $USER_REV_PRIMER"
else
  info "Reads to sample      : $N_SAMPLES per orientation"
  info "Max IUPAC mismatches : $MAX_MISMATCHES"
  info "Max 5'-offset        : $MAX_OFFSET bp (handles linkers / frameshift Ns)"
  info "Min match rate       : $MIN_MATCH_RATE"
  [[ -n "$EXTRA_PRIMERS" ]] && info "Extra primers TSV    : $EXTRA_PRIMERS"
  $REPORT_BEST && info "Mode                 : --report-best (tentative below threshold)"
  $DEBUG_CSV   && info "Debug CSV            : will be written to mbx_primer_debug_matrix.csv"
fi
$DRY_RUN     && warn "DRY-RUN: detection skipped."
sep

# ── Collect candidate directories ────────────────────────────────────────────
WORK_DIR="$(mktemp -d -t mbx_primer_XXXXXX 2>/dev/null \
              || mktemp -d "/tmp/mbx_primer_${$}_$(date +%s)")"
trap 'rm -rf "$WORK_DIR"' EXIT

DIRS_FILE="$WORK_DIR/dirs.txt"
SEEN_FILE="$WORK_DIR/seen.txt"
touch "$DIRS_FILE" "$SEEN_FILE"

_add_dir() {
  local d
  d="$(_abspath "$1" 2>/dev/null)" || return 0
  [[ -d "$d" ]] || return 0
  grep -qxF "$d" "$DIRS_FILE" 2>/dev/null || echo "$d" >> "$DIRS_FILE"
}

_add_dir "$FASTQ_DIR"
if $SEARCH_PARENT; then
  _P="$(dirname "$FASTQ_DIR")"
  [[ "$_P" != "$FASTQ_DIR" ]] && _add_dir "$_P"
fi
if $SEARCH_CHILD; then
  while IFS= read -r -d '' c; do
    [[ -d "$c" ]] && _add_dir "$c"
  done < <(find "$FASTQ_DIR" -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null || true)
fi

# ── Collect FASTQ files (gzipped or uncompressed) ─────────────────────────────
FASTQ_LIST="$WORK_DIR/files.txt"
touch "$FASTQ_LIST"

while IFS= read -r d; do
  while IFS= read -r -d '' f; do
    fp="$(_abspath "$f" 2>/dev/null)" || continue
    grep -qxF "$fp" "$SEEN_FILE" 2>/dev/null && continue
    echo "$fp" >> "$SEEN_FILE"
    echo "$fp" >> "$FASTQ_LIST"
  done < <(find "$d" -maxdepth 1 \( -iname "*.fastq.gz" -o -iname "*.fq.gz" \
                                    -o -iname "*.fastq"   -o -iname "*.fq"    \) \
                      -type f -print0 2>/dev/null || true)
done < "$DIRS_FILE"

TOTAL="$(wc -l < "$FASTQ_LIST" | tr -d ' ')"
[[ "$TOTAL" -eq 0 ]] && err "No FASTQ files (.fastq[.gz] / .fq[.gz]) found.
  Searched:
$(while IFS= read -r d; do echo "    $d"; done < "$DIRS_FILE")"

info "Found $TOTAL FASTQ file(s)"

# ── Robust R1/R2 detection ────────────────────────────────────────────────────
# Accepts:  _R1_, _R2_, _R1.fastq, _R2.fastq, _1.fastq, _2.fastq,
#           .R1., .R2.,  -R1-, -R2-, plus uppercase/lowercase.
# Falls back to "_1." / "_2." only if no _R1/_R2 found.
R1_LIST="$WORK_DIR/r1.txt"; R2_LIST="$WORK_DIR/r2.txt"
SE_LIST="$WORK_DIR/se.txt"
touch "$R1_LIST" "$R2_LIST" "$SE_LIST"

while IFS= read -r f; do
  nm="$(basename "$f")"
  if echo "$nm" | grep -qiE '(^|[._-])R1([._-]|$)'; then echo "$f" >> "$R1_LIST"
  elif echo "$nm" | grep -qiE '(^|[._-])R2([._-]|$)'; then echo "$f" >> "$R2_LIST"
  else echo "$f" >> "$SE_LIST"
  fi
done < "$FASTQ_LIST"

R1_COUNT="$(wc -l < "$R1_LIST" | tr -d ' ')"
R2_COUNT="$(wc -l < "$R2_LIST" | tr -d ' ')"
SE_COUNT="$(wc -l < "$SE_LIST" | tr -d ' ')"

# Fallback: try _1/_2 if no _R1/_R2 detected
if [[ "$R1_COUNT" -eq 0 && "$R2_COUNT" -eq 0 && "$SE_COUNT" -gt 0 ]]; then
  > "$R1_LIST"; > "$R2_LIST"; > "$SE_LIST"
  while IFS= read -r f; do
    nm="$(basename "$f")"
    if   echo "$nm" | grep -qiE '(^|[._-])1([._-]|$)'; then echo "$f" >> "$R1_LIST"
    elif echo "$nm" | grep -qiE '(^|[._-])2([._-]|$)'; then echo "$f" >> "$R2_LIST"
    else echo "$f" >> "$SE_LIST"
    fi
  done < "$FASTQ_LIST"
  R1_COUNT="$(wc -l < "$R1_LIST" | tr -d ' ')"
  R2_COUNT="$(wc -l < "$R2_LIST" | tr -d ' ')"
  SE_COUNT="$(wc -l < "$SE_LIST" | tr -d ' ')"
fi

# Final fallback: treat unclassified files as single-end R1
if [[ "$R1_COUNT" -eq 0 && "$SE_COUNT" -gt 0 ]]; then
  cat "$SE_LIST" >> "$R1_LIST"
  R1_COUNT="$(wc -l < "$R1_LIST" | tr -d ' ')"
  SE_COUNT=0
fi

if [[ "$R2_COUNT" -gt 0 ]]; then READ_TYPE="paired"; else READ_TYPE="single"; fi

info "Read type            : $READ_TYPE"
info "R1 files             : $R1_COUNT"
[[ "$READ_TYPE" == "paired" ]] && info "R2 files             : $R2_COUNT"
[[ "$SE_COUNT" -gt 0 ]] && warn "Unclassified files: $SE_COUNT (could not detect R1/R2 pattern)"
sep

[[ "$R1_COUNT" -eq 0 ]] && err "No R1 files found.
  File names should contain _R1_ / _R2_ (Illumina) or _1 / _2 (SRA-style)."

if $DRY_RUN; then warn "Dry-run: skipping primer detection."; exit 0; fi

# ── Fast-path: user supplied primers manually ─────────────────────────────────
# When the orchestrator (or expert user) passes --forward-primer / --reverse-primer
# we trust them: skip read-scanning entirely and write a USER_SUPPLIED info file.
if [[ -n "$USER_FWD_PRIMER" || -n "$USER_REV_PRIMER" ]]; then
  info "User-supplied primers — skipping read-based detection."
  python3 - "$OUTPUT_TXT" "$READ_TYPE" "$FASTQ_DIR" "$USER_FWD_PRIMER" "$USER_REV_PRIMER" <<'PYUSER'
import datetime, sys, unicodedata
out_txt, read_type, fastq_dir, fwd_raw, rev_raw = sys.argv[1:6]
ALLOWED = set("ACGTUMRWSYKBDHVN")
HOMO = {"\u041D":"N","\u043D":"N","\u0410":"A","\u0430":"A","\u0421":"C","\u0441":"C","\u0420":"R","\u0440":"R"}
def clean(s, label):
    s = unicodedata.normalize("NFKC", s or "")
    out = []
    for ch in s:
        ch = HOMO.get(ch, ch)
        u = ch.upper()
        if u in ALLOWED: out.append(u)
        elif u and not u.isspace():
            print(f"[WARN]  Non-IUPAC '{ch}' (U+{ord(ch):04X}) in {label} — dropped.", file=sys.stderr)
    return "".join(out)

fwd = clean(fwd_raw, "--forward-primer")
rev = clean(rev_raw, "--reverse-primer") if read_type == "paired" else ""
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

L = [
    "# mbx_primer_info.txt",
    "# Generated by mbx_primer_identifier.sh v2.0  (USER_SUPPLIED mode)",
    f"# Date              : {now}",
    f"# FASTQ directory   : {fastq_dir}",
    "# Detection skipped — primers provided manually via --forward-primer / --reverse-primer.",
    "#",
    f"READ_TYPE={read_type}",
    "DETECTION_STATUS=USER_SUPPLIED",
    "CONFIDENCE_LEVEL=USER",
    "INFERRED_REGION=None",
    "DETECTION_NOTE=Primers were provided manually by the user; no read-scan was performed.",
    "",
]
if fwd:
    L += [
        f"FORWARD_PRIMER_NAME=user-supplied",
        f"FORWARD_PRIMER_SEQUENCE={fwd}",
        f"FORWARD_PRIMER_LENGTH={len(fwd)}",
        f"FORWARD_PRIMER_TARGET_REGION=Unknown",
        f"FORWARD_PRIMER_MATCH_RATE=N/A",
        f"FORWARD_PRIMER_MATCHED_READS=N/A",
        f"FORWARD_PRIMER_FOUND_AT_OFFSET=N/A",
    ]
else:
    L += ["FORWARD_PRIMER_NAME=None","FORWARD_PRIMER_SEQUENCE=None",
          "FORWARD_PRIMER_LENGTH=None","FORWARD_PRIMER_TARGET_REGION=None",
          "FORWARD_PRIMER_MATCH_RATE=None","FORWARD_PRIMER_MATCHED_READS=None",
          "FORWARD_PRIMER_FOUND_AT_OFFSET=None"]
L.append("")
if read_type == "paired":
    if rev:
        L += [
            f"REVERSE_PRIMER_NAME=user-supplied",
            f"REVERSE_PRIMER_SEQUENCE={rev}",
            f"REVERSE_PRIMER_LENGTH={len(rev)}",
            f"REVERSE_PRIMER_TARGET_REGION=Unknown",
            f"REVERSE_PRIMER_MATCH_RATE=N/A",
            f"REVERSE_PRIMER_MATCHED_READS=N/A",
            f"REVERSE_PRIMER_FOUND_AT_OFFSET=N/A",
        ]
    else:
        L += ["REVERSE_PRIMER_NAME=None","REVERSE_PRIMER_SEQUENCE=None",
              "REVERSE_PRIMER_LENGTH=None","REVERSE_PRIMER_TARGET_REGION=None",
              "REVERSE_PRIMER_MATCH_RATE=None","REVERSE_PRIMER_MATCHED_READS=None",
              "REVERSE_PRIMER_FOUND_AT_OFFSET=None"]
else:
    L.append("REVERSE_PRIMER_NAME=N/A (single-end)")

with open(out_txt, "w") as f:
    f.write("\n".join(L) + "\n")
print(f"[OK]    User-supplied primer info written → {out_txt}")
PYUSER
  sep
  ok "User-supplied primers recorded successfully."
  echo ""
  cat "$OUTPUT_TXT" | grep -v '^#' | grep -v '^$' | sed 's/^/  /'
  echo ""
  exit 0
fi

# ── Python detection engine ───────────────────────────────────────────────────
export MBX_R1_LIST="$R1_LIST"
export MBX_R2_LIST="$R2_LIST"
export MBX_READ_TYPE="$READ_TYPE"
export MBX_N_SAMPLES="$N_SAMPLES"
export MBX_MAX_MM="$MAX_MISMATCHES"
export MBX_MIN_RATE="$MIN_MATCH_RATE"
export MBX_MAX_OFFSET="$MAX_OFFSET"
export MBX_OUTPUT_TXT="$OUTPUT_TXT"
export MBX_FASTQ_DIR="$FASTQ_DIR"
export MBX_TIMESTAMP="$TIMESTAMP"
export MBX_EXTRA_PRIMERS="$EXTRA_PRIMERS"
export MBX_REPORT_BEST="$( $REPORT_BEST && echo 1 || echo 0 )"
export MBX_DEBUG_CSV="$( $DEBUG_CSV && echo "$DEBUG_OUT" || echo "" )"
export MBX_STRICT="$( $STRICT && echo 1 || echo 0 )"

# Write the python engine to a PID-based temp file so heredoc quirks never bite.
PYTMP="${WORK_DIR}/mbx_primer_engine.py"

cat > "$PYTMP" <<'PYEOF'
"""
mbx_primer_identifier — Python detection engine v2.0
====================================================

Three-tier detection:

  TIER 1  Direct primer match (IUPAC sliding window, 3 orientations).
  TIER 2  V-region anchor motifs — recognises that primers were *already
          trimmed* and infers which V-region was amplified.
  TIER 3  Failure with rich diagnostics + concrete next-steps.

Treats N as wild on **both** read and primer sides (canonical IUPAC).
Sanitises the primer database to strict ASCII IUPAC (kills the well-known
Cyrillic-Н-pasted-from-Excel bug class).
"""

from __future__ import annotations
import csv
import gzip
import os
import re
import sys
import unicodedata
from collections import Counter
from pathlib import Path
from typing import Iterator


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Sequence sanitisation + IUPAC machinery
# ─────────────────────────────────────────────────────────────────────────────
ALLOWED_IUPAC = set("ACGTUMRWSYKBDHVN")

# Cyrillic / look-alike ↔ Latin DNA letters that show up when sequences are
# pasted from Excel / Word / certain web pages.
HOMOGLYPHS = {
    "\u0410": "A", "\u0430": "A",  # Cyrillic A
    "\u0412": "B", "\u0432": "B",  # Cyrillic V (looks like B)
    "\u0421": "C", "\u0441": "C",  # Cyrillic S (looks like C)
    "\u0415": "E", "\u0435": "E",  # Cyrillic Ie
    "\u041D": "N", "\u043D": "N",  # ★ Cyrillic En (the bug we hit)
    "\u041E": "O", "\u043E": "O",
    "\u0420": "R", "\u0440": "R",  # Cyrillic Er
    "\u0422": "T", "\u0442": "T",
    "\u0425": "X", "\u0445": "X",
    "\u04A9": "C",                  # Coptic
    "\u00DF": "B",                  # German ß sometimes
    "\u2014": "-", "\u2013": "-",   # em / en dash
    "\u00A0": "",                   # non-breaking space
}

def sanitise_seq(seq: str, name: str = "?", strict: bool = False) -> str:
    """
    Normalise to NFKC, replace common homoglyphs, uppercase, strip non-IUPAC.
    Returns a clean IUPAC-only string.  Raises if STRICT and a non-IUPAC char
    is encountered.
    """
    s = unicodedata.normalize("NFKC", seq)
    out_chars = []
    for ch in s:
        # First map known homoglyphs
        if ch in HOMOGLYPHS:
            ch = HOMOGLYPHS[ch]
        if not ch:
            continue
        u = ch.upper()
        if u in ALLOWED_IUPAC:
            out_chars.append(u)
        else:
            msg = (f"Non-IUPAC char {ch!r} (U+{ord(ch):04X}) "
                   f"in primer '{name}' — dropped.")
            if strict:
                raise ValueError(msg)
            print(f"[WARN]  {msg}", file=sys.stderr)
    return "".join(out_chars)


IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),    "C": frozenset("C"),
    "G": frozenset("G"),    "T": frozenset("T"),
    "U": frozenset("T"),
    "M": frozenset("AC"),   "R": frozenset("AG"),
    "W": frozenset("AT"),   "S": frozenset("CG"),
    "Y": frozenset("CT"),   "K": frozenset("GT"),
    "B": frozenset("CGT"),  "D": frozenset("AGT"),
    "H": frozenset("ACT"),  "V": frozenset("ACG"),
    "N": frozenset("ACGTN"),
}

RC_TABLE = str.maketrans("ACGTUMRWSYKBDHVNacgtumrwsykbdhvn",
                          "TGCAAKYWSRMVHDBNtgcaakywsrmvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def iupac_mismatches(read: str, primer: str) -> int:
    """
    Count positions where the read base is NOT in the IUPAC expansion of the
    primer base.  Treats N (read or primer) as wild — matches anything.
    """
    mm = 0
    for rb, pb in zip(read, primer):
        rb_u = rb.upper()
        pb_u = pb.upper()
        if rb_u == "N" or pb_u == "N":
            continue
        allowed = IUPAC.get(pb_u, frozenset(pb_u))
        if rb_u not in allowed:
            mm += 1
    return mm


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Primer database (all sequences ASCII IUPAC, sanitised at load)
# ─────────────────────────────────────────────────────────────────────────────
# Each entry: (name, seq, region, dir)   dir = "F" (forward) or "R" (reverse)
PRIMER_DB_RAW: list[dict] = [
    # ── V1-V2 ───────────────────────────────────────────────────────────────
    {"name": "27F",                "seq": "AGAGTTTGATCMTGGCTCAG",     "region": "V1-V2/V1-V3/V1-V9", "dir": "F"},
    {"name": "8F",                 "seq": "AGRGTTTGATYMTGGCTYAG",     "region": "V1-V2/V1-V3",       "dir": "F"},
    {"name": "338R",               "seq": "TGCTGCCTCCCGTAGGAGT",      "region": "V1-V2",             "dir": "R"},
    {"name": "341R",               "seq": "CTGCWGCCNCCCGTAGG",        "region": "V1-V2",             "dir": "R"},

    # ── V1-V3 ───────────────────────────────────────────────────────────────
    {"name": "519R",               "seq": "GWATTACCGCGGCKGCTG",       "region": "V1-V3",             "dir": "R"},
    {"name": "534R",               "seq": "ATTACCGCGGCTGCTGG",        "region": "V1-V3",             "dir": "R"},
    {"name": "534R (8F-pair)",     "seq": "TBACCGCGGCTGCTGGCAC",      "region": "V1-V3",             "dir": "R"},
    {"name": "518R",               "seq": "WTTACCGCGGCTGCTGG",        "region": "V1-V3",             "dir": "R"},
    {"name": "68F",                "seq": "TNANACATGCAAGTCGRRCG",     "region": "V1-V3",             "dir": "F"},

    # ── V1-V9 (full length) ─────────────────────────────────────────────────
    {"name": "1492R",              "seq": "TACGGYTACCTTGTTACGACTT",   "region": "V1-V9",             "dir": "R"},
    {"name": "1492R (short)",      "seq": "GGYTACCTTGTTACGACTT",      "region": "V1-V9",             "dir": "R"},

    # ── V3 only ─────────────────────────────────────────────────────────────
    {"name": "338F",               "seq": "ACWCCTACGGGNGGCWG",        "region": "V3",                "dir": "F"},
    {"name": "537R",               "seq": "GWNTACCGCGGCKGCT",         "region": "V3",                "dir": "R"},

    # ── V3-V4 (most common Illumina pair) ───────────────────────────────────
    # 341F (Klindworth 2013 / Bakt_341F).  785R and 805R are aliases.
    {"name": "341F",               "seq": "CCTACGGGNGGCWGCAG",        "region": "V3-V4",             "dir": "F"},
    {"name": "341F (Klindworth)",  "seq": "CCTACGGGNBGCASCAG",        "region": "V3-V4",             "dir": "F"},
    {"name": "Bakt_341F",          "seq": "CCTAYGGGRBGCASCAG",        "region": "V3-V4",             "dir": "F"},
    {"name": "785R",               "seq": "GACTACHVGGGTATCTAATCC",    "region": "V3-V4",             "dir": "R"},
    {"name": "805R",               "seq": "GACTACHVGGGTATCTAATCC",    "region": "V3-V4",             "dir": "R"},
    {"name": "Bakt_805R",          "seq": "GACTACNVGGGTATCTAATCC",    "region": "V3-V4",             "dir": "R"},
    # PRK degenerate universal prokaryote pair (Yu 2005)
    {"name": "PRK341F",            "seq": "CCTAYGGGRBGCASCAG",        "region": "V3-V4",             "dir": "F"},
    {"name": "PRK806R",            "seq": "GGACTACNNGGGTATCTAAT",     "region": "V3-V4",             "dir": "R"},

    # ── V3-V5 ───────────────────────────────────────────────────────────────
    {"name": "341F (older)",       "seq": "CCTACGGGAGGCAGCAG",        "region": "V3-V5",             "dir": "F"},
    {"name": "926R (V3-V5)",       "seq": "CCGTCAATTCMTTTGAGTTT",     "region": "V3-V5",             "dir": "R"},

    # ── V4 (EMP classic) ────────────────────────────────────────────────────
    {"name": "515F (Caporaso)",    "seq": "GTGCCAGCMGCCGCGGTAA",      "region": "V4",                "dir": "F"},
    {"name": "806R (Caporaso)",    "seq": "GGACTACHVGGGTWTCTAAT",     "region": "V4",                "dir": "R"},
    # EMP UPDATED — Parada 515F + Apprill 806R (current EMP)
    {"name": "515F (Parada)",      "seq": "GTGYCAGCMGCCGCGGTAA",      "region": "V4",                "dir": "F"},
    {"name": "806R (Apprill)",     "seq": "GGACTACNVGGGTWTCTAAT",     "region": "V4",                "dir": "R"},

    # ── V4-V5 ───────────────────────────────────────────────────────────────
    {"name": "926R (Quince V4-V5)","seq": "CCGYCAATTYMTTTRAGTTT",     "region": "V4-V5",             "dir": "R"},

    # ── V5-V7 (chloroplast-skipping, plant friendly) ────────────────────────
    {"name": "799F",               "seq": "AACMGGATTAGATACCCKG",      "region": "V5-V7",             "dir": "F"},
    {"name": "1193R",              "seq": "ACGTCATCCCCACCTTCC",       "region": "V5-V7",             "dir": "R"},
    {"name": "1391R",              "seq": "GACGGGCGGTGWGTRCA",        "region": "V5-V7",             "dir": "R"},

    # ── V6-V7 ───────────────────────────────────────────────────────────────
    {"name": "967F",               "seq": "CAACGCGAAGAACCTTACC",      "region": "V6-V7",             "dir": "F"},

    # ── V6-V8 ───────────────────────────────────────────────────────────────
    {"name": "926F",               "seq": "AAACTYAAAKGAATTGACGG",     "region": "V6-V8",             "dir": "F"},
    {"name": "1392R",              "seq": "ACGGGCGGTGTGTRC",          "region": "V6-V8",             "dir": "R"},
    {"name": "B969F",              "seq": "ACGCGHNRAACCTTACC",        "region": "V6-V8",             "dir": "F"},
    {"name": "BA1406R",            "seq": "CGACRRCATGCANCACCT",       "region": "V6-V8",             "dir": "R"},
    {"name": "BA1406R (v2)",       "seq": "CGACRRCCATGCANCACCT",      "region": "V6-V8",             "dir": "R"},
    {"name": "968F",               "seq": "AACGCGAAGAACCTTAC",        "region": "V6-V8",             "dir": "F"},
    {"name": "1401R",              "seq": "CGGTGTGTACAAGACCC",        "region": "V6-V8",             "dir": "R"},

    # ── Nadkarni 2002 universal (qPCR-style, V3 region) ─────────────────────
    {"name": "Nadkarni-F",         "seq": "TCCTACGGGAGGCAGCAGT",      "region": "V3",                "dir": "F"},
    {"name": "Nadkarni-R",         "seq": "GGACTACCAGGGTATCTAATCCTGTT","region": "V3-V4",            "dir": "R"},
]


def load_primer_db(extra_path: str = "", strict: bool = False) -> list[dict]:
    """Sanitise built-ins, then merge user-supplied TSV (name<TAB>seq<TAB>region<TAB>F|R)."""
    db: list[dict] = []
    seen_keys: set[tuple[str, str]] = set()
    for p in PRIMER_DB_RAW:
        clean = sanitise_seq(p["seq"], name=p["name"], strict=strict)
        if not clean:
            print(f"[WARN]  Primer {p['name']!r} has empty cleaned seq — skipping.", file=sys.stderr)
            continue
        key = (clean, p["dir"])
        if key in seen_keys:
            continue
        seen_keys.add(key)
        db.append({"name": p["name"], "seq": clean,
                   "region": p["region"], "dir": p["dir"], "source": "builtin"})
    if extra_path:
        try:
            with open(extra_path) as fh:
                for raw in fh:
                    raw = raw.rstrip("\n")
                    if not raw or raw.startswith("#"): continue
                    parts = raw.split("\t")
                    if len(parts) < 4:
                        print(f"[WARN]  --extra-primers row malformed (need 4 TSV cols): {raw!r}",
                              file=sys.stderr)
                        continue
                    name, seq, region, direction = parts[0], parts[1], parts[2], parts[3].strip().upper()
                    if direction not in ("F", "R"):
                        print(f"[WARN]  --extra-primers dir must be F or R: {raw!r}", file=sys.stderr)
                        continue
                    clean = sanitise_seq(seq, name=name, strict=strict)
                    if not clean:
                        continue
                    key = (clean, direction)
                    if key in seen_keys: continue
                    seen_keys.add(key)
                    db.append({"name": f"{name} (user)", "seq": clean,
                               "region": region, "dir": direction, "source": "user"})
        except Exception as e:
            print(f"[WARN]  Could not read --extra-primers file: {e}", file=sys.stderr)
    return db


# ─────────────────────────────────────────────────────────────────────────────
# 3.  V-region anchor motifs (for "primers were already trimmed" detection)
# ─────────────────────────────────────────────────────────────────────────────
# Each anchor is a short conserved IUPAC motif that should appear at offset
# 0–2 of *trimmed* reads.  When TIER-1 primer detection fails but an anchor
# matches strongly, the script reports DETECTION_STATUS=TRIMMED.
#
# Sources: Lane 1991, Pruesse 2007 SILVA seed alignment, manual inspection
# of E. coli K-12 16S (J01695) at primer-landing sites.
V_REGION_ANCHORS_R1: list[dict] = [
    {"name": "V4_5p_after_515F",
     "motif": "TACGNAGG",
     "region": "V4",
     "expected_primers": "515F-Parada / 806R-Apprill (most likely) or 515F-Caporaso / 806R-Caporaso"},
    {"name": "V4_5p_alt",
     "motif": "TACGTAGG",
     "region": "V4",
     "expected_primers": "515F / 806R family"},
    {"name": "V3_5p_after_341F",
     "motif": "NGGCWGCAG",
     "region": "V3-V4",
     "expected_primers": "341F / 805R (Bakt) family"},
    {"name": "V1_5p_after_27F",
     "motif": "GGCGGCRTGCYTAA",
     "region": "V1-V2 or V1-V3",
     "expected_primers": "27F / 338R or 27F / 519R"},
    {"name": "V5_5p_after_799F",
     "motif": "TAGATACCC",
     "region": "V5-V7",
     "expected_primers": "799F / 1193R (chloroplast-skipping)"},
]
# R2 reads after reverse primer trimming run 3'→5' along 16S; anchors below
# are the *first bases of R2* (i.e. immediately downstream of the rev primer
# landing site, on the rev-comp strand).
V_REGION_ANCHORS_R2: list[dict] = [
    {"name": "V4_3p_after_806R",
     "motif": "TAGATACCC",        # rev-comp of common bases 5'-of-806R
     "region": "V4",
     "expected_primers": "*-806R reverse primer"},
    {"name": "V3-V4_3p_after_805R",
     "motif": "TAGATACCCNGGTAGTCC",
     "region": "V3-V4",
     "expected_primers": "*-805R reverse primer"},
]


def anchor_match(read: str, motif: str, max_mm: int = 1) -> bool:
    """Anchored IUPAC match at offset 0 with up to max_mm mismatches."""
    if len(read) < len(motif):
        return False
    return iupac_mismatches(read[:len(motif)], motif) <= max_mm


def detect_trimmed_region(reads: list[str], anchors: list[dict],
                          min_rate: float = 0.30) -> dict | None:
    """Return best-matching anchor, or None."""
    best = None
    for anc in anchors:
        hits = sum(1 for r in reads if anchor_match(r, anc["motif"], max_mm=1))
        rate = hits / len(reads) if reads else 0.0
        if rate >= min_rate and (best is None or rate > best["rate"]):
            best = {**anc, "rate": rate, "matched": hits, "total": len(reads)}
    return best


# ─────────────────────────────────────────────────────────────────────────────
# 4.  FASTQ reading + matching engine
# ─────────────────────────────────────────────────────────────────────────────
def open_fastq(path: str):
    """Transparently open .fastq or .fastq.gz."""
    if path.endswith((".gz", ".GZ")):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "r", errors="replace")


def reads_from_file(path: str, n: int) -> Iterator[str]:
    """Yield up to n sequence strings from a (possibly gzipped) FASTQ file."""
    count = 0
    try:
        fh = open_fastq(path)
    except Exception as e:
        print(f"[WARN]  Could not open {path}: {e}", file=sys.stderr)
        return
    with fh:
        while count < n:
            header = fh.readline()
            if not header: break
            seq = fh.readline().strip()
            fh.readline(); fh.readline()
            if seq:
                yield seq
                count += 1


def load_reads(file_list_path: str, n: int) -> list[str]:
    """Load up to n reads spread across all listed files (round-robin-ish)."""
    files = [f for f in Path(file_list_path).read_text().splitlines() if f.strip()]
    if not files: return []
    per_file = max(1, n // len(files))
    out: list[str] = []
    for f in files:
        out.extend(reads_from_file(f, per_file))
        if len(out) >= n: break
    return out[:n]


def slide_match(read: str, primer: str, max_mm: int, max_offset: int) -> tuple[bool, int, int]:
    """
    Try every offset in [0..max_offset], return (matched, best_offset, best_mm).
    """
    plen = len(primer)
    best_mm = plen
    best_off = -1
    for offset in range(max_offset + 1):
        window = read[offset : offset + plen]
        if len(window) < plen: break
        mm = iupac_mismatches(window, primer)
        if mm < best_mm:
            best_mm = mm; best_off = offset
        if mm == 0: break  # short-circuit
    matched = best_mm <= max_mm
    return matched, (best_off if matched else -1), best_mm


def score_primers(reads: list[str], primers: list[dict],
                  max_mm: int, max_offset: int, as_rc: bool = False) -> dict[str, dict]:
    """
    For every primer (or its rev-comp if as_rc), count matches across reads.
    Returns dict keyed by composite name "<name>|<dir>".
    """
    results: dict[str, dict] = {}
    n = len(reads) if reads else 0
    for p in primers:
        seq = revcomp(p["seq"]) if as_rc else p["seq"]
        plen = len(seq)
        if plen == 0: continue
        matched = 0
        offsets: Counter = Counter()
        mm_sum = 0
        for read in reads:
            hit, off, mm = slide_match(read, seq, max_mm, max_offset)
            if hit:
                matched += 1
                offsets[off] += 1
                mm_sum += mm
        rate = matched / n if n > 0 else 0.0
        modal_offset = offsets.most_common(1)[0][0] if offsets else -1
        avg_mm = (mm_sum / matched) if matched > 0 else 0.0
        results[f"{p['name']}|{p['dir']}"] = {
            "name": p["name"], "seq": p["seq"], "region": p["region"],
            "dir": p["dir"], "source": p.get("source", "builtin"),
            "matched": matched, "total": n, "rate": rate,
            "modal_offset": modal_offset, "avg_mismatch": avg_mm,
            "as_rc": as_rc,
        }
    return results


def best_primer(scores: dict[str, dict], direction: str,
                min_rate: float) -> dict | None:
    """Highest-rate primer in `direction` whose rate >= min_rate, else None."""
    cands = [v for v in scores.values() if v["dir"] == direction and v["rate"] >= min_rate]
    if not cands: return None
    return max(cands, key=lambda x: x["rate"])


def best_primer_any(scores: dict[str, dict], direction: str) -> dict | None:
    """Highest-rate primer (no threshold) — used for --report-best."""
    cands = [v for v in scores.values() if v["dir"] == direction]
    if not cands: return None
    return max(cands, key=lambda x: x["rate"])


def consensus_bases(reads: list[str], length: int) -> str:
    out = []
    for i in range(length):
        bases = [r[i].upper() for r in reads if len(r) > i and r[i].upper() in "ACGTN"]
        if not bases: out.append("?"); continue
        out.append(Counter(bases).most_common(1)[0][0])
    return "".join(out)


def position_base_table(reads: list[str], length: int) -> str:
    """Pretty per-position base composition table (TIER-3 diagnostic)."""
    lines = ["pos   A%    C%    G%    T%    N%   most"]
    for i in range(length):
        bases = [r[i].upper() for r in reads if len(r) > i]
        n = len(bases) or 1
        c = Counter(bases)
        a, cn, g, t, nn = c.get("A", 0), c.get("C", 0), c.get("G", 0), c.get("T", 0), c.get("N", 0)
        top = c.most_common(1)[0][0] if bases else "?"
        lines.append(f"{i:3d}  {100*a/n:5.1f} {100*cn/n:5.1f} "
                     f"{100*g/n:5.1f} {100*t/n:5.1f} {100*nn/n:5.1f}   {top}")
    return "\n".join(lines)


def confidence_label(rate: float) -> str:
    if rate >= 0.80: return "HIGH"
    if rate >= 0.30: return "MEDIUM"
    if rate >= 0.05: return "LOW"
    return "NONE"


def write_debug_csv(path: str, all_scores: dict[str, dict[str, dict]]) -> None:
    """Dump every primer × orientation rate to CSV."""
    cols = ["orientation", "name", "dir", "region", "matched", "total",
            "rate", "modal_offset", "avg_mismatch", "as_rc"]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for ori, scores in all_scores.items():
            for v in scores.values():
                w.writerow([ori, v["name"], v["dir"], v["region"],
                            v["matched"], v["total"], f"{v['rate']:.4f}",
                            v["modal_offset"], f"{v['avg_mismatch']:.2f}",
                            int(v["as_rc"])])


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Top-level driver
# ─────────────────────────────────────────────────────────────────────────────
def banner(s: str) -> None:
    print(f"\n  ── {s} ──")


def fmt_top(scores: dict[str, dict], direction: str, n: int = 5) -> str:
    rows = sorted([v for v in scores.values() if v["dir"] == direction],
                  key=lambda x: x["rate"], reverse=True)[:n]
    if not rows: return "  (no candidates)"
    out = []
    for c in rows:
        bar = "█" * int(c["rate"] * 40)
        off = f"offset={c['modal_offset']}" if c["modal_offset"] >= 0 else ""
        out.append(f"  {c['name']:<28} {c['rate']*100:5.1f}%  {bar:<40}  {off}")
    return "\n".join(out)


def fmt_primer(label: str, best: dict | None, fallback_note: str) -> list[str]:
    out = []
    if best:
        out.append(f"{label}_NAME={best['name']}")
        out.append(f"{label}_SEQUENCE={best['seq']}")
        out.append(f"{label}_LENGTH={len(best['seq'])}")
        out.append(f"{label}_TARGET_REGION={best['region']}")
        out.append(f"{label}_MATCH_RATE={best['rate']:.4f}")
        out.append(f"{label}_MATCHED_READS={best['matched']}/{best['total']}")
        out.append(f"{label}_FOUND_AT_OFFSET={best['modal_offset']}")
    else:
        out.append(f"# {fallback_note}")
        out.append(f"{label}_NAME=None")
        out.append(f"{label}_SEQUENCE=None")
        out.append(f"{label}_LENGTH=None")
        out.append(f"{label}_TARGET_REGION=None")
        out.append(f"{label}_MATCH_RATE=None")
        out.append(f"{label}_MATCHED_READS=None")
        out.append(f"{label}_FOUND_AT_OFFSET=None")
    return out


def main() -> None:
    r1_list    = os.environ["MBX_R1_LIST"]
    r2_list    = os.environ["MBX_R2_LIST"]
    read_type  = os.environ["MBX_READ_TYPE"]
    n_samples  = int(os.environ["MBX_N_SAMPLES"])
    max_mm     = int(os.environ["MBX_MAX_MM"])
    min_rate   = float(os.environ["MBX_MIN_RATE"])
    max_offset = int(os.environ["MBX_MAX_OFFSET"])
    out_txt    = os.environ["MBX_OUTPUT_TXT"]
    fastq_dir  = os.environ["MBX_FASTQ_DIR"]
    timestamp  = os.environ["MBX_TIMESTAMP"]
    extra_path = os.environ.get("MBX_EXTRA_PRIMERS", "")
    report_best = os.environ.get("MBX_REPORT_BEST", "0") == "1"
    debug_csv  = os.environ.get("MBX_DEBUG_CSV", "")
    strict     = os.environ.get("MBX_STRICT", "0") == "1"

    # 5.1  Load + sanitise primer DB
    db = load_primer_db(extra_path, strict=strict)
    n_builtin = sum(1 for p in db if p["source"] == "builtin")
    n_user    = sum(1 for p in db if p["source"] == "user")
    min_plen = min(len(p["seq"]) for p in db) if db else 0
    max_plen = max(len(p["seq"]) for p in db) if db else 0
    print(f"[INFO]  Primer DB loaded: {len(db)} entries  ({n_builtin} built-in + {n_user} user)")
    print(f"[INFO]  Primer length range: {min_plen}–{max_plen} nt")
    print(f"[INFO]  Sliding window offset: 0..{max_offset} bp  |  Max mismatches: {max_mm}")

    # 5.2  Load reads
    print(f"[INFO]  Sampling up to {n_samples} reads per orientation...")
    r1_reads = load_reads(r1_list, n_samples)
    print(f"[INFO]  R1 reads sampled: {len(r1_reads)}")
    if not r1_reads:
        print("[ERROR] No reads could be read from R1 files.", file=sys.stderr)
        sys.exit(2)

    r2_reads: list[str] = []
    if read_type == "paired":
        r2_reads = load_reads(r2_list, n_samples)
        print(f"[INFO]  R2 reads sampled: {len(r2_reads)}")

    # Read length distribution (used in TIER-3)
    r1_lens = [len(r) for r in r1_reads]
    r2_lens = [len(r) for r in r2_reads]
    print(f"[INFO]  R1 length: median={sorted(r1_lens)[len(r1_lens)//2]} "
          f"min={min(r1_lens)} max={max(r1_lens)}")
    if r2_lens:
        print(f"[INFO]  R2 length: median={sorted(r2_lens)[len(r2_lens)//2]} "
              f"min={min(r2_lens)} max={max(r2_lens)}")

    # 5.3  Diagnostic consensus (search-window-wide)
    diag_len = max_plen + max_offset
    r1_cons = consensus_bases(r1_reads, diag_len)
    print(f"[DIAG]  R1 read start consensus ({diag_len}bp window): {r1_cons}")
    r2_cons = ""
    if r2_reads:
        r2_cons = consensus_bases(r2_reads, diag_len)
        print(f"[DIAG]  R2 read start consensus ({diag_len}bp window): {r2_cons}")

    # 5.4  TIER 1 — direct primer detection at all 3 orientations
    r1_fwd = score_primers(r1_reads, db, max_mm, max_offset, as_rc=False)
    r1_rc  = score_primers(r1_reads, db, max_mm, max_offset, as_rc=True)
    r2_fwd: dict[str, dict] = {}
    r2_rc:  dict[str, dict] = {}
    if r2_reads:
        r2_fwd = score_primers(r2_reads, db, max_mm, max_offset, as_rc=False)
        r2_rc  = score_primers(r2_reads, db, max_mm, max_offset, as_rc=True)

    if debug_csv:
        all_scores = {
            "R1_fwd": r1_fwd, "R1_rc": r1_rc,
            "R2_fwd": r2_fwd, "R2_rc": r2_rc,
        }
        write_debug_csv(debug_csv, all_scores)
        print(f"[DEBUG] Wrote debug matrix → {debug_csv}")

    # Three orientation strategies (same as before, but cleaner)
    fwd_A = best_primer(r1_fwd, "F", min_rate)
    rev_A = best_primer(r2_rc,  "R", min_rate) if r2_rc  else None
    rev_B = best_primer(r2_fwd, "R", min_rate) if r2_fwd else None
    fwd_C = best_primer(r2_fwd, "F", min_rate) if r2_fwd else None
    rev_C = best_primer(r1_rc,  "R", min_rate)

    rate = lambda p: p["rate"] if p else 0.0
    rate_A = rate(fwd_A) + rate(rev_A)
    rate_B = rate(fwd_A) + rate(rev_B)
    rate_C = rate(fwd_C) + rate(rev_C)

    if max(rate_A, rate_B, rate_C) <= 0:
        fwd_best, rev_best = fwd_A, rev_A
        orientation = "RC (R1=forward, R2=RC of reverse)"
    elif rate_B >= rate_A and rate_B >= rate_C:
        fwd_best, rev_best = fwd_A, rev_B
        orientation = "direct (R1=forward 5'→3', R2=reverse 5'→3' WITHOUT rev-comp)"
    elif rate_C > rate_A and rate_C > rate_B:
        fwd_best, rev_best = fwd_C, rev_C
        orientation = "swapped (R1 contains REVERSE primer, R2 contains FORWARD primer)"
    else:
        fwd_best, rev_best = fwd_A, rev_A
        orientation = "RC (R1=forward, R2=RC of reverse)"

    # 5.5  Decide DETECTION_STATUS
    fwd_rate = rate(fwd_best)
    rev_rate = rate(rev_best) if read_type == "paired" else None
    have_fwd = fwd_best is not None
    have_rev = (read_type != "paired") or (rev_best is not None)

    detection_status = "UNKNOWN"
    inferred_region  = ""
    detection_note   = ""

    if have_fwd and have_rev:
        detection_status = "DETECTED"
        confidence = confidence_label(min(fwd_rate, rev_rate if rev_rate is not None else fwd_rate))
        detection_note = f"Primer pair detected with {confidence.lower()} confidence."
    else:
        # 5.6  TIER 2 — V-region anchor detection
        anchor_r1 = detect_trimmed_region(r1_reads, V_REGION_ANCHORS_R1, min_rate=0.30)
        anchor_r2 = detect_trimmed_region(r2_reads, V_REGION_ANCHORS_R2, min_rate=0.30) if r2_reads else None
        if anchor_r1 or anchor_r2:
            detection_status = "TRIMMED"
            anc = anchor_r1 or anchor_r2
            inferred_region = anc["region"]
            confidence = confidence_label(anc["rate"])
            detection_note = (f"Primers appear ALREADY TRIMMED. "
                              f"Reads start with conserved 16S anchor "
                              f"'{anc['motif']}' ({anc['rate']*100:.1f}% of reads). "
                              f"Inferred region: {inferred_region}. "
                              f"Likely original primers: {anc['expected_primers']}.")
        else:
            detection_status = "UNKNOWN"
            confidence = "NONE"
            detection_note = ("No primers above threshold and no V-region anchor matched. "
                              "Reads may be from a non-16S library, an unsupported primer set, "
                              "or have heavy 5' contamination. See diagnostics below.")

    # ── Print human-readable summary ────────────────────────────────────────
    print()
    banner("Top 5 forward primer candidates (R1, forward)")
    print(fmt_top(r1_fwd, "F", 5))
    if read_type == "paired":
        print()
        if orientation.startswith("direct"):
            banner("Top 5 reverse primer candidates (R2, FORWARD orientation — direct)")
            print(fmt_top(r2_fwd, "R", 5))
        elif orientation.startswith("swapped"):
            banner("Top 5 forward primer candidates (R2)")
            print(fmt_top(r2_fwd, "F", 5))
            banner("Top 5 reverse primer candidates (R1, rev-comp)")
            print(fmt_top(r1_rc, "R", 5))
        else:
            banner("Top 5 reverse primer candidates (R2, rev-comp)")
            print(fmt_top(r2_rc, "R", 5))

    print()
    print(f"  Orientation tested:  {orientation}")
    print(f"  Detection status:    {detection_status}")
    print(f"  Confidence:          {confidence}")
    if inferred_region:
        print(f"  Inferred region:     {inferred_region}")
    print(f"  Note:                {detection_note}")
    print()

    # ── TIER 3 diagnostics — only when status != DETECTED ───────────────────
    diag_lines: list[str] = []
    if detection_status != "DETECTED":
        diag_lines.append("# ── TIER-3 DIAGNOSTICS (no high-confidence primer pair) ─────────")
        diag_lines.append(f"# R1 reads sampled: {len(r1_reads)}")
        diag_lines.append(f"# R1 length      : median={sorted(r1_lens)[len(r1_lens)//2]} "
                          f"min={min(r1_lens)} max={max(r1_lens)}")
        diag_lines.append(f"# R1 consensus   : {r1_cons}")
        if r2_reads:
            diag_lines.append(f"# R2 length      : median={sorted(r2_lens)[len(r2_lens)//2]} "
                              f"min={min(r2_lens)} max={max(r2_lens)}")
            diag_lines.append(f"# R2 consensus   : {r2_cons}")
        diag_lines.append("# Per-position base composition (R1, first 30 bp):")
        for ln in position_base_table(r1_reads, 30).splitlines():
            diag_lines.append(f"#   {ln}")
        # Top 3 most frequent 5'-prefixes (10 bp) — useful when reads are
        # primer-trimmed or contaminated by linkers.
        prefix_counter = Counter(r[:10] for r in r1_reads if len(r) >= 10)
        top_pref = prefix_counter.most_common(3)
        diag_lines.append("# Top 3 R1 5'-prefixes (10 bp):")
        for p, c in top_pref:
            diag_lines.append(f"#   {p}  {c}/{len(r1_reads)}  ({100*c/len(r1_reads):.1f}%)")
        diag_lines.append("")

    if diag_lines:
        for ln in diag_lines:
            print(ln)

    # ── Build output ────────────────────────────────────────────────────────
    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        f"# mbx_primer_info.txt",
        f"# Generated by mbx_primer_identifier.sh v2.0",
        f"# Date              : {now}",
        f"# FASTQ directory   : {fastq_dir}",
        f"# Reads sampled     : R1={len(r1_reads)}" + (f"  R2={len(r2_reads)}" if r2_reads else ""),
        f"# Max mismatches    : {max_mm}",
        f"# Max 5'-offset     : {max_offset}",
        f"# Min match rate    : {min_rate}",
        f"# Orientation       : {orientation}",
        f"# R1 consensus      : {r1_cons}",
    ]
    if r2_cons:
        lines.append(f"# R2 consensus      : {r2_cons}")
    lines += [
        f"#",
        f"# CONSUMERS:",
        f"#   mbx_dada2_parameter_finder.sh  (reads FORWARD/REVERSE_PRIMER_SEQUENCE)",
        f"#   mbx_classifier_arranger.sh     (reads FORWARD/REVERSE_PRIMER_SEQUENCE)",
        f"#   mbx_final_report.sh            (reads everything for the report)",
        f"",
        f"READ_TYPE={read_type}",
        f"DETECTION_STATUS={detection_status}",
        f"CONFIDENCE_LEVEL={confidence}",
        f"INFERRED_REGION={inferred_region or 'None'}",
        f"DETECTION_NOTE={detection_note}",
        f"",
    ]

    if detection_status == "DETECTED":
        lines += fmt_primer("FORWARD_PRIMER", fwd_best, "Forward primer detected.")
        lines.append("")
        if read_type == "paired":
            lines += fmt_primer("REVERSE_PRIMER", rev_best, "Reverse primer detected.")
        else:
            lines.append("# Single-end data — no reverse primer.")
            lines.append("REVERSE_PRIMER_NAME=N/A (single-end)")

    elif detection_status == "TRIMMED":
        # Output the most likely original primers for traceability, but flag
        # that they are NOT in the reads any more.
        lines.append("# Primers were already trimmed off by the sequencing facility.")
        lines.append("# The likely ORIGINAL primers (NOT present in reads) are listed below.")
        lines.append("# DOWNSTREAM HINT: pass trim-left=0 / --no-trim-primers to DADA2.")
        # If --report-best, fill in the highest-rate (likely <5%) primer; else None.
        rb_fwd = best_primer_any(r1_fwd, "F") if report_best else None
        rb_rev = (best_primer_any(r2_fwd, "R") if orientation.startswith("direct")
                  else best_primer_any(r2_rc, "R")) if (read_type == "paired" and report_best) else None
        lines += fmt_primer("FORWARD_PRIMER", rb_fwd,
                            "No forward primer in reads (already trimmed).")
        lines.append("")
        if read_type == "paired":
            lines += fmt_primer("REVERSE_PRIMER", rb_rev,
                                "No reverse primer in reads (already trimmed).")
        else:
            lines.append("# Single-end data.")
            lines.append("REVERSE_PRIMER_NAME=N/A (single-end)")

    else:  # UNKNOWN
        lines.append("# Primers could NOT be identified above the match-rate threshold,")
        lines.append("# and no V-region anchor matched.  This usually means one of:")
        lines.append("#   (1) Non-16S library (ITS, 18S, or shotgun) — not yet supported.")
        lines.append("#   (2) Unusual primer set not in the built-in DB.")
        lines.append("#   (3) Heavy 5'-end contamination (e.g. very long custom linker).")
        lines.append("# Suggested next steps:")
        lines.append("#   - Re-run with: --mismatches 4 --offset 35 --min-rate 0.02")
        lines.append("#   - Provide custom primers via:    --extra-primers my_primers.tsv")
        lines.append("#   - Or skip detection and supply primers directly to:")
        lines.append("#       mbx_dada2_parameter_finder.sh --forward-primer ... --reverse-primer ...")
        rb_fwd = best_primer_any(r1_fwd, "F") if report_best else None
        rb_rev = (best_primer_any(r2_fwd, "R") if orientation.startswith("direct")
                  else best_primer_any(r2_rc, "R")) if (read_type == "paired" and report_best) else None
        lines += fmt_primer("FORWARD_PRIMER", rb_fwd, "No forward primer detected.")
        lines.append("")
        if read_type == "paired":
            lines += fmt_primer("REVERSE_PRIMER", rb_rev, "No reverse primer detected.")
        else:
            lines.append("# Single-end data.")
            lines.append("REVERSE_PRIMER_NAME=N/A (single-end)")

    if diag_lines:
        lines.append("")
        lines.extend(diag_lines)

    Path(out_txt).write_text("\n".join(lines) + "\n")
    print(f"[OK]    Primer info written → {out_txt}")


if __name__ == "__main__":
    main()
PYEOF

info "Running 3-tier IUPAC-aware primer detection engine..."
python3 "$PYTMP"

# ── Print result summary ──────────────────────────────────────────────────────
sep
ok "Detection complete!"
echo ""
echo "  ── Result file: $OUTPUT_TXT ──"
echo ""
grep -v '^#' "$OUTPUT_TXT" | grep -v '^$' | while IFS= read -r line; do
  echo "  $line"
done
echo ""
sep
echo ""

DET_STATUS="$(grep '^DETECTION_STATUS=' "$OUTPUT_TXT" | cut -d= -f2)"
INF_REGION="$(grep '^INFERRED_REGION='  "$OUTPUT_TXT" | cut -d= -f2)"
FWD_SEQ="$(grep '^FORWARD_PRIMER_SEQUENCE=' "$OUTPUT_TXT" | cut -d= -f2)"
REV_SEQ="$(grep '^REVERSE_PRIMER_SEQUENCE=' "$OUTPUT_TXT" | cut -d= -f2)"
FWD_NAME="$(grep '^FORWARD_PRIMER_NAME=' "$OUTPUT_TXT" | cut -d= -f2)"
REV_NAME="$(grep '^REVERSE_PRIMER_NAME=' "$OUTPUT_TXT" | cut -d= -f2)"

echo "  ── How to use this result in the pipeline ──────────────────────"
echo ""
case "$DET_STATUS" in
  DETECTED)
    if [[ "$FWD_SEQ" != "None" && "$REV_SEQ" != "None" && "$REV_SEQ" != "N/A (single-end)" ]]; then
      echo "  # Pass primers to mbx_dada2_parameter_finder.sh:"
      printf '  mbx_dada2_parameter_finder.sh \\\n'
      printf '    /path/to/2_first_artifact_file/Paired_End_artifact.qza \\\n'
      printf '    --forward-primer %s \\\n' "$FWD_SEQ"
      printf '    --reverse-primer %s\n'   "$REV_SEQ"
    elif [[ "$FWD_SEQ" != "None" ]]; then
      echo "  # Single-end — pass forward primer only:"
      printf '  mbx_dada2_parameter_finder.sh \\\n'
      printf '    /path/to/artifact.qza \\\n'
      printf '    --forward-primer %s\n' "$FWD_SEQ"
    fi
    ;;
  TRIMMED)
    echo "  Primers were already TRIMMED OFF before delivery."
    echo "  Inferred amplicon region: ${INF_REGION:-unknown}"
    echo ""
    echo "  → mbx_dada2_parameter_finder.sh will use trim-left=0 automatically."
    echo "  → mbx_classifier_arranger.sh will use the FULL backbone reads"
    echo "    (no extract-reads step needed since primers aren't in the data)."
    ;;
  UNKNOWN)
    echo "  Detection failed — no primer above threshold AND no V-region anchor."
    echo ""
    echo "  Try one of:"
    echo "    1) Re-run with relaxed thresholds:"
    echo "         mbx_primer_identifier.sh '$FASTQ_DIR' \\"
    echo "             --mismatches 4 --offset 35 --min-rate 0.02 --report-best"
    echo ""
    echo "    2) Supply your primers via --extra-primers (TSV: name<TAB>seq<TAB>region<TAB>F|R)"
    echo ""
    echo "    3) Skip detection and pass primers directly to DADA2 step:"
    echo "         mbx_dada2_parameter_finder.sh /path/to/artifact.qza \\"
    echo "             --forward-primer GTGYCAGCMGCCGCGGTAA \\"
    echo "             --reverse-primer GGACTACNVGGGTWTCTAAT"
    ;;
esac

echo ""
echo "  Output tree:"
echo "  $(dirname "$PRIMER_DIR")/"
echo "  └── 0_primer_handling/"
echo "      └── mbx_primer_info.txt"
[[ -f "$DEBUG_OUT" ]] && echo "      └── mbx_primer_debug_matrix.csv"
echo ""
