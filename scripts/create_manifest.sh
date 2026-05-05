#!/usr/bin/env bash
# =============================================================================
#  create_manifest.sh
#  Generate a QIIME2-compatible manifest TSV from a directory of .fastq.gz files
#
#  Compatible with bash 3.2+ (macOS default shell — no brew/conda bash needed)
#
#  OUTPUT STRUCTURE (sibling of the given FASTQ directory):
#    <parent_of_fastq_dir>/mbX_pro_outputs_<YYYYMMDD_HHMMSS>/
#    └── 1_manifest_file/
#        └── manifest.txt
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
warn()     { echo "[WARN]  $*" >&2; }
info()     { echo "[INFO]  $*"; }
sep()      { echo "────────────────────────────────────────────────────────────────"; }
to_upper() { echo "$1" | tr '[:lower:]' '[:upper:]'; }

# bash 3.2-safe absolute path (no realpath needed)
_abspath() {
  if [[ -d "$1" ]]; then
    cd "$1" && pwd
  elif [[ -f "$1" ]]; then
    echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  else
    return 1
  fi
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

create_manifest.sh — Generate a QIIME2-compatible manifest from .fastq.gz files

USAGE:
  bash create_manifest.sh <fastq_dir> [OPTIONS]
  create_manifest.sh <fastq_dir> [OPTIONS]       (if installed in ~/bin)

DESCRIPTION:
  Scans <fastq_dir> (and optionally its parent and child directories) for
  .fastq.gz files, infers sample IDs and read orientation (R1/R2), and writes
  a QIIME2-compatible tab-separated manifest file.

  Paired-end columns : sample-id | forward-absolute-filepath | reverse-absolute-filepath
  Single-end columns : sample-id | absolute-filepath

OUTPUT (placed ALONGSIDE the FASTQ directory, not inside it):
  .../rough_2/FASTQ/                              <- your input
  .../rough_2/mbX_pro_outputs_YYYYMMDD_HHMMSS/   <- output lands here
              └── 1_manifest_file/
                  └── manifest.txt

OPTIONS:
  --no-search-parent   Skip searching the parent directory
  --no-search-child    Skip searching child (sub-)directories
  -h, --help           Show this help message and exit

EXAMPLES:
  bash create_manifest.sh /path/to/FASTQ
  bash create_manifest.sh /path/to/FASTQ --no-search-parent
  bash create_manifest.sh /path/to/FASTQ --no-search-parent --no-search-child

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
FASTQ_DIR=""
SEARCH_PARENT=true
SEARCH_CHILD=true

[[ $# -eq 0 ]] && usage

for arg in "$@"; do
  case "$arg" in
    -h|--help)          usage ;;
    --no-search-parent) SEARCH_PARENT=false ;;
    --no-search-child)  SEARCH_CHILD=false ;;
    -*)  err "Unknown option: '${arg}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$FASTQ_DIR" ]]; then
        FASTQ_DIR="$arg"
      else
        err "Unexpected extra argument: '${arg}'  —  only one positional argument (<fastq_dir>) is accepted."
      fi
      ;;
  esac
done

[[ -z "$FASTQ_DIR" ]] && err "No FASTQ directory provided.  Run with --help for usage."
[[ -d "$FASTQ_DIR" ]] || err "Directory does not exist: '${FASTQ_DIR}'"

FASTQ_DIR="$(_abspath "$FASTQ_DIR")"

# ── Create output directories ─────────────────────────────────────────────────
# Output lands ALONGSIDE (sibling of) the FASTQ directory by default.
# When invoked by the mbXPro orchestrator, MBX_OUT_DIR is exported and we
# reuse that single shared directory so all 18 steps land in ONE place.
if [[ -n "${MBX_OUT_DIR:-}" ]]; then
  OUT_ROOT="$MBX_OUT_DIR"
  TIMESTAMP="$(basename "$OUT_ROOT" | sed -E 's/^mbX_pro_outputs_//')"
  [[ -z "$TIMESTAMP" ]] && TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
else
  TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
  OUT_ROOT="$(dirname "$FASTQ_DIR")/mbX_pro_outputs_${TIMESTAMP}"
fi
MANIFEST_DIR="${OUT_ROOT}/1_manifest_file"
MANIFEST_FILE="${MANIFEST_DIR}/manifest.txt"

mkdir -p "$MANIFEST_DIR" \
  || err "Could not create output directory: '${MANIFEST_DIR}'  —  check write permissions."

sep
info "Run timestamp      : $TIMESTAMP"
info "Output root        : $OUT_ROOT"
info "FASTQ directory    : $FASTQ_DIR"
info "Search parent dir  : $SEARCH_PARENT"
info "Search child dirs  : $SEARCH_CHILD"
sep

# ── Temp workspace (auto-cleaned on exit) ─────────────────────────────────────
# R1/ and R2/ subdirs act as our "associative array" — bash 3.2-safe
# Each file is named after the sample ID and its content is the fastq filepath
WORK_DIR="$(mktemp -d)"
mkdir -p "$WORK_DIR/R1" "$WORK_DIR/R2"
DIRS_FILE="$WORK_DIR/candidate_dirs.txt"
SEEN_FILE="$WORK_DIR/seen_files.txt"
FASTQ_LIST="$WORK_DIR/fastq_files.txt"
touch "$DIRS_FILE" "$SEEN_FILE" "$FASTQ_LIST"

trap 'rm -rf "$WORK_DIR"' EXIT

# ── Collect candidate directories ─────────────────────────────────────────────
_add_dir() {
  local d
  d="$(_abspath "$1" 2>/dev/null)" || return 0
  [[ -d "$d" ]] || return 0
  grep -qxF "$d" "$DIRS_FILE" 2>/dev/null || echo "$d" >> "$DIRS_FILE"
}

_add_dir "$FASTQ_DIR"

if $SEARCH_PARENT; then
  _PARENT="$(dirname "$FASTQ_DIR")"
  [[ "$_PARENT" != "$FASTQ_DIR" ]] && _add_dir "$_PARENT"
fi

if $SEARCH_CHILD; then
  while IFS= read -r -d '' child; do
    [[ -d "$child" ]] && _add_dir "$child"
  done < <(find "$FASTQ_DIR" -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null || true)
fi

info "Directories to search ($(wc -l < "$DIRS_FILE" | tr -d ' ')):"
while IFS= read -r d; do info "  └─ $d"; done < "$DIRS_FILE"
sep

# ── Discover .fastq.gz files (deduplicated by absolute path) ──────────────────
while IFS= read -r d; do
  while IFS= read -r -d '' f; do
    fp="$(_abspath "$f" 2>/dev/null)" || continue
    grep -qxF "$fp" "$SEEN_FILE" 2>/dev/null && continue
    echo "$fp" >> "$SEEN_FILE"
    echo "$fp" >> "$FASTQ_LIST"
  done < <(find "$d" -maxdepth 1 -iname "*.fastq.gz" -type f -print0 2>/dev/null || true)
done < "$DIRS_FILE"

TOTAL_FILES=$(wc -l < "$FASTQ_LIST" | tr -d ' ')
[[ "$TOTAL_FILES" -eq 0 ]] && err "No .fastq.gz files were found in the provided directory, one level up, or one level down.
  Searched:
$(while IFS= read -r d; do echo "    $d"; done < "$DIRS_FILE")"

info "Total .fastq.gz files discovered: $TOTAL_FILES"
sep

# ── Helper: detect read direction ─────────────────────────────────────────────
_detect_read() {
  local nm
  nm="$(basename "$1")"
  if echo "$nm" | grep -qiE 'r1'; then echo "R1"
  elif echo "$nm" | grep -qiE 'r2'; then echo "R2"
  else echo "NA"
  fi
}

# ── Helper: extract sample ID ─────────────────────────────────────────────────
_extract_sample() {
  local nm id
  nm="$(basename "$1")"

  # Strategy 1: match s(ample)?[-_]?[0-9]+  (mirrors R function regex)
  id="$(echo "$nm" | grep -oiE 's(ample)?[-_]?[0-9]+' | head -1)"

  # Strategy 2: strip extension & read tag, trim surrounding punctuation
  if [[ -z "$id" ]]; then
    id="$(printf '%s' "$nm" | perl -pe '
      s/\.fastq\.gz$//i;
      s/r[12]//gi;
      s/[._-]+$//;
      s/^[._-]+//;
      chomp;
    ')"
  fi

  to_upper "$id"
}

# ── Process each file into WORK_DIR/R1/ or WORK_DIR/R2/ ──────────────────────
# Each slot file is named after the sample_id and contains the fastq filepath.
# Duplicate detection: if the slot file already exists → two files for same sample+read.

PROCESSED=0
UNDETECTED=0

while IFS= read -r f; do
  _read="$(_detect_read "$f")"
  if [[ "$_read" == "NA" ]]; then
    warn "No R1/R2 tag detected in: $(basename "$f")  —  skipping."
    UNDETECTED=$(( UNDETECTED + 1 ))
    continue
  fi

  _sid="$(_extract_sample "$f")"
  SLOT="$WORK_DIR/$_read/$_sid"

  if [[ -f "$SLOT" ]]; then
    EXISTING="$(cat "$SLOT")"
    err "Duplicate $_read file found for sample '$_sid':
  File 1: $EXISTING
  File 2: $f
  Each sample must map to exactly one file per read direction."
  fi

  printf '%s' "$f" > "$SLOT"
  PROCESSED=$(( PROCESSED + 1 ))
done < "$FASTQ_LIST"

[[ "$PROCESSED" -eq 0 ]] && \
  err "Could not detect R1 or R2 in any filename.
  Ensure filenames contain 'R1' or 'R2' (case-insensitive) to indicate read direction."

if [[ "$UNDETECTED" -gt 0 ]]; then
  warn "$UNDETECTED file(s) were skipped because no R1/R2 tag was found."
fi

# ── Determine paired-end vs single-end ───────────────────────────────────────
R2_COUNT=$(find "$WORK_DIR/R2" -type f 2>/dev/null | wc -l | tr -d ' ')
[[ "$R2_COUNT" -gt 0 ]] && READ_TYPE="paired" || READ_TYPE="single"
info "Inferred read type: $READ_TYPE"

# ── Collect all unique sample IDs (sorted alphabetically) ────────────────────
SAMPLE_IDS=()
while IFS= read -r sid; do
  SAMPLE_IDS+=("$sid")
done < <(
  find "$WORK_DIR/R1" "$WORK_DIR/R2" -type f 2>/dev/null \
    | while IFS= read -r p; do basename "$p"; done \
    | sort -u
)

# ── Build & validate manifest — PAIRED-END ────────────────────────────────────
if [[ "$READ_TYPE" == "paired" ]]; then

  MISSING=()
  for sid in "${SAMPLE_IDS[@]}"; do
    if [[ ! -f "$WORK_DIR/R1/$sid" || ! -f "$WORK_DIR/R2/$sid" ]]; then
      MISSING+=("$sid")
    fi
  done

  if [[ ${#MISSING[@]} -gt 0 ]]; then
    err "The following samples are missing either R1 or R2:
$(printf '  - %s\n' "${MISSING[@]}")
Both R1 and R2 must be present for every sample in paired-end mode."
  fi

  printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n' \
    > "$MANIFEST_FILE"

  for sid in "${SAMPLE_IDS[@]}"; do
    r1="$(cat "$WORK_DIR/R1/$sid")"
    r2="$(cat "$WORK_DIR/R2/$sid")"
    printf '%s\t%s\t%s\n' "$sid" "$r1" "$r2" >> "$MANIFEST_FILE"
  done

# ── Build & validate manifest — SINGLE-END ───────────────────────────────────
else

  R1_COUNT=$(find "$WORK_DIR/R1" -type f 2>/dev/null | wc -l | tr -d ' ')
  [[ "$R1_COUNT" -eq 0 ]] && err "Single-end mode inferred, but no R1 files were found."

  printf 'sample-id\tabsolute-filepath\n' > "$MANIFEST_FILE"

  for sid in "${SAMPLE_IDS[@]}"; do
    [[ -f "$WORK_DIR/R1/$sid" ]] || continue
    fp="$(cat "$WORK_DIR/R1/$sid")"
    printf '%s\t%s\n' "$sid" "$fp" >> "$MANIFEST_FILE"
  done

fi

N_SAMPLES="${#SAMPLE_IDS[@]}"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
info "✔  Read type          : $READ_TYPE"
info "✔  Files processed    : $PROCESSED"
info "✔  Samples in manifest: $N_SAMPLES"
info "✔  Manifest written   : $MANIFEST_FILE"
sep
echo ""
echo "  Next step — import into QIIME2:"
if [[ "$READ_TYPE" == "paired" ]]; then
  printf '  qiime tools import \\\n'
  printf "    --type 'SampleData[PairedEndSequencesWithQuality]' \\\n"
  printf '    --input-path  %s \\\n' "$MANIFEST_FILE"
  printf '    --output-path paired-end-demux.qza \\\n'
  printf "    --input-format PairedEndFastqManifestPhred33V2\n"
else
  printf '  qiime tools import \\\n'
  printf "    --type 'SampleData[SequencesWithQuality]' \\\n"
  printf '    --input-path  %s \\\n' "$MANIFEST_FILE"
  printf '    --output-path single-end-demux.qza \\\n'
  printf "    --input-format SingleEndFastqManifestPhred33V2\n"
fi
echo ""
