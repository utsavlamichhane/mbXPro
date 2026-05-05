#!/usr/bin/env bash
# =============================================================================
#  mbx_dada2_run.sh
#  Run QIIME2 DADA2 denoising using parameters from mbx_dada2_parameter_finder
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  USAGE:
#    mbx_dada2_run.sh <dada2_parameters.txt> <metadata.txt> [OPTIONS]
#
#  Full output structure:
#    mbX_pro_outputs_<timestamp>/
#    ├── 1_manifest_file/
#    ├── 2_first_artifact_file/
#    ├── 3_dada2_parameters/
#    │   └── dada2_parameters.txt      <- you provide this path
#    └── 4_dada2_outputs/
#        ├── feature_table.qza
#        ├── feature_table.qzv
#        ├── representative_sequences.qza
#        ├── representative_sequences.qzv
#        ├── dada2_stats.qza
#        └── dada2_stats.qzv
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
step() {
  echo ""
  echo "┌─────────────────────────────────────────────────────────────────"
  echo "│  Step $*"
  echo "└─────────────────────────────────────────────────────────────────"
}
sep()  { echo "────────────────────────────────────────────────────────────────"; }

# bash 3.2-safe absolute path
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

mbx_dada2_run.sh — Run QIIME2 DADA2 denoising + generate all summary artifacts

USAGE:
  mbx_dada2_run.sh <dada2_parameters.txt> <metadata.txt> [OPTIONS]

DESCRIPTION:
  Reads DADA2 parameters from the .txt file produced by mbx_dada2_parameter_finder.sh,
  validates the metadata file, then runs the full DADA2 denoising pipeline:

    Step 1  Validate metadata file (format, headers, sample IDs)
    Step 2  Parse DADA2 parameters + locate .qza artifact
    Step 3  qiime dada2 denoise-paired
    Step 4  qiime metadata tabulate       → dada2_stats.qzv
    Step 5  qiime feature-table summarize → feature_table.qzv
    Step 6  qiime feature-table tabulate-seqs → representative_sequences.qzv

  All outputs go into 4_dada2_outputs/ inside the same mbX_pro_outputs_* directory.

OUTPUT:
  mbX_pro_outputs_<timestamp>/
  └── 4_dada2_outputs/
      ├── feature_table.qza
      ├── feature_table.qzv
      ├── representative_sequences.qza
      ├── representative_sequences.qzv
      ├── dada2_stats.qza
      └── dada2_stats.qzv

OPTIONS:
  --threads <N>    Number of CPU threads for DADA2
                   Default: auto-detected from available cores
                   Override with e.g. --threads 4
  --dry-run        Print all commands without executing them
  -h, --help       Show this help message and exit

EXAMPLES:
  mbx_dada2_run.sh \
    /path/to/3_dada2_parameters/dada2_parameters.txt \
    /path/to/metadata.txt

  mbx_dada2_run.sh \
    /path/to/3_dada2_parameters/dada2_parameters.txt \
    /path/to/metadata.txt \
    --threads 8

METADATA REQUIREMENTS (from mbX R package validation):
  • Format  : tab-separated (.txt/.tsv) or comma-separated (.csv)
  • First column header must be one of (case-insensitive):
      id, sampleid, sample id, sample-id,
      featureid, feature id, feature-id
  • Sample IDs : no duplicates, no empty values
  • Sample IDs : no leading/trailing whitespace
  • QIIME2 note: if a #q2:types directive row is present it is allowed

COMMON ERRORS AND FIXES:
  "qiime: command not found"
    → conda activate qiime2-amplicon-2025.4

  "dada2_parameters.txt: no trunc-len values found"
    → Run mbx_dada2_parameter_finder.sh first to generate the parameters file.

  "Artifact file not found"
    → The parameters file records the .qza path. Ensure the artifact still
      exists at the path listed in dada2_parameters.txt.

  "metadata first column header is invalid"
    → The first column must be named one of: sample-id, id, sampleid, etc.
      Open your metadata file and rename the first column header accordingly.

  "duplicate sample IDs found"
    → Each row in metadata must have a unique sample ID.

  "DADA2 denoising failed"
    → Common causes: too few reads, very low quality data, or truncation
      lengths shorter than primers. Re-run mbx_dada2_parameter_finder.sh
      or manually adjust --p-trunc-len-f / --p-trunc-len-r.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
PARAM_TXT=""
METADATA_PATH=""
THREADS="auto"
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)    usage ;;
    --dry-run)    DRY_RUN=true; shift ;;
    --threads)
      THREADS="${2:?Missing value for --threads}"; shift 2 ;;
    -*)
      err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$PARAM_TXT" ]]; then
        PARAM_TXT="$1"; shift
      elif [[ -z "$METADATA_PATH" ]]; then
        METADATA_PATH="$1"; shift
      else
        err "Unexpected extra argument: '${1}'
  Usage: mbx_dada2_run.sh <dada2_parameters.txt> <metadata.txt> [OPTIONS]"
      fi
      ;;
  esac
done

[[ -z "$PARAM_TXT" ]]     && err "No dada2_parameters.txt provided.  Run with --help for usage."
[[ -z "$METADATA_PATH" ]] && err "No metadata file provided.  Run with --help for usage."

[[ -f "$PARAM_TXT" ]]     || err "Parameters file does not exist: '${PARAM_TXT}'
  → Run mbx_dada2_parameter_finder.sh first."
[[ -f "$METADATA_PATH" ]] || err "Metadata file does not exist: '${METADATA_PATH}'
  → Check the path for typos."

PARAM_TXT="$(_abspath "$PARAM_TXT")"
METADATA_PATH="$(_abspath "$METADATA_PATH")"

# Auto-detect cores unless user passed an explicit number
if [[ "$THREADS" == "auto" || "$THREADS" == "0" ]]; then
  if command -v nproc &>/dev/null; then
    THREADS="$(nproc)"
  elif command -v sysctl &>/dev/null; then
    THREADS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
  else
    THREADS=1
    warn "Could not auto-detect core count. Falling back to 1 thread."
  fi
  info "Auto-detected $THREADS CPU thread(s)."
fi

# ── Verify QIIME2 ─────────────────────────────────────────────────────────────
command -v qiime &>/dev/null || err "qiime command not found.
  → conda activate qiime2-amplicon-2025.4"

sep
info "Parameters file    : $PARAM_TXT"
info "Metadata file      : $METADATA_PATH"
info "Threads            : $THREADS"
$DRY_RUN && warn "DRY-RUN mode — commands will be printed but NOT executed."
sep

# ─────────────────────────────────────────────────────────────────────────────
step "1/6 — Validate metadata file"
# ─────────────────────────────────────────────────────────────────────────────
# Mirrors the validation logic from mbX_functions.R:
#   - File extension must be txt/tsv/csv
#   - First column header must be a QIIME2-valid sample-id alias
#   - No empty sample IDs
#   - No duplicate sample IDs
#   - No leading/trailing whitespace in sample IDs

META_EXT="${METADATA_PATH##*.}"
META_EXT_LOWER="$(echo "$META_EXT" | tr '[:upper:]' '[:lower:]')"

case "$META_EXT_LOWER" in
  txt|tsv)
    META_DELIM=$'\t'
    ;;
  csv)
    META_DELIM=","
    ;;
  xls|xlsx)
    err "Excel metadata files (.xls/.xlsx) are not directly readable in bash.
  → Export your metadata to a tab-separated .txt file and re-run.
  → In Excel: File → Save As → Tab Delimited Text (.txt)"
    ;;
  *)
    err "Unsupported metadata file format: '.${META_EXT}'
  → Accepted formats: .txt (tab-separated), .tsv, .csv
  → This matches the format validation in the mbX R package."
    ;;
esac

info "Metadata format    : $META_EXT_LOWER (delimiter: $([ "$META_DELIM" = $'\t' ] && echo 'TAB' || echo 'COMMA'))"

# ── Check first column header ─────────────────────────────────────────────────
# Read first line, extract first field (handles both tab and comma)
FIRST_HEADER="$(head -1 "$METADATA_PATH" | cut -d"$(printf '%s' "$META_DELIM")" -f1)"
# Trim whitespace and lowercase (mirrors R's tolower(trimws(...)))
FIRST_HEADER_CLEAN="$(printf '%s' "$FIRST_HEADER" | tr '[:upper:]' '[:lower:]' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"

# Valid headers from mbX_functions.R (lines 60-61 and 747-748)
VALID_HEADERS="id|sampleid|sample id|sample-id|featureid|feature id|feature-id"

# Remove BOM character if present (common in Excel-exported files)
FIRST_HEADER_CLEAN="$(printf '%s' "$FIRST_HEADER_CLEAN" | sed 's/^\xef\xbb\xbf//')"

if ! echo "$FIRST_HEADER_CLEAN" | grep -qiE "^(id|sampleid|sample id|sample-id|featureid|feature id|feature-id)$"; then
  err "Metadata first column header is invalid.
  Found    : '${FIRST_HEADER}'
  Expected : one of — id, sampleid, sample id, sample-id,
                       featureid, feature id, feature-id
             (case-insensitive, as validated by the mbX R package)
  → Open your metadata file and rename the first column header.
  → QIIME2 standard is 'sample-id'."
fi
ok "First column header is valid: '${FIRST_HEADER}'"

# ── Skip #q2:types directive row if present (QIIME2-specific) ─────────────────
SECOND_LINE="$(sed -n '2p' "$METADATA_PATH" | cut -d"$(printf '%s' "$META_DELIM")" -f1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
if [[ "$SECOND_LINE" == "#q2:types" ]]; then
  DATA_START=3
  info "QIIME2 #q2:types directive row detected — will skip for validation."
else
  DATA_START=2
fi

# ── Extract all sample IDs (first column, skipping header [+ directive]) ──────
WORK_DIR_META="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR_META"' EXIT

tail -n "+${DATA_START}" "$METADATA_PATH" \
  | cut -d"$(printf '%s' "$META_DELIM")" -f1 \
  | grep -v '^[[:space:]]*$' \
  > "$WORK_DIR_META/sample_ids.txt" || true

TOTAL_SAMPLE_IDS="$(wc -l < "$WORK_DIR_META/sample_ids.txt" | tr -d ' ')"
[[ "$TOTAL_SAMPLE_IDS" -eq 0 ]] && err "Metadata has no data rows after the header.
  → Ensure your metadata file contains at least one sample row."

info "Sample rows found  : $TOTAL_SAMPLE_IDS"

# ── Check for empty sample IDs ────────────────────────────────────────────────
EMPTY_COUNT=0
while IFS= read -r sid; do
  trimmed="$(printf '%s' "$sid" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  [[ -z "$trimmed" ]] && EMPTY_COUNT=$(( EMPTY_COUNT + 1 ))
done < "$WORK_DIR_META/sample_ids.txt"

[[ "$EMPTY_COUNT" -gt 0 ]] && err "$EMPTY_COUNT empty sample ID(s) found in metadata.
  → Every row must have a value in the first (sample-id) column.
  → Open the metadata file and fill in or remove the empty rows."
ok "No empty sample IDs found."

# ── Check for duplicate sample IDs ───────────────────────────────────────────
sort "$WORK_DIR_META/sample_ids.txt" > "$WORK_DIR_META/sorted_ids.txt"
sort -u "$WORK_DIR_META/sample_ids.txt" > "$WORK_DIR_META/unique_ids.txt"

TOTAL_SORTED="$(wc -l < "$WORK_DIR_META/sorted_ids.txt" | tr -d ' ')"
TOTAL_UNIQUE="$(wc -l < "$WORK_DIR_META/unique_ids.txt" | tr -d ' ')"

if [[ "$TOTAL_SORTED" -ne "$TOTAL_UNIQUE" ]]; then
  DUPS="$(comm -23 "$WORK_DIR_META/sorted_ids.txt" "$WORK_DIR_META/unique_ids.txt" | head -10)"
  err "Duplicate sample IDs found in metadata:
$(printf '  - %s\n' $DUPS)
  → Each sample must appear exactly once in the metadata file.
  → Check for accidental duplicate rows and remove them."
fi
ok "No duplicate sample IDs found."

# ── Check for leading/trailing whitespace in sample IDs ───────────────────────
WS_COUNT=0
while IFS= read -r sid; do
  trimmed="$(printf '%s' "$sid" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  [[ "$sid" != "$trimmed" ]] && WS_COUNT=$(( WS_COUNT + 1 ))
done < "$WORK_DIR_META/sample_ids.txt"

if [[ "$WS_COUNT" -gt 0 ]]; then
  err "$WS_COUNT sample ID(s) have leading or trailing whitespace.
  → QIIME2 will fail to match these IDs to the artifact.
  → Open the metadata and remove all extra spaces from sample ID values."
fi
ok "No whitespace issues in sample IDs."

# ── Check for QIIME2-incompatible characters in sample IDs ───────────────────
BAD_CHAR_COUNT=0
BAD_SAMPLES=""
while IFS= read -r sid; do
  if echo "$sid" | grep -qE '[^A-Za-z0-9._\-]'; then
    BAD_CHAR_COUNT=$(( BAD_CHAR_COUNT + 1 ))
    BAD_SAMPLES="${BAD_SAMPLES}  - ${sid}\n"
  fi
done < "$WORK_DIR_META/sample_ids.txt"

if [[ "$BAD_CHAR_COUNT" -gt 0 ]]; then
  warn "$BAD_CHAR_COUNT sample ID(s) contain characters that QIIME2 may not accept.
  QIIME2 recommends only: letters, numbers, periods (.), hyphens (-), underscores (_)
  Flagged IDs:
$(printf '%b' "$BAD_SAMPLES")
  → Consider renaming these samples to avoid potential import errors."
fi

ok "Metadata validation passed (${TOTAL_SAMPLE_IDS} samples)."

# ─────────────────────────────────────────────────────────────────────────────
step "2/6 — Parse DADA2 parameters and locate artifact"
# ─────────────────────────────────────────────────────────────────────────────

# Extract the four core parameters from the .txt file
# The file format (from create_dada2_parameters_txt.sh) has lines like:
#   --p-trunc-len-f 248
#   --p-trunc-len-r 233
#   --p-trim-left-f 18
#   --p-trim-left-r 22

_parse_param() {
  local key="$1"
  local val
  val="$(grep -m1 "^${key} " "$PARAM_TXT" | awk '{print $2}' || true)"
  echo "$val"
}

TRUNC_F="$(_parse_param "--p-trunc-len-f")"
TRUNC_R="$(_parse_param "--p-trunc-len-r")"
TRIM_F="$(_parse_param  "--p-trim-left-f")"
TRIM_R="$(_parse_param  "--p-trim-left-r")"

_check_param() {
  local name="$1" val="$2"
  [[ -z "$val" ]] && err "Could not parse '${name}' from: $PARAM_TXT
  → Ensure the file was generated by mbx_dada2_parameter_finder.sh.
  → Expected lines like:  --p-trunc-len-f 248"
  echo "$val" | grep -qE '^[0-9]+$' || err "Parameter '${name}' is not a valid integer: '${val}'
  → Open $PARAM_TXT and verify the value manually."
}

_check_param "--p-trunc-len-f" "$TRUNC_F"
_check_param "--p-trunc-len-r" "$TRUNC_R"
_check_param "--p-trim-left-f" "$TRIM_F"
_check_param "--p-trim-left-r" "$TRIM_R"

info "trunc-len-f        : $TRUNC_F"
info "trunc-len-r        : $TRUNC_R"
info "trim-left-f        : $TRIM_F bp"
info "trim-left-r        : $TRIM_R bp"

# Sanity check: trunc must be > trim
if [[ "$TRUNC_F" -le "$TRIM_F" ]]; then
  err "trunc-len-f ($TRUNC_F) must be greater than trim-left-f ($TRIM_F).
  → Re-run mbx_dada2_parameter_finder.sh — the current values are inconsistent."
fi
if [[ "$TRUNC_R" -le "$TRIM_R" ]]; then
  err "trunc-len-r ($TRUNC_R) must be greater than trim-left-r ($TRIM_R).
  → Re-run mbx_dada2_parameter_finder.sh — the current values are inconsistent."
fi

# ── Locate the .qza artifact ──────────────────────────────────────────────────
# The parameters file has a line: "Input artifact: /path/to/Paired_End_artifact.qza"
QZA_FROM_PARAMS="$(grep -m1 "^Input artifact:" "$PARAM_TXT" | sed 's/^Input artifact:[[:space:]]*//' || true)"

if [[ -z "$QZA_FROM_PARAMS" ]]; then
  err "Could not find 'Input artifact:' line in: $PARAM_TXT
  → The parameters file may be corrupted or hand-edited.
  → Re-run mbx_dada2_parameter_finder.sh to regenerate it."
fi

[[ -f "$QZA_FROM_PARAMS" ]] || err "Artifact file listed in parameters file does not exist:
  Listed   : $QZA_FROM_PARAMS
  → Ensure the artifact file has not been moved or deleted.
  → If it was moved, update the 'Input artifact:' line in $PARAM_TXT"

INPUT_QZA="$QZA_FROM_PARAMS"
info "Input artifact     : $INPUT_QZA"

# ── Resolve output directory ───────────────────────────────────────────────────
# params live in: mbX_pro_outputs_*/3_dada2_parameters/dada2_parameters.txt
PARAM_DIR="$(dirname "$PARAM_TXT")"
PARAM_DIR_NAME="$(basename "$PARAM_DIR")"
MBX_OUT_DIR="$(dirname "$PARAM_DIR")"

if [[ "$PARAM_DIR_NAME" != "3_dada2_parameters" ]]; then
  warn "Parameters file parent folder is '${PARAM_DIR_NAME}', not '3_dada2_parameters'."
  warn "Placing 4_dada2_outputs/ alongside the parameters file's parent folder anyway."
fi

DADA2_OUT_DIR="${MBX_OUT_DIR}/4_dada2_outputs"
mkdir -p "$DADA2_OUT_DIR" \
  || err "Could not create output directory: '${DADA2_OUT_DIR}'  —  check write permissions."

info "Output directory   : $DADA2_OUT_DIR"

# ── Define output paths ────────────────────────────────────────────────────────
FEATURE_TABLE_QZA="${DADA2_OUT_DIR}/feature_table.qza"
FEATURE_TABLE_QZV="${DADA2_OUT_DIR}/feature_table.qzv"
REP_SEQS_QZA="${DADA2_OUT_DIR}/representative_sequences.qza"
REP_SEQS_QZV="${DADA2_OUT_DIR}/representative_sequences.qzv"
DADA2_STATS_QZA="${DADA2_OUT_DIR}/dada2_stats.qza"
DADA2_STATS_QZV="${DADA2_OUT_DIR}/dada2_stats.qzv"

# ── Dry-run helper ─────────────────────────────────────────────────────────────
run_cmd() {
  echo ""
  echo "  \$ $*"
  if ! $DRY_RUN; then
    "$@" || err "Command failed: $1
  → Check the QIIME2 output above for details."
  fi
}

# ── Print full plan ────────────────────────────────────────────────────────────
sep
echo ""
echo "  ── DADA2 parameters to be used ───────────────────────────────"
echo "    --p-trunc-len-f  $TRUNC_F"
echo "    --p-trunc-len-r  $TRUNC_R"
echo "    --p-trim-left-f  $TRIM_F"
echo "    --p-trim-left-r  $TRIM_R"
echo "    --p-n-threads    $THREADS"
echo ""
echo "  ── Input files ────────────────────────────────────────────────"
echo "    Artifact  : $INPUT_QZA"
echo "    Metadata  : $METADATA_PATH"
echo ""
echo "  ── Output files ───────────────────────────────────────────────"
echo "    $FEATURE_TABLE_QZA"
echo "    $FEATURE_TABLE_QZV"
echo "    $REP_SEQS_QZA"
echo "    $REP_SEQS_QZV"
echo "    $DADA2_STATS_QZA"
echo "    $DADA2_STATS_QZV"
sep

# ─────────────────────────────────────────────────────────────────────────────
step "3/6 — DADA2 denoise-paired  (this is the longest step)"
# ─────────────────────────────────────────────────────────────────────────────

run_cmd qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$INPUT_QZA" \
  --p-trunc-len-f "$TRUNC_F" \
  --p-trunc-len-r "$TRUNC_R" \
  --p-trim-left-f "$TRIM_F" \
  --p-trim-left-r "$TRIM_R" \
  --p-n-threads   "$THREADS" \
  --o-table                    "$FEATURE_TABLE_QZA" \
  --o-representative-sequences "$REP_SEQS_QZA" \
  --o-denoising-stats          "$DADA2_STATS_QZA"

$DRY_RUN || ok "DADA2 denoising complete."

# ─────────────────────────────────────────────────────────────────────────────
step "4/6 — Tabulate DADA2 stats → dada2_stats.qzv"
# ─────────────────────────────────────────────────────────────────────────────

run_cmd qiime metadata tabulate \
  --m-input-file "$DADA2_STATS_QZA" \
  --o-visualization "$DADA2_STATS_QZV"

$DRY_RUN || ok "Stats visualization written: $DADA2_STATS_QZV"

# ─────────────────────────────────────────────────────────────────────────────
step "5/6 — Summarize feature table → feature_table.qzv"
# ─────────────────────────────────────────────────────────────────────────────

run_cmd qiime feature-table summarize \
  --i-table "$FEATURE_TABLE_QZA" \
  --m-sample-metadata-file "$METADATA_PATH" \
  --o-visualization "$FEATURE_TABLE_QZV"

$DRY_RUN || ok "Feature table visualization written: $FEATURE_TABLE_QZV"

# ─────────────────────────────────────────────────────────────────────────────
step "6/6 — Tabulate representative sequences → representative_sequences.qzv"
# ─────────────────────────────────────────────────────────────────────────────

run_cmd qiime feature-table tabulate-seqs \
  --i-data "$REP_SEQS_QZA" \
  --o-visualization "$REP_SEQS_QZV"

$DRY_RUN || ok "Representative sequences visualization written: $REP_SEQS_QZV"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete — no commands were executed."
else
  ok "All 6 steps completed successfully!"
fi
sep
echo ""
echo "  Output structure:"
echo "  $MBX_OUT_DIR/"
echo "  ├── 1_manifest_file/"
echo "  ├── 2_first_artifact_file/"
echo "  ├── 3_dada2_parameters/"
echo "  └── 4_dada2_outputs/"
echo "      ├── feature_table.qza"
echo "      ├── feature_table.qzv            ← feature counts per sample"
echo "      ├── representative_sequences.qza"
echo "      ├── representative_sequences.qzv  ← BLAST-able ASV sequences"
echo "      ├── dada2_stats.qza"
echo "      └── dada2_stats.qzv              ← reads in/out per step"
echo ""
echo "  View all .qzv files at: https://view.qiime2.org"
echo ""
echo "  ── Suggested next steps ─────────────────────────────────────────"
echo "  Taxonomic classification:"
printf '    qiime feature-classifier classify-sklearn \\\n'
printf '      --i-classifier  /path/to/classifier.qza \\\n'
printf '      --i-reads       %s \\\n' "$REP_SEQS_QZA"
printf '      --o-classification %s/taxonomy.qza\n' "$DADA2_OUT_DIR"
echo ""
echo "  Phylogenetic tree:"
printf '    qiime phylogeny align-to-tree-mafft-fasttree \\\n'
printf '      --i-sequences   %s \\\n' "$REP_SEQS_QZA"
printf '      --o-alignment   %s/aligned_seqs.qza \\\n' "$DADA2_OUT_DIR"
printf '      --o-masked-alignment %s/masked_aligned_seqs.qza \\\n' "$DADA2_OUT_DIR"
printf '      --o-tree        %s/unrooted_tree.qza \\\n' "$DADA2_OUT_DIR"
printf '      --o-rooted-tree %s/rooted_tree.qza\n' "$DADA2_OUT_DIR"
echo ""
