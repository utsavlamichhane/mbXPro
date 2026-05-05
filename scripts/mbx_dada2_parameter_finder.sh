#!/usr/bin/env bash
# =============================================================================
#  mbx_dada2_parameter_finder.sh
#  Run QIIME2 demux QC steps + auto-recommend DADA2 parameters
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  INPUT  : Paired_End_artifact.qza  (from 2_first_artifact_file/)
#  OUTPUT : All files placed in 3_dada2_parameters/ alongside the other dirs
#
#  Full output structure:
#    mbX_pro_outputs_<timestamp>/
#    ├── 1_manifest_file/
#    │   └── manifest.txt
#    ├── 2_first_artifact_file/
#    │   ├── Paired_End_artifact.qza      <- you provide this path
#    │   └── Paired_End_artifact.qzv
#    └── 3_dada2_parameters/
#        ├── exported_demux_summary/      qiime tools export  (qzv)
#        ├── extracted_qza/               qiime tools extract (qza)
#        ├── extracted_qzv/               qiime tools extract (qzv)
#        ├── demux_read_counts.qza        per-sample read counts artifact
#        ├── demux_read_counts.qzv        read counts visualization
#        ├── exported_demux_read_counts/  qiime tools export  (counts qzv)
#        ├── extracted_demux_read_counts/ qiime tools extract (counts qza)
#        └── dada2_parameters.txt         recommended DADA2 parameters
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

mbx_dada2_parameter_finder.sh — QIIME2 demux QC + DADA2 parameter recommendation

USAGE:
  mbx_dada2_parameter_finder.sh <artifact.qza> [OPTIONS]

DESCRIPTION:
  Takes a QIIME2 SampleData[PairedEndSequencesWithQuality] .qza artifact
  produced by artifact_creator.sh and:

    Step 1  qiime tools peek            — verify artifact type
    Step 2  qiime tools export/extract  — unpack demux summary (.qzv)
    Step 3  qiime tools export/extract  — unpack artifact (.qza)
    Step 4  qiime demux tabulate-read-counts  — per-sample read counts
    Step 5  qiime metadata tabulate     — make read counts viewable (.qzv)
    Step 6  export/extract read counts
    Step 7  create_dada2_parameters_txt.sh  — recommend DADA2 parameters

  All outputs land in 3_dada2_parameters/ alongside 1_manifest_file/ and
  2_first_artifact_file/ inside the same mbX_pro_outputs_* directory.

OPTIONS:
  --forward-primer <SEQ>    Forward primer sequence (e.g. GTGYCAGCMGCCGCGGTAA)
                            Overrides any auto-detected primer.  Length sets --p-trim-left-f.
  --reverse-primer <SEQ>    Reverse primer sequence (e.g. GGACTACNVGGGTWTCTAAT)
                            Overrides any auto-detected primer.  Length sets --p-trim-left-r.
  --assume-primer-length N  Force --p-trim-left-f and --p-trim-left-r to N bp.
                            Overrides auto-detection.  Use 0 if you are sure
                            primers were already trimmed before delivery.
  --primer-info <FILE>      Explicit path to mbx_primer_info.txt (optional;
                            auto-discovered from <artifact>/../0_primer_handling/
                            when omitted).
  --amplicon-length <N>     Expected amplicon insert length (optional but helpful)
  --min-overlap <N>         Minimum merge overlap in bp (default: 12)
  --max-ee-f <N>            Max expected errors forward (default: 2.0)
  --max-ee-r <N>            Max expected errors reverse (default: 2.0)
  --dry-run                 Print all commands without executing them
  -h, --help                Show this help message and exit

PRIMER / TRIM-LEFT BEHAVIOUR:
  This script follows the verdict written by mbx_primer_identifier.sh
  into 0_primer_handling/mbx_primer_info.txt (auto-discovered):

    DETECTION_STATUS=DETECTED       primers found in the reads.
                                    Use FORWARD_PRIMER_SEQUENCE / REVERSE_PRIMER_SEQUENCE
                                    -> trim-left = primer length on each side.

    DETECTION_STATUS=USER_SUPPLIED  user supplied primers via --forward-primer / --reverse-primer.
                                    Treat exactly like DETECTED.

    DETECTION_STATUS=TRIMMED        sequencing facility already removed primers.
                                    -> trim-left = 0 on BOTH sides (NEVER 20).

    DETECTION_STATUS=UNKNOWN        detection failed (we genuinely don't know).
                                    -> trim-left = 20 on both sides (defensive default).

    no info file present            same as UNKNOWN -> trim-left = 20.

  Explicit flags ALWAYS win over auto-detection:
    --forward-primer / --reverse-primer    set those primers' lengths
    --assume-primer-length 0               force trim-left = 0
    --assume-primer-length 20              force trim-left = 20

EXAMPLES:
  # No primers — trim-left defaults to 20 on each side
  mbx_dada2_parameter_finder.sh /path/to/2_first_artifact_file/Paired_End_artifact.qza

  # 515F / 806R primers — trim-left set to primer lengths (19 and 20 bp)
  mbx_dada2_parameter_finder.sh Paired_End_artifact.qza \
    --forward-primer GTGYCAGCMGCCGCGGTAA \
    --reverse-primer GGACTACNVGGGTWTCTAAT

  # With amplicon length hint for better truncation estimation
  mbx_dada2_parameter_finder.sh Paired_End_artifact.qza \
    --forward-primer GTGYCAGCMGCCGCGGTAA \
    --reverse-primer GGACTACNVGGGTWTCTAAT \
    --amplicon-length 253

COMMON ERRORS AND FIXES:
  "qiime: command not found"
    → Activate your QIIME2 environment first:
      conda activate qiime2-amplicon-2025.4

  "create_dada2_parameters_txt.sh: command not found"
    → Install it to ~/bin:
      cp create_dada2_parameters_txt.sh ~/bin/
      chmod +x ~/bin/create_dada2_parameters_txt.sh

  "artifact does not exist"
    → Run artifact_creator.sh first to generate the .qza file.

  "not SampleData[PairedEndSequencesWithQuality]"
    → This script currently supports paired-end data only.
      Check the artifact type with: qiime tools peek <file.qza>

  "No candidate truncation lengths satisfy the constraints"
    → Data quality may be low. Try relaxing thresholds:
      --max-ee-f 5.0 --max-ee-r 5.0 --min-overlap 8

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
QZA_PATH=""
FORWARD_PRIMER=""
REVERSE_PRIMER=""
ASSUME_PRIMER_LEN=""              # empty = unset; will be auto-decided below
ASSUME_PRIMER_LEN_USER_SET=false  # tracks whether the flag was passed explicitly
PRIMER_INFO_OVERRIDE=""
AMPLICON_LENGTH=""
MIN_OVERLAP="12"
MAX_EE_F="2.0"
MAX_EE_R="2.0"
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)              usage ;;
    --dry-run)              DRY_RUN=true; shift ;;
    --forward-primer|-f)
      FORWARD_PRIMER="${2:?Missing value for $1}"; shift 2 ;;
    --reverse-primer|-r)
      REVERSE_PRIMER="${2:?Missing value for $1}"; shift 2 ;;
    --assume-primer-length)
      ASSUME_PRIMER_LEN="${2:?Missing value for $1}"
      ASSUME_PRIMER_LEN_USER_SET=true
      shift 2 ;;
    --primer-info)
      PRIMER_INFO_OVERRIDE="${2:?Missing value for $1}"; shift 2 ;;
    --amplicon-length)
      AMPLICON_LENGTH="${2:?Missing value for $1}"; shift 2 ;;
    --min-overlap)
      MIN_OVERLAP="${2:?Missing value for $1}"; shift 2 ;;
    --max-ee-f)
      MAX_EE_F="${2:?Missing value for $1}"; shift 2 ;;
    --max-ee-r)
      MAX_EE_R="${2:?Missing value for $1}"; shift 2 ;;
    -*)
      err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$QZA_PATH" ]]; then
        QZA_PATH="$1"; shift
      else
        err "Unexpected extra argument: '${1}'  —  only one positional argument (<artifact.qza>) is accepted."
      fi
      ;;
  esac
done

# ── Validate input ────────────────────────────────────────────────────────────
[[ -z "$QZA_PATH" ]] && err "No .qza artifact provided.  Run with --help for usage."
[[ -f "$QZA_PATH" ]] || err "Artifact file does not exist: '${QZA_PATH}'
  → Run artifact_creator.sh first to generate the .qza file."

QZA_PATH="$(_abspath "$QZA_PATH")"
QZV_PATH="${QZA_PATH%.qza}.qzv"

[[ -f "$QZV_PATH" ]] || \
  warn "Matching .qzv not found alongside .qza: $QZV_PATH
  → The export/extract steps for the .qzv will be skipped if this is missing."

# ── Verify QIIME2 is available ────────────────────────────────────────────────
command -v qiime &>/dev/null || err "qiime command not found.
  → Activate your QIIME2 conda environment:
    conda activate qiime2-amplicon-2025.4"

# ── Verify create_dada2_parameters_txt.sh is available ───────────────────────
PARAM_FINDER_CMD=""
if command -v create_dada2_parameters_txt.sh &>/dev/null; then
  PARAM_FINDER_CMD="create_dada2_parameters_txt.sh"
elif [[ -x "$HOME/bin/create_dada2_parameters_txt.sh" ]]; then
  PARAM_FINDER_CMD="$HOME/bin/create_dada2_parameters_txt.sh"
elif [[ -f "$(dirname "$0")/create_dada2_parameters_txt.sh" ]]; then
  PARAM_FINDER_CMD="$(dirname "$0")/create_dada2_parameters_txt.sh"
else
  err "create_dada2_parameters_txt.sh not found.
  → Install it to ~/bin/:
    cp create_dada2_parameters_txt.sh ~/bin/
    chmod +x ~/bin/create_dada2_parameters_txt.sh"
fi

# ── Resolve mbX_pro_outputs root and build output dir ────────────────────────
# QZA lives in: .../mbX_pro_outputs_*/2_first_artifact_file/file.qza
# We walk up two levels to find the root.

ARTIFACT_FILE_DIR="$(dirname "$QZA_PATH")"         # 2_first_artifact_file/
ARTIFACT_DIR_NAME="$(basename "$ARTIFACT_FILE_DIR")" # should be 2_first_artifact_file
MBX_OUT_DIR="$(dirname "$ARTIFACT_FILE_DIR")"       # mbX_pro_outputs_*/

if [[ "$ARTIFACT_DIR_NAME" != "2_first_artifact_file" ]]; then
  warn "The artifact's parent folder is '${ARTIFACT_DIR_NAME}', not '2_first_artifact_file'."
  warn "Placing 3_dada2_parameters/ alongside the artifact's parent folder anyway."
fi

DADA2_DIR="${MBX_OUT_DIR}/3_dada2_parameters"
mkdir -p "$DADA2_DIR" \
  || err "Could not create output directory: '${DADA2_DIR}'  —  check write permissions."

# ── Auto-discover mbx_primer_info.txt and apply DETECTION_STATUS rule ─────────
# This implements the contract documented in the help text:
#   DETECTED / USER_SUPPLIED -> use FORWARD_PRIMER_SEQUENCE / REVERSE_PRIMER_SEQUENCE
#   TRIMMED                  -> trim-left = 0
#   UNKNOWN / no file         -> trim-left = 20  (the python helper's default)
#
# Explicit flags (--forward-primer, --reverse-primer, --assume-primer-length)
# always override the auto-detection.

PRIMER_INFO_TXT=""
if [[ -n "$PRIMER_INFO_OVERRIDE" ]]; then
  PRIMER_INFO_TXT="$PRIMER_INFO_OVERRIDE"
elif [[ -f "${MBX_OUT_DIR}/0_primer_handling/mbx_primer_info.txt" ]]; then
  PRIMER_INFO_TXT="${MBX_OUT_DIR}/0_primer_handling/mbx_primer_info.txt"
fi

DETECTION_STATUS=""
DET_FWD_SEQ=""
DET_REV_SEQ=""
if [[ -n "$PRIMER_INFO_TXT" && -f "$PRIMER_INFO_TXT" ]]; then
  DETECTION_STATUS="$(grep '^DETECTION_STATUS=' "$PRIMER_INFO_TXT" | cut -d= -f2 | tr -d '[:space:]')"
  DET_FWD_SEQ="$(grep '^FORWARD_PRIMER_SEQUENCE=' "$PRIMER_INFO_TXT" | cut -d= -f2- | tr -d '[:space:]')"
  DET_REV_SEQ="$(grep '^REVERSE_PRIMER_SEQUENCE=' "$PRIMER_INFO_TXT" | cut -d= -f2- | tr -d '[:space:]')"
  [[ "$DET_FWD_SEQ" == "None" || "$DET_FWD_SEQ" == "N/A" ]] && DET_FWD_SEQ=""
  [[ "$DET_REV_SEQ" == "None" || "$DET_REV_SEQ" == "N/A" || "$DET_REV_SEQ" == "N/A(single-end)" ]] && DET_REV_SEQ=""
fi

# Apply auto-rules (only when the user did not pass an explicit override)
PRIMER_RULE_APPLIED="(none — defaults)"
if [[ -z "$FORWARD_PRIMER" && -z "$REVERSE_PRIMER" ]] && ! $ASSUME_PRIMER_LEN_USER_SET; then
  case "$DETECTION_STATUS" in
    DETECTED|USER_SUPPLIED)
      [[ -n "$DET_FWD_SEQ" ]] && FORWARD_PRIMER="$DET_FWD_SEQ"
      [[ -n "$DET_REV_SEQ" ]] && REVERSE_PRIMER="$DET_REV_SEQ"
      PRIMER_RULE_APPLIED="DETECTION_STATUS=$DETECTION_STATUS -> primers loaded from mbx_primer_info.txt"
      ;;
    TRIMMED)
      ASSUME_PRIMER_LEN="0"
      PRIMER_RULE_APPLIED="DETECTION_STATUS=TRIMMED -> --p-trim-left-f 0  --p-trim-left-r 0"
      ;;
    UNKNOWN|"")
      PRIMER_RULE_APPLIED="DETECTION_STATUS=${DETECTION_STATUS:-(no info file)} -> --p-trim-left-f 20  --p-trim-left-r 20 (defensive default)"
      ;;
    *)
      PRIMER_RULE_APPLIED="DETECTION_STATUS=$DETECTION_STATUS -> defaults (trim-left=20)"
      ;;
  esac
fi

# ── Summary header ────────────────────────────────────────────────────────────
sep
info "Input artifact     : $QZA_PATH"
info "Input summary      : $QZV_PATH"
info "Output directory   : $DADA2_DIR"
if [[ -n "$PRIMER_INFO_TXT" ]]; then
  info "Primer info file   : $PRIMER_INFO_TXT"
  info "Primer rule        : $PRIMER_RULE_APPLIED"
else
  info "Primer info file   : (not found — using defaults)"
fi
[[ -n "$FORWARD_PRIMER" ]] \
  && info "Forward primer     : $FORWARD_PRIMER  (length: ${#FORWARD_PRIMER} bp → trim-left-f)" \
  || info "Forward primer     : not provided  (trim-left-f will use ${ASSUME_PRIMER_LEN:-20} bp)"
[[ -n "$REVERSE_PRIMER" ]] \
  && info "Reverse primer     : $REVERSE_PRIMER  (length: ${#REVERSE_PRIMER} bp → trim-left-r)" \
  || info "Reverse primer     : not provided  (trim-left-r will use ${ASSUME_PRIMER_LEN:-20} bp)"
[[ -n "$AMPLICON_LENGTH" ]] \
  && info "Amplicon length    : $AMPLICON_LENGTH bp" \
  || info "Amplicon length    : not provided  (will be estimated from data)"
info "Min overlap        : $MIN_OVERLAP bp"
info "Max EE (F/R)       : $MAX_EE_F / $MAX_EE_R"
$DRY_RUN && warn "DRY-RUN mode — commands will be printed but NOT executed."
sep

# ── Dry-run helper — run or just print ───────────────────────────────────────
run_cmd() {
  echo "  \$ $*"
  if ! $DRY_RUN; then
    "$@" || err "Command failed: $*
  → Check the output above for details from QIIME2."
  fi
}

# ─────────────────────────────────────────────────────────────────────────────
step "1/7 — Peek at artifact to verify type"
# ─────────────────────────────────────────────────────────────────────────────

run_cmd qiime tools peek "$QZA_PATH"
[[ -f "$QZV_PATH" ]] && run_cmd qiime tools peek "$QZV_PATH"

if ! $DRY_RUN; then
  PEEK_OUT="$(qiime tools peek "$QZA_PATH" 2>&1 || true)"
  if echo "$PEEK_OUT" | grep -q "SampleData\[PairedEndSequencesWithQuality\]"; then
    ok "Artifact type confirmed: SampleData[PairedEndSequencesWithQuality]"
  else
    err "Artifact is not SampleData[PairedEndSequencesWithQuality].
  Peek output:
$PEEK_OUT
  → This script currently supports paired-end data only.
  → Check artifact type with: qiime tools peek $QZA_PATH"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "2/7 — Export demux summary (.qzv)"
# ─────────────────────────────────────────────────────────────────────────────

EXPORTED_SUMMARY_DIR="${DADA2_DIR}/exported_demux_summary"
mkdir -p "$EXPORTED_SUMMARY_DIR"

if [[ -f "$QZV_PATH" ]]; then
  run_cmd qiime tools export \
    --input-path  "$QZV_PATH" \
    --output-path "$EXPORTED_SUMMARY_DIR"
  ok "Exported demux summary → $EXPORTED_SUMMARY_DIR"
else
  warn "Skipping .qzv export — file not found: $QZV_PATH"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "3/7 — Extract artifact and summary (qza + qzv)"
# ─────────────────────────────────────────────────────────────────────────────

EXTRACTED_QZA_DIR="${DADA2_DIR}/extracted_qza"
EXTRACTED_QZV_DIR="${DADA2_DIR}/extracted_qzv"
mkdir -p "$EXTRACTED_QZA_DIR" "$EXTRACTED_QZV_DIR"

run_cmd qiime tools extract \
  --input-path  "$QZA_PATH" \
  --output-path "$EXTRACTED_QZA_DIR"
ok "Extracted .qza → $EXTRACTED_QZA_DIR"

if [[ -f "$QZV_PATH" ]]; then
  run_cmd qiime tools extract \
    --input-path  "$QZV_PATH" \
    --output-path "$EXTRACTED_QZV_DIR"
  ok "Extracted .qzv → $EXTRACTED_QZV_DIR"
else
  warn "Skipping .qzv extract — file not found: $QZV_PATH"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "4/7 — Generate per-sample read counts"
# ─────────────────────────────────────────────────────────────────────────────

READ_COUNTS_QZA="${DADA2_DIR}/demux_read_counts.qza"

run_cmd qiime demux tabulate-read-counts \
  --i-sequences "$QZA_PATH" \
  --o-counts    "$READ_COUNTS_QZA"
ok "Read counts artifact → $READ_COUNTS_QZA"

# ─────────────────────────────────────────────────────────────────────────────
step "5/7 — Make read counts viewable (.qzv)"
# ─────────────────────────────────────────────────────────────────────────────

READ_COUNTS_QZV="${DADA2_DIR}/demux_read_counts.qzv"

run_cmd qiime metadata tabulate \
  --m-input-file "$READ_COUNTS_QZA" \
  --o-visualization "$READ_COUNTS_QZV"
ok "Read counts visualization → $READ_COUNTS_QZV"

# ─────────────────────────────────────────────────────────────────────────────
step "6/7 — Export and extract read counts"
# ─────────────────────────────────────────────────────────────────────────────

EXPORTED_COUNTS_DIR="${DADA2_DIR}/exported_demux_read_counts"
EXTRACTED_COUNTS_DIR="${DADA2_DIR}/extracted_demux_read_counts"
mkdir -p "$EXPORTED_COUNTS_DIR" "$EXTRACTED_COUNTS_DIR"

run_cmd qiime tools export \
  --input-path  "$READ_COUNTS_QZV" \
  --output-path "$EXPORTED_COUNTS_DIR"
ok "Exported read counts → $EXPORTED_COUNTS_DIR"

run_cmd qiime tools extract \
  --input-path  "$READ_COUNTS_QZA" \
  --output-path "$EXTRACTED_COUNTS_DIR"
ok "Extracted read counts → $EXTRACTED_COUNTS_DIR"

# ─────────────────────────────────────────────────────────────────────────────
step "7/7 — Recommend DADA2 parameters"
# ─────────────────────────────────────────────────────────────────────────────

PARAM_TXT="${DADA2_DIR}/dada2_parameters.txt"

# Build the create_dada2_parameters_txt.sh argument list
PARAM_ARGS=(-i "$QZA_PATH" -o "$PARAM_TXT")

[[ -n "$FORWARD_PRIMER"   ]] && PARAM_ARGS+=(--forward-primer       "$FORWARD_PRIMER")
[[ -n "$REVERSE_PRIMER"   ]] && PARAM_ARGS+=(--reverse-primer       "$REVERSE_PRIMER")
[[ -n "$ASSUME_PRIMER_LEN" ]] && PARAM_ARGS+=(--assume-primer-length "$ASSUME_PRIMER_LEN")
[[ -n "$AMPLICON_LENGTH"  ]] && PARAM_ARGS+=(--amplicon-length      "$AMPLICON_LENGTH")
PARAM_ARGS+=(--min-overlap "$MIN_OVERLAP")
PARAM_ARGS+=(--max-ee-f    "$MAX_EE_F")
PARAM_ARGS+=(--max-ee-r    "$MAX_EE_R")

if $DRY_RUN; then
  echo "  \$ $PARAM_FINDER_CMD ${PARAM_ARGS[*]}"
else
  info "Running parameter recommender — this may take 1-2 minutes..."
  if ! "$PARAM_FINDER_CMD" "${PARAM_ARGS[@]}"; then
    err "create_dada2_parameters_txt.sh failed.
  Common causes:
    1. The .qza artifact is corrupted or not a paired-end demux artifact.
    2. Too few reads in the artifact (need at least a few hundred per sample).
    3. Data quality is extremely low — try: --max-ee-f 5.0 --max-ee-r 5.0
    4. No candidate truncation lengths found — try: --min-overlap 8
  → Parameters file would have been written to: $PARAM_TXT"
  fi
  ok "DADA2 parameters written → $PARAM_TXT"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Final summary
# ─────────────────────────────────────────────────────────────────────────────
sep
ok "All 7 steps completed successfully!"
sep
echo ""
echo "  Output structure:"
echo "  $MBX_OUT_DIR/"
echo "  ├── 1_manifest_file/"
echo "  │   └── manifest.txt"
echo "  ├── 2_first_artifact_file/"
echo "  │   ├── Paired_End_artifact.qza"
echo "  │   └── Paired_End_artifact.qzv"
echo "  └── 3_dada2_parameters/"
echo "      ├── exported_demux_summary/"
echo "      ├── extracted_qza/"
echo "      ├── extracted_qzv/"
echo "      ├── demux_read_counts.qza"
echo "      ├── demux_read_counts.qzv"
echo "      ├── exported_demux_read_counts/"
echo "      ├── extracted_demux_read_counts/"
echo "      └── dada2_parameters.txt     ← your recommended parameters"
echo ""

if ! $DRY_RUN && [[ -f "$PARAM_TXT" ]]; then
  echo "  ── Recommended DADA2 parameters ──────────────────────────────"
  grep "^--p-" "$PARAM_TXT" | while IFS= read -r line; do
    echo "  $line"
  done
  echo ""
  echo "  ── Ready-to-run QIIME2 command ────────────────────────────────"
  grep -A 10 "^qiime dada2" "$PARAM_TXT" | while IFS= read -r line; do
    echo "  $line"
  done
  echo ""
fi

echo "  View .qzv files at: https://view.qiime2.org"
echo ""
