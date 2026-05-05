#!/usr/bin/env bash
# =============================================================================
#  artifact_creator.sh
#  Import a QIIME2 manifest into a .qza artifact
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  Reads the manifest header to auto-detect paired-end vs single-end.
#  Output is placed in the SAME mbX_pro_outputs_* directory that
#  create_manifest.sh created, under a new subfolder: 2_first_artifact_file/
#
#  USAGE:
#    artifact_creator.sh <path/to/manifest.txt>
#    artifact_creator.sh <path/to/manifest.txt> --dry-run
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
  echo "  Need help? Common fixes are listed at the bottom of --help" >&2
  echo "" >&2
  exit 1
}
warn() { echo "[WARN]  $*" >&2; }
info() { echo "[INFO]  $*"; }
ok()   { echo "[OK]    $*"; }
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

artifact_creator.sh — Import a QIIME2 manifest into a .qza artifact

USAGE:
  artifact_creator.sh <manifest.txt> [OPTIONS]

DESCRIPTION:
  Reads the header of the manifest file to automatically detect whether the
  data is paired-end or single-end, then runs qiime tools import with the
  correct parameters.

  The output .qza file is placed INSIDE the same mbX_pro_outputs_* directory
  that create_manifest.sh created, under a new subfolder:

    mbX_pro_outputs_<timestamp>/
    ├── 1_manifest_file/
    │   └── manifest.txt          <- you provide this path
    └── 2_first_artifact_file/
        └── Paired_End_artifact.qza   (or Single_End_artifact.qza)

OPTIONS:
  --dry-run    Print the qiime command that would run, but do not execute it
  -h, --help   Show this help message and exit

EXAMPLES:
  artifact_creator.sh /path/to/mbX_pro_outputs_20250416_143022/1_manifest_file/manifest.txt
  artifact_creator.sh /path/to/manifest.txt --dry-run

COMMON ERRORS AND FIXES:
  "qiime: command not found"
    → Activate your QIIME2 conda environment first:
      conda activate qiime2-amplicon-2025.4

  "manifest file does not exist"
    → Check the path you passed; run create_manifest.sh first if needed.

  "could not detect read type from manifest header"
    → The manifest header must contain 'forward-absolute-filepath' (paired)
      or 'absolute-filepath' (single). Re-run create_manifest.sh to regenerate.

  "parent folder is not named 1_manifest_file"
    → artifact_creator.sh expects the manifest to live inside a
      1_manifest_file/ folder (as created by create_manifest.sh).
      Use --manifest-dir to override the output location if needed.

  "No such file or directory" inside the qiime import step
    → One or more filepaths inside manifest.txt do not exist on this machine.
      Open the manifest and verify the paths are correct and accessible.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MANIFEST_PATH=""
DRY_RUN=false

[[ $# -eq 0 ]] && usage

for arg in "$@"; do
  case "$arg" in
    -h|--help)  usage ;;
    --dry-run)  DRY_RUN=true ;;
    -*)  err "Unknown option: '${arg}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MANIFEST_PATH" ]]; then
        MANIFEST_PATH="$arg"
      else
        err "Unexpected extra argument: '${arg}'  —  only one positional argument (manifest path) is accepted."
      fi
      ;;
  esac
done

[[ -z "$MANIFEST_PATH" ]] && err "No manifest file provided.  Run with --help for usage."
[[ -f "$MANIFEST_PATH" ]] || err "Manifest file does not exist: '${MANIFEST_PATH}'
  → Did you run create_manifest.sh first?
  → Check the path for typos."

MANIFEST_PATH="$(_abspath "$MANIFEST_PATH")"

sep
info "Manifest file      : $MANIFEST_PATH"

# ── Verify QIIME2 is available ────────────────────────────────────────────────
if ! command -v qiime &>/dev/null; then
  err "qiime command not found.
  → Activate your QIIME2 conda environment before running this script:
    conda activate qiime2-amplicon-2025.4"
fi

QIIME_VERSION="$(qiime info 2>/dev/null | grep 'QIIME 2 release' | awk '{print $NF}' || echo "unknown")"
info "QIIME2 version     : $QIIME_VERSION"

# ── Detect read type from manifest header ────────────────────────────────────
HEADER="$(head -1 "$MANIFEST_PATH")"

if echo "$HEADER" | grep -q "forward-absolute-filepath"; then
  READ_TYPE="paired"
elif echo "$HEADER" | grep -q "absolute-filepath"; then
  READ_TYPE="single"
else
  err "Could not detect read type from the manifest header.
  Header found : $HEADER
  Expected     : 'sample-id<TAB>forward-absolute-filepath<TAB>reverse-absolute-filepath'  (paired-end)
              or 'sample-id<TAB>absolute-filepath'                                         (single-end)
  → Re-run create_manifest.sh to regenerate a valid manifest."
fi

info "Detected read type : $READ_TYPE"

# ── Validate manifest has data rows ───────────────────────────────────────────
DATA_ROWS=$(tail -n +2 "$MANIFEST_PATH" | grep -c '[^[:space:]]' || true)
[[ "$DATA_ROWS" -eq 0 ]] && err "The manifest file has a header but no data rows.
  → Re-run create_manifest.sh — no samples were written to the manifest."
info "Sample rows found  : $DATA_ROWS"

# ── Validate filepaths inside the manifest ────────────────────────────────────
info "Validating filepaths inside manifest..."
MISSING_COUNT=0
while IFS=$'\t' read -r _ col2 col3; do
  for fp in "$col2" "${col3:-}"; do
    [[ -z "$fp" ]] && continue
    if [[ ! -f "$fp" ]]; then
      warn "File not found: $fp"
      MISSING_COUNT=$(( MISSING_COUNT + 1 ))
    fi
  done
done < <(tail -n +2 "$MANIFEST_PATH")

if [[ "$MISSING_COUNT" -gt 0 ]]; then
  err "$MISSING_COUNT filepath(s) in the manifest do not exist on disk.
  → Open the manifest and verify all paths are correct and the files are accessible.
  → If files are on an external drive, make sure it is mounted."
fi
ok "All filepaths in manifest are accessible."

# ── Resolve output directory ──────────────────────────────────────────────────
# Expected structure from create_manifest.sh:
#   mbX_pro_outputs_<timestamp>/
#   └── 1_manifest_file/
#       └── manifest.txt     <- MANIFEST_PATH
#
# We walk up two levels to find mbX_pro_outputs_* and place output alongside.

MANIFEST_PARENT_DIR="$(dirname "$MANIFEST_PATH")"          # 1_manifest_file/
MANIFEST_PARENT_NAME="$(basename "$MANIFEST_PARENT_DIR")"  # should be "1_manifest_file"
MBX_OUT_DIR="$(dirname "$MANIFEST_PARENT_DIR")"            # mbX_pro_outputs_*/

if [[ "$MANIFEST_PARENT_NAME" != "1_manifest_file" ]]; then
  warn "The manifest's parent folder is '${MANIFEST_PARENT_NAME}', not '1_manifest_file'."
  warn "Expected the structure created by create_manifest.sh."
  warn "Output will still be placed alongside the manifest's parent folder."
fi

if [[ ! -d "$MBX_OUT_DIR" ]]; then
  err "Could not locate the mbX_pro_outputs directory at: $MBX_OUT_DIR
  → Make sure you pass the manifest.txt path exactly as produced by create_manifest.sh."
fi

ARTIFACT_DIR="${MBX_OUT_DIR}/2_first_artifact_file"
mkdir -p "$ARTIFACT_DIR" \
  || err "Could not create output directory: '${ARTIFACT_DIR}'  —  check write permissions."

info "Output directory   : $ARTIFACT_DIR"

# ── Set QIIME2 import parameters based on read type ──────────────────────────
if [[ "$READ_TYPE" == "paired" ]]; then
  QIIME_TYPE="SampleData[PairedEndSequencesWithQuality]"
  QIIME_FORMAT="PairedEndFastqManifestPhred33V2"
  OUTPUT_QZA="${ARTIFACT_DIR}/Paired_End_artifact.qza"
  OUTPUT_QZV="${ARTIFACT_DIR}/Paired_End_artifact.qzv"
else
  QIIME_TYPE="SampleData[SequencesWithQuality]"
  QIIME_FORMAT="SingleEndFastqManifestPhred33V2"
  OUTPUT_QZA="${ARTIFACT_DIR}/Single_End_artifact.qza"
  OUTPUT_QZV="${ARTIFACT_DIR}/Single_End_artifact.qzv"
fi

info "QIIME2 type        : $QIIME_TYPE"
info "QIIME2 format      : $QIIME_FORMAT"
info "Output artifact    : $OUTPUT_QZA"
info "Output summary     : $OUTPUT_QZV"
sep

# ── Build the qiime command ───────────────────────────────────────────────────
QIIME_CMD="qiime tools import \
  --type \"${QIIME_TYPE}\" \
  --input-format ${QIIME_FORMAT} \
  --input-path  ${MANIFEST_PATH} \
  --output-path ${OUTPUT_QZA}"

echo ""
echo "  Commands to run:"
echo "  ┌─────────────────────────────────────────────────────────────"
echo "  │  Step 1 — Import FASTQ files into a QIIME2 artifact"
printf '  │  qiime tools import \\\n'
printf '  │    --type "%s" \\\n'       "$QIIME_TYPE"
printf '  │    --input-format %s \\\n' "$QIIME_FORMAT"
printf '  │    --input-path  %s \\\n'  "$MANIFEST_PATH"
printf '  │    --output-path %s\n'     "$OUTPUT_QZA"
echo "  │"
echo "  │  Step 2 — Summarize demultiplexing into a visualization"
printf '  │  qiime demux summarize \\\n'
printf '  │    --i-data %s \\\n'           "$OUTPUT_QZA"
printf '  │    --o-visualization %s\n'     "$OUTPUT_QZV"
echo "  └─────────────────────────────────────────────────────────────"
echo ""

# ── Dry run exit ──────────────────────────────────────────────────────────────
if $DRY_RUN; then
  warn "Dry-run mode: command was NOT executed."
  sep
  exit 0
fi

# ── Run QIIME2 import ─────────────────────────────────────────────────────────
info "Step 1/2 — Running QIIME2 import — this may take a few minutes..."
sep

if ! qiime tools import \
      --type "$QIIME_TYPE" \
      --input-format "$QIIME_FORMAT" \
      --input-path  "$MANIFEST_PATH" \
      --output-path "$OUTPUT_QZA"; then

  err "qiime tools import failed.
  Common causes:
    1. One or more FASTQ files listed in the manifest are missing or unreadable.
    2. The FASTQ files are corrupted or not properly gzipped.
    3. The Phred encoding does not match (Phred33 is assumed — check your data).
    4. Insufficient disk space for the output .qza file.
  → Re-run with --dry-run to inspect the exact command being used.
  → Run: qiime tools import --help   for full parameter documentation."
fi

ok "Import complete: $OUTPUT_QZA"
sep

# ── Run QIIME2 demux summarize ────────────────────────────────────────────────
info "Step 2/2 — Running demux summarize to generate visualization..."
sep

if ! qiime demux summarize \
      --i-data "$OUTPUT_QZA" \
      --o-visualization "$OUTPUT_QZV"; then

  err "qiime demux summarize failed.
  Common causes:
    1. The .qza artifact was not created correctly in Step 1.
    2. Insufficient disk space for the .qzv file.
    3. Corrupted or truncated FASTQ data inside the artifact.
  → The .qza file was successfully created at: $OUTPUT_QZA
  → You can retry the summarize step manually:
    qiime demux summarize \\
      --i-data $OUTPUT_QZA \\
      --o-visualization $OUTPUT_QZV"
fi

# ── Final summary ─────────────────────────────────────────────────────────────
sep
ok "All steps completed successfully!"
ok "Read type           : $READ_TYPE"
ok "Samples imported    : $DATA_ROWS"
ok "Artifact  (.qza)    : $OUTPUT_QZA"
ok "Summary   (.qzv)    : $OUTPUT_QZV"
sep
echo ""
echo "  Output structure:"
echo "  $(dirname "$ARTIFACT_DIR")/"
echo "  ├── 1_manifest_file/"
echo "  │   └── manifest.txt"
echo "  └── 2_first_artifact_file/"
if [[ "$READ_TYPE" == "paired" ]]; then
  echo "      ├── Paired_End_artifact.qza"
  echo "      └── Paired_End_artifact.qzv"
else
  echo "      ├── Single_End_artifact.qza"
  echo "      └── Single_End_artifact.qzv"
fi
echo ""
echo "  To view the summary, drag the .qzv file into:"
echo "  https://view.qiime2.org"
echo ""
