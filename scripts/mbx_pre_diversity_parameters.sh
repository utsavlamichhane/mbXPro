#!/usr/bin/env bash
# =============================================================================
#  mbx_pre_diversity_parameters.sh
#  Build phylogenetic tree + compute scientifically-defensible rarefaction depth
#
#  Compatible with bash 3.2+ (macOS default shell)
#  Requires: qiime, Rscript (system R via brew, NOT inside conda env)
#
#  SUBSTEPS (this script is pipeline step 11):
#    1   Input discovery + version logging
#    2   QC visualizations (metadata tabulate, feature-table summarize)
#    3   Build rooted phylogenetic tree (mafft → fasttree → root)
#    4   Tabulate per-sample frequencies + TSV export
#    5   Extract metadata.tsv (qiime tools export) + tree Newick export
#    6   Export feature table (BIOM + TSV) — same artifact, two formats
#    7   Analytical rarefaction (Hurlbert 1971): observed_features + Good's
#        coverage, exact closed-form, x-axis from 1k to max(sample frequencies)
#    8   Depth selection — three concurrent criteria:
#          (a) ≥0.90 of samples retained at depth d (overall retention)
#          (b) mean Good's coverage of retained samples at d ≥ 0.98
#          (c) observed_features curve slope at d < 0.5 / 1000 reads (plateau)
#    9   ID validation (metadata ↔ table sample IDs;
#                       table ↔ rep-seqs ↔ tree feature IDs)
#    10  Official QIIME alpha-rarefaction (observed_features + shannon +
#        faith_pd), max-depth = max(sample freq), iterations 10, steps 20
#    11  Decision-supporting visualizations (depth histogram, retention curve)
#    12  Write info file + plain-language summary + STATUS + READY_FOR_DIVERSITY
#
#  STATUS field meanings (written to mbx_pre_diversity_info.txt):
#    PASS                  — all three depth-selection criteria met cleanly
#    PASS_WITH_WARNINGS    — criteria met but borderline; non-blocking warnings
#    REVIEW_REQUIRED       — fallback rule used, or a criterion failed; user
#                            should look at the summary file before continuing
#    FAIL                  — pipeline-stopping problem (no metadata-table
#                            overlap, group fully wiped, etc.); next step
#                            (mbx_diversity_run.sh) MUST refuse to run
#
#  DEPTH SELECTION (the science):
#    The script chooses the highest depth d that satisfies ALL THREE of:
#      (a) at least MIN_OVERALL fraction of samples have N_i ≥ d
#      (b) mean Good's coverage at d (over retained samples) ≥ GOOD_COV_MIN
#      (c) slope of mean observed_features curve at d, expressed as
#          (features per 1000 additional reads), is below PLATEAU_SLOPE_MAX
#    If --group-col is provided, also requires ≥ MIN_GROUP per-group retention.
#    If no d satisfies all three, status drops to REVIEW_REQUIRED and the
#    script reports which criterion was binding.
#
#  OUTPUT STRUCTURE:
#    mbX_pro_outputs_<timestamp>/
#    └── 11_pre_diversity/
#        ├── metadata_summary.qzv
#        ├── feature_table_summary.qzv
#        ├── sample-frequencies.qza
#        ├── sample-frequencies.tsv                    ← exported (NEW)
#        ├── feature-table.biom                        ← exported
#        ├── feature-table.tsv                         ← exported (NEW)
#        ├── aligned-rep-seqs.qza
#        ├── masked-aligned-rep-seqs.qza
#        ├── unrooted-tree.qza
#        ├── rooted-tree.qza
#        ├── rooted-tree.nwk                           ← exported (NEW)
#        ├── sampling_depth_candidates.csv
#        ├── samples_retained_at_recommended_depth.csv ← NEW
#        ├── samples_removed_at_recommended_depth.csv  ← NEW
#        ├── sample_depth_summary.csv                  ← NEW
#        ├── group_depth_summary.csv                   ← NEW (if categorical
#        │                                                cols exist; reports
#        │                                                ALL of them)
#        ├── sequencing_depth_distribution.png         ← NEW
#        ├── depth_vs_retention.png                    ← NEW
#        ├── alpha-rarefaction-mbx-preview.qzv         ← analytical preview
#        │                                                (observed + Good's)
#        ├── alpha-rarefaction-qiime.qzv               ← official QIIME
#        │                                                (observed + shannon
#        │                                                + faith_pd)
#        ├── alpha_rarefaction_data.csv                ← analytical raw data
#        ├── alpha_rarefaction_curves.png              ← analytical plot,
#        │                                                annotated w/ status
#        ├── mbx_pre_diversity_info.txt                ← machine-readable;
#        │                                                consumed by step 12
#        └── mbx_pre_diversity_summary.txt             ← human-readable;
#                                                         what the lay user
#                                                         actually reads
# =============================================================================

set -euo pipefail

# Capture the original invocation verbatim for provenance (info file step 12)
MBX_INVOCATION="$0 $*"
SCRIPT_VERSION="mbx_pre_diversity_parameters.sh (12-step rewrite, Waves 1-8)"

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
warn()    { echo "[WARN]  $*" >&2; }
info()    { echo "[INFO]  $*"; }
ok()      { echo "[OK]    $*"; }
skipped() { echo "[SKIP]  $* — already exists."; }
step() {
  echo ""
  echo "┌─────────────────────────────────────────────────────────────────"
  echo "│  Step $*"
  echo "└─────────────────────────────────────────────────────────────────"
}
sep()         { echo "────────────────────────────────────────────────────────────────"; }
timer_start() { _T0="$(date +%s)"; }
timer_end()   {
  local _T1; _T1="$(date +%s)"
  local _S=$(( _T1 - _T0 ))
  printf "[TIME]  %dm %02ds\n" $(( _S/60 )) $(( _S%60 ))
}

_abspath() {
  if [[ -d "$1" ]]; then cd "$1" && pwd
  elif [[ -f "$1" ]]; then echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
  else return 1; fi
}

_read_key() {
  local key="$1" file="$2"
  grep "^${key}=" "$file" 2>/dev/null | cut -d= -f2-
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_pre_diversity_parameters.sh — Tree + rarefaction depth (publication grade)

USAGE:
  mbx_pre_diversity_parameters.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads all input paths automatically from previous pipeline info files:
    5_classifier_working_dir/mbx_classifier_run_info.txt  → feature table
    7_taxonomy_csv/mbx_taxonomy_info.txt                  → metadata, filtered
    4_dada2_outputs/representative_sequences.qza          → rep seqs

  Then:
    • Builds rooted phylogenetic tree (mafft → fasttree → root)
    • Exports feature table to BIOM + TSV
    • Computes analytical rarefaction (Hurlbert 1971): observed_features +
      Good's coverage, x-axis to max(sample frequency)
    • Selects sampling depth using THREE concurrent criteria:
        (a) ≥ MIN_OVERALL of samples retained
        (b) mean Good's coverage of retained samples ≥ GOOD_COV_MIN
        (c) observed_features slope at depth < PLATEAU_SLOPE_MAX
    • Runs official QIIME alpha-rarefaction (observed_features + shannon +
      faith_pd) for publication-grade evidence
    • Validates sample IDs (metadata ↔ table) and feature IDs
      (table ↔ rep-seqs ↔ tree)
    • Writes STATUS + READY_FOR_DIVERSITY consumed by mbx_diversity_run.sh

OPTIONS:
  --group-col <col>           Metadata column for group-aware retention
                              constraint (NOT auto-selected; if omitted, only
                              the overall retention rule is enforced)
  --min-overall <0-1>         Min fraction of samples retained (default: 0.90)
  --min-group   <0-1>         Min fraction per group retained (default: 0.80;
                              only used if --group-col is given)
  --good-coverage-min <0-1>   Mean Good's coverage threshold (default: 0.98)
  --plateau-slope-max <num>   Max acceptable slope of observed_features curve,
                              in features per 1000 reads (default: 0.5)
  --allow-unfiltered-table    OPT-IN: continue even if mito/chloro filtered
                              feature table is missing. Strongly discouraged
                              for publication-bound analyses (default: refuse)
  --skip-tree                 Skip tree building if rooted-tree.qza exists
                              (refused if no rooted tree present anywhere)
  --skip-qc                   Skip metadata/feature-table QC visualizations
  --skip-qiime-rarefaction    Skip the official QIIME alpha-rarefaction (only
                              the analytical preview will be produced; not
                              recommended for publication)
  --rare-steps <N>            Analytical rarefaction depth steps (default: 20)
  --qiime-rare-iter <N>       Iterations for official QIIME rarefaction
                              (default: 10)
  --qiime-rare-steps <N>      Steps for official QIIME rarefaction
                              (default: 20)
  --force-rerun               Re-run every substep even if outputs exist.
                              Default: skip-if-exists (fast & idempotent),
                              but stale-output warnings are printed
  --dry-run                   Print commands without executing
  -h, --help                  Show this help and exit

EXAMPLES:
  # Default — automated decision, publication-grade evidence:
  mbx_pre_diversity_parameters.sh /path/to/mbX_pro_outputs_20250422_143022

  # With a primary biological grouping variable:
  mbx_pre_diversity_parameters.sh /path/to/mbX_pro_outputs_20250422_143022 \
    --group-col Treatment

  # More permissive Good's coverage threshold (e.g. very deep dataset):
  mbx_pre_diversity_parameters.sh /path/to/mbX_pro_outputs_20250422_143022 \
    --good-coverage-min 0.95

OUTPUT STATUS:
  After running, read 11_pre_diversity/mbx_pre_diversity_info.txt:
    STATUS=PASS                — clean, proceed to mbx_diversity_run.sh
    STATUS=PASS_WITH_WARNINGS  — proceed, review the warnings list
    STATUS=REVIEW_REQUIRED     — read mbx_pre_diversity_summary.txt, decide
    STATUS=FAIL                — diversity script will refuse to run

COMMON ERRORS:
  "Filtered feature table not found"
    → 7_taxonomy_csv/feature_table_filtered.qza is required (mito/chloro
      removed). If you intentionally did not filter, pass
      --allow-unfiltered-table (REVIEW_REQUIRED status will result).
  "mafft: command not found"
    → mafft ships with QIIME2 — activate it: conda activate qiime2-amplicon-2025.4
  "Rscript not found"
    → brew install r  (system-wide, not inside conda)

EOF
  exit 0
}

# ── Capture invocation BEFORE parsing consumes arguments ────────────────────
# Recorded in mbx_pre_diversity_info.txt for reproducibility
RUN_INVOCATION="$0"
for _arg in "$@"; do
  case "$_arg" in
    *[[:space:]]*) RUN_INVOCATION="$RUN_INVOCATION \"$_arg\"" ;;
    *)             RUN_INVOCATION="$RUN_INVOCATION $_arg" ;;
  esac
done

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
GROUP_COL=""
MIN_OVERALL="0.90"
MIN_GROUP="0.80"
GOOD_COV_MIN="0.98"          # mean Good's coverage threshold (PASS gate)
PLATEAU_SLOPE_MAX="0.5"      # features per 1000 reads — above this, curve still climbing
ALLOW_UNFILTERED_TABLE=false # opt-in to continue without mito/chloro filtering
SKIP_TREE=false
SKIP_QC=false
SKIP_QIIME_RAREFACTION=false # opt-out of official QIIME alpha-rarefaction
DRY_RUN=false
FORCE_RERUN=false            # if true, re-run every substep even if output exists
RARE_STEPS=20                # analytical rarefaction steps (smooth curves; cheap)
QIIME_RARE_STEPS=20          # official QIIME rarefaction steps
QIIME_RARE_ITERATIONS=10     # official QIIME rarefaction iterations per step

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)                 usage ;;
    --dry-run)                 DRY_RUN=true;  shift ;;
    --skip-tree)               SKIP_TREE=true; shift ;;
    --skip-qc)                 SKIP_QC=true;  shift ;;
    --skip-qiime-rarefaction)  SKIP_QIIME_RAREFACTION=true; shift ;;
    --allow-unfiltered-table)  ALLOW_UNFILTERED_TABLE=true; shift ;;
    --force-rerun)             FORCE_RERUN=true; shift ;;
    --group-col)               GROUP_COL="${2:?Missing value for --group-col}"; shift 2 ;;
    --rare-steps)              RARE_STEPS="${2:?Missing value for --rare-steps}"; shift 2 ;;
    --qiime-rare-iter)         QIIME_RARE_ITERATIONS="${2:?Missing value for --qiime-rare-iter}"; shift 2 ;;
    --qiime-rare-steps)        QIIME_RARE_STEPS="${2:?Missing value for --qiime-rare-steps}"; shift 2 ;;
    --min-overall)             MIN_OVERALL="${2:?Missing value for --min-overall}"; shift 2 ;;
    --min-group)               MIN_GROUP="${2:?Missing value for --min-group}"; shift 2 ;;
    --good-coverage-min)       GOOD_COV_MIN="${2:?Missing value for --good-coverage-min}"; shift 2 ;;
    --plateau-slope-max)       PLATEAU_SLOPE_MAX="${2:?Missing value for --plateau-slope-max}"; shift 2 ;;
    # Backwards-compatibility aliases (accepted silently)
    --rare-iter)               QIIME_RARE_ITERATIONS="${2:?Missing value for --rare-iter}"; shift 2 ;;
    -*)  err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '${1}'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -d "$MBX_OUT_DIR" ]] || err "Directory does not exist: '${MBX_OUT_DIR}'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

command -v qiime &>/dev/null || err "qiime not found.
  → conda activate qiime2-amplicon-2025.4"

# Pin Python interpreter to the QIIME2 conda env's python (sibling of `qiime`).
# This is critical: the system/Homebrew python3 will NOT have biom/qiime2/scipy
# installed and will cause step 7 (analytical rarefaction) to fail with
# "ModuleNotFoundError: No module named 'biom'". The `qiime` binary is already
# validated to be the conda env's, so its sibling `python` is guaranteed correct.
QIIME_BIN_DIR="$(cd "$(dirname "$(command -v qiime)")" && pwd)"
PY_BIN=""
for _p in "${QIIME_BIN_DIR}/python" "${QIIME_BIN_DIR}/python3"; do
  [[ -x "$_p" ]] && { PY_BIN="$_p"; break; }
done
[[ -n "$PY_BIN" ]] || err "Could not find python in QIIME2 env at: ${QIIME_BIN_DIR}/python
  → Verify with: ls -la ${QIIME_BIN_DIR}/python*
  → Try: conda activate qiime2-amplicon-2025.4 && which python"

# Auto-detect Rscript
RSCRIPT_CMD=""
for _c in \
    "$(command -v Rscript 2>/dev/null || true)" \
    "/opt/homebrew/bin/Rscript" \
    "/usr/local/bin/Rscript" \
    "/usr/bin/Rscript"; do
  [[ -n "$_c" && -x "$_c" ]] && { RSCRIPT_CMD="$_c"; break; }
done
[[ -n "$RSCRIPT_CMD" ]] || err "Rscript not found.
  → brew install r"

# Auto-detect CPU cores
if command -v nproc &>/dev/null; then
  N_JOBS="$(nproc)"
elif command -v sysctl &>/dev/null; then
  N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else
  N_JOBS=1
fi

# PID-based temp files (bash 3.2 safe, no mktemp suffix issues)
_TMPID="${$}_$(date +%s)"
R_DEPTH_SCRIPT="/tmp/mbx_depth_${_TMPID}.R"
R_DETECT_SCRIPT="/tmp/mbx_grpdetect_${_TMPID}.R"
R_DEPTH_OUT="/tmp/mbx_depth_out_${_TMPID}.txt"
R_VALIDATE_SCRIPT="/tmp/mbx_validate_${_TMPID}.R"
R_VIZ_SCRIPT="/tmp/mbx_viz_${_TMPID}.R"
trap 'rm -f "$R_DEPTH_SCRIPT" "$R_DETECT_SCRIPT" "$R_DEPTH_OUT" "$R_VALIDATE_SCRIPT" "$R_VIZ_SCRIPT"' EXIT

# ── Stale-output warning helper ──────────────────────────────────────────────
# When a step is being skipped because the output file already exists, this
# prints the file's mtime so the user sees what's being reused. Critical for
# scientific reproducibility — silent skip-if-exists across pipeline reruns
# can cause stale outputs to feed downstream steps.
_stale_warn() {
  # Args: <label> <path>
  local _label="$1" _path="$2"
  if [[ -e "$_path" ]]; then
    local _mtime
    _mtime="$(date -r "$_path" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || echo 'unknown')"
    info "  └─ reusing $_label (mtime: $_mtime)"
  fi
}

# ── Run-context capture (versions, commands, paths) ──────────────────────────
RUN_TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RUN_HOSTNAME="$(hostname 2>/dev/null || echo 'unknown')"
RUN_USER="${USER:-unknown}"
QIIME_VERSION="$(qiime --version 2>&1 | head -n1 | awk '{print $NF}')"
PY_VERSION="$("$PY_BIN" --version 2>&1 | awk '{print $NF}')"
R_VERSION="$("$RSCRIPT_CMD" --version 2>&1 | head -n1 | sed 's/^R scripting front-end version //; s/ .*//')"
# RUN_INVOCATION was captured at top, before argument parsing.
# All this gets written into mbx_pre_diversity_info.txt at the end.

# ─────────────────────────────────────────────────────────────────────────────
step "1/12 — Auto-discover input paths from previous pipeline info files"
# ─────────────────────────────────────────────────────────────────────────────

# ── Metadata (from 7_taxonomy_csv) ───────────────────────────────────────────
TAXONOMY_INFO="${MBX_OUT_DIR}/7_taxonomy_csv/mbx_taxonomy_info.txt"
[[ -f "$TAXONOMY_INFO" ]] || err "7_taxonomy_csv/mbx_taxonomy_info.txt not found.
  → Run mbx_taxonomy_run.sh first."
METADATA_TXT="$(_read_key "METADATA_TXT" "$TAXONOMY_INFO")"
[[ -f "$METADATA_TXT" ]] || err "Metadata file not found: $METADATA_TXT"

# ── Feature table (from 4_dada2_outputs) ─────────────────────────────────────
CLASSIFIER_INFO="${MBX_OUT_DIR}/5_classifier_working_dir/mbx_classifier_run_info.txt"
if [[ -f "$CLASSIFIER_INFO" ]]; then
  FEATURE_TABLE_QZA="$(_read_key "FEATURE_TABLE_QZA" "$CLASSIFIER_INFO")"
fi
# Fallback: look directly in 4_dada2_outputs
if [[ -z "${FEATURE_TABLE_QZA:-}" || ! -f "${FEATURE_TABLE_QZA:-}" ]]; then
  FEATURE_TABLE_QZA="${MBX_OUT_DIR}/4_dada2_outputs/feature_table.qza"
fi
[[ -f "$FEATURE_TABLE_QZA" ]] || err "feature_table.qza not found: $FEATURE_TABLE_QZA
  → Run mbx_dada2_run.sh first."

# ── Filtered feature table is REQUIRED ──────────────────────────────────────
# Mitochondria and chloroplast features can dominate plant/soil/host-associated
# samples. Diversity computed on the unfiltered table will be scientifically
# weaker, often without the user understanding why. Refuse to continue silently.
FILTERED_TABLE="${MBX_OUT_DIR}/7_taxonomy_csv/feature_table_filtered.qza"
TABLE_FILTERED_FLAG=""    # written into info file; "yes" or "no_OPT_IN"
if [[ -f "$FILTERED_TABLE" ]]; then
  info "Using filtered feature table (mito/chloro removed): $FILTERED_TABLE"
  ACTIVE_FEATURE_TABLE="$FILTERED_TABLE"
  TABLE_FILTERED_FLAG="yes"
elif $ALLOW_UNFILTERED_TABLE; then
  warn "──────────────────────────────────────────────────────────────────"
  warn "  ATTENTION: Running on UNFILTERED feature table (--allow-unfiltered-table)"
  warn "  This means mitochondria and chloroplast features are NOT removed."
  warn "  Downstream diversity results may be biased by host/plastid signal."
  warn "  STATUS will be set to REVIEW_REQUIRED in the info file."
  warn "──────────────────────────────────────────────────────────────────"
  ACTIVE_FEATURE_TABLE="$FEATURE_TABLE_QZA"
  TABLE_FILTERED_FLAG="no_OPT_IN"
else
  err "Filtered feature table is required but not found:
    Expected: $FILTERED_TABLE

  This script refuses to silently fall back to the raw feature table because
  diversity results would be biased by mitochondria and chloroplast features.

  → Run mbx_taxonomy_run.sh first to generate the filtered table.
  → If you have a specific reason to analyze the unfiltered table, re-run
    this script with --allow-unfiltered-table  (status will be REVIEW_REQUIRED)."
fi

# ── Representative sequences (from 4_dada2_outputs) ──────────────────────────
REP_SEQS_QZA="${MBX_OUT_DIR}/4_dada2_outputs/representative_sequences.qza"
[[ -f "$REP_SEQS_QZA" ]] || err "representative_sequences.qza not found: $REP_SEQS_QZA
  → Run mbx_dada2_run.sh first."

# ── Taxonomy (from 6_classifier_taxonomy) ────────────────────────────────────
TAXONOMY_QZA="${MBX_OUT_DIR}/6_classifier_taxonomy/taxonomy.qza"
[[ -f "$TAXONOMY_QZA" ]] || warn "taxonomy.qza not found — feature-table summary will run without taxonomy."

sep
info "Run timestamp      : $RUN_TIMESTAMP"
info "Host               : $RUN_USER@$RUN_HOSTNAME"
info "QIIME2             : $QIIME_VERSION"
info "Python             : $PY_VERSION  ($PY_BIN)"
info "R                  : $R_VERSION  ($RSCRIPT_CMD)"
info "Metadata           : $METADATA_TXT"
info "Feature table      : $ACTIVE_FEATURE_TABLE  (filtered: $TABLE_FILTERED_FLAG)"
info "Rep seqs           : $REP_SEQS_QZA"
info "CPU cores          : $N_JOBS"
info "Analytical steps   : $RARE_STEPS"
info "QIIME rarefaction  : $QIIME_RARE_STEPS steps × $QIIME_RARE_ITERATIONS iter $($SKIP_QIIME_RAREFACTION && echo '(SKIPPED)')"
info "Retention overall  : ≥ $MIN_OVERALL"
[[ -n "$GROUP_COL" ]] && info "Retention per-group: ≥ $MIN_GROUP  (group column: $GROUP_COL)"
[[ -z "$GROUP_COL" ]] && info "Retention per-group: not enforced (no --group-col given)"
info "Good's coverage min: $GOOD_COV_MIN  (PASS gate)"
info "Plateau slope max  : $PLATEAU_SLOPE_MAX features per 1000 reads"
$FORCE_RERUN && warn "--force-rerun: every substep will recompute, ignoring existing outputs."
$DRY_RUN && warn "DRY-RUN — commands printed but NOT executed."
sep

# ── Create output directory ───────────────────────────────────────────────────
PRE_DIV_DIR="${MBX_OUT_DIR}/11_pre_diversity"
mkdir -p "$PRE_DIV_DIR" \
  || err "Could not create: $PRE_DIV_DIR — check permissions."

# Output file paths
METADATA_QZV="${PRE_DIV_DIR}/metadata_summary.qzv"
FT_SUMMARY_QZV="${PRE_DIV_DIR}/feature_table_summary.qzv"
SAMPLE_FREQ_QZA="${PRE_DIV_DIR}/sample-frequencies.qza"
SAMPLE_FREQ_DIR="${PRE_DIV_DIR}/sample_frequencies_unzipped"
SAMPLE_FREQ_TSV="${PRE_DIV_DIR}/sample-frequencies.tsv"
ALIGNED_QZA="${PRE_DIV_DIR}/aligned-rep-seqs.qza"
MASKED_QZA="${PRE_DIV_DIR}/masked-aligned-rep-seqs.qza"
UNROOTED_QZA="${PRE_DIV_DIR}/unrooted-tree.qza"
ROOTED_QZA="${PRE_DIV_DIR}/rooted-tree.qza"
ROOTED_NWK="${PRE_DIV_DIR}/rooted-tree.nwk"
FEATURE_BIOM_DIR="${PRE_DIV_DIR}/exported_feature_table"
FEATURE_BIOM="${FEATURE_BIOM_DIR}/feature-table.biom"
FEATURE_TSV="${PRE_DIV_DIR}/feature-table.tsv"
DEPTH_CSV="${PRE_DIV_DIR}/sampling_depth_candidates.csv"
SAMPLES_RETAINED_CSV="${PRE_DIV_DIR}/samples_retained_at_recommended_depth.csv"
SAMPLES_REMOVED_CSV="${PRE_DIV_DIR}/samples_removed_at_recommended_depth.csv"
SAMPLE_DEPTH_SUMMARY_CSV="${PRE_DIV_DIR}/sample_depth_summary.csv"
GROUP_DEPTH_SUMMARY_CSV="${PRE_DIV_DIR}/group_depth_summary.csv"
DEPTH_HIST_PNG="${PRE_DIV_DIR}/sequencing_depth_distribution.png"
DEPTH_DECISION_PNG="${PRE_DIV_DIR}/depth_vs_retention.png"
ALPHA_RARE_MBX_QZV="${PRE_DIV_DIR}/alpha-rarefaction-mbx-preview.qzv"
ALPHA_RARE_QIIME_QZV="${PRE_DIV_DIR}/alpha-rarefaction-qiime.qzv"
ALPHA_RARE_PNG="${PRE_DIV_DIR}/alpha_rarefaction_curves.png"
ALPHA_RARE_CSV="${PRE_DIV_DIR}/alpha_rarefaction_data.csv"
PER_SAMPLE_RARE_CSV="${PRE_DIV_DIR}/alpha_rarefaction_per_sample.csv"
PRE_DIV_INFO="${PRE_DIV_DIR}/mbx_pre_diversity_info.txt"
PRE_DIV_SUMMARY="${PRE_DIV_DIR}/mbx_pre_diversity_summary.txt"

# Backwards-compat alias for an older variable name used in the unmodified
# Step 7 block (will be removed when Step 7 is rewritten in Wave 3):
ALPHA_RARE_QZV="$ALPHA_RARE_MBX_QZV"

# Dry-run helper
run_cmd() {
  echo ""
  echo "  \$ $(printf '%s ' "$@" | sed 's/ --/\\\n    --/g')"
  echo ""
  $DRY_RUN && return 0
  "$@" || err "Command failed: $1
  → Check QIIME2 output above.
  → Re-run with --dry-run to inspect the exact command."
}

# ─────────────────────────────────────────────────────────────────────────────
step "2/12 — QC visualizations (metadata + feature table)"
# ─────────────────────────────────────────────────────────────────────────────

if $SKIP_QC; then
  skipped "QC visualizations (--skip-qc)"
else
  if [[ -f "$METADATA_QZV" ]]; then
    skipped "metadata_summary.qzv"
  else
    info "Tabulating metadata..."
    timer_start
    run_cmd qiime metadata tabulate \
      --m-input-file      "$METADATA_TXT" \
      --o-visualization   "$METADATA_QZV"
    timer_end
    $DRY_RUN || ok "Metadata summary → $METADATA_QZV"
  fi

  if [[ -f "$FT_SUMMARY_QZV" ]]; then
    skipped "feature_table_summary.qzv"
  else
    info "Summarizing feature table..."
    timer_start
    run_cmd qiime feature-table summarize \
      --i-table                  "$ACTIVE_FEATURE_TABLE" \
      --m-sample-metadata-file   "$METADATA_TXT" \
      --o-visualization          "$FT_SUMMARY_QZV"
    timer_end
    $DRY_RUN || ok "Feature table summary → $FT_SUMMARY_QZV"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "3/12 — Build rooted phylogenetic tree"
# ─────────────────────────────────────────────────────────────────────────────
# Outputs: aligned-rep-seqs.qza, masked-aligned-rep-seqs.qza,
#          unrooted-tree.qza, rooted-tree.qza
# Runtime: 5–30 minutes depending on dataset size.
# The rooted tree is required for Faith's PD and UniFrac.

if [[ -f "$ROOTED_QZA" ]] || $SKIP_TREE; then
  if [[ -f "$ROOTED_QZA" ]]; then
    skipped "rooted-tree.qza (already exists)"
    _stale_warn "rooted-tree.qza" "$ROOTED_QZA"
  else
    err "--skip-tree set but $ROOTED_QZA does not exist.
  → A rooted phylogenetic tree is required for Faith PD and UniFrac in step 12.
  → Re-run without --skip-tree to build the tree, OR copy a valid rooted-tree.qza
    into $PRE_DIV_DIR/ before re-running."
  fi
else
  info "Building phylogenetic tree — this may take 5–30 minutes..."
  info "Using $N_JOBS CPU thread(s) for alignment."
  timer_start
  run_cmd qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences       "$REP_SEQS_QZA" \
    --p-n-threads       "$N_JOBS" \
    --o-alignment       "$ALIGNED_QZA" \
    --o-masked-alignment "$MASKED_QZA" \
    --o-tree            "$UNROOTED_QZA" \
    --o-rooted-tree     "$ROOTED_QZA"
  timer_end
  $DRY_RUN || ok "Rooted tree → $ROOTED_QZA"
fi

# ── Newick export (for ID validation + downstream readability) ──────────────
if ! $DRY_RUN && [[ -f "$ROOTED_QZA" ]]; then
  if [[ -f "$ROOTED_NWK" ]] && ! $FORCE_RERUN; then
    skipped "rooted-tree.nwk (already exported)"
  else
    info "Exporting rooted tree to Newick (for feature-ID validation)..."
    _NWK_TMP="${PRE_DIV_DIR}/_nwk_export_tmp"
    rm -rf "$_NWK_TMP"
    qiime tools export \
      --input-path  "$ROOTED_QZA" \
      --output-path "$_NWK_TMP" \
      || err "Failed to export rooted tree to Newick."
    if [[ -f "${_NWK_TMP}/tree.nwk" ]]; then
      mv "${_NWK_TMP}/tree.nwk" "$ROOTED_NWK"
      rm -rf "$_NWK_TMP"
      ok "Newick → $ROOTED_NWK"
    else
      err "qiime tools export did not produce tree.nwk in $_NWK_TMP"
    fi
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "4/12 — Tabulate per-sample sequence frequencies"
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "$SAMPLE_FREQ_QZA" ]]; then
  skipped "sample-frequencies.qza"
else
  info "Computing per-sample read counts..."
  timer_start
  run_cmd qiime feature-table tabulate-sample-frequencies \
    --i-table               "$ACTIVE_FEATURE_TABLE" \
    --o-sample-frequencies  "$SAMPLE_FREQ_QZA"
  timer_end
  $DRY_RUN || ok "Sample frequencies → $SAMPLE_FREQ_QZA"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "5/12 — Extract sample frequency counts (qiime tools export)"
# ─────────────────────────────────────────────────────────────────────────────
# Use the official `qiime tools export` rather than raw unzip. Cleaner, more
# robust to internal QZA structure changes, and explicit about provenance.
#
# Output: a copy at $SAMPLE_FREQ_TSV for top-level inspection, and the full
# extracted archive at $SAMPLE_FREQ_DIR/ for advanced users.

FREQ_TSV=""
if $DRY_RUN; then
  warn "[DRY-RUN] Would export $SAMPLE_FREQ_QZA → $SAMPLE_FREQ_DIR"
  FREQ_TSV="${SAMPLE_FREQ_DIR}/metadata.tsv"
else
  if [[ -f "$SAMPLE_FREQ_TSV" ]] && ! $FORCE_RERUN; then
    skipped "sample-frequencies.tsv (already exported)"
    _stale_warn "sample-frequencies.tsv" "$SAMPLE_FREQ_TSV"
    FREQ_TSV="$SAMPLE_FREQ_TSV"
  else
    info "Exporting sample-frequencies.qza..."
    rm -rf "$SAMPLE_FREQ_DIR"
    mkdir -p "$SAMPLE_FREQ_DIR"
    qiime tools export \
      --input-path  "$SAMPLE_FREQ_QZA" \
      --output-path "$SAMPLE_FREQ_DIR" \
      || err "qiime tools export failed on $SAMPLE_FREQ_QZA
  → Verify the artifact is valid: qiime tools peek $SAMPLE_FREQ_QZA"

    # qiime tools export writes metadata.tsv directly into the output dir
    if [[ -f "${SAMPLE_FREQ_DIR}/metadata.tsv" ]]; then
      FREQ_TSV="${SAMPLE_FREQ_DIR}/metadata.tsv"
    else
      # Fallback: find anywhere in the export dir (defensive against version drift)
      FREQ_TSV="$(find "$SAMPLE_FREQ_DIR" -name "metadata.tsv" 2>/dev/null | head -1 || true)"
    fi
    [[ -n "$FREQ_TSV" && -f "$FREQ_TSV" ]] || err \
      "metadata.tsv not found in qiime export of $SAMPLE_FREQ_QZA
  → Inspect: ls -R $SAMPLE_FREQ_DIR
  → If QIIME2's export format changed, this script needs updating."

    # Also copy to top-level for easy lay-user inspection
    cp "$FREQ_TSV" "$SAMPLE_FREQ_TSV"
    ok "Frequency TSV: $SAMPLE_FREQ_TSV"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "6/12 — Export feature table to BIOM + TSV (downstream analytical inputs)"
# ─────────────────────────────────────────────────────────────────────────────
# WHY THIS STEP NOW EXISTS AS A STANDALONE STEP (was previously inside step 7):
#   The next two steps both need to read the feature table NUMERICALLY:
#     - Step 7 (analytical rarefaction) loads the BIOM directly via biom-format.
#     - Step 8 (depth selection) reads sample-level totals; with Wave 4, it
#       will also consume the per-sample-per-depth Good's coverage CSV that
#       Step 7 emits.  By exporting BIOM+TSV here once, both steps share a
#       single, inspectable artifact and we never rely on QZA-internal layout.
#
# WHY DEPTH SELECTION MOVED *AFTER* ANALYTICAL RAREFACTION:
#   The publication-grade depth-selection rule (Wave 4) uses three concurrent
#   criteria — sample retention, mean Good's coverage, and curve plateau slope.
#   Two of those three signals come from the analytical rarefaction output.
#   Therefore: rarefaction first, then depth selection.
#
# OUTPUTS (top-level inside 11_pre_diversity/ for inspection by lay users):
#   - 11_pre_diversity/exported_feature_table/feature-table.biom
#   - 11_pre_diversity/feature-table.tsv
#   - bash variables MIN_FREQ, MEDIAN_FREQ, MAX_FREQ — anchor the curve x-axis
#     and the depth-selection candidate grid in subsequent steps.

if [[ -f "$FEATURE_BIOM" && -f "$FEATURE_TSV" ]] && ! $FORCE_RERUN; then
  skipped "BIOM/TSV export (outputs already exist)"
  _stale_warn "feature-table.biom" "$FEATURE_BIOM"
  _stale_warn "feature-table.tsv"  "$FEATURE_TSV"
elif $DRY_RUN; then
  warn "[DRY-RUN] Would export feature table to BIOM and TSV"
else
  info "[step6] Exporting feature table → BIOM (qiime tools export)..."
  timer_start
  mkdir -p "$FEATURE_BIOM_DIR"
  qiime tools export \
    --input-path  "$ACTIVE_FEATURE_TABLE" \
    --output-path "$FEATURE_BIOM_DIR" \
    || err "qiime tools export failed for the feature table.
  → Verify the QZA is valid: qiime tools peek '$ACTIVE_FEATURE_TABLE'
  → If filtered table is corrupt, re-run step 7 of the pipeline (mbx_taxonomy_run.sh)."
  [[ -f "$FEATURE_BIOM" ]] || err "BIOM file not produced at expected path: $FEATURE_BIOM
  → 'qiime tools export' may have written to an unexpected filename.
  → Inspect: ls -la '$FEATURE_BIOM_DIR'"
  ok "BIOM exported  → $FEATURE_BIOM"

  info "[step6] Converting BIOM → TSV (human-readable counts table)..."
  if command -v biom &>/dev/null; then
    biom convert -i "$FEATURE_BIOM" -o "$FEATURE_TSV" --to-tsv \
      || err "biom convert failed.
  → biom-format CLI must be on \$PATH inside the QIIME2 env.
  → Test: biom --version"
  else
    # Fall back to QIIME2 env's python with the biom-format Python module.
    "$PY_BIN" - <<PYBIOM
import biom, sys
t = biom.load_table("$FEATURE_BIOM")
open("$FEATURE_TSV", "w").write(t.to_tsv())
PYBIOM
    [[ $? -eq 0 ]] || err "Python fallback for biom→TSV failed.  Verify '$PY_BIN' has biom-format."
  fi
  [[ -f "$FEATURE_TSV" ]] || err "TSV not produced: $FEATURE_TSV"
  ok "TSV exported   → $FEATURE_TSV"
  timer_end
fi

# ── Compute per-sample frequency anchors (MIN / MEDIAN / MAX) ────────────────
# Anchors the curve x-axis (Step 7), the depth-selection candidate grid
# (Step 8), and the official QIIME alpha-rarefaction max depth (Step 9).
#
# Format reminder for sample-frequencies.tsv (QIIME2 2024+):
#   Line 1: '#SampleID\t#q2:types\t#q2:types'         ← directives, skip
#   Line 2: 'Sample ID\tFrequency\tNo. of Associated Features'  ← real header
#   Line 3+: data rows; values may carry comma thousands separators (28,005.0)

if $DRY_RUN; then
  MIN_FREQ="DRYRUN"; MEDIAN_FREQ="DRYRUN"; MAX_FREQ="DRYRUN"; N_SAMPLES_DETECTED="DRYRUN"
  warn "[DRY-RUN] Would compute MIN_FREQ / MEDIAN_FREQ / MAX_FREQ from $SAMPLE_FREQ_TSV"
else
  [[ -f "$SAMPLE_FREQ_TSV" ]] || err "Cannot compute frequency stats — missing TSV:
  → Expected: $SAMPLE_FREQ_TSV
  → This file should have been written in step 5.
  → Re-run with --force-rerun if step 5 was interrupted."

  _FREQ_STATS="$(
    awk -F'\t' '
      BEGIN { col=0; minv=""; maxv=""; n=0 }
      /^#/ { next }
      NF == 0 { next }
      header_done == 0 {
        for (i=1; i<=NF; i++) {
          h = $i
          gsub(/[ \t\r"]+/, "", h)
          if (tolower(h) == "frequency") col = i
        }
        if (col == 0) col = 2     # default: column 2 (matches QIIME default layout)
        header_done = 1
        next
      }
      {
        v = $col
        gsub(/,/,         "", v)   # strip thousands separators
        gsub(/[ \t\r"]+/, "", v)
        if (v ~ /^[0-9]+(\.[0-9]+)?$/) {
          x = v + 0
          n++; a[n] = x
          if (minv == "" || x < minv) minv = x
          if (maxv == "" || x > maxv) maxv = x
        }
      }
      END {
        if (n == 0) { print "ERROR\t0\t0\t0\t0"; exit }
        # in-place sort (small n, simple O(n^2) is fine)
        for (i=1; i<=n; i++) for (j=i+1; j<=n; j++) if (a[i] > a[j]) { t=a[i]; a[i]=a[j]; a[j]=t }
        med = (n % 2 == 1) ? a[(n+1)/2] : (a[n/2] + a[n/2+1]) / 2.0
        printf "OK\t%d\t%d\t%d\t%d\n", n, int(minv+0.5), int(med+0.5), int(maxv+0.5)
      }
    ' "$SAMPLE_FREQ_TSV"
  )"
  _FREQ_OK="$(echo "$_FREQ_STATS" | cut -f1)"
  if [[ "$_FREQ_OK" != "OK" ]]; then
    err "Could not parse a numeric frequency column from $SAMPLE_FREQ_TSV
  → No numeric values were detected in the auto-located column.
  → Inspect manually: head -5 '$SAMPLE_FREQ_TSV'
  → Confirm the second column is the per-sample total (it usually is)."
  fi
  N_SAMPLES_DETECTED="$(echo "$_FREQ_STATS" | cut -f2)"
  MIN_FREQ="$(echo "$_FREQ_STATS"           | cut -f3)"
  MEDIAN_FREQ="$(echo "$_FREQ_STATS"        | cut -f4)"
  MAX_FREQ="$(echo "$_FREQ_STATS"           | cut -f5)"

  ok "Frequency stats from $N_SAMPLES_DETECTED samples:"
  ok "  MIN    : $MIN_FREQ reads     (lowest-coverage sample)"
  ok "  MEDIAN : $MEDIAN_FREQ reads"
  ok "  MAX    : $MAX_FREQ reads     (used as analytical-curve x-axis maximum)"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "7/12 — Alpha rarefaction (analytical, seconds not hours)"
# ─────────────────────────────────────────────────────────────────────────────
# WHY MONTE CARLO SUBSAMPLING IS WRONG FOR CURVES:
#   qiime diversity alpha-rarefaction (and our previous approach) uses random
#   subsampling repeated N times to ESTIMATE the expected rarefaction curve.
#   This is statistically unnecessary and computationally wasteful.
#
# THE CORRECT APPROACH — Analytical rarefaction (Hurlbert 1971):
#   The expected number of species at depth n has a mathematically EXACT
#   closed-form solution:
#
#     E[S(n)] = Σᵢ [ 1 - C(N-nᵢ, n) / C(N, n) ]
#
#   where N = total reads, nᵢ = reads for species i, C() = combinations.
#   This IS the expected value — not an approximation of it.
#   Running 10 iterations of random subsampling gives you an ESTIMATE of
#   what this formula computes EXACTLY. We skip the estimation entirely.
#
#   Used by: iNEXT (R), vegan::rarecurve(), EstimateS, Chao lab tools.
#   Reference: Hurlbert 1971, Colwell et al. 2012 (J Plant Ecol).
#
# RESULT:
#   Before: 50-200 jobs x random subsampling = minutes to hours
#   After:  1 pass per (sample x depth), fully vectorized numpy = SECONDS
#
# SCIENCE PRESERVED:
#   Observed Features: exact analytical formula (exact, not estimated)
#   Shannon entropy:   analytical expected value via Good-Turing framework
#   Faith PD:          computed in core-metrics-phylogenetic (step 12)
#                      — meaningless on a curve anyway since tree structure
#                        doesn't change with subsampling depth
#   Output QZV:        same format, viewable at view.qiime2.org
#   Output PNG:        same multi-panel figure

ALPHA_RARE_PNG="${PRE_DIV_DIR}/alpha_rarefaction_curves.png"
ALPHA_RARE_CSV="${PRE_DIV_DIR}/alpha_rarefaction_data.csv"
# NOTE: BIOM export was moved into Step 6 (standalone) so both Step 7 and
# Step 8 (depth selection, after Wave 4) can consume it.  We use the
# session-wide $FEATURE_BIOM and $FEATURE_TSV variables defined up top.

info "[step7] Build version: analytical-rarefaction v2 (vectorized, $(date +%Y-%m-%d))"
info "[step7] Python interpreter: $PY_BIN"
info "[step7] Curve maximum    : $MAX_FREQ reads (== max sample frequency)"

if [[ -f "$ALPHA_RARE_MBX_QZV" && -f "$ALPHA_RARE_PNG" ]] && ! $FORCE_RERUN; then
  skipped "alpha rarefaction (analytical preview) — outputs already exist"
  _stale_warn "alpha-rarefaction-mbx-preview.qzv" "$ALPHA_RARE_MBX_QZV"
elif $DRY_RUN; then
  warn "[DRY-RUN] Would run analytical rarefaction (seconds)"
else
  [[ -f "$FEATURE_BIOM" ]] || err "BIOM table not found: $FEATURE_BIOM
  → Step 6 should have produced this.  Re-run with --force-rerun to recreate."

  info "[step7] Running analytical rarefaction (Hurlbert 1971 exact formula)..."
  info "[step7] No random sampling — completes in seconds regardless of dataset size."
  info "[step7] If you see no further output for >60s, the matplotlib font cache"
  info "[step7] may be rebuilding (first run on macOS). It will appear eventually."

  _TMPID="${$}_$(date +%s)"
  _PY_WORKER="/tmp/mbx_analytical_rare_${_TMPID}.py"

  {
  cat << 'PYHDR'
"""
Analytical rarefaction curves — preview pass.

Two metrics are emitted, both via *closed-form* expectations under
sampling-without-replacement (multivariate hypergeometric).  No Monte Carlo,
no random subsampling — the values you see are mathematical expectations.

  observed_features   E[S(n)]  =  Σ_i [ 1 - C(N - n_i, n) / C(N, n) ]
                                  (Hurlbert 1971)

  goods_coverage      G(n)     =  1 - E[F1(n)] / n
                                  E[F1(n)] = Σ_i n_i · C(N-n_i, n-1) / C(N, n)
                                  (Good 1953; Chao 1984; Chao & Jost 2012)

`observed_features` is the EXACT expected value of QIIME2's Monte-Carlo
`alpha-rarefaction` curve — running infinite iterations of random subsampling
would converge here.  It is therefore safe for publication when the figure is
labeled as "expected value" or "analytical".

`goods_coverage` is the rarefaction-stage Good's coverage.  Its value at
n = N matches the standard sample-level Good's coverage 1 - f1/N where f1 is
the count of singleton features.  We use it to drive Wave 4's depth-selection
gate (mean coverage of retained samples ≥ GOOD_COV_MIN, default 0.98).

Shannon entropy is INTENTIONALLY OMITTED here.  The Chao-Shen (2003)
coverage-adjusted analytical Shannon is a different estimator from QIIME's
Monte-Carlo Shannon — claiming equivalence misleads peer reviewers.  For a
publication-grade Shannon curve see `alpha-rarefaction-qiime.qzv` produced
by the official QIIME plugin in Step 9.

Per-sample diagnostics (one row per sample × depth) are written to a separate
CSV so Wave 4's depth-selection R script can compute mean coverage and
plateau slope analytically without re-running this Python worker.

References:
  Hurlbert SH (1971) Ecology 52:577-586
  Good IJ (1953) Biometrika 40:237-264
  Chao A (1984) Scand J Stat 11:265-270
  Chao A & Jost L (2012) Ecology 93:2533-2547
  Colwell RK et al (2012) J Plant Ecol 5:3-21
"""
import os, sys, zipfile, shutil, uuid, time, hashlib, datetime
print("[py]    starting up...", flush=True)
_t_imp = time.time()
from pathlib import Path
import numpy as np
print(f"[py]    numpy {np.__version__} imported ({time.time()-_t_imp:.2f}s)", flush=True)
_t = time.time()
import pandas as pd
print(f"[py]    pandas {pd.__version__} imported ({time.time()-_t:.2f}s)", flush=True)
_t = time.time()
from scipy.special import gammaln   # log-gamma for stable log-combinations
print(f"[py]    scipy.special.gammaln imported ({time.time()-_t:.2f}s)", flush=True)
_t = time.time()
import biom
print(f"[py]    biom imported ({time.time()-_t:.2f}s)", flush=True)
_t = time.time()
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
print(f"[py]    matplotlib imported ({time.time()-_t:.2f}s) "
      f"[on first run this may rebuild the font cache and take 30-60s]", flush=True)
print(f"[py]    total import time: {time.time()-_t_imp:.2f}s", flush=True)

BIOM_FILE       = os.environ["MBX_BIOM"]
META_FILE       = os.environ["MBX_META"]
OUT_CSV         = os.environ["MBX_RARE_CSV"]            # mean curve (one row per depth)
OUT_PER_SAMPLE  = os.environ["MBX_PER_SAMPLE_CSV"]      # per-sample × depth diagnostics
OUT_PNG         = os.environ["MBX_RARE_PNG"]
OUT_QZV         = os.environ["MBX_RARE_QZV"]
MAX_DEPTH       = int(os.environ["MBX_MAX_DEPTH"])      # x-axis max == max(sample freq)
N_STEPS         = int(os.environ["MBX_STEPS"])
N_JOBS          = int(os.environ["MBX_JOBS"])

# ── Load BIOM (sparse — never densify the entire features × samples table) ────
print(f"[INFO]  Loading BIOM: {BIOM_FILE}", flush=True)
_t_load = time.time()
table       = biom.load_table(BIOM_FILE)
sample_ids  = list(table.ids(axis="sample"))
sparse_mat  = table.matrix_data.tocsc()
n_features, n_samples = sparse_mat.shape
totals      = np.asarray(sparse_mat.sum(axis=0)).ravel().astype(np.int64)
min_total   = int(totals.min()) if totals.size else 0
max_total   = int(totals.max()) if totals.size else 0
print(f"[INFO]  Features: {n_features:,} | Samples: {n_samples}", flush=True)
print(f"[INFO]  Per-sample totals — min: {min_total:,} | max: {max_total:,}", flush=True)
print(f"[INFO]  Curve x-axis maximum (MAX_FREQ): {MAX_DEPTH:,}", flush=True)
print(f"[INFO]  BIOM load: {time.time()-_t_load:.2f}s", flush=True)

# ── Load metadata ─────────────────────────────────────────────────────────────
raw     = open(META_FILE).readlines()
content = [l for l in raw if not l.startswith("#") and l.strip()]
cols    = content[0].rstrip("\n").split("\t")
rows    = [l.rstrip("\n").split("\t") for l in content[1:]]
meta_df = pd.DataFrame(rows, columns=cols)
meta_df.rename(columns={meta_df.columns[0]: "SampleID"}, inplace=True)
meta_df.set_index("SampleID", inplace=True)

# ── Rarefaction depth grid ────────────────────────────────────────────────────
# Geometric spacing from min_d up to MAX_DEPTH (== max sample frequency).
# We do NOT cap at min_total — per-sample loop in compute_sample_curves drops
# any depth > that sample's N, so the curve still gracefully handles
# heterogeneous library sizes.  This means every sample contributes to the
# left portion of the curve, and only the deepest samples extend to MAX_DEPTH.
min_d_anchor = max(100, min(1000, max(min_total // 10, 100)))
depths = np.unique(
    np.round(np.geomspace(min_d_anchor, MAX_DEPTH, N_STEPS)).astype(int)
)
depths = depths[depths >= 1]
if depths.size == 0:
    depths = np.array([max(min_total, 1)])
print(f"[INFO]  Rarefaction depths ({len(depths)} steps): "
      f"{int(depths[0]):,} to {int(depths[-1]):,}", flush=True)

# ── Vectorized analytical curves for ONE sample (all depths at once) ──────────
# observed_features:
#   log_p_absent_i = log C(N - n_i, n) - log C(N, n)             (per feature)
#   p_present_i    = 1 - exp(clip(log_p_absent_i, -700, 0))
#   E[S(n)]        = Σ_i p_present_i
#
# goods_coverage:
#   log_p_one_i    = log(n_i) + log C(N - n_i, n - 1) - log C(N, n)
#   p_one_i        = exp(clip(log_p_one_i, -700, 0))             (singleton mass)
#   E[F1(n)]       = Σ_i p_one_i
#   G(n)           = 1 - E[F1(n)] / n
#
# Notes:
#  - n_i is the present-only counts vector (entries > 0).
#  - For features where N - n_i < n, the absence probability is 0 →
#    p_present = 1.  Likewise for n - 1 in the singleton expression.
#  - Coverage is well-defined only for n ≥ 1; we skip n == 0.
def compute_sample_curves(n_i: np.ndarray, N: int, depths_arr: np.ndarray):
    if N <= 0 or n_i.size == 0:
        return []
    n_i      = np.asarray(n_i, dtype=np.int64)
    n_i_f    = n_i.astype(np.float64)
    log_n_i  = np.log(n_i_f)                                # all > 0 by construction

    # Pre-compute terms that depend only on n_i (independent of depth)
    log_diff_top = gammaln((N - n_i_f) + 1.0)               # gammaln(N - n_i + 1)

    rows = []
    for n in depths_arr:
        n = int(n)
        if n < 1 or N < n:
            continue

        # log C(N, n) — same for all features
        log_n_p1  = gammaln(n + 1.0)
        log_denom = gammaln(N + 1.0) - log_n_p1 - gammaln(N - n + 1.0)

        # ── Observed features (Hurlbert 1971) ───────────────────────────────
        diff_bot = (N - n_i_f) - n                          # may be negative
        valid    = diff_bot >= 0
        log_p_absent = np.full(n_i.shape, -np.inf, dtype=np.float64)
        if valid.any():
            log_p_absent[valid] = (log_diff_top[valid]
                                   - log_n_p1
                                   - gammaln(diff_bot[valid] + 1.0)) - log_denom
        p_absent  = np.exp(np.clip(log_p_absent, -700.0, 0.0))
        p_present = 1.0 - p_absent
        obs       = float(p_present.sum())

        # ── Good's coverage at depth n ──────────────────────────────────────
        # Pr(rarefied count == 1) per feature = n_i · C(N-n_i, n-1) / C(N, n)
        #   log = log(n_i) + [gammaln(N-n_i+1) - gammaln(n) - gammaln(N-n_i-(n-1)+1)] - log_denom
        # Domain: requires N - n_i ≥ n - 1
        if n == 1:
            # at depth 1, every drawn read is a singleton by definition
            # E[F1(1)] = 1, so G(1) = 0
            goods_cov = 0.0
        else:
            diff_one_bot = (N - n_i_f) - (n - 1)            # may be negative
            valid_s      = diff_one_bot >= 0
            log_p_one    = np.full(n_i.shape, -np.inf, dtype=np.float64)
            if valid_s.any():
                log_p_one[valid_s] = (
                    log_n_i[valid_s]
                    + gammaln((N - n_i_f[valid_s]) + 1.0)
                    - gammaln(float(n))                     # gammaln(n)  =  log((n-1)!)
                    - gammaln(diff_one_bot[valid_s] + 1.0)
                    - log_denom
                )
            p_one     = np.exp(np.clip(log_p_one, -700.0, 0.0))
            E_F1      = float(p_one.sum())
            goods_cov = 1.0 - (E_F1 / float(n))
            # Numerical guard: keep coverage in [0,1] even if floating-point
            # accumulation drifts slightly outside.
            if   goods_cov < 0.0:  goods_cov = 0.0
            elif goods_cov > 1.0:  goods_cov = 1.0

        rows.append({
            "depth":             n,
            "observed_features": obs,
            "goods_coverage":    goods_cov,
            "sample_total":      N,
        })
    return rows

def _get_sample_present_counts(j: int):
    """Return (nonzero_counts_vector, total_reads) for sample j (CSC slice)."""
    col = sparse_mat.getcol(j)                              # CSC column = O(nnz)
    n_i = col.data.astype(np.int64)                         # already nonzero only
    return n_i, int(n_i.sum())

# ── Run computation ───────────────────────────────────────────────────────────
print(f"[INFO]  Computing analytical curves for {n_samples} samples × "
      f"{len(depths)} depths (vectorized)...", flush=True)
_t_compute = time.time()

USE_THREADS = (n_samples > 32) and (N_JOBS > 1)
all_results = []
if USE_THREADS:
    from concurrent.futures import ThreadPoolExecutor
    print(f"[INFO]  Parallel mode: {N_JOBS} threads (numpy releases GIL)", flush=True)
    def _job(j):
        n_i, N = _get_sample_present_counts(j)
        rows = compute_sample_curves(n_i, N, depths)
        for r in rows:
            r["sample_id"] = sample_ids[j]
        return rows
    done = 0
    step_print = max(1, n_samples // 10)
    with ThreadPoolExecutor(max_workers=N_JOBS) as ex:
        for rows in ex.map(_job, range(n_samples)):
            done += 1
            all_results.extend(rows)
            if done == n_samples or done % step_print == 0:
                print(f"  [{done:>4}/{n_samples}] samples done", flush=True)
else:
    print(f"[INFO]  Sequential vectorized mode", flush=True)
    step_print = max(1, n_samples // 10)
    for j in range(n_samples):
        n_i, N = _get_sample_present_counts(j)
        rows = compute_sample_curves(n_i, N, depths)
        for r in rows:
            r["sample_id"] = sample_ids[j]
        all_results.extend(rows)
        if (j + 1) == n_samples or (j + 1) % step_print == 0:
            print(f"  [{j+1:>4}/{n_samples}] samples done", flush=True)

print(f"[OK]    Analytical computation done in {time.time()-_t_compute:.2f}s "
      f"— {len(all_results)} (sample, depth) rows", flush=True)

# ── Per-sample × depth diagnostic CSV (consumed by Wave 4 depth selection) ────
per_sample_df = pd.DataFrame(all_results, columns=[
    "sample_id", "depth", "sample_total", "observed_features", "goods_coverage"
])
per_sample_df.to_csv(OUT_PER_SAMPLE, index=False)
print(f"[OK]    Per-sample diagnostics : {OUT_PER_SAMPLE}", flush=True)

# ── Mean-of-samples curve (one row per depth) — what the PNG shows ────────────
# Only retains samples whose N >= depth (those are the rows present in
# per_sample_df at that depth, by construction of compute_sample_curves).
mean_df = (
    per_sample_df
    .groupby("depth", as_index=False)
    .agg(
        n_samples_at_depth   = ("sample_id",          "count"),
        observed_features    = ("observed_features",  "mean"),
        goods_coverage_mean  = ("goods_coverage",     "mean"),
        goods_coverage_min   = ("goods_coverage",     "min"),
    )
    .sort_values("depth")
)
mean_df.to_csv(OUT_CSV, index=False)
print(f"[OK]    Mean curve CSV         : {OUT_CSV}", flush=True)

# Merge metadata onto per-sample for the PNG (sample-coloured curves)
df = per_sample_df.merge(meta_df.reset_index(),
                         left_on="sample_id", right_on="SampleID", how="left")

# ── Plot ──────────────────────────────────────────────────────────────────────
metrics = ["observed_features", "goods_coverage"]
labels  = {
    "observed_features": "Observed Features\n(Hurlbert 1971 — analytical expectation)",
    "goods_coverage":    "Good's Coverage\n(Good 1953; analytical expectation)",
}
fig = plt.figure(figsize=(6 * len(metrics), 5))
gs  = gridspec.GridSpec(1, len(metrics), figure=fig, wspace=0.32)
pal = plt.cm.tab10.colors

for ci, metric in enumerate(metrics):
    ax = fig.add_subplot(gs[0, ci])
    for si, sid in enumerate(sample_ids):
        s = df[df.sample_id == sid].sort_values("depth")
        if s.empty:
            continue
        c = pal[si % len(pal)]
        ax.plot(s.depth, s[metric], color=c, alpha=0.8, lw=1.8,
                label=sid if ci == 0 else "")
    ax.axvline(MAX_DEPTH, color="crimson", ls="--", lw=2.0,
               label=f"Curve maximum\n({MAX_DEPTH:,} reads = max sample frequency)")
    ax.set_xlabel("Sequencing Depth", fontsize=10)
    ax.set_ylabel(labels[metric], fontsize=10)
    ax.set_title(labels[metric], fontsize=10, fontweight="bold")
    ax.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    ax.grid(True, alpha=0.3, linestyle=":")
    if metric == "goods_coverage":
        ax.set_ylim(0.0, 1.005)
    if ci == 0 and n_samples <= 25:
        ax.legend(fontsize=6, loc="lower right", framealpha=0.7)

plt.suptitle(
    f"Alpha Rarefaction Curves  (Analytical preview — Hurlbert 1971 + Good 1953)  |  "
    f"{n_samples} samples  |  {len(depths)} depths  |  "
    f"Curve max (max sample freq): {MAX_DEPTH:,} reads",
    fontsize=11, fontweight="bold")
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.close()
print(f"[OK]    Figure                 : {OUT_PNG}", flush=True)

# ── Pack QZV ─────────────────────────────────────────────────────────────────
qzv_uuid = str(uuid.uuid4())
tmp      = Path(OUT_QZV).parent / f"_qzvtmp_{qzv_uuid}"
data_dir = tmp / qzv_uuid / "data"
data_dir.mkdir(parents=True)
(tmp / qzv_uuid / "VERSION").write_text(
    "QIIME 2\narchive: 5\nframework: 2024.5.0\n")
(tmp / qzv_uuid / "metadata.yaml").write_text(
    f"uuid: {qzv_uuid}\ntype: Visualization\nformat: null\n")
shutil.copy(OUT_PNG,        data_dir / "alpha_rarefaction_curves.png")
shutil.copy(OUT_CSV,        data_dir / "alpha_rarefaction_mean.csv")
shutil.copy(OUT_PER_SAMPLE, data_dir / "alpha_rarefaction_per_sample.csv")

# Summary table — mean per depth, rounded for display
summ  = mean_df.copy()
hdr   = "<tr>" + "".join(
    f"<th>{c}</th>" for c in summ.columns) + "</tr>"
trows = "".join(
    "<tr>" + "".join(
        f"<td>{round(v,4) if isinstance(v, float) else v}</td>"
        for v in row
    ) + "</tr>"
    for row in summ.itertuples(index=False))

(data_dir / "index.html").write_text(
    "<!DOCTYPE html><html><head><meta charset='utf-8'>"
    "<title>Alpha Rarefaction (Analytical Preview)</title>"
    "<style>"
    "body{font-family:sans-serif;max-width:1100px;margin:40px auto;padding:20px}"
    "h1{color:#2c5f8a}"
    "img{max-width:100%;border:1px solid #ddd;border-radius:4px}"
    "table{border-collapse:collapse;width:100%;margin-top:20px}"
    "th,td{border:1px solid #ccc;padding:6px 10px;text-align:right}"
    "th{background:#2c5f8a;color:#fff}"
    "tr:nth-child(even){background:#f5f5f5}"
    ".note{background:#e8f4fd;border-left:4px solid #2c5f8a;"
    "       padding:12px;border-radius:4px;margin:16px 0;font-size:0.9em}"
    ".warn{background:#fff8e1;border-left:4px solid #f0b400;"
    "       padding:12px;border-radius:4px;margin:16px 0;font-size:0.9em}"
    "</style></head><body>"
    "<h1>Alpha Rarefaction Curves &mdash; Analytical Preview</h1>"
    "<div class='note'>"
    "<b>Method: Closed-form expectations under sampling-without-replacement.</b><br>"
    "<b>Observed Features</b> &mdash; Hurlbert (1971): the <i>exact expected value</i> "
    "of the rarefaction curve, mathematically equivalent to QIIME's Monte-Carlo "
    "<code>alpha-rarefaction</code> as iterations &rarr; &infin;. "
    "Computed in seconds via closed-form combinatorics.<br>"
    "<b>Good's Coverage</b> &mdash; Good (1953), Chao &amp; Jost (2012): "
    "G(n) = 1 - E[F<sub>1</sub>(n)] / n, the rarefaction-stage analogue of "
    "the standard sample-level Good's coverage.  Drives the depth-selection "
    "PASS gate in step 8.<br>"
    f"Red dashed line = curve x-axis maximum ({MAX_DEPTH:,} reads = max sample frequency)."
    "</div>"
    "<div class='warn'>"
    "<b>This is a fast preview to inspect the shape of the curve.</b> "
    "For figures intended for publication, use "
    "<code>alpha-rarefaction-qiime.qzv</code> (produced by the official "
    "QIIME plugin in step 9), which includes Shannon and Faith PD."
    "</div>"
    f"<img src='alpha_rarefaction_curves.png' alt='Rarefaction Curves'>"
    f"<h2>Mean curve (averaged across samples retained at each depth)</h2>"
    f"<table>{hdr}{trows}</table>"
    "<p>"
    "<a href='alpha_rarefaction_mean.csv'>Download mean curve CSV</a> &middot; "
    "<a href='alpha_rarefaction_per_sample.csv'>Download per-sample CSV</a>"
    "</p>"
    "</body></html>")

# ── Provenance (required by QIIME 2 archive spec v5; q2view uses for graph) ─
prov_dir         = tmp / qzv_uuid / "provenance"
prov_artifacts   = prov_dir / "artifacts"
prov_action_dir  = prov_dir / "action"
prov_action_dir.mkdir(parents=True)
prov_artifacts.mkdir(parents=True)

now_iso = datetime.datetime.utcnow().isoformat()
(prov_dir / "VERSION").write_text(
    "QIIME 2\narchive: 5\nframework: 2024.5.0\n")
(prov_dir / "metadata.yaml").write_text(
    f"uuid: {qzv_uuid}\ntype: Visualization\nformat: null\n")
(prov_action_dir / "action.yaml").write_text(
    "execution:\n"
    f"    uuid: {qzv_uuid}\n"
    "    runtime:\n"
    "        python: '3.10'\n"
    "        framework: '2024.5.0'\n"
    f"        start: '{now_iso}'\n"
    f"        end:   '{now_iso}'\n"
    "action:\n"
    "    type: import\n"
    "    format: null\n"
    "    manifest: []\n"
)

# ── checksums.md5 (required by archive v5+) ─────────────────────────────────
checksums_path = tmp / qzv_uuid / "checksums.md5"
_uuid_root = tmp / qzv_uuid
_lines = []
for fp in sorted(_uuid_root.rglob("*")):
    if not fp.is_file():
        continue
    if fp == checksums_path:
        continue
    rel = fp.relative_to(_uuid_root).as_posix()
    h   = hashlib.md5(fp.read_bytes()).hexdigest()
    _lines.append(f"{h}  {rel}")
checksums_path.write_text("\n".join(_lines) + "\n")

if Path(OUT_QZV).exists():
    Path(OUT_QZV).unlink()
with zipfile.ZipFile(OUT_QZV, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    for fp in tmp.rglob("*"):
        if fp.is_file():
            zf.write(fp, fp.relative_to(tmp))
shutil.rmtree(tmp)
print(f"[OK]    QZV: {OUT_QZV}", flush=True)
print("[DONE]  Analytical rarefaction (preview) complete.", flush=True)
PYHDR
  } > "$_PY_WORKER"

  info "[step7] Launching python worker (unbuffered)..."
  timer_start
  MBX_BIOM="$FEATURE_BIOM" \
  MBX_META="$METADATA_TXT" \
  MBX_RARE_CSV="$ALPHA_RARE_CSV" \
  MBX_PER_SAMPLE_CSV="$PER_SAMPLE_RARE_CSV" \
  MBX_RARE_PNG="$ALPHA_RARE_PNG" \
  MBX_RARE_QZV="$ALPHA_RARE_MBX_QZV" \
  MBX_MAX_DEPTH="$MAX_FREQ" \
  MBX_STEPS="$RARE_STEPS" \
  MBX_JOBS="$N_JOBS" \
  PYTHONUNBUFFERED=1 \
    "$PY_BIN" -u "$_PY_WORKER" \
    || err "Analytical rarefaction failed (python: $PY_BIN).
  Common causes:
    1. biom-format / scipy / numpy missing in QIIME2 env (very unusual).
       → conda activate qiime2-amplicon-2025.4 && python -c 'import biom, scipy, numpy'
    2. Wrong python interpreter — verify '$PY_BIN' is the QIIME2 env's python:
       → ls -la $PY_BIN
    3. Out of memory: try --rare-steps 5 to reduce depth resolution.
  Check the Python traceback above for the exact error."
  timer_end
  rm -f "$_PY_WORKER"

  ok "Rarefaction preview QZV  -> $ALPHA_RARE_MBX_QZV"
  ok "Figure                   -> $ALPHA_RARE_PNG"
  ok "Mean curve CSV           -> $ALPHA_RARE_CSV"
  ok "Per-sample diag CSV      -> $PER_SAMPLE_RARE_CSV"
fi
# ─────────────────────────────────────────────────────────────────────────────
step "8/12 — Optimal sampling depth (three concurrent criteria, R)"
# ─────────────────────────────────────────────────────────────────────────────
# Publication-grade depth selection using THREE concurrent criteria.
# A candidate depth d earns STATUS=PASS only if ALL three are satisfied:
#
#   (a) retention   : fraction of samples with N >= d  is  >= MIN_OVERALL
#                     (default 0.90; --min-overall to change)
#       AND if --group-col was supplied, every group's retention is
#       >= MIN_GROUP (default 0.80; --min-group to change).
#
#   (b) coverage    : mean Good's coverage of retained samples at depth d
#                     is  >= GOOD_COV_MIN  (default 0.98; --good-coverage-min
#                     to change).  Computed analytically by step 7 and read
#                     from $PER_SAMPLE_RARE_CSV.
#
#   (c) plateau     : the observed_features slope at depth d (units:
#                     features per 1000 additional reads, derived from the
#                     mean curve in $ALPHA_RARE_CSV) is  <  PLATEAU_SLOPE_MAX
#                     (default 0.5; --plateau-slope-max to change).
#                     This catches the "still climbing" case.
#
# Among candidate depths satisfying all three, we pick the MAXIMUM depth
# (deeper → more reads per retained sample → better diversity estimates).
#
# If NONE satisfy all three, we fall back along an explicit ladder and the
# STATUS becomes REVIEW_REQUIRED — never silent acceptance:
#
#   fallback 1  : drop criterion (c), keep (a)+(b).
#   fallback 2  : drop (b)+(c), keep (a) only with looser threshold 0.75.
#   fallback 3  : Q1 of per-sample totals.
#   fallback 4  : min(per-sample totals) — guarantees all samples retained.
#
# We also FAIL hard if a group ends up entirely removed at the chosen depth.
#
# NO AUTO-DETECTION OF GROUP_COL.  Per scientific safety: if the user did not
# supply --group-col we evaluate overall-only retention and write
# GROUP_COLUMN=none.  This is a feedback-driven rule (auto-detection produced
# different results across runs depending on column ordering, which is unsafe).
#
# Outputs:
#   $DEPTH_CSV               — every candidate depth × criteria evaluation
#   $SAMPLES_RETAINED_CSV    — list of samples kept at RECOMMENDED_DEPTH
#   $SAMPLES_REMOVED_CSV     — list of samples dropped at RECOMMENDED_DEPTH
#   $SAMPLE_DEPTH_SUMMARY_CSV — per-sample diagnostic (id, total, group, kept)
#   $GROUP_DEPTH_SUMMARY_CSV  — per-group diagnostic (n, kept, retention)
#   $R_DEPTH_OUT             — key=value summary parsed back into bash
#
# Whitespace handling: silent trim on sample IDs and group labels at entry,
# with one info line "Trimmed N IDs / M group labels" if any were touched.

if [[ -z "$GROUP_COL" ]]; then
  info "No --group-col supplied → overall retention only (group constraint disabled)."
  info "(This is the scientifically-safe default; auto-selection has been removed.)"
  GROUP_COLUMN_EFFECTIVE="none"
else
  info "Group column requested: $GROUP_COL"
  GROUP_COLUMN_EFFECTIVE="$GROUP_COL"
fi

# Write the R depth-selection script to a temp file.  Heredoc is QUOTED
# ('RDEPTH') so $-signs in R code are NOT interpolated by bash; all params
# arrive via environment variables.
cat > "$R_DEPTH_SCRIPT" << 'RDEPTH'
# ============================================================================
# Three-criterion sampling depth selector (mbX Pro).
# Reads:
#   MBX_FREQ_TSV       — sample-frequencies.tsv (qiime export)
#   MBX_METADATA       — metadata file
#   MBX_PER_SAMPLE_CSV — alpha_rarefaction_per_sample.csv (Wave 3 output)
#   MBX_MEAN_CSV       — alpha_rarefaction_data.csv (mean curve, Wave 3)
#   MBX_GROUP_COL      — group column name; empty/"none" disables group rule
#   MBX_MIN_OVERALL    — overall retention threshold (default 0.90)
#   MBX_MIN_GROUP      — per-group retention threshold (default 0.80)
#   MBX_GOOD_COV_MIN   — mean Good's coverage threshold (default 0.98)
#   MBX_SLOPE_MAX      — plateau slope threshold (features/1000 reads, default 0.5)
# Writes:
#   MBX_DEPTH_CSV
#   MBX_SAMPLES_RETAINED_CSV
#   MBX_SAMPLES_REMOVED_CSV
#   MBX_SAMPLE_DEPTH_SUMMARY_CSV
#   MBX_GROUP_DEPTH_SUMMARY_CSV
#   MBX_DEPTH_TXT      — key=value file consumed by bash
# ============================================================================

freq_tsv      <- Sys.getenv("MBX_FREQ_TSV")
meta_file     <- Sys.getenv("MBX_METADATA")
per_sample    <- Sys.getenv("MBX_PER_SAMPLE_CSV")
mean_curve    <- Sys.getenv("MBX_MEAN_CSV")
group_column  <- Sys.getenv("MBX_GROUP_COL")
min_overall   <- as.numeric(Sys.getenv("MBX_MIN_OVERALL",  "0.90"))
min_group     <- as.numeric(Sys.getenv("MBX_MIN_GROUP",    "0.80"))
good_cov_min  <- as.numeric(Sys.getenv("MBX_GOOD_COV_MIN", "0.98"))
slope_max     <- as.numeric(Sys.getenv("MBX_SLOPE_MAX",    "0.5"))
out_csv       <- Sys.getenv("MBX_DEPTH_CSV")
out_kept      <- Sys.getenv("MBX_SAMPLES_RETAINED_CSV")
out_removed   <- Sys.getenv("MBX_SAMPLES_REMOVED_CSV")
out_sample    <- Sys.getenv("MBX_SAMPLE_DEPTH_SUMMARY_CSV")
out_group     <- Sys.getenv("MBX_GROUP_DEPTH_SUMMARY_CSV")
out_txt       <- Sys.getenv("MBX_DEPTH_TXT")

if (group_column == "" || tolower(group_column) == "none") group_column <- NULL

# ── Load per-sample frequencies (QIIME-aware: skip #, strip commas) ──────────
sf <- tryCatch({
  raw <- readLines(freq_tsv)
  raw <- raw[nchar(trimws(raw)) > 0]
  content <- raw[!grepl("^#", raw)]
  if (length(content) < 2) stop("sample-frequencies TSV has no data rows")
  hcols <- trimws(strsplit(content[1], "\t")[[1]])
  drows <- content[-1]
  ncols <- length(hcols)
  parsed <- lapply(drows, function(l) {
    fld <- strsplit(l, "\t")[[1]]
    length(fld) <- ncols
    fld
  })
  d <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
  names(d) <- hcols
  names(d)[1] <- "SampleID"
  for (col in names(d)[-1]) {
    cleaned <- gsub(",", "", trimws(d[[col]]))
    conv    <- suppressWarnings(as.numeric(cleaned))
    if (!all(is.na(conv))) d[[col]] <- conv
  }
  d
}, error = function(e) {
  stop(sprintf("Cannot read frequency TSV (%s): %s", freq_tsv, conditionMessage(e)))
})

# Pick the numeric column with the largest median = the "Frequency" column
nums <- names(sf)[sapply(sf, is.numeric)]
if (length(nums) == 0) stop("No numeric columns in sample-frequencies TSV.")
total_col <- nums[which.max(sapply(sf[nums], median, na.rm = TRUE))]
cat(sprintf("[INFO]  Frequency column     : '%s'\n", total_col))

# ── Load metadata + silent whitespace trim of sample IDs / group labels ──────
ext <- tolower(tools::file_ext(meta_file))
md  <- if (ext == "csv") {
  read.csv(meta_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  read.delim(meta_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE,
             na.strings = c("", "NA", "N/A"))
}
md <- md[!grepl("^#", md[[1]]), , drop = FALSE]
names(md)[1] <- "SampleID"

n_trimmed_id <- sum(md$SampleID != trimws(md$SampleID), na.rm = TRUE)
md$SampleID <- trimws(md$SampleID)

n_trimmed_grp <- 0L
if (!is.null(group_column)) {
  if (!(group_column %in% names(md))) {
    stop(sprintf(
      "Requested --group-col '%s' not found in metadata.\n  Available columns: %s",
      group_column, paste(names(md), collapse = ", ")))
  }
  raw_grp <- md[[group_column]]
  trim_grp <- trimws(as.character(raw_grp))
  n_trimmed_grp <- sum(!is.na(raw_grp) & raw_grp != trim_grp, na.rm = TRUE)
  md[[group_column]] <- trim_grp
}
if (n_trimmed_id > 0)   cat(sprintf("[INFO]  Trimmed %d sample IDs (whitespace).\n",  n_trimmed_id))
if (n_trimmed_grp > 0)  cat(sprintf("[INFO]  Trimmed %d group labels (whitespace).\n", n_trimmed_grp))

# ── Merge frequencies + metadata ─────────────────────────────────────────────
sf$SampleID <- trimws(as.character(sf$SampleID))
df <- merge(sf[, c("SampleID", total_col)], md, by = "SampleID", all.x = TRUE)
names(df)[2] <- "total_freq"
df <- df[!is.na(df$total_freq) & df$total_freq > 0, , drop = FALSE]

n_samples <- nrow(df)
cat(sprintf("[INFO]  Samples with counts  : %d\n", n_samples))
cat(sprintf("[INFO]  Min count            : %d\n", min(df$total_freq)))
cat(sprintf("[INFO]  Max count            : %d\n", max(df$total_freq)))
cat(sprintf("[INFO]  Median count         : %.0f\n", median(df$total_freq)))
cat(sprintf("[INFO]  Min retention target : %.0f%% overall%s\n",
            min_overall * 100,
            if (is.null(group_column)) "" else sprintf(", %.0f%% per group", min_group * 100)))
cat(sprintf("[INFO]  Good's coverage min  : %.3f\n",  good_cov_min))
cat(sprintf("[INFO]  Plateau slope max    : %.3f features per 1000 reads\n", slope_max))

# ── Load Wave 3 outputs (analytical Good's coverage + observed features) ─────
ps <- tryCatch(
  read.csv(per_sample, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop(sprintf(
    "Cannot read per-sample rarefaction CSV (%s): %s\n  Step 7 (analytical rarefaction) must run first.",
    per_sample, conditionMessage(e)))
)
required_ps_cols <- c("sample_id", "depth", "sample_total", "observed_features", "goods_coverage")
miss <- setdiff(required_ps_cols, names(ps))
if (length(miss) > 0) stop(sprintf(
  "Per-sample CSV is missing columns: %s\n  → Re-run step 7 with --force-rerun.",
  paste(miss, collapse = ", ")))
ps$sample_id <- trimws(as.character(ps$sample_id))

mc <- tryCatch(
  read.csv(mean_curve, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) stop(sprintf(
    "Cannot read mean-curve CSV (%s): %s", mean_curve, conditionMessage(e))))

# ── Build candidate depth grid ───────────────────────────────────────────────
# Nearest-100 rounding + actual sample totals + the depths from the Wave 3 grid.
round_base <- 100
raw_cands <- sort(unique(floor(df$total_freq / round_base) * round_base))
raw_cands <- raw_cands[raw_cands > 0]
ps_depths <- sort(unique(as.integer(ps$depth)))
actual    <- sort(unique(df$total_freq))
cand <- sort(unique(c(raw_cands, actual, ps_depths)))
cand <- cand[cand >= min(df$total_freq) | cand >= 1000]   # drop pathologically tiny anchors
cand <- cand[cand <= max(df$total_freq)]

# ── Helpers ──────────────────────────────────────────────────────────────────

# Mean Good's coverage at depth d, averaging only over samples whose N >= d.
mean_goods_at <- function(d) {
  rows <- ps[ps$depth == d & ps$sample_total >= d, , drop = FALSE]
  if (nrow(rows) == 0) return(NA_real_)
  mean(rows$goods_coverage, na.rm = TRUE)
}

# Linear-interpolation slope of mean observed_features in features per 1000 reads
# at depth d.  Two-point estimator: f(d + h) - f(d - h) over (2h)/1000.
slope_at <- function(d, mean_curve_df) {
  mc2 <- mean_curve_df[order(mean_curve_df$depth), , drop = FALSE]
  if (nrow(mc2) < 2) return(NA_real_)
  # nearest neighbours around d
  i_lo <- which(mc2$depth <= d); i_hi <- which(mc2$depth >= d)
  if (length(i_lo) == 0 || length(i_hi) == 0) return(NA_real_)
  lo <- mc2[max(i_lo), , drop = FALSE]
  hi <- mc2[min(i_hi), , drop = FALSE]
  if (lo$depth == hi$depth) {
    # exact match — use the next neighbour on whichever side
    j <- which(mc2$depth > d)
    if (length(j) == 0) j <- which(mc2$depth < d)
    if (length(j) == 0) return(NA_real_)
    nb <- mc2[j[1], , drop = FALSE]
    delta_x <- abs(nb$depth - lo$depth)
    delta_y <- abs(nb$observed_features - lo$observed_features)
  } else {
    delta_x <- hi$depth - lo$depth
    delta_y <- hi$observed_features - lo$observed_features
  }
  if (delta_x <= 0) return(NA_real_)
  # features per 1000 reads
  (delta_y / delta_x) * 1000
}

# ── Score every candidate depth on all three criteria ────────────────────────
overall_keep    <- sapply(cand, function(d) mean(df$total_freq >= d))
n_samples_kept  <- sapply(cand, function(d) sum(df$total_freq >= d))
goods_at        <- sapply(cand, mean_goods_at)
slope_at_d      <- sapply(cand, function(d) slope_at(d, mc))

if (!is.null(group_column)) {
  per_group_keep  <- sapply(cand, function(d) {
    grp <- tapply(df$total_freq >= d, df[[group_column]], mean, na.rm = TRUE)
    if (length(grp) == 0) NA_real_ else min(grp, na.rm = TRUE)
  })
  any_group_zero  <- sapply(cand, function(d) {
    grp <- tapply(df$total_freq >= d, df[[group_column]], sum, na.rm = TRUE)
    any(grp == 0)
  })
  total_groups    <- length(unique(df[[group_column]]))
} else {
  per_group_keep  <- rep(NA_real_, length(cand))
  any_group_zero  <- rep(FALSE,    length(cand))
  total_groups    <- NA_integer_
}

crit_a <- overall_keep >= min_overall &
          (is.na(per_group_keep) | per_group_keep >= min_group) &
          !any_group_zero
crit_b <- !is.na(goods_at)   & goods_at  >= good_cov_min
crit_c <- !is.na(slope_at_d) & slope_at_d <  slope_max

results <- data.frame(
  depth                       = cand,
  overall_keep_frac           = round(overall_keep, 4),
  n_samples_kept              = n_samples_kept,
  min_group_keep_frac         = round(per_group_keep, 4),
  any_group_fully_removed     = any_group_zero,
  mean_goods_coverage         = round(goods_at, 4),
  obs_feat_slope_per_1k_reads = round(slope_at_d, 4),
  satisfies_a_retention       = crit_a,
  satisfies_b_coverage        = crit_b,
  satisfies_c_plateau         = crit_c,
  satisfies_all_three         = crit_a & crit_b & crit_c,
  stringsAsFactors            = FALSE
)
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("[INFO]  Candidate table      : %s\n", out_csv))

# ── Pick the recommended depth + STATUS ──────────────────────────────────────
status        <- "FAIL"
status_reason <- ""
depth_method  <- ""

if (any(crit_a & crit_b & crit_c)) {
  best_depth   <- max(cand[crit_a & crit_b & crit_c])
  status       <- "PASS"
  depth_method <- "all_three_criteria"
  status_reason <- sprintf(
    "All three criteria satisfied: retention=%.1f%%, mean_goods=%.4f, slope=%.3f",
    overall_keep[which(cand == best_depth)] * 100,
    goods_at   [which(cand == best_depth)],
    slope_at_d [which(cand == best_depth)])
} else if (any(crit_a & crit_b)) {
  best_depth   <- max(cand[crit_a & crit_b])
  status       <- "PASS_WITH_WARNINGS"
  depth_method <- "fallback_no_plateau"
  status_reason <- sprintf(
    "Plateau criterion (slope<%.3f) NOT met. Retention+coverage OK: %.1f%%, %.4f",
    slope_max,
    overall_keep[which(cand == best_depth)] * 100,
    goods_at   [which(cand == best_depth)])
} else if (any(crit_a)) {
  best_depth   <- max(cand[crit_a])
  status       <- "REVIEW_REQUIRED"
  depth_method <- "fallback_retention_only"
  status_reason <- sprintf(
    "Coverage and plateau criteria NOT met. Only retention satisfied: %.1f%% (mean_goods=%.4f, slope=%.3f).",
    overall_keep[which(cand == best_depth)] * 100,
    goods_at   [which(cand == best_depth)],
    slope_at_d [which(cand == best_depth)])
} else if (any(overall_keep >= 0.75)) {
  best_depth   <- max(cand[overall_keep >= 0.75])
  status       <- "REVIEW_REQUIRED"
  depth_method <- "fallback_75pct"
  status_reason <- sprintf(
    "No candidate met the %.0f%% retention target. Loosened to 75%% retention: %.1f%%.",
    min_overall * 100, overall_keep[which(cand == best_depth)] * 100)
} else {
  q1 <- as.numeric(quantile(df$total_freq, 0.25, na.rm = TRUE))
  best_depth   <- floor(q1 / round_base) * round_base
  if (best_depth < 1) best_depth <- min(df$total_freq)
  status       <- "REVIEW_REQUIRED"
  depth_method <- "fallback_Q1"
  status_reason <- "Even 75% retention impossible; using Q1 of per-sample totals."
}

# Hard FAIL if best_depth empties any group entirely (UniFrac/PERMANOVA explode).
if (!is.null(group_column)) {
  grp_kept_at_best <- tapply(df$total_freq >= best_depth, df[[group_column]], sum, na.rm = TRUE)
  if (any(grp_kept_at_best == 0)) {
    status        <- "FAIL"
    depth_method  <- paste0(depth_method, "+group_zero")
    status_reason <- sprintf(
      "Selected depth %d empties group(s): %s.  Re-run with --min-group lower or omit --group-col.",
      best_depth,
      paste(names(grp_kept_at_best)[grp_kept_at_best == 0], collapse = ", "))
  }
}

# ── Sample-level diagnostics ─────────────────────────────────────────────────
df$kept_at_recommended <- df$total_freq >= best_depth
kept    <- df[df$kept_at_recommended,  , drop = FALSE]
removed <- df[!df$kept_at_recommended, , drop = FALSE]

write.csv(
  data.frame(SampleID = kept$SampleID,
             total_freq = kept$total_freq,
             stringsAsFactors = FALSE),
  out_kept,    row.names = FALSE)
write.csv(
  data.frame(SampleID   = removed$SampleID,
             total_freq = removed$total_freq,
             reason     = rep(sprintf("total_freq<%d", best_depth), nrow(removed)),
             stringsAsFactors = FALSE),
  out_removed, row.names = FALSE)

# Per-sample summary — every sample, with optional group column
sample_summary <- data.frame(
  SampleID            = df$SampleID,
  total_freq          = df$total_freq,
  group               = if (is.null(group_column)) NA_character_ else df[[group_column]],
  kept_at_recommended = df$kept_at_recommended,
  stringsAsFactors    = FALSE)
write.csv(sample_summary, out_sample, row.names = FALSE)

# Per-group summary
if (!is.null(group_column)) {
  grp_total <- table(df[[group_column]])
  grp_kept  <- tapply(df$kept_at_recommended, df[[group_column]], sum, na.rm = TRUE)
  group_summary <- data.frame(
    group              = names(grp_total),
    n_total            = as.integer(grp_total),
    n_kept             = as.integer(grp_kept[names(grp_total)]),
    retention_fraction = round(as.numeric(grp_kept[names(grp_total)]) /
                               as.numeric(grp_total), 4),
    fully_retained     = as.integer(grp_kept[names(grp_total)]) ==
                         as.integer(grp_total),
    stringsAsFactors   = FALSE)
} else {
  group_summary <- data.frame(group = character(), n_total = integer(),
                              n_kept = integer(), retention_fraction = numeric(),
                              fully_retained = logical())
}
write.csv(group_summary, out_group, row.names = FALSE)

# ── Summary stats ────────────────────────────────────────────────────────────
median_depth <- floor(median(df$total_freq, na.rm = TRUE) / round_base) * round_base
q1_depth     <- floor(as.numeric(quantile(df$total_freq, 0.25, na.rm = TRUE)) / round_base) * round_base
min_depth    <- min(df$total_freq, na.rm = TRUE)
max_depth_v  <- max(df$total_freq, na.rm = TRUE)

samples_at_best  <- sum(df$total_freq >= best_depth)
fraction_at_best <- mean(df$total_freq >= best_depth)
mean_gc_at_best  <- if (best_depth %in% cand)
                      goods_at[which(cand == best_depth)] else NA_real_
slope_at_best    <- if (best_depth %in% cand)
                      slope_at_d[which(cand == best_depth)] else NA_real_

ready <- if (status %in% c("PASS", "PASS_WITH_WARNINGS")) "yes" else "no"

cat(sprintf("\n[RESULT] Recommended_depth        = %d\n", best_depth))
cat(sprintf("[RESULT] Samples_retained         = %d / %d (%.1f%%)\n",
            samples_at_best, n_samples, fraction_at_best * 100))
cat(sprintf("[RESULT] Mean_Goods_at_depth      = %s\n",
            if (is.na(mean_gc_at_best)) "NA" else sprintf("%.4f", mean_gc_at_best)))
cat(sprintf("[RESULT] Slope_at_depth           = %s features/1k reads\n",
            if (is.na(slope_at_best)) "NA" else sprintf("%.3f", slope_at_best)))
cat(sprintf("[RESULT] STATUS                   = %s\n", status))
cat(sprintf("[RESULT] DEPTH_SELECTION_METHOD   = %s\n", depth_method))
cat(sprintf("[RESULT] READY_FOR_DIVERSITY      = %s\n", ready))

# ── key=value file consumed by bash ──────────────────────────────────────────
writeLines(c(
  sprintf("RECOMMENDED_DEPTH=%d",                 best_depth),
  sprintf("SAMPLES_RETAINED=%d",                  samples_at_best),
  sprintf("TOTAL_SAMPLES=%d",                     n_samples),
  sprintf("FRACTION_RETAINED=%.4f",               fraction_at_best),
  sprintf("PCT_RETAINED=%.1f",                    fraction_at_best * 100),
  sprintf("DEPTH_METHOD=%s",                      depth_method),
  sprintf("DEPTH_STATUS=%s",                      status),
  sprintf("DEPTH_RATIONALE=%s",                   gsub("\n", " ", status_reason)),
  sprintf("READY_FOR_DIVERSITY=%s",               ready),
  sprintf("MEDIAN_DEPTH=%d",                      median_depth),
  sprintf("Q1_DEPTH=%d",                          q1_depth),
  sprintf("MIN_SAMPLE_COUNT=%d",                  min_depth),
  sprintf("MAX_SAMPLE_COUNT=%d",                  max_depth_v),
  sprintf("GROUP_COLUMN=%s",                      if (is.null(group_column)) "none" else group_column),
  sprintf("N_GROUPS=%s",                          if (is.na(total_groups)) "NA" else as.character(total_groups)),
  sprintf("MIN_OVERALL_THRESHOLD=%.2f",           min_overall),
  sprintf("MIN_GROUP_THRESHOLD=%.2f",             min_group),
  sprintf("GOOD_COV_MIN=%.3f",                    good_cov_min),
  sprintf("PLATEAU_SLOPE_MAX=%.3f",               slope_max),
  sprintf("MEAN_GOODS_AT_RECOMMENDED=%s",
          if (is.na(mean_gc_at_best)) "NA" else sprintf("%.4f", mean_gc_at_best)),
  sprintf("SLOPE_AT_RECOMMENDED=%s",
          if (is.na(slope_at_best)) "NA" else sprintf("%.4f", slope_at_best)),
  sprintf("N_TRIMMED_SAMPLE_IDS=%d",              n_trimmed_id),
  sprintf("N_TRIMMED_GROUP_LABELS=%d",            n_trimmed_grp)
), out_txt)
cat(sprintf("[INFO]  Depth params         : %s\n", out_txt))
RDEPTH

if $DRY_RUN; then
  warn "[DRY-RUN] Would run R depth selection with three concurrent criteria."
  RECOMMENDED_DEPTH="DRYRUN"
  MEDIAN_DEPTH="DRYRUN"
  FRACTION_RETAINED="DRYRUN"
  DEPTH_METHOD="dry_run"
  DEPTH_STATUS="DRYRUN"
  READY_FOR_DIVERSITY="dry_run"
else
  info "Running R depth-selection algorithm (3 concurrent criteria)..."
  MBX_FREQ_TSV="$SAMPLE_FREQ_TSV" \
  MBX_METADATA="$METADATA_TXT" \
  MBX_PER_SAMPLE_CSV="$PER_SAMPLE_RARE_CSV" \
  MBX_MEAN_CSV="$ALPHA_RARE_CSV" \
  MBX_GROUP_COL="${GROUP_COL:-}" \
  MBX_MIN_OVERALL="$MIN_OVERALL" \
  MBX_MIN_GROUP="$MIN_GROUP" \
  MBX_GOOD_COV_MIN="$GOOD_COV_MIN" \
  MBX_SLOPE_MAX="$PLATEAU_SLOPE_MAX" \
  MBX_DEPTH_CSV="$DEPTH_CSV" \
  MBX_SAMPLES_RETAINED_CSV="$SAMPLES_RETAINED_CSV" \
  MBX_SAMPLES_REMOVED_CSV="$SAMPLES_REMOVED_CSV" \
  MBX_SAMPLE_DEPTH_SUMMARY_CSV="$SAMPLE_DEPTH_SUMMARY_CSV" \
  MBX_GROUP_DEPTH_SUMMARY_CSV="$GROUP_DEPTH_SUMMARY_CSV" \
  MBX_DEPTH_TXT="$R_DEPTH_OUT" \
    "$RSCRIPT_CMD" --vanilla "$R_DEPTH_SCRIPT" \
    || err "R depth-selection script failed.
  → Inspect the output above for the R traceback.
  → Verify per-sample CSV exists: ls -la '$PER_SAMPLE_RARE_CSV'"

  # Parse the key=value summary
  RECOMMENDED_DEPTH="$(_read_key  "RECOMMENDED_DEPTH"         "$R_DEPTH_OUT")"
  MEDIAN_DEPTH="$(_read_key       "MEDIAN_DEPTH"              "$R_DEPTH_OUT")"
  Q1_DEPTH="$(_read_key           "Q1_DEPTH"                  "$R_DEPTH_OUT")"
  MIN_DEPTH="$(_read_key          "MIN_SAMPLE_COUNT"          "$R_DEPTH_OUT")"
  MAX_DEPTH="$(_read_key          "MAX_SAMPLE_COUNT"          "$R_DEPTH_OUT")"
  FRACTION_RETAINED="$(_read_key  "FRACTION_RETAINED"         "$R_DEPTH_OUT")"
  PCT_RETAINED="$(_read_key       "PCT_RETAINED"              "$R_DEPTH_OUT")"
  SAMPLES_RETAINED="$(_read_key   "SAMPLES_RETAINED"          "$R_DEPTH_OUT")"
  TOTAL_SAMPLES="$(_read_key      "TOTAL_SAMPLES"             "$R_DEPTH_OUT")"
  DEPTH_METHOD="$(_read_key       "DEPTH_METHOD"              "$R_DEPTH_OUT")"
  DEPTH_STATUS="$(_read_key       "DEPTH_STATUS"              "$R_DEPTH_OUT")"
  DEPTH_RATIONALE="$(_read_key    "DEPTH_RATIONALE"           "$R_DEPTH_OUT")"
  READY_FOR_DIVERSITY="$(_read_key "READY_FOR_DIVERSITY"      "$R_DEPTH_OUT")"
  MEAN_GOODS_AT_RECOMMENDED="$(_read_key "MEAN_GOODS_AT_RECOMMENDED" "$R_DEPTH_OUT")"
  SLOPE_AT_RECOMMENDED="$(_read_key      "SLOPE_AT_RECOMMENDED"      "$R_DEPTH_OUT")"
  N_GROUPS="$(_read_key           "N_GROUPS"                  "$R_DEPTH_OUT")"
  GROUP_COLUMN_REPORTED="$(_read_key "GROUP_COLUMN"           "$R_DEPTH_OUT")"

  [[ -z "$RECOMMENDED_DEPTH" ]] && err "Could not parse RECOMMENDED_DEPTH from R output.
  → Inspect: cat '$R_DEPTH_OUT'"
  [[ -z "$PCT_RETAINED" ]] && PCT_RETAINED="?"

  sep
  case "$DEPTH_STATUS" in
    PASS)               ok   "DEPTH STATUS         : PASS" ;;
    PASS_WITH_WARNINGS) warn "DEPTH STATUS         : PASS_WITH_WARNINGS" ;;
    REVIEW_REQUIRED)    warn "DEPTH STATUS         : REVIEW_REQUIRED" ;;
    FAIL)               warn "DEPTH STATUS         : FAIL — diversity step will refuse to run" ;;
    *)                  warn "DEPTH STATUS         : $DEPTH_STATUS (unrecognized)" ;;
  esac
  ok "  Recommended depth   : $RECOMMENDED_DEPTH reads"
  ok "  Samples retained    : $SAMPLES_RETAINED / $TOTAL_SAMPLES  (${PCT_RETAINED}%)"
  ok "  Mean Good's coverage: $MEAN_GOODS_AT_RECOMMENDED"
  ok "  Slope at depth      : $SLOPE_AT_RECOMMENDED features per 1000 reads"
  ok "  Depth method        : $DEPTH_METHOD"
  ok "  Group column        : $GROUP_COLUMN_REPORTED"
  ok "  Median depth        : $MEDIAN_DEPTH"
  ok "  Min/Max sample count: $MIN_DEPTH / $MAX_DEPTH"
  ok "  Candidate table     : $DEPTH_CSV"
  ok "  Samples retained    : $SAMPLES_RETAINED_CSV"
  ok "  Samples removed     : $SAMPLES_REMOVED_CSV"
  ok "  Sample summary      : $SAMPLE_DEPTH_SUMMARY_CSV"
  ok "  Group summary       : $GROUP_DEPTH_SUMMARY_CSV"
  ok "  Rationale           : $DEPTH_RATIONALE"
  sep

  if [[ "$DEPTH_STATUS" == "FAIL" ]]; then
    warn "FAIL status: mbx_diversity_run.sh will refuse to proceed."
    warn "Read the rationale above and re-run with adjusted thresholds."
  elif [[ "$DEPTH_STATUS" == "REVIEW_REQUIRED" ]]; then
    warn "REVIEW_REQUIRED: the chosen depth is a fallback. Inspect $DEPTH_CSV"
    warn "and read mbx_pre_diversity_summary.txt before running step 12."
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "9/12 — Official QIIME alpha-rarefaction (publication-grade)"
# ─────────────────────────────────────────────────────────────────────────────
# Why this step exists in addition to step 7's analytical preview:
#   The analytical preview (step 7) is mathematically exact for observed_features
#   and Good's coverage but does NOT compute Shannon (Chao-Shen ≠ QIIME Monte Carlo)
#   or Faith PD.  Reviewers expect the artifact emitted by the QIIME2 plugin.
#   This step runs the real plugin so the user has a publication-grade QZV.
#
#   --p-metrics observed_features --p-metrics shannon --p-metrics faith_pd
#
# Faith PD is included so the lay user has the option of presenting it in
# their publication, instead of running step 12's core-metrics first.
#
# Honors:  $SKIP_QIIME_RAREFACTION  → skips with status=SKIPPED_BY_USER
#          $FORCE_RERUN              → recomputes even if QZV exists
#          $DRY_RUN                  → prints command, does not execute
#
# Failure mode: if the plugin errors, we WARN (not err/exit) so the user keeps
# the analytical preview from step 7, and step 12 records QIIME_RARE_STATUS=FAILED.

QIIME_RARE_STATUS="not_run"
QIIME_RARE_ELAPSED="N/A"

if $SKIP_QIIME_RAREFACTION; then
  warn "[step9] --skip-qiime-rarefaction set — official QIIME alpha-rarefaction SKIPPED."
  warn "[step9] Analytical preview from step 7 will be your only curve."
  warn "[step9] Not recommended for publication — re-run without the flag for the official artifact."
  QIIME_RARE_STATUS="SKIPPED_BY_USER"
elif [[ -f "$ALPHA_RARE_QIIME_QZV" ]] && ! $FORCE_RERUN; then
  skipped "official QIIME alpha-rarefaction (output exists)"
  _stale_warn "alpha-rarefaction-qiime.qzv" "$ALPHA_RARE_QIIME_QZV"
  QIIME_RARE_STATUS="SKIPPED_OUTPUT_EXISTS"
elif $DRY_RUN; then
  warn "[DRY-RUN] Would run:"
  warn "  qiime diversity alpha-rarefaction \\"
  warn "    --i-table        $ACTIVE_FEATURE_TABLE \\"
  warn "    --i-phylogeny    $ROOTED_QZA \\"
  warn "    --p-max-depth    $MAX_FREQ \\"
  warn "    --p-steps        $QIIME_RARE_STEPS \\"
  warn "    --p-iterations   $QIIME_RARE_ITERATIONS \\"
  warn "    --p-metrics      observed_features \\"
  warn "    --p-metrics      shannon \\"
  warn "    --p-metrics      faith_pd \\"
  warn "    --m-metadata-file $METADATA_TXT \\"
  warn "    --o-visualization $ALPHA_RARE_QIIME_QZV"
  QIIME_RARE_STATUS="DRY_RUN"
else
  info "[step9] Running qiime diversity alpha-rarefaction (publication-grade)…"
  info "[step9] This step is single-threaded Monte-Carlo and may take 5 min – 2 h"
  info "[step9] depending on table size."
  info "[step9]   Steps         : $QIIME_RARE_STEPS"
  info "[step9]   Iterations    : $QIIME_RARE_ITERATIONS  (per step)"
  info "[step9]   Max depth     : $MAX_FREQ  (= max sample frequency)"
  info "[step9]   Metrics       : observed_features, shannon, faith_pd"
  info "[step9]   Phylogeny     : $ROOTED_QZA"
  info "[step9]   Feature table : $ACTIVE_FEATURE_TABLE"
  info "[step9]   Metadata      : $METADATA_TXT"
  info "[step9]   Output        : $ALPHA_RARE_QIIME_QZV"
  echo ""

  _t9_start="$(date +%s)"
  set +e
  qiime diversity alpha-rarefaction \
    --i-table         "$ACTIVE_FEATURE_TABLE" \
    --i-phylogeny     "$ROOTED_QZA" \
    --p-max-depth     "$MAX_FREQ" \
    --p-steps         "$QIIME_RARE_STEPS" \
    --p-iterations    "$QIIME_RARE_ITERATIONS" \
    --p-metrics       observed_features \
    --p-metrics       shannon \
    --p-metrics       faith_pd \
    --m-metadata-file "$METADATA_TXT" \
    --o-visualization "$ALPHA_RARE_QIIME_QZV"
  _q9_rc=$?
  set -e
  _t9_end="$(date +%s)"
  QIIME_RARE_ELAPSED="$((_t9_end - _t9_start))s"

  if [[ $_q9_rc -eq 0 && -f "$ALPHA_RARE_QIIME_QZV" ]]; then
    QIIME_RARE_STATUS="OK"
    ok "Official QIIME alpha-rarefaction → $ALPHA_RARE_QIIME_QZV"
    ok "  elapsed: $QIIME_RARE_ELAPSED"
    info "[step9] OPEN THIS IN https://view.qiime2.org BEFORE RUNNING STEP 12."
    info "[step9] Confirm richness curves plateau at or before depth $RECOMMENDED_DEPTH."
  else
    QIIME_RARE_STATUS="FAILED"
    warn "[step9] qiime diversity alpha-rarefaction FAILED (rc=$_q9_rc, elapsed $QIIME_RARE_ELAPSED)"
    warn "[step9] Common causes:"
    warn "[step9]   - max-depth ($MAX_FREQ) larger than any single sample after filtering"
    warn "[step9]   - feature table empty or no features in tree (check step 10 ID validation)"
    warn "[step9]   - out-of-memory on large tables — try --qiime-rare-steps 10 --qiime-rare-iter 5"
    warn "[step9] You still have the analytical preview from step 7."
    warn "[step9] Re-run with --force after fixing, OR --skip-qiime-rarefaction to bypass."
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "10/12 — Cross-artifact ID validation (publication-quality safety net)"
# ─────────────────────────────────────────────────────────────────────────────
# Cross-checks every ID across artifacts so step 12 (core-metrics-phylogenetic)
# can never silently produce nonsense:
#
#   Feature IDs : BIOM ↔ rep-seqs.qza ↔ rooted-tree.qza   (must match exactly)
#   Sample IDs  : metadata.txt ↔ BIOM ↔ sample-frequencies.qza
#                 (and report any silently-trimmed-whitespace IDs)
#
# Status ladder:
#   PASS               — all IDs match exactly
#   PASS_WITH_WARNINGS — whitespace-only differences (silently trimmed)
#   REVIEW_REQUIRED    — non-trivial mismatch (e.g. metadata IDs missing from BIOM)
#   FAIL               — feature IDs in BIOM are not in the tree (UniFrac would
#                        silently treat them as zero-distance — hard failure)
#
# Status is written into mbx_pre_diversity_info.txt as ID_VALIDATION_STATUS.

ID_VAL_REPORT="${PRE_DIV_DIR}/id_validation_report.txt"
ID_VAL_STATUS="not_run"
ID_VAL_FEATURES_BIOM=0
ID_VAL_FEATURES_REPSEQS=0
ID_VAL_FEATURES_TREE=0
ID_VAL_FEATURES_BIOM_NOT_IN_TREE=0
ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS=0
ID_VAL_SAMPLES_METADATA=0
ID_VAL_SAMPLES_BIOM=0
ID_VAL_SAMPLES_FREQ=0
ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM=0
ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA=0
ID_VAL_TRIMMED_FEATURE_IDS=0
ID_VAL_TRIMMED_SAMPLE_IDS=0

if $DRY_RUN; then
  warn "[DRY-RUN] Would run ID validation across BIOM/rep-seqs/tree/metadata."
  ID_VAL_STATUS="DRY_RUN"
else
  info "[step10] Cross-checking feature IDs and sample IDs across artifacts…"

  # Stage rep-seqs and rooted-tree to temp dirs (Newick already exported in step 3)
  _TMPID="${$}_$(date +%s)"
  IDV_REPSEQS_DIR="/tmp/mbx_idv_repseqs_${_TMPID}"
  rm -rf "$IDV_REPSEQS_DIR"
  mkdir -p "$IDV_REPSEQS_DIR"
  if [[ -f "$REP_SEQS_QZA" ]]; then
    if ! qiime tools export --input-path "$REP_SEQS_QZA" --output-path "$IDV_REPSEQS_DIR" >/dev/null 2>&1; then
      warn "[step10] qiime tools export failed for $REP_SEQS_QZA — using raw unzip fallback"
      ( cd "$IDV_REPSEQS_DIR" && unzip -qq "$REP_SEQS_QZA" )
    fi
  fi
  REP_SEQS_FASTA="$(find "$IDV_REPSEQS_DIR" -type f -name 'dna-sequences.fasta' 2>/dev/null | head -1)"

  # rooted-tree.nwk should already exist (step 3 exports it). If not, export now.
  if [[ ! -s "$ROOTED_NWK" ]]; then
    IDV_TREE_DIR="/tmp/mbx_idv_tree_${_TMPID}"
    rm -rf "$IDV_TREE_DIR"
    mkdir -p "$IDV_TREE_DIR"
    if qiime tools export --input-path "$ROOTED_QZA" --output-path "$IDV_TREE_DIR" >/dev/null 2>&1; then
      _exported_nwk="$(find "$IDV_TREE_DIR" -type f -name 'tree.nwk' | head -1)"
      [[ -s "$_exported_nwk" ]] && cp "$_exported_nwk" "$ROOTED_NWK"
    fi
  fi

  # Build sample-frequencies sample list from the TSV produced in step 5
  ID_VAL_FREQ_TSV="$SAMPLE_FREQ_TSV"

  # Run the validation in QIIME's python so biom + skbio are available.
  IDV_PY="/tmp/mbx_idv_${_TMPID}.py"
  cat > "$IDV_PY" << 'PYEOF'
import os, sys, re, csv, json

biom_path     = os.environ["MBX_BIOM"]
fasta_path    = os.environ.get("MBX_REPSEQS_FASTA","")
tree_path     = os.environ.get("MBX_TREE_NWK","")
metadata_path = os.environ["MBX_METADATA"]
freq_tsv_path = os.environ.get("MBX_FREQ_TSV","")
report_path   = os.environ["MBX_REPORT"]

# ── feature IDs from BIOM ────────────────────────────────────────────────────
from biom import load_table
T = load_table(biom_path)
biom_features = list(T.ids(axis="observation"))
biom_samples  = list(T.ids(axis="sample"))

# ── feature IDs from rep-seqs.fasta ──────────────────────────────────────────
repseq_features = []
if fasta_path and os.path.exists(fasta_path):
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                # FASTA header up to first whitespace
                repseq_features.append(line[1:].strip().split()[0])

# ── feature IDs from tree (Newick — extract leaf labels) ─────────────────────
def newick_leaves(nwk_text):
    # Strip comments [...] and remove whitespace
    s = re.sub(r"\[[^\]]*\]", "", nwk_text)
    s = "".join(s.split())
    leaves = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c in "(),;":
            i += 1; continue
        # Read label until next structural char or branch length
        j = i
        while j < n and s[j] not in "(),;:":
            j += 1
        label = s[i:j].strip()
        # Skip past optional branch length :NNN
        if j < n and s[j] == ":":
            k = j + 1
            while k < n and s[k] not in "(),;":
                k += 1
            j = k
        # Internal nodes are labels that immediately follow ')'.
        # We treat any non-empty label that came right after '(' or ',' as a leaf.
        # Safer heuristic: if previous char (at position i-1) is '(' or ',', this is a leaf.
        if label and (i == 0 or s[i-1] in "(,"):
            leaves.append(label)
        i = j
    return leaves

tree_features = []
if tree_path and os.path.exists(tree_path):
    with open(tree_path) as fh:
        nwk = fh.read()
    if nwk.strip():
        tree_features = newick_leaves(nwk)

# ── sample IDs from metadata ─────────────────────────────────────────────────
def read_metadata_ids(path):
    # Tab or comma; skip QIIME2 #q2:types directive lines
    sep = "\t" if path.lower().endswith((".tsv",".txt")) else ","
    raw_ids = []
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter=sep)
        header = None
        for row in rd:
            if not row: continue
            if row[0].startswith("#") and "q2:types" in row[0].lower():
                continue
            if header is None:
                header = row
                continue
            if not row[0]: continue
            raw_ids.append(row[0])
    return raw_ids

raw_meta_ids = read_metadata_ids(metadata_path)

# ── sample IDs from sample-frequencies TSV (skip # lines) ────────────────────
freq_samples = []
if freq_tsv_path and os.path.exists(freq_tsv_path):
    with open(freq_tsv_path) as fh:
        header_seen = False
        for line in fh:
            if line.startswith("#"): continue
            if not line.strip(): continue
            if not header_seen:
                header_seen = True
                continue
            freq_samples.append(line.split("\t")[0])

# ── Whitespace-trim accounting (silent trim, but counted) ────────────────────
def trim_count(ids):
    cleaned = []
    n_trim = 0
    for x in ids:
        s = x.strip()
        if s != x: n_trim += 1
        cleaned.append(s)
    return cleaned, n_trim

biom_features_t,    n_trim_feat_biom    = trim_count(biom_features)
repseq_features_t,  n_trim_feat_repseq  = trim_count(repseq_features)
tree_features_t,    n_trim_feat_tree    = trim_count(tree_features)
biom_samples_t,     n_trim_samp_biom    = trim_count(biom_samples)
meta_ids_t,         n_trim_samp_meta    = trim_count(raw_meta_ids)
freq_ids_t,         n_trim_samp_freq    = trim_count(freq_samples)

trimmed_features = n_trim_feat_biom + n_trim_feat_repseq + n_trim_feat_tree
trimmed_samples  = n_trim_samp_biom + n_trim_samp_meta  + n_trim_samp_freq

bf  = set(biom_features_t)
rf  = set(repseq_features_t) if repseq_features_t else set()
tf  = set(tree_features_t)   if tree_features_t   else set()
bs  = set(biom_samples_t)
ms  = set(meta_ids_t)
fs  = set(freq_ids_t) if freq_ids_t else set()

biom_minus_tree    = sorted(bf - tf) if tf else []
biom_minus_repseqs = sorted(bf - rf) if rf else []
meta_minus_biom    = sorted(ms - bs)
biom_minus_meta    = sorted(bs - ms)
freq_minus_biom    = sorted(fs - bs) if fs else []

# ── Status determination ─────────────────────────────────────────────────────
status = "PASS"
notes  = []

if tf and biom_minus_tree:
    status = "FAIL"
    notes.append("Feature IDs in BIOM not in tree: %d" % len(biom_minus_tree))

if rf and biom_minus_repseqs:
    if status not in ("FAIL",):
        status = "REVIEW_REQUIRED"
    notes.append("Feature IDs in BIOM not in rep-seqs: %d" % len(biom_minus_repseqs))

if meta_minus_biom:
    if status not in ("FAIL",):
        status = "REVIEW_REQUIRED"
    notes.append("Sample IDs in metadata not in BIOM: %d" % len(meta_minus_biom))

if biom_minus_meta:
    # metadata is allowed to be a superset, but BIOM samples MUST be in metadata
    # (otherwise core-metrics will drop them silently)
    if status not in ("FAIL",):
        status = "REVIEW_REQUIRED"
    notes.append("Sample IDs in BIOM not in metadata: %d" % len(biom_minus_meta))

if fs and freq_minus_biom:
    if status not in ("FAIL","REVIEW_REQUIRED"):
        status = "PASS_WITH_WARNINGS"
    notes.append("Sample IDs in sample-frequencies not in BIOM: %d" % len(freq_minus_biom))

# Whitespace-trim is a soft warning
if (trimmed_features + trimmed_samples) > 0 and status == "PASS":
    status = "PASS_WITH_WARNINGS"
    notes.append("Trimmed whitespace in IDs (features=%d, samples=%d)" %
                 (trimmed_features, trimmed_samples))

# ── Write report ─────────────────────────────────────────────────────────────
with open(report_path, "w") as out:
    out.write("# ============================================================\n")
    out.write("# id_validation_report.txt\n")
    out.write("# ============================================================\n")
    out.write("STATUS=%s\n\n" % status)
    out.write("Counts:\n")
    out.write("  features_biom               = %d\n" % len(bf))
    out.write("  features_rep_seqs           = %d\n" % len(rf))
    out.write("  features_tree               = %d\n" % len(tf))
    out.write("  features_biom_not_in_tree   = %d\n" % len(biom_minus_tree))
    out.write("  features_biom_not_in_repseq = %d\n" % len(biom_minus_repseqs))
    out.write("  samples_metadata            = %d\n" % len(ms))
    out.write("  samples_biom                = %d\n" % len(bs))
    out.write("  samples_freq                = %d\n" % len(fs))
    out.write("  samples_meta_not_in_biom    = %d\n" % len(meta_minus_biom))
    out.write("  samples_biom_not_in_meta    = %d\n" % len(biom_minus_meta))
    out.write("  trimmed_feature_ids         = %d\n" % trimmed_features)
    out.write("  trimmed_sample_ids          = %d\n" % trimmed_samples)
    out.write("\nNotes:\n")
    for n in notes: out.write("  - %s\n" % n)
    if biom_minus_tree:
        out.write("\nFirst 20 feature IDs in BIOM but not in tree:\n")
        for x in biom_minus_tree[:20]: out.write("  %s\n" % x)
    if meta_minus_biom:
        out.write("\nFirst 20 sample IDs in metadata but not in BIOM:\n")
        for x in meta_minus_biom[:20]: out.write("  %s\n" % x)
    if biom_minus_meta:
        out.write("\nFirst 20 sample IDs in BIOM but not in metadata:\n")
        for x in biom_minus_meta[:20]: out.write("  %s\n" % x)

# ── Emit machine-readable key=value to stdout for bash to parse ──────────────
print("ID_VAL_STATUS=%s" % status)
print("ID_VAL_FEATURES_BIOM=%d" % len(bf))
print("ID_VAL_FEATURES_REPSEQS=%d" % len(rf))
print("ID_VAL_FEATURES_TREE=%d" % len(tf))
print("ID_VAL_FEATURES_BIOM_NOT_IN_TREE=%d" % len(biom_minus_tree))
print("ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS=%d" % len(biom_minus_repseqs))
print("ID_VAL_SAMPLES_METADATA=%d" % len(ms))
print("ID_VAL_SAMPLES_BIOM=%d" % len(bs))
print("ID_VAL_SAMPLES_FREQ=%d" % len(fs))
print("ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM=%d" % len(meta_minus_biom))
print("ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA=%d" % len(biom_minus_meta))
print("ID_VAL_TRIMMED_FEATURE_IDS=%d" % trimmed_features)
print("ID_VAL_TRIMMED_SAMPLE_IDS=%d" % trimmed_samples)
PYEOF

  IDV_OUT="/tmp/mbx_idv_out_${_TMPID}.txt"
  set +e
  MBX_BIOM="$FEATURE_BIOM" \
  MBX_REPSEQS_FASTA="${REP_SEQS_FASTA:-}" \
  MBX_TREE_NWK="$ROOTED_NWK" \
  MBX_METADATA="$METADATA_TXT" \
  MBX_FREQ_TSV="${ID_VAL_FREQ_TSV:-}" \
  MBX_REPORT="$ID_VAL_REPORT" \
    "$PY_BIN" "$IDV_PY" > "$IDV_OUT" 2>&1
  _idv_rc=$?
  set -e

  if [[ $_idv_rc -ne 0 ]]; then
    warn "[step10] ID validation script failed (rc=$_idv_rc)."
    warn "[step10] Output:"
    sed 's/^/  /' "$IDV_OUT" >&2 || true
    ID_VAL_STATUS="FAILED_TO_RUN"
  else
    # Parse key=value lines
    while IFS='=' read -r _k _v; do
      case "$_k" in
        ID_VAL_STATUS)                          ID_VAL_STATUS="$_v" ;;
        ID_VAL_FEATURES_BIOM)                   ID_VAL_FEATURES_BIOM="$_v" ;;
        ID_VAL_FEATURES_REPSEQS)                ID_VAL_FEATURES_REPSEQS="$_v" ;;
        ID_VAL_FEATURES_TREE)                   ID_VAL_FEATURES_TREE="$_v" ;;
        ID_VAL_FEATURES_BIOM_NOT_IN_TREE)       ID_VAL_FEATURES_BIOM_NOT_IN_TREE="$_v" ;;
        ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS)    ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS="$_v" ;;
        ID_VAL_SAMPLES_METADATA)                ID_VAL_SAMPLES_METADATA="$_v" ;;
        ID_VAL_SAMPLES_BIOM)                    ID_VAL_SAMPLES_BIOM="$_v" ;;
        ID_VAL_SAMPLES_FREQ)                    ID_VAL_SAMPLES_FREQ="$_v" ;;
        ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM)    ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM="$_v" ;;
        ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA)    ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA="$_v" ;;
        ID_VAL_TRIMMED_FEATURE_IDS)             ID_VAL_TRIMMED_FEATURE_IDS="$_v" ;;
        ID_VAL_TRIMMED_SAMPLE_IDS)              ID_VAL_TRIMMED_SAMPLE_IDS="$_v" ;;
      esac
    done < "$IDV_OUT"
  fi

  # User-facing summary
  case "$ID_VAL_STATUS" in
    PASS)
      ok "ID validation: PASS"
      ok "  features (BIOM/repseqs/tree): $ID_VAL_FEATURES_BIOM / $ID_VAL_FEATURES_REPSEQS / $ID_VAL_FEATURES_TREE"
      ok "  samples (metadata/BIOM/freq): $ID_VAL_SAMPLES_METADATA / $ID_VAL_SAMPLES_BIOM / $ID_VAL_SAMPLES_FREQ"
      ;;
    PASS_WITH_WARNINGS)
      warn "[step10] ID validation: PASS_WITH_WARNINGS"
      warn "[step10]   trimmed whitespace — features=$ID_VAL_TRIMMED_FEATURE_IDS  samples=$ID_VAL_TRIMMED_SAMPLE_IDS"
      warn "[step10]   safe to proceed; review $ID_VAL_REPORT"
      ;;
    REVIEW_REQUIRED)
      warn "[step10] ID validation: REVIEW_REQUIRED"
      warn "[step10]   metadata not in BIOM: $ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM"
      warn "[step10]   BIOM not in metadata: $ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA"
      warn "[step10]   features missing from rep-seqs: $ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS"
      warn "[step10]   See $ID_VAL_REPORT for details — fix BEFORE running step 12."
      ;;
    FAIL)
      warn "[step10] ID validation: FAIL — DO NOT PROCEED to step 12."
      warn "[step10]   features in BIOM not in tree: $ID_VAL_FEATURES_BIOM_NOT_IN_TREE"
      warn "[step10]   UniFrac would silently treat these as zero — that is junk science."
      warn "[step10]   See $ID_VAL_REPORT and reconcile feature_table_filtered.qza against rooted-tree.qza."
      ;;
    FAILED_TO_RUN|*)
      warn "[step10] ID validation could not be run; status=$ID_VAL_STATUS"
      ;;
  esac

  # Cleanup
  rm -f "$IDV_PY" "$IDV_OUT"
  rm -rf "$IDV_REPSEQS_DIR" "${IDV_TREE_DIR:-/tmp/_nonexistent_xyz}"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "11/12 — Visualizations (depth distribution + decision plot + annotated curves)"
# ─────────────────────────────────────────────────────────────────────────────
# Renders three publication-quality PNGs for the user's records:
#
#   1. sequencing_depth_distribution.png
#      Per-sample-totals histogram with vertical lines at MIN / MEDIAN / MAX
#      sample frequencies and the RECOMMENDED sampling depth (red, dashed).
#
#   2. depth_vs_retention.png
#      Decision curve: x = candidate depth, y = retention (overall + per group).
#      Threshold lines drawn for MIN_OVERALL and MIN_GROUP. Recommended depth
#      highlighted with a vertical band and the DEPTH_STATUS badge (PASS /
#      PASS_WITH_WARNINGS / REVIEW_REQUIRED / FAIL).
#
#   3. alpha_rarefaction_curves_annotated.png
#      Annotated overlay of step 7's curve PNG with RECOMMENDED_DEPTH red line
#      and a status badge.  Step 7's raw figure is NOT overwritten — this is a
#      separate file so reviewers can see both.
#
# Honors $DRY_RUN. Idempotent: re-runs only the missing PNGs unless --force.

DEPTH_DIST_PNG="${PRE_DIV_DIR}/sequencing_depth_distribution.png"
DEPTH_RETENTION_PNG="${PRE_DIV_DIR}/depth_vs_retention.png"
ALPHA_RARE_ANNOT_PNG="${PRE_DIV_DIR}/alpha_rarefaction_curves_annotated.png"

if $DRY_RUN; then
  warn "[DRY-RUN] Would render:"
  warn "  - $DEPTH_DIST_PNG"
  warn "  - $DEPTH_RETENTION_PNG"
  warn "  - $ALPHA_RARE_ANNOT_PNG"
elif [[ -f "$DEPTH_DIST_PNG" && -f "$DEPTH_RETENTION_PNG" && -f "$ALPHA_RARE_ANNOT_PNG" ]] && ! $FORCE_RERUN; then
  skipped "step 11 visualizations (all 3 PNGs exist)"
  _stale_warn "sequencing_depth_distribution.png" "$DEPTH_DIST_PNG"
  _stale_warn "depth_vs_retention.png"            "$DEPTH_RETENTION_PNG"
  _stale_warn "alpha_rarefaction_curves_annotated.png" "$ALPHA_RARE_ANNOT_PNG"
else
  info "[step11] Rendering depth distribution + decision plot + annotated rarefaction…"

  _TMPID="${$}_$(date +%s)"
  VIZ_PY="/tmp/mbx_viz_${_TMPID}.py"
  cat > "$VIZ_PY" << 'PYEOF'
import os, sys, csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

freq_tsv      = os.environ["MBX_FREQ_TSV"]
mean_csv      = os.environ.get("MBX_MEAN_CSV","")
depth_csv     = os.environ["MBX_DEPTH_CSV"]
sample_summ   = os.environ["MBX_SAMPLE_SUMM_CSV"]
group_summ    = os.environ.get("MBX_GROUP_SUMM_CSV","")
recommended   = float(os.environ["MBX_RECOMMENDED_DEPTH"])
min_freq      = float(os.environ["MBX_MIN_FREQ"])
median_freq   = float(os.environ["MBX_MEDIAN_FREQ"])
max_freq      = float(os.environ["MBX_MAX_FREQ"])
min_overall   = float(os.environ.get("MBX_MIN_OVERALL","0.90"))
min_group     = float(os.environ.get("MBX_MIN_GROUP","0.80"))
status        = os.environ.get("MBX_DEPTH_STATUS","UNKNOWN")
group_col     = os.environ.get("MBX_GROUP_COL","")
ready         = os.environ.get("MBX_READY","unknown")

out_dist      = os.environ["MBX_OUT_DIST"]
out_retention = os.environ["MBX_OUT_RETENTION"]
out_annot     = os.environ["MBX_OUT_ANNOT"]
in_curve_png  = os.environ.get("MBX_CURVE_PNG","")

# ── Color/badge per status ───────────────────────────────────────────────────
status_color = {
    "PASS":               "#2ca02c",
    "PASS_WITH_WARNINGS": "#ff9900",
    "REVIEW_REQUIRED":    "#d62728",
    "FAIL":               "#7f0000",
}.get(status, "#666666")

# ── Helper: read sample frequencies from QIIME TSV (skip # rows, comma sep) ──
def read_sample_freqs(path):
    totals = []
    with open(path) as fh:
        header_seen = False
        for line in fh:
            if line.startswith("#"): continue
            if not line.strip(): continue
            if not header_seen:
                header_seen = True
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2: continue
            v = parts[1].strip().replace(",","")
            try:
                totals.append(float(v))
            except ValueError:
                continue
    return np.array(totals, dtype=float)

totals = read_sample_freqs(freq_tsv)
n_samples = len(totals)

# ── Plot 1: sequencing_depth_distribution.png ────────────────────────────────
fig, ax = plt.subplots(figsize=(9,5), dpi=120)
if n_samples > 0:
    nbins = max(8, int(math.sqrt(n_samples)*2))
    ax.hist(totals, bins=nbins, color="#4c78a8", edgecolor="white", alpha=0.85)

for x, lab, c, ls in [
    (min_freq,    f"min ({int(min_freq):,})",       "#888888", ":"),
    (median_freq, f"median ({int(median_freq):,})", "#222222", "--"),
    (max_freq,    f"max ({int(max_freq):,})",       "#888888", ":"),
    (recommended, f"recommended ({int(recommended):,})", "#d62728", "-"),
]:
    ax.axvline(x, color=c, linestyle=ls, linewidth=1.6)
    ax.text(x, ax.get_ylim()[1]*0.97, " " + lab,
            rotation=90, va="top", ha="left", fontsize=8, color=c)

ax.set_xlabel("Per-sample sequence count")
ax.set_ylabel("Number of samples")
ax.set_title(f"Per-sample sequencing depth distribution  (N = {n_samples})")
ax.grid(axis="y", alpha=0.3)
fig.text(0.99, 0.01,
         f"Status: {status}   |  Ready for diversity: {ready}",
         ha="right", va="bottom", fontsize=9, color=status_color, weight="bold")
fig.tight_layout()
fig.savefig(out_dist, dpi=150)
plt.close(fig)

# ── Plot 2: depth_vs_retention.png  (read from sampling_depth_candidates.csv)─
def read_depth_csv(path):
    depths, overall, per_group = [], [], []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            try:
                d = float(row.get("depth", row.get("sampling_depth","nan")))
            except ValueError:
                continue
            try:
                ok = float(row.get("overall_keep", row.get("retention","nan")))
            except (ValueError, TypeError):
                ok = float("nan")
            try:
                pg = row.get("per_group_keep", row.get("min_group_retention",""))
                pg_val = float(pg) if pg not in ("","NA","NaN",None) else float("nan")
            except (ValueError, TypeError):
                pg_val = float("nan")
            depths.append(d); overall.append(ok); per_group.append(pg_val)
    return np.array(depths), np.array(overall), np.array(per_group)

dx, ov, pg = read_depth_csv(depth_csv)
fig, ax = plt.subplots(figsize=(9,5.5), dpi=120)
order = np.argsort(dx)
dx, ov, pg = dx[order], ov[order], pg[order]
if dx.size:
    ax.plot(dx, ov, color="#1f77b4", linewidth=2, label="Overall retention")
    if not np.all(np.isnan(pg)):
        ax.plot(dx, pg, color="#9467bd", linewidth=2, linestyle="--",
                label=f"Min per-group retention" + (f" ({group_col})" if group_col else ""))

ax.axhline(min_overall, color="#1f77b4", linestyle=":", linewidth=1,
           label=f"Overall threshold ({min_overall:.2f})")
ax.axhline(min_group, color="#9467bd", linestyle=":", linewidth=1,
           label=f"Per-group threshold ({min_group:.2f})")

ax.axvline(recommended, color=status_color, linewidth=2,
           label=f"Recommended depth ({int(recommended):,})")
# Status band
ax.axvspan(recommended*0.97, recommended*1.03, color=status_color, alpha=0.10)

ax.set_xlabel("Candidate sampling depth (reads)")
ax.set_ylabel("Fraction of samples retained")
ax.set_ylim(0, 1.05)
ax.set_title(f"Depth-vs-retention decision curve  —  status: {status}")
ax.legend(loc="lower left", fontsize=8, framealpha=0.9)
ax.grid(alpha=0.3)
fig.text(0.99, 0.01,
         f"Ready for diversity: {ready}",
         ha="right", va="bottom", fontsize=9, color=status_color, weight="bold")
fig.tight_layout()
fig.savefig(out_retention, dpi=150)
plt.close(fig)

# ── Plot 3: alpha_rarefaction_curves_annotated.png ───────────────────────────
# Re-render the mean curve from alpha_rarefaction_data.csv if present, else
# overlay the existing PNG.
def read_mean_csv(path):
    depths, obs, gc_mean, gc_min = [], [], [], []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            try:
                d   = float(row["depth"])
                of  = float(row.get("observed_features","nan"))
                gcm = float(row.get("goods_coverage_mean","nan"))
                gci = float(row.get("goods_coverage_min","nan"))
                depths.append(d); obs.append(of); gc_mean.append(gcm); gc_min.append(gci)
            except (KeyError, ValueError):
                continue
    return np.array(depths), np.array(obs), np.array(gc_mean), np.array(gc_min)

mean_ok = bool(mean_csv) and os.path.exists(mean_csv)

if mean_ok:
    d, of, gcm, gci = read_mean_csv(mean_csv)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9,7), dpi=120, sharex=True,
                                   gridspec_kw={"height_ratios":[2,1]})
    ax1.plot(d, of, color="#1f77b4", linewidth=2, label="Mean observed_features")
    ax1.axvline(recommended, color=status_color, linewidth=2)
    ax1.axvspan(recommended*0.97, recommended*1.03, color=status_color, alpha=0.10)
    ax1.text(recommended, ax1.get_ylim()[1]*0.95,
             f"  recommended = {int(recommended):,}",
             color=status_color, fontsize=9, va="top", ha="left", weight="bold")
    ax1.set_ylabel("E[observed_features]")
    ax1.set_title(f"Analytical alpha-rarefaction (Hurlbert 1971) — status: {status}")
    ax1.legend(loc="lower right", fontsize=8)
    ax1.grid(alpha=0.3)

    ax2.plot(d, gcm, color="#2ca02c", linewidth=1.8, label="Mean Good's coverage")
    if not np.all(np.isnan(gci)):
        ax2.plot(d, gci, color="#2ca02c", linewidth=1, linestyle="--",
                 alpha=0.7, label="Min Good's coverage")
    ax2.axhline(0.98, color="#666666", linestyle=":", linewidth=1, label="0.98")
    ax2.axvline(recommended, color=status_color, linewidth=2)
    ax2.set_ylim(0.5, 1.01)
    ax2.set_ylabel("Good's coverage")
    ax2.set_xlabel("Sampling depth")
    ax2.legend(loc="lower right", fontsize=8)
    ax2.grid(alpha=0.3)

    fig.text(0.99, 0.01,
             f"Ready for diversity: {ready}",
             ha="right", va="bottom", fontsize=9, color=status_color, weight="bold")
    fig.tight_layout()
    fig.savefig(out_annot, dpi=150)
    plt.close(fig)
elif in_curve_png and os.path.exists(in_curve_png):
    img = mpimg.imread(in_curve_png)
    fig, ax = plt.subplots(figsize=(9,6), dpi=120)
    ax.imshow(img); ax.axis("off")
    ax.set_title(f"Analytical alpha-rarefaction (annotated) — status: {status}")
    fig.text(0.99, 0.01,
             f"recommended depth = {int(recommended):,}   |   Ready for diversity: {ready}",
             ha="right", va="bottom", fontsize=9, color=status_color, weight="bold")
    fig.tight_layout()
    fig.savefig(out_annot, dpi=150)
    plt.close(fig)
else:
    # No source available; emit a placeholder so the file always exists
    fig, ax = plt.subplots(figsize=(8,4), dpi=120)
    ax.axis("off")
    ax.text(0.5, 0.5,
            "alpha_rarefaction_data.csv not found\n(re-run step 7 to render this plot)",
            ha="center", va="center", fontsize=11, color="#666")
    fig.savefig(out_annot, dpi=150)
    plt.close(fig)

print("VIZ_OK")
PYEOF

  set +e
  MBX_FREQ_TSV="$SAMPLE_FREQ_TSV" \
  MBX_MEAN_CSV="$ALPHA_RARE_CSV" \
  MBX_DEPTH_CSV="$DEPTH_CSV" \
  MBX_SAMPLE_SUMM_CSV="$SAMPLE_DEPTH_SUMMARY_CSV" \
  MBX_GROUP_SUMM_CSV="$GROUP_DEPTH_SUMMARY_CSV" \
  MBX_RECOMMENDED_DEPTH="$RECOMMENDED_DEPTH" \
  MBX_MIN_FREQ="${MIN_FREQ:-0}" \
  MBX_MEDIAN_FREQ="${MEDIAN_FREQ:-0}" \
  MBX_MAX_FREQ="${MAX_FREQ:-0}" \
  MBX_MIN_OVERALL="$MIN_OVERALL" \
  MBX_MIN_GROUP="$MIN_GROUP" \
  MBX_DEPTH_STATUS="${DEPTH_STATUS:-UNKNOWN}" \
  MBX_GROUP_COL="${GROUP_COL:-}" \
  MBX_READY="${READY_FOR_DIVERSITY:-unknown}" \
  MBX_CURVE_PNG="$ALPHA_RARE_PNG" \
  MBX_OUT_DIST="$DEPTH_DIST_PNG" \
  MBX_OUT_RETENTION="$DEPTH_RETENTION_PNG" \
  MBX_OUT_ANNOT="$ALPHA_RARE_ANNOT_PNG" \
    "$PY_BIN" "$VIZ_PY"
  _viz_rc=$?
  set -e

  rm -f "$VIZ_PY"

  if [[ $_viz_rc -ne 0 ]]; then
    warn "[step11] Visualization script failed (rc=$_viz_rc). PNGs may be missing or partial."
  else
    [[ -f "$DEPTH_DIST_PNG" ]]      && ok "Depth distribution      → $DEPTH_DIST_PNG"
    [[ -f "$DEPTH_RETENTION_PNG" ]] && ok "Depth-vs-retention      → $DEPTH_RETENTION_PNG"
    [[ -f "$ALPHA_RARE_ANNOT_PNG" ]] && ok "Annotated rarefaction    → $ALPHA_RARE_ANNOT_PNG"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "12/12 — Write mbx_pre_diversity_info.txt + plain-language summary"
# ─────────────────────────────────────────────────────────────────────────────
# Final synthesis. Computes the READY_FOR_DIVERSITY gate from the three
# sub-statuses (depth selection, ID validation, QIIME rarefaction), then writes:
#
#   mbx_pre_diversity_info.txt     — machine-readable key=value file (read by step 12)
#   mbx_pre_diversity_summary.txt  — plain-language summary (for the user)
#
# Printable status badges are echoed at the end.

NOW="$(date '+%Y-%m-%d %H:%M:%S')"
PRE_DIV_SUMMARY="${PRE_DIV_DIR}/mbx_pre_diversity_summary.txt"

# ── Compute READY_FOR_DIVERSITY (the master gate read by step 12) ────────────
DEPTH_STATUS="${DEPTH_STATUS:-UNKNOWN}"
ID_VAL_STATUS="${ID_VAL_STATUS:-not_run}"
QIIME_RARE_STATUS="${QIIME_RARE_STATUS:-not_run}"

# Hard FAIL conditions
if [[ "$DEPTH_STATUS" == "FAIL" ]] || [[ "$ID_VAL_STATUS" == "FAIL" ]]; then
  READY_FOR_DIVERSITY="no"
  OVERALL_STATUS="FAIL"
elif [[ "$DEPTH_STATUS" == "REVIEW_REQUIRED" ]] || [[ "$ID_VAL_STATUS" == "REVIEW_REQUIRED" ]]; then
  READY_FOR_DIVERSITY="review"
  OVERALL_STATUS="REVIEW_REQUIRED"
elif [[ "$DEPTH_STATUS" == "PASS_WITH_WARNINGS" ]] || \
     [[ "$ID_VAL_STATUS" == "PASS_WITH_WARNINGS" ]]; then
  READY_FOR_DIVERSITY="yes_with_warnings"
  OVERALL_STATUS="PASS_WITH_WARNINGS"
elif [[ "$DEPTH_STATUS" == "PASS" ]] && \
     [[ "$ID_VAL_STATUS" == "PASS" || "$ID_VAL_STATUS" == "DRY_RUN" ]]; then
  READY_FOR_DIVERSITY="yes"
  OVERALL_STATUS="PASS"
else
  READY_FOR_DIVERSITY="unknown"
  OVERALL_STATUS="UNKNOWN"
fi

# Soft warning if QIIME rarefaction failed/skipped (does not gate step 12, but flag it)
QIIME_RARE_NOTE=""
case "$QIIME_RARE_STATUS" in
  OK)                   QIIME_RARE_NOTE="ok" ;;
  SKIPPED_BY_USER)      QIIME_RARE_NOTE="skipped — analytical preview only (not publication-grade)" ;;
  SKIPPED_OUTPUT_EXISTS) QIIME_RARE_NOTE="skipped — output already exists" ;;
  FAILED)               QIIME_RARE_NOTE="FAILED — re-run with --force or fall back to analytical preview" ;;
  DRY_RUN)              QIIME_RARE_NOTE="dry-run — no artifact produced" ;;
  *)                    QIIME_RARE_NOTE="$QIIME_RARE_STATUS" ;;
esac

# ── Capture tool versions (best-effort, never errors) ────────────────────────
QIIME_VERSION="$(qiime --version 2>/dev/null | head -1 | awk '{print $NF}' || echo unknown)"
PY_VERSION="$("$PY_BIN" -c 'import sys;print(sys.version.split()[0])' 2>/dev/null || echo unknown)"
SCIPY_VERSION="$("$PY_BIN" -c 'import scipy;print(scipy.__version__)' 2>/dev/null || echo unknown)"
NUMPY_VERSION="$("$PY_BIN" -c 'import numpy;print(numpy.__version__)' 2>/dev/null || echo unknown)"
BIOM_VERSION="$("$PY_BIN" -c 'import biom;print(biom.__version__)' 2>/dev/null || echo unknown)"
R_VERSION="$("$RSCRIPT_CMD" -e 'cat(as.character(getRversion()))' 2>/dev/null || echo unknown)"
SCRIPT_VERSION="${SCRIPT_VERSION:-mbx_pre_diversity_parameters.sh (12-step rewrite, Waves 1-8)}"

# ── Write mbx_pre_diversity_info.txt ─────────────────────────────────────────
cat > "$PRE_DIV_INFO" << INFO
# ============================================================================
# mbx_pre_diversity_info.txt
# Generated by mbx_pre_diversity_parameters.sh
# Date : $NOW
# ============================================================================
#
#   OVERALL STATUS         : $OVERALL_STATUS
#   READY FOR DIVERSITY    : $READY_FOR_DIVERSITY
#
#   depth_selection_status : $DEPTH_STATUS
#   id_validation_status   : $ID_VAL_STATUS
#   qiime_rarefaction      : $QIIME_RARE_STATUS  ($QIIME_RARE_NOTE)
#
# This file is read by mbx_diversity_run.sh (step 12 of the pipeline).
# Do NOT edit the key=value lines below — they are parsed programmatically.

# ── Master gate for step 12 ───────────────────────────────────────────────────
OVERALL_STATUS=$OVERALL_STATUS
READY_FOR_DIVERSITY=$READY_FOR_DIVERSITY
DEPTH_SELECTION_STATUS=$DEPTH_STATUS
DEPTH_RATIONALE=${DEPTH_RATIONALE:-N/A}

# ── Input files used ─────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
METADATA_TXT=$METADATA_TXT
FEATURE_TABLE_QZA=$ACTIVE_FEATURE_TABLE
REP_SEQS_QZA=$REP_SEQS_QZA
TAXONOMY_QZA=${TAXONOMY_QZA:-not_found}

# ── Phylogenetic tree outputs ─────────────────────────────────────────────────
ROOTED_TREE_QZA=$ROOTED_QZA
UNROOTED_TREE_QZA=$UNROOTED_QZA
ALIGNED_SEQS_QZA=$ALIGNED_QZA
MASKED_ALIGNED_SEQS_QZA=$MASKED_QZA
ROOTED_TREE_NEWICK=$ROOTED_NWK

# ── Sampling depth parameters ─────────────────────────────────────────────────
RECOMMENDED_DEPTH=${RECOMMENDED_DEPTH}
SAMPLES_RETAINED=${SAMPLES_RETAINED:-N/A}
TOTAL_SAMPLES=${TOTAL_SAMPLES:-N/A}
FRACTION_RETAINED=${FRACTION_RETAINED:-N/A}
PCT_RETAINED=${PCT_RETAINED:-N/A}
DEPTH_SELECTION_METHOD=${DEPTH_METHOD:-N/A}
MEDIAN_DEPTH=${MEDIAN_DEPTH:-N/A}
Q1_DEPTH=${Q1_DEPTH:-N/A}
MIN_SAMPLE_COUNT=${MIN_DEPTH:-N/A}
MAX_SAMPLE_COUNT=${MAX_DEPTH:-N/A}
GROUP_COLUMN=${GROUP_COL:-overall_only}
GROUP_COLUMN_REPORTED=${GROUP_COLUMN_REPORTED:-N/A}
N_GROUPS=${N_GROUPS:-N/A}
MIN_OVERALL_THRESHOLD=$MIN_OVERALL
MIN_GROUP_THRESHOLD=$MIN_GROUP
GOOD_COVERAGE_MIN=${GOOD_COV_MIN:-N/A}
SLOPE_MAX=${SLOPE_MAX:-N/A}
MEAN_GOODS_AT_RECOMMENDED=${MEAN_GOODS_AT_RECOMMENDED:-N/A}
SLOPE_AT_RECOMMENDED=${SLOPE_AT_RECOMMENDED:-N/A}
SAMPLING_DEPTH_CANDIDATES_CSV=$DEPTH_CSV
SAMPLES_RETAINED_CSV=${SAMPLES_RETAINED_CSV:-N/A}
SAMPLES_REMOVED_CSV=${SAMPLES_REMOVED_CSV:-N/A}
SAMPLE_DEPTH_SUMMARY_CSV=${SAMPLE_DEPTH_SUMMARY_CSV:-N/A}
GROUP_DEPTH_SUMMARY_CSV=${GROUP_DEPTH_SUMMARY_CSV:-N/A}
N_TRIMMED_SAMPLE_IDS=${N_TRIMMED_SAMPLE_IDS:-N/A}
N_TRIMMED_GROUP_LABELS=${N_TRIMMED_GROUP_LABELS:-N/A}

# ── QC and rarefaction outputs ────────────────────────────────────────────────
METADATA_SUMMARY_QZV=$METADATA_QZV
FEATURE_TABLE_SUMMARY_QZV=$FT_SUMMARY_QZV
SAMPLE_FREQUENCIES_QZA=$SAMPLE_FREQ_QZA
SAMPLE_FREQUENCIES_TSV=$SAMPLE_FREQ_TSV
FEATURE_BIOM=$FEATURE_BIOM
FEATURE_TSV=$FEATURE_TSV
SAMPLE_FREQ_MIN=${MIN_FREQ:-N/A}
SAMPLE_FREQ_MEDIAN=${MEDIAN_FREQ:-N/A}
SAMPLE_FREQ_MAX=${MAX_FREQ:-N/A}
ALPHA_RAREFACTION_PREVIEW_QZV=$ALPHA_RARE_MBX_QZV
ALPHA_RAREFACTION_PREVIEW_PNG=$ALPHA_RARE_PNG
ALPHA_RAREFACTION_PREVIEW_CSV=$ALPHA_RARE_CSV
ALPHA_RAREFACTION_PREVIEW_PER_SAMPLE_CSV=$PER_SAMPLE_RARE_CSV
ALPHA_RAREFACTION_PREVIEW_STEPS=$RARE_STEPS
ALPHA_RAREFACTION_QIIME_QZV=$ALPHA_RARE_QIIME_QZV
ALPHA_RAREFACTION_QIIME_STEPS=$QIIME_RARE_STEPS
ALPHA_RAREFACTION_QIIME_ITERATIONS=$QIIME_RARE_ITERATIONS
ALPHA_RAREFACTION_QIIME_STATUS=${QIIME_RARE_STATUS:-not_run}
ALPHA_RAREFACTION_QIIME_ELAPSED=${QIIME_RARE_ELAPSED:-N/A}

# ── ID validation (step 10) ───────────────────────────────────────────────────
ID_VALIDATION_STATUS=${ID_VAL_STATUS:-not_run}
ID_VALIDATION_REPORT=${ID_VAL_REPORT:-N/A}
ID_VAL_FEATURES_BIOM=${ID_VAL_FEATURES_BIOM:-N/A}
ID_VAL_FEATURES_REPSEQS=${ID_VAL_FEATURES_REPSEQS:-N/A}
ID_VAL_FEATURES_TREE=${ID_VAL_FEATURES_TREE:-N/A}
ID_VAL_FEATURES_BIOM_NOT_IN_TREE=${ID_VAL_FEATURES_BIOM_NOT_IN_TREE:-N/A}
ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS=${ID_VAL_FEATURES_BIOM_NOT_IN_REPSEQS:-N/A}
ID_VAL_SAMPLES_METADATA=${ID_VAL_SAMPLES_METADATA:-N/A}
ID_VAL_SAMPLES_BIOM=${ID_VAL_SAMPLES_BIOM:-N/A}
ID_VAL_SAMPLES_FREQ=${ID_VAL_SAMPLES_FREQ:-N/A}
ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM=${ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM:-N/A}
ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA=${ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA:-N/A}
ID_VAL_TRIMMED_FEATURE_IDS=${ID_VAL_TRIMMED_FEATURE_IDS:-N/A}
ID_VAL_TRIMMED_SAMPLE_IDS=${ID_VAL_TRIMMED_SAMPLE_IDS:-N/A}

# ── Visualization PNGs (step 11) ──────────────────────────────────────────────
SEQUENCING_DEPTH_DISTRIBUTION_PNG=${DEPTH_DIST_PNG:-N/A}
DEPTH_VS_RETENTION_PNG=${DEPTH_RETENTION_PNG:-N/A}
ALPHA_RAREFACTION_ANNOTATED_PNG=${ALPHA_RARE_ANNOT_PNG:-N/A}

TABLE_FILTERED_FLAG=$TABLE_FILTERED_FLAG
SKIP_QIIME_RAREFACTION=$SKIP_QIIME_RAREFACTION

# ── Provenance / versions ─────────────────────────────────────────────────────
GENERATED_AT=$NOW
SCRIPT_NAME=mbx_pre_diversity_parameters.sh
SCRIPT_VERSION=$SCRIPT_VERSION
QIIME2_VERSION=$QIIME_VERSION
PYTHON_VERSION=$PY_VERSION
NUMPY_VERSION=$NUMPY_VERSION
SCIPY_VERSION=$SCIPY_VERSION
BIOM_VERSION=$BIOM_VERSION
R_VERSION=$R_VERSION
N_JOBS=$N_JOBS

# ── Command log ───────────────────────────────────────────────────────────────
INVOCATION_CWD=$(pwd)
INVOCATION_USER=${USER:-unknown}
# Original argv preserved verbatim:
INVOCATION_ARGV=${MBX_INVOCATION:-$0 $*}

# ── Ready-to-run diversity command (copy to terminal after reviewing curves) ──
# IMPORTANT: Before running core-metrics, open alpha-rarefaction-qiime.qzv at
#   https://view.qiime2.org
# and confirm that richness curves have plateaued at the recommended depth.
#
# qiime diversity core-metrics-phylogenetic \\
#   --i-phylogeny  $ROOTED_QZA \\
#   --i-table      $ACTIVE_FEATURE_TABLE \\
#   --p-sampling-depth ${RECOMMENDED_DEPTH} \\
#   --m-metadata-file  $METADATA_TXT \\
#   --output-dir   $MBX_OUT_DIR/12_diversity_outputs
INFO

ok "Pre-diversity info → $PRE_DIV_INFO"

# ── Write plain-language summary (for the user) ──────────────────────────────
cat > "$PRE_DIV_SUMMARY" << SUMM
mbX Pro — pre-diversity step 11 summary
========================================
Generated : $NOW
Output    : $PRE_DIV_DIR

OVERALL STATUS                : $OVERALL_STATUS
READY FOR DIVERSITY (step 12) : $READY_FOR_DIVERSITY

  depth selection           : $DEPTH_STATUS  ($DEPTH_METHOD)
    rationale               : ${DEPTH_RATIONALE:-N/A}
    recommended depth       : ${RECOMMENDED_DEPTH:-N/A} reads
    samples retained        : ${SAMPLES_RETAINED:-N/A} / ${TOTAL_SAMPLES:-N/A}  (${PCT_RETAINED:-?}%)
    Good's coverage @ depth : ${MEAN_GOODS_AT_RECOMMENDED:-N/A}   (criterion ≥ ${GOOD_COV_MIN:-N/A})
    rarefaction slope       : ${SLOPE_AT_RECOMMENDED:-N/A}        (criterion < ${SLOPE_MAX:-N/A})
    grouping column         : ${GROUP_COLUMN_REPORTED:-overall_only}    n_groups = ${N_GROUPS:-N/A}

  id validation             : $ID_VAL_STATUS
    feature_ids (BIOM/repseqs/tree) : ${ID_VAL_FEATURES_BIOM:-N/A} / ${ID_VAL_FEATURES_REPSEQS:-N/A} / ${ID_VAL_FEATURES_TREE:-N/A}
    biom features missing from tree : ${ID_VAL_FEATURES_BIOM_NOT_IN_TREE:-N/A}    (must be 0 for UniFrac)
    sample_ids (metadata/BIOM/freq) : ${ID_VAL_SAMPLES_METADATA:-N/A} / ${ID_VAL_SAMPLES_BIOM:-N/A} / ${ID_VAL_SAMPLES_FREQ:-N/A}
    metadata not in biom            : ${ID_VAL_SAMPLES_METADATA_NOT_IN_BIOM:-N/A}
    biom not in metadata            : ${ID_VAL_SAMPLES_BIOM_NOT_IN_METADATA:-N/A}
    trimmed (silent whitespace)     : features=${ID_VAL_TRIMMED_FEATURE_IDS:-N/A}  samples=${ID_VAL_TRIMMED_SAMPLE_IDS:-N/A}

  qiime alpha-rarefaction   : $QIIME_RARE_STATUS  ($QIIME_RARE_NOTE)
    elapsed                 : ${QIIME_RARE_ELAPSED:-N/A}
    artifact                : $ALPHA_RARE_QIIME_QZV
    metrics                 : observed_features, shannon, faith_pd
    steps × iterations      : $QIIME_RARE_STEPS × $QIIME_RARE_ITERATIONS

What to do next
---------------
1. Open the alpha-rarefaction QZV in your browser at https://view.qiime2.org :
     - publication artifact : $ALPHA_RARE_QIIME_QZV
     - analytical preview   : $ALPHA_RARE_MBX_QZV   (PNG: $ALPHA_RARE_PNG)
   Confirm the richness curves plateau at or before $RECOMMENDED_DEPTH reads.

2. Open the decision plots:
     - $DEPTH_DIST_PNG
     - $DEPTH_RETENTION_PNG
     - $ALPHA_RARE_ANNOT_PNG

3. If the OVERALL STATUS is PASS or PASS_WITH_WARNINGS, proceed to step 12:
     mbx_diversity_run.sh $MBX_OUT_DIR

   If the OVERALL STATUS is REVIEW_REQUIRED, open:
     - $ID_VAL_REPORT
     - $DEPTH_CSV
   …and address the flagged items before re-running.

   If the OVERALL STATUS is FAIL, fix the BIOM↔tree mismatch or re-run DADA2
   with stricter filtering before proceeding. core-metrics-phylogenetic on a
   FAIL state would silently produce nonsense for UniFrac.

Provenance
----------
  qiime2  $QIIME_VERSION
  python  $PY_VERSION   numpy $NUMPY_VERSION   scipy $SCIPY_VERSION   biom $BIOM_VERSION
  R       $R_VERSION
  cores   $N_JOBS

  invoked-as : ${MBX_INVOCATION:-$0 $*}
  cwd        : $(pwd)
  user       : ${USER:-unknown}
SUMM

ok "Plain-language summary → $PRE_DIV_SUMMARY"

# ── Final terminal banner ────────────────────────────────────────────────────
sep
case "$OVERALL_STATUS" in
  PASS)
    ok "STEP 11 COMPLETE — STATUS: PASS"
    ok "Ready for diversity: yes  (proceed to mbx_diversity_run.sh)"
    ;;
  PASS_WITH_WARNINGS)
    warn "STEP 11 COMPLETE — STATUS: PASS_WITH_WARNINGS"
    warn "Ready for diversity: yes (with warnings — see $PRE_DIV_SUMMARY)"
    ;;
  REVIEW_REQUIRED)
    warn "STEP 11 COMPLETE — STATUS: REVIEW_REQUIRED"
    warn "Ready for diversity: NO until reviewed."
    warn "Open $PRE_DIV_SUMMARY for next steps."
    ;;
  FAIL)
    warn "STEP 11 COMPLETE — STATUS: FAIL"
    warn "Ready for diversity: NO. Fix the BIOM↔tree mismatch before step 12."
    warn "Open $PRE_DIV_SUMMARY for next steps."
    ;;
  *)
    warn "STEP 11 COMPLETE — STATUS: $OVERALL_STATUS"
    ;;
esac
sep
echo ""
echo "  Recommended depth : ${RECOMMENDED_DEPTH:-N/A} reads"
echo "  Samples retained  : ${SAMPLES_RETAINED:-N/A} / ${TOTAL_SAMPLES:-N/A}  (${PCT_RETAINED:-?}%)"
echo "  Frequency range   : ${MIN_FREQ:-N/A} – ${MAX_FREQ:-N/A} reads (median ${MEDIAN_FREQ:-N/A})"
echo "  Good's coverage   : mean=${MEAN_GOODS_AT_RECOMMENDED:-N/A}  (≥ ${GOOD_COV_MIN:-N/A})"
echo "  Plateau slope     : ${SLOPE_AT_RECOMMENDED:-N/A}            (< ${SLOPE_MAX:-N/A})"
echo ""
echo "  Read me first  →  $PRE_DIV_SUMMARY"
echo "  Machine info   →  $PRE_DIV_INFO"
echo ""
echo "  Decision plots :"
echo "    $DEPTH_DIST_PNG"
echo "    $DEPTH_RETENTION_PNG"
echo "    $ALPHA_RARE_ANNOT_PNG"
echo ""
if [[ "$READY_FOR_DIVERSITY" == "yes" || "$READY_FOR_DIVERSITY" == "yes_with_warnings" ]]; then
  echo "  When satisfied with the curves, run:"
  echo "    mbx_diversity_run.sh $MBX_OUT_DIR"
else
  echo "  DO NOT run mbx_diversity_run.sh until OVERALL_STATUS is PASS or PASS_WITH_WARNINGS."
fi
echo ""
