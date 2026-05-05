#!/usr/bin/env bash
# =============================================================================
#  mbx_alpha_diversity_run.sh   —   Step 12 of the mbX Pro pipeline
# =============================================================================
#  PURPOSE
#  -------
#  Compute alpha diversity (5 metrics), build a single tidy
#  alpha_diversity.xlsx, then run group statistics (Kruskal-Wallis +
#  pairwise Dunn + CLD) and boxplots for every categorical metadata column.
#
#  All inputs are read from the previous step's machine-readable info file:
#      <mbX_pro_outputs_*>/11_pre_diversity/mbx_pre_diversity_info.txt
#  The user provides only the mbX_pro_outputs directory path.
#
#  WHAT THIS SCRIPT PRODUCES
#  -------------------------
#  <mbX_pro_outputs_*>/
#  └── 12_alpha_diversity_results/
#      ├── all_alpha_outputs/                       ← raw QIIME2 artifacts
#      │   ├── rarefied_table.qza
#      │   ├── observed_features_vector.qza
#      │   ├── shannon_vector.qza
#      │   ├── simpson_vector.qza
#      │   ├── pielou_evenness_vector.qza
#      │   ├── faith_pd_vector.qza
#      │   └── exported/
#      │       ├── observed_features/alpha-diversity.tsv
#      │       ├── shannon/alpha-diversity.tsv
#      │       ├── simpson/alpha-diversity.tsv
#      │       ├── pielou_evenness/alpha-diversity.tsv
#      │       └── faith_pd/alpha-diversity.tsv
#      ├── alpha_diversity.xlsx                     ← THE main table
#      │      Columns: sample-id, ASVs_or_Features, Shannon_Index,
#      │               Simpson_Diversity, Faith_Phylogenetic_Diversity,
#      │               Pielou_Evenness  (+ all metadata columns)
#      ├── stats_for_alpha_diversity/
#      │   └── <Variable>/
#      │       ├── KW_<metric>_by_<Variable>.xlsx
#      │       ├── Pairwise_<metric>_by_<Variable>.xlsx
#      │       ├── CLD_Summary_<metric>_by_<Variable>.xlsx
#      │       └── Summary_KW_all_metrics_by_<Variable>.xlsx
#      ├── boxplots_for_alpha_diversity/
#      │   └── <Variable>/
#      │       ├── boxplot_<metric>_by_<Variable>.png
#      │       └── boxplot_panel_all_metrics_by_<Variable>.png
#      └── mbx_alpha_diversity_info.txt
#
#  GATING
#  ------
#  This script will refuse to run if step 11 reported
#       OVERALL_STATUS=FAIL or READY_FOR_DIVERSITY=no
#  unless the user passes --force.  PASS_WITH_WARNINGS is allowed.
#
#  Compatible with bash 3.2+ (macOS default shell).
#  Requires:  QIIME2 conda env + Rscript (system-wide, outside conda).
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
warn()    { echo "[WARN]  $*" >&2; }
info()    { echo "[INFO]  $*"; }
ok()      { echo "[OK]    $*"; }
skipped() { echo "[SKIP]  $*"; }
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

# Sanitize a string for use as a directory / filename:
#   spaces → underscores, strip parens/brackets/slashes/etc.
_sanitize() {
  printf '%s' "$1" \
    | tr ' ' '_' \
    | tr -d '()[]{}/<>|\\:*?"' \
    | sed 's/^[_.-]*//;s/[_.-]*$//'
}

# Pretty-print a QIIME2 command before running it
cmd_show() {
  echo ""
  echo "  \$ $1"
  shift
  for arg in "$@"; do
    printf '    %s\n' "$arg"
  done
  echo ""
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_alpha_diversity_run.sh — Compute alpha diversity, run stats & boxplots
                              for every categorical metadata variable.

USAGE:
  mbx_alpha_diversity_run.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads <mbX_pro_outputs_dir>/11_pre_diversity/mbx_pre_diversity_info.txt
  to discover:
      • Filtered feature table (mito/chloro removed)
      • Rooted phylogenetic tree
      • Sampling depth recommended by step 11
      • Metadata file path
  Then:
      1. Rarefies the feature table at the recommended depth.
      2. Computes 5 alpha diversity metrics:
            ASVs_or_Features              (observed_features)
            Shannon_Index                 (shannon)
            Simpson_Diversity             (simpson;  1 - D)
            Pielou_Evenness               (pielou_e)
            Faith_Phylogenetic_Diversity  (faith_pd, requires tree)
      3. Builds alpha_diversity.xlsx (one row per sample, all 5 metrics + metadata).
      4. Auto-detects categorical metadata columns.
      5. For every (categorical variable × metric) combination:
            • Kruskal-Wallis     → KW_<metric>_by_<var>.xlsx
            • Pairwise Dunn (BH) → Pairwise_<metric>_by_<var>.xlsx
            • Compact letters    → CLD_Summary_<metric>_by_<var>.xlsx
            • Boxplot with CLD   → boxplot_<metric>_by_<var>.png
      6. Writes one summary xlsx per variable (KW p-values for all 5 metrics).
      7. Writes one 5-panel "all metrics" PNG per variable.

OPTIONS:
  --depth N         Override sampling depth (default: read from step 11)
  --skip-stats      Skip statistics block (only build alpha_diversity.xlsx)
  --skip-boxplots   Skip boxplot generation
  --force           Run even if step 11 OVERALL_STATUS != PASS / PASS_WITH_WARNINGS
  --force-rerun     Recompute everything (ignore existing outputs)
  --dry-run         Print commands; do not execute
  -h, --help        Show this help

EXAMPLES:
  mbx_alpha_diversity_run.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_alpha_diversity_run.sh /path/to/mbX_pro_outputs_20260417_121431 --depth 25000

PROVENANCE:
  All paths are read from 11_pre_diversity/mbx_pre_diversity_info.txt.
  Re-runs are safe — completed steps are skipped automatically.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
USER_DEPTH=""
SKIP_STATS=false
SKIP_BOXPLOTS=false
FORCE_GATE=false
FORCE_RERUN=false
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)        usage ;;
    --depth)          USER_DEPTH="${2:-}"; [[ -z "$USER_DEPTH" ]] && err "--depth requires a value"; shift 2 ;;
    --skip-stats)     SKIP_STATS=true; shift ;;
    --skip-boxplots)  SKIP_BOXPLOTS=true; shift ;;
    --force)          FORCE_GATE=true; shift ;;
    --force-rerun)    FORCE_RERUN=true; shift ;;
    --dry-run)        DRY_RUN=true; shift ;;
    -*)               err "Unknown option: '$1'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '$1'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -d "$MBX_OUT_DIR" ]] || err "Directory does not exist: '$MBX_OUT_DIR'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# Validate sanity: the path should look like an mbX_pro_outputs_* dir
case "$(basename "$MBX_OUT_DIR")" in
  mbX_pro_outputs_*|mbx_pro_outputs_*) : ;;
  *) warn "Directory name does not match 'mbX_pro_outputs_*': $(basename "$MBX_OUT_DIR")
        Continuing anyway — but make sure this is the right output root." ;;
esac

# ── CPU detection ─────────────────────────────────────────────────────────────
if   command -v nproc   &>/dev/null; then N_JOBS="$(nproc)"
elif command -v sysctl  &>/dev/null; then N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else                                       N_JOBS=1
fi

# ── Locate a system-wide Rscript with the required packages ──────────────────
# We deliberately AVOID conda-env's Rscript per the project rule (see CLAUDE.md):
# conda's R has Python ABI conflicts with QIIME2's Python.
#
# CRITICAL: when this script runs inside a QIIME2 conda env, conda exports
#   R_LIBS_USER=<conda_env>/lib/R/library/
# This poisons system R: it looks for compiled extensions (e.g. Rcpp.so) in
# the conda R library path, where they don't exist (or are ABI-incompatible).
#
# We solve this by ALWAYS invoking system Rscript via `_R` wrapper, which
# strips every R_* env var so the system R uses only its own libPaths.
_R() {
  env -u R_HOME \
      -u R_LIBS \
      -u R_LIBS_USER \
      -u R_LIBS_SITE \
      -u R_PROFILE \
      -u R_PROFILE_USER \
      -u R_ENVIRON \
      -u R_ENVIRON_USER \
      -u R_PAPERSIZE \
      -u R_INCLUDE_DIR \
      -u R_DOC_DIR \
      -u R_SHARE_DIR \
      "$RSCRIPT_CMD" "$@"
}

# Pick the first Rscript whose libPaths can load openxlsx (the one used by the
# rest of the mbX pipeline, e.g. installed by mbx_ezclean_all_levels.sh).
RSCRIPT_CMD=""
RPROBE='if (!requireNamespace("openxlsx", quietly = TRUE)) quit(status=1); cat("ok")'
_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}
for _c in \
    "/opt/homebrew/bin/Rscript" \
    "/usr/local/bin/Rscript" \
    "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
    "/usr/bin/Rscript" \
    "$(command -v Rscript 2>/dev/null || true)"; do
  [[ -n "$_c" && -x "$_c" ]] || continue
  if _strip_env "$_c" -e "$RPROBE" >/dev/null 2>&1; then
    RSCRIPT_CMD="$_c"; break
  fi
done

# Fallback: pick any executable Rscript and let auto-install handle missing pkgs.
if [[ -z "$RSCRIPT_CMD" ]]; then
  for _c in \
      "/opt/homebrew/bin/Rscript" \
      "/usr/local/bin/Rscript" \
      "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
      "/usr/bin/Rscript" \
      "$(command -v Rscript 2>/dev/null || true)"; do
    [[ -n "$_c" && -x "$_c" ]] && { RSCRIPT_CMD="$_c"; break; }
  done
fi

[[ -n "$RSCRIPT_CMD" ]] || err "Rscript not found.
  → Install R (system-wide):  brew install r
  → Then run: mbx_ezclean_all_levels.sh <any output dir>
    (this installs all required R packages: mbX + openxlsx + ggplot2 + ...)"

R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"

# ── PID-based temp file IDs (avoids mktemp suffix issues on macOS bash 3.2) ───
_TMPID="${$}_$(date +%s)"
trap 'rm -f /tmp/mbx_alpha_${_TMPID}*.R /tmp/mbx_alpha_${_TMPID}*.txt' EXIT

# =============================================================================
step "1/9 — Read mbx_pre_diversity_info.txt"
# =============================================================================
PRE_DIV_DIR="${MBX_OUT_DIR}/11_pre_diversity"
PRE_DIV_INFO="${PRE_DIV_DIR}/mbx_pre_diversity_info.txt"

[[ -f "$PRE_DIV_INFO" ]] || err "Pre-diversity info file not found:
    $PRE_DIV_INFO
  → Run mbx_pre_diversity_parameters.sh first (step 11)."

OVERALL_STATUS="$(_read_key OVERALL_STATUS         "$PRE_DIV_INFO")"
READY="$(        _read_key READY_FOR_DIVERSITY    "$PRE_DIV_INFO")"
METADATA_TXT="$( _read_key METADATA_TXT           "$PRE_DIV_INFO")"
FT_QZA="$(       _read_key FEATURE_TABLE_QZA      "$PRE_DIV_INFO")"
TREE_QZA="$(     _read_key ROOTED_TREE_QZA        "$PRE_DIV_INFO")"
REC_DEPTH="$(    _read_key RECOMMENDED_DEPTH      "$PRE_DIV_INFO")"
SAMPLES_RETAINED="$(_read_key SAMPLES_RETAINED    "$PRE_DIV_INFO")"
TOTAL_SAMPLES="$(   _read_key TOTAL_SAMPLES       "$PRE_DIV_INFO")"
PCT_RETAINED="$(    _read_key PCT_RETAINED        "$PRE_DIV_INFO")"

# Allow user to override the depth
if [[ -n "$USER_DEPTH" ]]; then
  case "$USER_DEPTH" in
    ''|*[!0-9]*) err "--depth must be a positive integer (got: '$USER_DEPTH')" ;;
  esac
  warn "Overriding sampling depth: ${REC_DEPTH} → ${USER_DEPTH}"
  REC_DEPTH="$USER_DEPTH"
fi

# Sanity checks
[[ -n "$METADATA_TXT" && -f "$METADATA_TXT" ]] || err "METADATA_TXT not usable: '$METADATA_TXT'
  → Re-run step 11 (mbx_pre_diversity_parameters.sh)."
[[ -n "$FT_QZA"       && -f "$FT_QZA"       ]] || err "FEATURE_TABLE_QZA not usable: '$FT_QZA'"
[[ -n "$TREE_QZA"     && -f "$TREE_QZA"     ]] || err "ROOTED_TREE_QZA not usable: '$TREE_QZA'"
[[ -n "$REC_DEPTH" ]]                          || err "RECOMMENDED_DEPTH missing in pre-diversity info."
case "$REC_DEPTH" in
  ''|*[!0-9]*) err "RECOMMENDED_DEPTH is not a positive integer: '$REC_DEPTH'" ;;
esac

# ── Gating ────────────────────────────────────────────────────────────────────
case "$OVERALL_STATUS" in
  PASS|PASS_WITH_WARNINGS) : ;;
  *)
    if $FORCE_GATE; then
      warn "Step 11 OVERALL_STATUS = '$OVERALL_STATUS' — proceeding because --force was given."
    else
      err "Step 11 OVERALL_STATUS = '$OVERALL_STATUS'.
  → Diversity analysis on a non-PASS state may produce misleading results.
  → Read $PRE_DIV_DIR/pre_diversity_summary.txt and resolve the issues, OR
  → re-run with --force to bypass this gate."
    fi ;;
esac
case "$READY" in
  yes|yes_with_warnings) : ;;
  *)
    if ! $FORCE_GATE; then
      err "Step 11 READY_FOR_DIVERSITY = '$READY'.
  → Re-run with --force to bypass." ; fi ;;
esac

NOW="$(date '+%Y%m%d_%H%M%S')"
NOW_PRETTY="$(date '+%Y-%m-%d %H:%M:%S')"

info "Run timestamp        : $NOW"
info "Host                 : $(whoami)@$(hostname)"
info "QIIME2 env           : ${CONDA_DEFAULT_ENV:-unknown}"
info "Rscript              : $RSCRIPT_CMD  (R $R_VERSION)"
info "Output root          : $MBX_OUT_DIR"
info "Metadata             : $METADATA_TXT"
info "Feature table        : $FT_QZA"
info "Rooted tree          : $TREE_QZA"
info "Sampling depth       : $REC_DEPTH${USER_DEPTH:+ (user override)}"
info "Samples (step 11)    : ${SAMPLES_RETAINED:-?} / ${TOTAL_SAMPLES:-?}  (${PCT_RETAINED:-?}%)"
info "CPU cores            : $N_JOBS"
info "Step-11 status       : $OVERALL_STATUS  /  ready=$READY"
sep

# ── Set up output directories ─────────────────────────────────────────────────
ALPHA_DIR="${MBX_OUT_DIR}/12_alpha_diversity_results"
QZA_DIR="${ALPHA_DIR}/all_alpha_outputs"
EXPORT_DIR="${QZA_DIR}/exported"
STATS_DIR="${ALPHA_DIR}/stats_for_alpha_diversity"
BOXPLOTS_DIR="${ALPHA_DIR}/boxplots_for_alpha_diversity"
ALPHA_XLSX="${ALPHA_DIR}/alpha_diversity.xlsx"
INFO_FILE="${ALPHA_DIR}/mbx_alpha_diversity_info.txt"

mkdir -p "$ALPHA_DIR" "$QZA_DIR" "$EXPORT_DIR" "$STATS_DIR" "$BOXPLOTS_DIR"

# Metric tables (parallel arrays — bash 3.2 compatible)
QIIME_METRICS=( "observed_features" "shannon" "simpson" "pielou_e" "faith_pd" )
FRIENDLY=(     "ASVs_or_Features"   "Shannon_Index" "Simpson_Diversity" "Pielou_Evenness" "Faith_Phylogenetic_Diversity" )
# Output filenames for each .qza (matches QIIME2 convention)
QZA_NAMES=(    "observed_features_vector.qza" "shannon_vector.qza" "simpson_vector.qza" "pielou_evenness_vector.qza" "faith_pd_vector.qza" )

N_METRICS=${#QIIME_METRICS[@]}

# =============================================================================
step "2/9 — Rarefy feature table at depth ${REC_DEPTH}"
# =============================================================================
RARE_QZA="${QZA_DIR}/rarefied_table.qza"

if [[ -f "$RARE_QZA" && "$FORCE_RERUN" == false ]]; then
  skipped "Rarefied table already exists: $(basename "$RARE_QZA")"
else
  cmd_show "qiime feature-table rarefy" \
    "--i-table $FT_QZA" \
    "--p-sampling-depth $REC_DEPTH" \
    "--o-rarefied-table $RARE_QZA"

  if ! $DRY_RUN; then
    timer_start
    qiime feature-table rarefy \
      --i-table             "$FT_QZA" \
      --p-sampling-depth    "$REC_DEPTH" \
      --o-rarefied-table    "$RARE_QZA" \
      || err "qiime feature-table rarefy failed.
  → Verify the depth ($REC_DEPTH) is reasonable for your data.
  → Min sample count from step 11: $(_read_key MIN_SAMPLE_COUNT "$PRE_DIV_INFO")"
    timer_end
  fi
  ok "Rarefied table → $RARE_QZA"
fi
sep

# =============================================================================
step "3/9 — Compute 5 alpha diversity metrics (QIIME2)"
# =============================================================================
i=0
while [[ $i -lt $N_METRICS ]]; do
  METRIC="${QIIME_METRICS[$i]}"
  FRIENDLY_NAME="${FRIENDLY[$i]}"
  OUT_QZA="${QZA_DIR}/${QZA_NAMES[$i]}"

  if [[ -f "$OUT_QZA" && "$FORCE_RERUN" == false ]]; then
    skipped "[$((i+1))/$N_METRICS] $METRIC ($FRIENDLY_NAME) — already computed"
    i=$(( i + 1 )); continue
  fi

  info "[$((i+1))/$N_METRICS] Computing $METRIC ($FRIENDLY_NAME)..."

  if [[ "$METRIC" == "faith_pd" ]]; then
    cmd_show "qiime diversity alpha-phylogenetic" \
      "--i-table $RARE_QZA" \
      "--i-phylogeny $TREE_QZA" \
      "--p-metric faith_pd" \
      "--o-alpha-diversity $OUT_QZA"

    if ! $DRY_RUN; then
      timer_start
      qiime diversity alpha-phylogenetic \
        --i-table          "$RARE_QZA" \
        --i-phylogeny      "$TREE_QZA" \
        --p-metric         faith_pd \
        --o-alpha-diversity "$OUT_QZA" \
        || err "qiime diversity alpha-phylogenetic faith_pd failed.
  → Make sure the rooted tree contains every feature in the rarefied table.
  → If feature IDs were trimmed at step 11, re-run step 11 to regenerate the tree."
      timer_end
    fi
  else
    cmd_show "qiime diversity alpha" \
      "--i-table $RARE_QZA" \
      "--p-metric $METRIC" \
      "--o-alpha-diversity $OUT_QZA"

    if ! $DRY_RUN; then
      timer_start
      qiime diversity alpha \
        --i-table          "$RARE_QZA" \
        --p-metric         "$METRIC" \
        --o-alpha-diversity "$OUT_QZA" \
        || err "qiime diversity alpha (--p-metric $METRIC) failed."
      timer_end
    fi
  fi
  ok "  → $(basename "$OUT_QZA")"
  i=$(( i + 1 ))
done
sep

# =============================================================================
step "4/9 — Export alpha vectors to TSV"
# =============================================================================
i=0
while [[ $i -lt $N_METRICS ]]; do
  METRIC="${QIIME_METRICS[$i]}"
  IN_QZA="${QZA_DIR}/${QZA_NAMES[$i]}"
  OUT_SUBDIR="${EXPORT_DIR}/${METRIC}"
  TSV_FILE="${OUT_SUBDIR}/alpha-diversity.tsv"

  if [[ -f "$TSV_FILE" && "$FORCE_RERUN" == false ]]; then
    skipped "TSV already exported: $METRIC"
    i=$(( i + 1 )); continue
  fi

  info "Exporting $METRIC..."
  rm -rf "$OUT_SUBDIR"
  mkdir -p "$OUT_SUBDIR"

  if ! $DRY_RUN; then
    qiime tools export \
      --input-path  "$IN_QZA" \
      --output-path "$OUT_SUBDIR" \
      || err "Failed to export $IN_QZA"
    [[ -f "$TSV_FILE" ]] || err "Expected TSV not found after export: $TSV_FILE"
  fi
  i=$(( i + 1 ))
done
ok "All 5 metrics exported to: $EXPORT_DIR"
sep

# =============================================================================
step "5/9 — Build merged alpha_diversity.xlsx"
# =============================================================================
if [[ -f "$ALPHA_XLSX" && "$FORCE_RERUN" == false ]]; then
  skipped "alpha_diversity.xlsx already exists (use --force-rerun to regenerate)"
else
  R_BUILD="/tmp/mbx_alpha_${_TMPID}_build.R"

  cat > "$R_BUILD" << 'RBUILD'
suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE))
    install.packages("openxlsx", repos = "https://cloud.r-project.org", quiet = TRUE)
  library(openxlsx)
})

args        <- commandArgs(trailingOnly = TRUE)
metadata    <- args[1]
out_xlsx    <- args[2]
export_dir  <- args[3]

# QIIME2 metric folder name → friendly column name
mapping <- list(
  "observed_features" = "ASVs_or_Features",
  "shannon"           = "Shannon_Index",
  "simpson"           = "Simpson_Diversity",
  "pielou_e"          = "Pielou_Evenness",
  "faith_pd"          = "Faith_Phylogenetic_Diversity"
)

read_alpha_tsv <- function(metric_dir) {
  tsv <- file.path(export_dir, metric_dir, "alpha-diversity.tsv")
  if (!file.exists(tsv)) stop(sprintf("Missing TSV: %s", tsv))
  d <- read.delim(tsv, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, comment.char = "")
  if (ncol(d) < 2) stop(sprintf("Unexpected TSV format: %s", tsv))
  names(d)[1] <- "sample-id"
  d[, 2] <- suppressWarnings(as.numeric(d[, 2]))
  d <- d[, 1:2, drop = FALSE]
  names(d)[2] <- mapping[[metric_dir]]
  d
}

cat("[INFO]  Reading 5 alpha-diversity TSVs...\n")
metrics <- names(mapping)
tabs <- lapply(metrics, read_alpha_tsv)
names(tabs) <- metrics

# Merge by sample-id (full outer join)
alpha <- Reduce(function(a, b) merge(a, b, by = "sample-id", all = TRUE), tabs)

# Reorder columns:  sample-id, ASVs_or_Features, Shannon_Index,
# Simpson_Diversity, Faith_Phylogenetic_Diversity, Pielou_Evenness
desired_order <- c("sample-id",
                   "ASVs_or_Features",
                   "Shannon_Index",
                   "Simpson_Diversity",
                   "Faith_Phylogenetic_Diversity",
                   "Pielou_Evenness")
alpha <- alpha[, desired_order, drop = FALSE]
cat(sprintf("[INFO]  Alpha matrix: %d samples × %d metrics.\n",
            nrow(alpha), ncol(alpha) - 1))

# Load metadata, drop QIIME #q2:types row, attach for stats step
ext <- tolower(tools::file_ext(metadata))
md <- if (ext == "csv") {
  read.csv(metadata, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  read.delim(metadata, header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE,
             na.strings = c("", "NA", "N/A", "na", "n/a"))
}
q2_rows <- grepl("^#", md[[1]])
if (any(q2_rows)) {
  md <- md[!q2_rows, , drop = FALSE]
  cat("[INFO]  Dropped QIIME2 #q2:types row(s) from metadata.\n")
}
names(md)[1] <- "sample-id"
md[["sample-id"]] <- trimws(as.character(md[["sample-id"]]))
alpha[["sample-id"]] <- trimws(as.character(alpha[["sample-id"]]))

# Report missing/dropped samples
in_alpha_only <- setdiff(alpha[["sample-id"]], md[["sample-id"]])
in_meta_only  <- setdiff(md[["sample-id"]], alpha[["sample-id"]])
if (length(in_alpha_only) > 0)
  cat(sprintf("[WARN]  %d sample(s) in alpha-diversity but NOT in metadata: %s\n",
              length(in_alpha_only), paste(head(in_alpha_only, 5), collapse=", ")))
if (length(in_meta_only) > 0)
  cat(sprintf("[INFO]  %d sample(s) in metadata but NOT in alpha-diversity (rarefaction dropped): %s\n",
              length(in_meta_only), paste(head(in_meta_only, 5), collapse=", ")))

merged <- merge(alpha, md, by = "sample-id", all.x = TRUE, all.y = FALSE)

# Round metric columns to 4 decimals for readability
metric_cols <- c("ASVs_or_Features", "Shannon_Index", "Simpson_Diversity",
                 "Faith_Phylogenetic_Diversity", "Pielou_Evenness")
for (mc in metric_cols) {
  if (mc %in% names(merged))
    merged[[mc]] <- round(merged[[mc]], 4)
}

wb <- createWorkbook()
addWorksheet(wb, "alpha_diversity")
writeData(wb, "alpha_diversity", merged, withFilter = TRUE, headerStyle = createStyle(
  textDecoration = "bold", border = "Bottom"))
freezePane(wb, "alpha_diversity", firstActiveRow = 2, firstActiveCol = 2)
setColWidths(wb, "alpha_diversity",
             cols = 1:ncol(merged),
             widths = "auto")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat(sprintf("[OK]    Wrote: %s  (%d rows × %d cols)\n",
            out_xlsx, nrow(merged), ncol(merged)))
RBUILD

  if $DRY_RUN; then
    warn "[DRY-RUN] Would build alpha_diversity.xlsx via Rscript."
  else
    _R --vanilla "$R_BUILD" "$METADATA_TXT" "$ALPHA_XLSX" "$EXPORT_DIR" \
      || err "Failed to build alpha_diversity.xlsx — see R output above."
  fi
  rm -f "$R_BUILD"
fi
ok "alpha_diversity.xlsx → $ALPHA_XLSX"
sep

# =============================================================================
step "6/9 — Detect categorical metadata columns"
# =============================================================================
CATS_FILE="/tmp/mbx_alpha_${_TMPID}_cats.txt"
R_DETECT="/tmp/mbx_alpha_${_TMPID}_detect.R"

cat > "$R_DETECT" << 'RDETECT'
args     <- commandArgs(trailingOnly = TRUE)
metadata <- args[1]

ext <- tolower(tools::file_ext(metadata))
md <- if (ext == "csv") {
  read.csv(metadata, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  read.delim(metadata, header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE,
             na.strings = c("", "NA", "N/A", "na", "n/a"))
}
q2_rows <- grepl("^#", md[[1]])
if (any(q2_rows)) {
  md <- md[!q2_rows, , drop = FALSE]
  cat("[INFO]  Removed QIIME2 #q2:types row(s).\n", file = stderr())
}
n_samples <- nrow(md)
cat(sprintf("[INFO]  Metadata: %d samples, %d columns.\n",
            n_samples, ncol(md)), file = stderr())

# Skip first column (sample-id) and any column that LOOKS like an ID column
# (sampleid, farmid, runid, etc.)  ezviz/ezstat already skip the first column;
# we additionally skip any column whose values are ALL UNIQUE since that's an
# implicit ID.
col_names <- names(md)[-1]

categorical          <- character(0)
skipped_numeric      <- character(0)
skipped_alluniq      <- character(0)
skipped_noval        <- character(0)
skipped_mingroupsize <- character(0)

for (col in col_names) {
  vals_raw <- md[[col]]
  vals     <- vals_raw[!is.na(vals_raw) & nchar(trimws(as.character(vals_raw))) > 0]
  if (length(vals) == 0)                             { skipped_noval    <- c(skipped_noval, col);   next }
  num_test <- suppressWarnings(as.numeric(vals))
  if (!any(is.na(num_test)))                         { skipped_numeric  <- c(skipped_numeric, col); next }
  n_unique <- length(unique(trimws(as.character(vals))))
  if (n_unique <= 1)                                 { skipped_noval    <- c(skipped_noval, col);   next }
  if (n_unique == n_samples)                         { skipped_alluniq  <- c(skipped_alluniq, col); next }
  group_counts <- table(trimws(as.character(vals)))
  if (max(group_counts) < 2)                         { skipped_mingroupsize <- c(skipped_mingroupsize, col); next }
  categorical <- c(categorical, col)
}

if (length(skipped_numeric)      > 0) cat(sprintf("[INFO]  Skipped numeric       : %s\n", paste(skipped_numeric, collapse=", ")),      file=stderr())
if (length(skipped_alluniq)      > 0) cat(sprintf("[INFO]  Skipped all-unique    : %s\n", paste(skipped_alluniq, collapse=", ")),      file=stderr())
if (length(skipped_noval)        > 0) cat(sprintf("[INFO]  Skipped empty/const   : %s\n", paste(skipped_noval, collapse=", ")),        file=stderr())
if (length(skipped_mingroupsize) > 0) cat(sprintf("[INFO]  Skipped singleton grp : %s\n", paste(skipped_mingroupsize, collapse=", ")), file=stderr())

if (length(categorical) == 0) {
  cat("[ERROR] No categorical columns found in metadata.\n", file = stderr())
  quit(status = 1)
}

cat(sprintf("[INFO]  Categorical columns : %d → %s\n",
            length(categorical), paste(categorical, collapse=", ")), file=stderr())
cat(paste(categorical, collapse="\n"), "\n", sep="")
RDETECT

if $DRY_RUN; then
  warn "[DRY-RUN] Would detect categorical columns from $METADATA_TXT"
  CATEGORICAL_COLS=( "Treatment" "SampleType" )
else
  _R --vanilla "$R_DETECT" "$METADATA_TXT" > "$CATS_FILE" 2>&1 \
    || err "Categorical column detection failed.
  → Verify metadata file: $METADATA_TXT
  → It must have at least one non-numeric grouping column."
  grep "^\[" "$CATS_FILE" || true
  CATEGORICAL_COLS=()
  while IFS= read -r _col; do
    [[ -n "$_col" ]] && CATEGORICAL_COLS+=("$_col")
  done < <(grep -v "^\[" "$CATS_FILE" | grep -v '^$' || true)
  [[ ${#CATEGORICAL_COLS[@]} -eq 0 ]] && err "No categorical columns detected.
  → Verify your metadata has ≥1 grouping column with string values."
fi
rm -f "$R_DETECT"
ok "Categorical columns: ${CATEGORICAL_COLS[*]}"

# Build sanitized name list (parallel to CATEGORICAL_COLS)
SANITIZED_COLS=()
for col in "${CATEGORICAL_COLS[@]}"; do
  SANITIZED_COLS+=( "$(_sanitize "$col")" )
done
sep

# =============================================================================
step "7/9 — Statistics  (Kruskal-Wallis + pairwise Dunn + CLD)"
# =============================================================================
if $SKIP_STATS; then
  warn "Skipping statistics block (--skip-stats)."
else
  R_STATS="/tmp/mbx_alpha_${_TMPID}_stats.R"

  # Build the variables list as a comma-separated quoted string for R
  VARS_R=""
  for col in "${CATEGORICAL_COLS[@]}"; do
    _esc="${col//\"/\\\"}"
    if [[ -z "$VARS_R" ]]; then VARS_R="\"$_esc\""
    else                        VARS_R="$VARS_R, \"$_esc\""; fi
  done

  # Convert bash true/false → R TRUE/FALSE for heredoc injection
  if $SKIP_BOXPLOTS; then SKIP_BOXPLOTS_R="TRUE"; else SKIP_BOXPLOTS_R="FALSE"; fi

  cat > "$R_STATS" << RSTATS
# =============================================================================
# alpha-diversity stats + boxplots
# =============================================================================
suppressPackageStartupMessages({
  required_pkgs <- c("openxlsx", "ggplot2", "dunn.test", "multcompView")
  to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install) > 0) {
    cat(sprintf("[INFO]  Installing R package(s): %s\n", paste(to_install, collapse = ", ")))
    install.packages(to_install, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(openxlsx)
  library(ggplot2)
  library(dunn.test)
  library(multcompView)
})

ALPHA_XLSX     <- "${ALPHA_XLSX}"
STATS_DIR      <- "${STATS_DIR}"
BOXPLOTS_DIR   <- "${BOXPLOTS_DIR}"
SKIP_BOXPLOTS  <- ${SKIP_BOXPLOTS_R}
VARS           <- c(${VARS_R})

METRICS <- c("ASVs_or_Features",
             "Shannon_Index",
             "Simpson_Diversity",
             "Faith_Phylogenetic_Diversity",
             "Pielou_Evenness")

# Pretty labels for plots
METRIC_LABELS <- c(
  ASVs_or_Features              = "Observed ASVs / Features",
  Shannon_Index                 = "Shannon Index (H')",
  Simpson_Diversity             = "Simpson Diversity (1 - D)",
  Faith_Phylogenetic_Diversity  = "Faith's Phylogenetic Diversity",
  Pielou_Evenness               = "Pielou's Evenness (J')"
)

cat(sprintf("[INFO]  Reading alpha xlsx : %s\n", ALPHA_XLSX))
df <- read.xlsx(ALPHA_XLSX, sheet = "alpha_diversity", check.names = FALSE)
names(df)[1] <- "sample-id"
cat(sprintf("[INFO]  %d samples × %d columns loaded.\n", nrow(df), ncol(df)))

# ── Compute Kruskal-Wallis + pairwise Dunn + CLD ──────────────────────────────
# Returns a list with:  kw_df, pairwise_df, cld_df, plot (ggplot)
compute_one <- function(d, var, metric) {
  d <- d[!is.na(d[[var]]) & !is.na(d[[metric]]) &
         nchar(trimws(as.character(d[[var]]))) > 0, , drop = FALSE]
  d[[var]] <- factor(trimws(as.character(d[[var]])))
  if (length(unique(d[[var]])) < 2) {
    return(list(error = sprintf("Variable '%s' has <2 groups after dropping NAs.", var)))
  }

  # Drop singleton groups (n < 2) — KW + Dunn require ≥ 2 per group.
  group_n <- table(d[[var]])
  bad     <- names(group_n)[group_n < 2]
  if (length(bad) > 0) {
    d <- d[!d[[var]] %in% bad, , drop = FALSE]
    d[[var]] <- droplevels(d[[var]])
    if (length(unique(d[[var]])) < 2)
      return(list(error = sprintf(
        "After dropping singleton groups (%s), <2 groups remain.",
        paste(bad, collapse = ", "))))
  }
  small <- names(group_n)[group_n >= 2 & group_n < 3]

  formula_str <- sprintf("\`%s\` ~ \`%s\`", metric, var)
  kw <- kruskal.test(as.formula(formula_str), data = d)

  # Effect size: epsilon-squared  (Tomczak & Tomczak 2014)
  H <- as.numeric(kw\$statistic)
  N <- nrow(d)
  k <- length(unique(d[[var]]))
  eps_sq <- (H - k + 1) / (N - k)
  if (!is.finite(eps_sq)) eps_sq <- NA_real_
  eps_sq <- max(0, eps_sq)

  kw_df <- data.frame(
    metric            = metric,
    variable          = var,
    n_samples         = N,
    n_groups          = k,
    smallest_group_n  = min(group_n[group_n >= 2]),
    KW_chi_squared    = round(H, 4),
    KW_df             = as.integer(kw\$parameter),
    KW_p_value        = signif(kw\$p.value, 6),
    epsilon_squared   = round(eps_sq, 4),
    significant_at_05 = kw\$p.value < 0.05,
    note              = if (length(small) > 0)
                          sprintf("groups n<3 (results unstable): %s", paste(small, collapse=", "))
                        else "",
    stringsAsFactors  = FALSE
  )

  # ── Pairwise post-hoc test (BH-corrected for k>=3) ──────────────────────────
  groups_kept <- levels(d[[var]])
  pw_df <- NULL
  if (k == 2) {
    # KW with 2 groups is mathematically equivalent to Mann-Whitney U / Wilcoxon
    # rank-sum.  We use wilcox.test here because dunn.test has a known
    # 2-group indexing bug ("incorrect number of dimensions" via Psort[1,i]).
    wt <- tryCatch(
      wilcox.test(as.formula(formula_str), data = d, exact = FALSE),
      error = function(e) NULL)
    pw_df <- data.frame(
      comparison       = paste(groups_kept, collapse = " - "),
      test             = "Wilcoxon rank-sum (Mann-Whitney U)",
      statistic        = if (!is.null(wt)) round(as.numeric(wt\$statistic), 4) else NA_real_,
      p_unadjusted     = if (!is.null(wt)) signif(wt\$p.value, 6)               else NA_real_,
      p_BH_adjusted    = if (!is.null(wt)) signif(wt\$p.value, 6)               else NA_real_,
      significant_BH05 = if (!is.null(wt)) !is.na(wt\$p.value) && wt\$p.value < 0.05 else NA,
      stringsAsFactors = FALSE
    )
  } else if (k >= 3) {
    invisible(capture.output(
      dt <- dunn.test::dunn.test(x = d[[metric]], g = d[[var]],
                                 method = "bh",
                                 kw = FALSE, label = TRUE, table = FALSE,
                                 list = FALSE, alpha = 0.05)
    ))
    pw_df <- data.frame(
      comparison       = dt\$comparisons,
      test             = "Dunn (BH-adjusted)",
      statistic        = round(dt\$Z, 4),
      p_unadjusted     = signif(dt\$P, 6),
      p_BH_adjusted    = signif(dt\$P.adjusted, 6),
      significant_BH05 = dt\$P.adjusted < 0.05,
      stringsAsFactors = FALSE
    )
  }

  # ── Compact letter display (CLD) ────────────────────────────────────────────
  cld_letters <- setNames(rep("a", length(groups_kept)), groups_kept)
  kw_sig <- !is.na(kw\$p.value) && kw\$p.value < 0.05
  if (!is.null(pw_df) && nrow(pw_df) > 0 && kw_sig) {
    pmat_vec <- pw_df\$p_BH_adjusted
    names(pmat_vec) <- gsub("\\\\s*-\\\\s*", "-", pw_df\$comparison)
    cld_obj <- tryCatch(
      multcompView::multcompLetters(pmat_vec, threshold = 0.05),
      error = function(e) NULL)
    if (!is.null(cld_obj) && !is.null(cld_obj\$Letters))
      cld_letters[names(cld_obj\$Letters)] <- cld_obj\$Letters
  }

  cld_df <- data.frame(
    group   = groups_kept,
    n       = as.integer(table(d[[var]])[groups_kept]),
    median  = round(tapply(d[[metric]], d[[var]], median, na.rm = TRUE)[groups_kept], 4),
    mean    = round(tapply(d[[metric]], d[[var]], mean,   na.rm = TRUE)[groups_kept], 4),
    sd      = round(tapply(d[[metric]], d[[var]], sd,     na.rm = TRUE)[groups_kept], 4),
    IQR     = round(tapply(d[[metric]], d[[var]], IQR,    na.rm = TRUE)[groups_kept], 4),
    cld     = cld_letters[groups_kept],
    stringsAsFactors = FALSE
  )
  cld_df <- cld_df[order(-cld_df\$median), , drop = FALSE]

  # ── Boxplot (CLD letters above each box) ────────────────────────────────────
  ymax <- max(d[[metric]], na.rm = TRUE)
  ymin <- min(d[[metric]], na.rm = TRUE)
  yrange <- if (ymax > ymin) ymax - ymin else 1
  letter_y <- ymax + 0.06 * yrange

  cld_plot <- data.frame(group = factor(names(cld_letters), levels = groups_kept),
                         cld   = unname(cld_letters),
                         y     = letter_y, stringsAsFactors = FALSE)

  subtitle <- sprintf("KW p = %s   |   epsilon^2 = %s   |   n = %d in %d groups",
                      format.pval(kw\$p.value, digits = 3, eps = 1e-6),
                      ifelse(is.na(eps_sq), "NA", sprintf("%.3f", eps_sq)),
                      N, k)

  p <- ggplot(d, aes(x = .data[[var]], y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA, fill = "#cfe2f3", colour = "#1f2d3d",
                 width = 0.6, alpha = 0.85) +
    geom_jitter(width = 0.15, height = 0, size = 1.6,
                alpha = 0.7, colour = "#1f2d3d") +
    geom_text(data = cld_plot,
              aes(x = group, y = y, label = cld),
              inherit.aes = FALSE,
              size = 5, fontface = "bold", colour = "#1f2d3d") +
    labs(title    = sprintf("%s by %s", METRIC_LABELS[metric], var),
         subtitle = subtitle,
         x = var, y = METRIC_LABELS[metric]) +
    theme_classic(base_size = 13) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 11),
          axis.text.x   = element_text(angle = 25, hjust = 1))

  list(kw_df = kw_df, pw_df = pw_df, cld_df = cld_df, plot = p)
}

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^[._-]+|[._-]+\$", "", x)
  x
}

n_pass <- 0L
n_fail <- 0L
fail_log <- character(0)
overall_summary <- list()

for (var in VARS) {
  safe_var   <- sanitize(var)
  stat_subdir <- file.path(STATS_DIR,    safe_var)
  box_subdir  <- file.path(BOXPLOTS_DIR, safe_var)
  dir.create(stat_subdir, recursive = TRUE, showWarnings = FALSE)
  if (!SKIP_BOXPLOTS) dir.create(box_subdir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n  ── %s  →  %s/  &  %s/\n", var, stat_subdir, box_subdir))

  per_var_summary <- data.frame()
  panel_plots <- list()

  for (metric in METRICS) {
    if (!(metric %in% names(df))) {
      cat(sprintf("[WARN]  Metric '%s' not in alpha_diversity.xlsx — skipping.\n", metric))
      next
    }

    res <- tryCatch(compute_one(df, var, metric),
                    error = function(e) list(error = conditionMessage(e)))

    if (!is.null(res\$error)) {
      cat(sprintf("[FAIL]  %s × %s   →   %s\n", var, metric, res\$error))
      n_fail <- n_fail + 1L
      fail_log <- c(fail_log, sprintf("%s × %s : %s", var, metric, res\$error))
      next
    }

    # Write 3 xlsx files
    kw_xlsx <- file.path(stat_subdir,
      sprintf("KW_%s_by_%s.xlsx", metric, safe_var))
    pw_xlsx <- file.path(stat_subdir,
      sprintf("Pairwise_%s_by_%s.xlsx", metric, safe_var))
    cld_xlsx <- file.path(stat_subdir,
      sprintf("CLD_Summary_%s_by_%s.xlsx", metric, safe_var))

    write.xlsx(res\$kw_df,  kw_xlsx,  overwrite = TRUE,
               headerStyle = createStyle(textDecoration = "bold"))
    if (!is.null(res\$pw_df))
      write.xlsx(res\$pw_df, pw_xlsx, overwrite = TRUE,
                 headerStyle = createStyle(textDecoration = "bold"))
    write.xlsx(res\$cld_df, cld_xlsx, overwrite = TRUE,
               headerStyle = createStyle(textDecoration = "bold"))

    cat(sprintf("[OK]    %s × %s   p=%-10s eps2=%-6s\n",
                var, metric,
                format.pval(res\$kw_df\$KW_p_value, digits = 3, eps = 1e-6),
                ifelse(is.na(res\$kw_df\$epsilon_squared), "NA",
                       sprintf("%.3f", res\$kw_df\$epsilon_squared))))

    # Boxplot
    if (!SKIP_BOXPLOTS) {
      png_path <- file.path(box_subdir,
        sprintf("boxplot_%s_by_%s.png", metric, safe_var))
      ggsave(png_path, res\$plot,
             width = 7, height = 5, dpi = 300, bg = "white")
      panel_plots[[metric]] <- res\$plot
    }

    per_var_summary <- rbind(per_var_summary, res\$kw_df)
    n_pass <- n_pass + 1L
  }

  # Per-variable summary xlsx (one row per metric)
  if (nrow(per_var_summary) > 0) {
    summary_xlsx <- file.path(stat_subdir,
      sprintf("Summary_KW_all_metrics_by_%s.xlsx", safe_var))
    write.xlsx(per_var_summary, summary_xlsx, overwrite = TRUE,
               headerStyle = createStyle(textDecoration = "bold"))
    cat(sprintf("[OK]    Variable summary  → %s\n", basename(summary_xlsx)))
    overall_summary[[var]] <- per_var_summary
  }

  # 5-panel "all metrics" PNG  (uses patchwork if available, else cowplot)
  if (!SKIP_BOXPLOTS && length(panel_plots) > 0) {
    panel_path <- file.path(box_subdir,
      sprintf("boxplot_panel_all_metrics_by_%s.png", safe_var))

    pw_ok <- requireNamespace("patchwork", quietly = TRUE)
    if (!pw_ok) {
      try(install.packages("patchwork", repos = "https://cloud.r-project.org",
                           quiet = TRUE), silent = TRUE)
      pw_ok <- requireNamespace("patchwork", quietly = TRUE)
    }
    if (pw_ok) {
      library(patchwork)
      combined <- Reduce(\`+\`, panel_plots) +
        plot_layout(ncol = 2) +
        plot_annotation(title = sprintf("Alpha diversity by %s", var))
      ggsave(panel_path, combined,
             width = 14, height = 10, dpi = 300, bg = "white")
      cat(sprintf("[OK]    Panel PNG (5-metric)  → %s\n", basename(panel_path)))
    } else {
      cat(sprintf("[WARN]  'patchwork' unavailable — skipped panel PNG for %s\n", var))
    }
  }
}

# ── Global summary across ALL variables (one matrix-style xlsx) ──────────────
if (length(overall_summary) > 0) {
  global <- do.call(rbind, overall_summary)
  rownames(global) <- NULL
  global_xlsx <- file.path(STATS_DIR, "Summary_KW_all_variables_all_metrics.xlsx")
  write.xlsx(global, global_xlsx, overwrite = TRUE,
             headerStyle = createStyle(textDecoration = "bold"))
  cat(sprintf("\n[OK]    Global KW summary  → %s\n", global_xlsx))
}

cat(sprintf("\n[SUMMARY]  PASSED=%d   FAILED=%d\n", n_pass, n_fail))
if (n_fail > 0) {
  cat("[SUMMARY]  Failures:\n")
  for (f in fail_log) cat(sprintf("  - %s\n", f))
}
RSTATS

  if $DRY_RUN; then
    warn "[DRY-RUN] Would run stats + boxplots via Rscript."
  else
    _R --vanilla "$R_STATS" \
      || err "Stats / boxplots R script failed — see output above."
  fi
  rm -f "$R_STATS"
fi
sep

# =============================================================================
step "8/9 — Boxplots status"
# =============================================================================
# Boxplots are produced inside the same R script as stats (above); this step
# exists only to give a clean section break + status print.
if $SKIP_BOXPLOTS; then
  warn "Boxplots were SKIPPED (--skip-boxplots)."
else
  ok "Boxplots written to: $BOXPLOTS_DIR/<Variable>/"
  ok "Per-variable 5-panel PNGs available alongside individual metric PNGs."
fi
sep

# =============================================================================
step "9/9 — Write mbx_alpha_diversity_info.txt"
# =============================================================================
QIIME_VER="$(qiime --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo unknown)"
PY_BIN="$(command -v python || true)"
PY_VER="$($PY_BIN -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>/dev/null || echo unknown)"

# Build VARIABLE= lines for the info file
VAR_LINES=""
for col in "${CATEGORICAL_COLS[@]}"; do
  VAR_LINES="${VAR_LINES}VARIABLE=${col}
"
done

cat > "$INFO_FILE" << INFO
# ============================================================================
# mbx_alpha_diversity_info.txt
# Generated by mbx_alpha_diversity_run.sh
# Date : $NOW_PRETTY
# ============================================================================
#
# This file is read by the next pipeline step (beta diversity).
# Do NOT edit the key=value lines below — they are parsed programmatically.

# ── Inputs ───────────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
PRE_DIVERSITY_INFO=$PRE_DIV_INFO
METADATA_TXT=$METADATA_TXT
FEATURE_TABLE_QZA=$FT_QZA
ROOTED_TREE_QZA=$TREE_QZA
RECOMMENDED_DEPTH=$REC_DEPTH
USER_OVERRIDE_DEPTH=${USER_DEPTH:-no}

# ── Output locations ─────────────────────────────────────────────────────────
ALPHA_OUT_DIR=$ALPHA_DIR
ALPHA_QZA_DIR=$QZA_DIR
ALPHA_EXPORT_DIR=$EXPORT_DIR
ALPHA_DIVERSITY_XLSX=$ALPHA_XLSX
STATS_DIR=$STATS_DIR
BOXPLOTS_DIR=$BOXPLOTS_DIR

# ── Alpha vector files ───────────────────────────────────────────────────────
RAREFIED_TABLE_QZA=$RARE_QZA
OBSERVED_FEATURES_QZA=${QZA_DIR}/observed_features_vector.qza
SHANNON_QZA=${QZA_DIR}/shannon_vector.qza
SIMPSON_QZA=${QZA_DIR}/simpson_vector.qza
PIELOU_EVENNESS_QZA=${QZA_DIR}/pielou_evenness_vector.qza
FAITH_PD_QZA=${QZA_DIR}/faith_pd_vector.qza

# ── Categorical variables analysed ───────────────────────────────────────────
N_VARIABLES=${#CATEGORICAL_COLS[@]}
${VAR_LINES}

# ── Provenance ───────────────────────────────────────────────────────────────
GENERATED_AT=$NOW_PRETTY
SCRIPT_NAME=mbx_alpha_diversity_run.sh
QIIME2_VERSION=$QIIME_VER
PYTHON_VERSION=$PY_VER
R_VERSION=$R_VERSION
N_JOBS=$N_JOBS
SKIPPED_STATS=$SKIP_STATS
SKIPPED_BOXPLOTS=$SKIP_BOXPLOTS
INVOCATION_USER=${USER:-unknown}
INVOCATION_CWD=$(pwd)
INVOCATION_ARGV=$0 $*
INFO

ok "Info file → $INFO_FILE"
sep

# ── Final banner ──────────────────────────────────────────────────────────────
echo ""
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║  STEP 12 COMPLETE  —  ALPHA DIVERSITY                       ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "  Main table        : $ALPHA_XLSX"
echo "  QIIME2 artifacts  : $QZA_DIR/"
echo "  Stats             : $STATS_DIR/<Variable>/"
echo "  Boxplots          : $BOXPLOTS_DIR/<Variable>/"
echo "  Variables         : ${CATEGORICAL_COLS[*]}"
echo ""
echo "  Open the main table in Excel / Numbers:"
echo "    open '$ALPHA_XLSX'"
echo ""
echo "  Each variable folder contains:"
echo "    • KW_<metric>_by_<var>.xlsx         (Kruskal-Wallis chi-sq, p, epsilon^2)"
echo "    • Pairwise_<metric>_by_<var>.xlsx   (Dunn test, BH-corrected)"
echo "    • CLD_Summary_<metric>_by_<var>.xlsx (compact letter display)"
echo "    • Summary_KW_all_metrics_by_<var>.xlsx (all 5 metrics in one sheet)"
echo "    • boxplot_<metric>_by_<var>.png     (1 PNG per metric)"
echo "    • boxplot_panel_all_metrics_by_<var>.png  (5-panel)"
echo ""
echo "  Provenance (machine-readable) → $INFO_FILE"
echo ""
