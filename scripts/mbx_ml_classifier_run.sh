#!/usr/bin/env bash
# =============================================================================
#  mbx_ml_classifier_run.sh  (step 16)
#  Random Forest biomarker classifier — predict each categorical metadata
#  variable from the microbiome composition (all 7 taxonomic levels).
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  PURPOSE
#    For every (taxonomic_level x categorical_variable) combination, this
#    script trains a Random Forest classifier and reports:
#
#      - Cross-validated accuracy + AUC + F1 + sensitivity/specificity
#      - Confusion matrix (heatmap)
#      - ROC curves (one-vs-rest for multi-class)
#      - Permutation feature importance (top-20 barplot + full xlsx)
#      - SHAP-style per-sample feature contributions (Saabas-equivalent)
#      - Per-sample predicted-vs-actual table
#      - Per-variable summary comparing accuracy/AUC across all 7 levels
#
#  WHY add Random Forest after ANCOMBC2?
#    ANCOMBC2 answers "which taxa differ between groups?" (descriptive).
#    Random Forest answers "can we PREDICT the group from the taxa?"
#    (predictive).  Reviewers commonly ask BOTH questions.  A taxon that
#    is statistically different but doesn't help the classifier may be a
#    false-positive; a taxon that helps the classifier but isn't flagged
#    by ANCOMBC2 may be a non-linear / interaction effect.
#
#  GATING (reads from previous steps' info files)
#    7_taxonomy_csv/mbx_taxonomy_info.txt    -> METADATA_TXT
#    8_cleaned_files/mbx_ezclean_info.txt    -> 7 cleaned xlsx paths
#
#  OUTPUT STRUCTURE
#    16_ml_biomarkers/
#      working_dir_ml/                 (intermediate fold-level files)
#      <variable>/                     (one dir per categorical column)
#         RF_<level>_by_<variable>/
#            model_metrics.xlsx
#            confusion_matrix.{png,pdf}
#            roc_curves.{png,pdf}
#            feature_importance.xlsx
#            top20_importance.{png,pdf}
#            shap_per_sample.{png,pdf}
#            predicted_vs_actual.xlsx
#            model.rds
#         Summary_RF_<variable>.xlsx       (one row per level)
#      mbx_ml_biomarkers_info.txt
#
#  RUNTIME
#    ~5-30 s per (level x variable) combination on N=20 samples.
#    Total: 7 levels x N_categorical_variables ~ 2-15 minutes.
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

# Sanitize a string for safe use as a directory name
_sanitize_dirname() {
  printf '%s' "$1" \
    | tr ' ' '_' \
    | tr -d '()[]{}/<>|\\:*?"' \
    | sed 's/^[_.-]*//;s/[_.-]*$//'
}

# ── R env wrapper — strip conda's R_LIBS_USER pollution ──────────────────────
# Same pattern as steps 12-15.  When run from inside a conda env, conda
# exports R_LIBS_USER / R_HOME pointing at the conda env's R library tree.
# System Rscript (Homebrew/CRAN) then tries to load packages from those paths
# and fails with errors like:
#     Error: shared object 'Rcpp.so' not found
#     Error in library.dynam: ... 'methods.dylib' not found
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}
_R() { _strip_env "$RSCRIPT_CMD" "$@"; }

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_ml_classifier_run.sh — Random Forest biomarker classifier (step 16)

USAGE:
  mbx_ml_classifier_run.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Trains a Random Forest classifier per (taxonomic_level x categorical_variable)
  combination using the cleaned tables from step 8 and the metadata from
  step 7.  Produces cross-validated metrics, ROC curves, permutation feature
  importance, and SHAP-style per-sample contributions for every combination.

OPTIONS:
  --levels <list>       Comma-separated subset of {d,p,c,o,f,g,s}
                        (default: d,p,c,o,f,g,s — all 7 levels)
  --num-trees <N>       Random Forest trees per model (default: 500)
  --seed <N>            Random seed for reproducibility (default: 42)
  --skip-shap           Skip SHAP-style per-sample contribution plots
  --skip-install        Don't auto-install missing R packages
  --force-rerun         Recompute even if model.rds already exists
  --dry-run             Print the planned R calls without running them
  -h, --help            Show this help and exit

CROSS-VALIDATION
  Auto-selected per (level x variable) combination:
    N >= 20 samples : 5-fold stratified cross-validation
    N <  20 samples : leave-one-out cross-validation (LOOCV)

OUTPUTS
  16_ml_biomarkers/<variable>/RF_<level_name>_by_<variable>/
    model_metrics.xlsx          accuracy, AUC, F1, sens/spec, OOB
    confusion_matrix.{png,pdf}  per-class confusion heatmap
    roc_curves.{png,pdf}        one-vs-rest ROC + AUC labels
    feature_importance.xlsx     full permutation importance table
    top20_importance.{png,pdf}  ranked horizontal barplot
    shap_per_sample.{png,pdf}   per-sample local importance heatmap
    predicted_vs_actual.xlsx    predicted class + per-class probabilities
    model.rds                   ranger fitted model (reusable in R)
  16_ml_biomarkers/<variable>/Summary_RF_<variable>.xlsx
                                  one row per level: accuracy, AUC, top taxa
  16_ml_biomarkers/mbx_ml_biomarkers_info.txt

EXAMPLE
  mbx_ml_classifier_run.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_ml_classifier_run.sh /path/to/mbX_pro_outputs_20260417_121431 --levels g,f
EOF
  exit 0
}

# ── Argument parsing ──────────────────────────────────────────────────────────
MBX_OUT_DIR=""
LEVELS_ARG=""
NUM_TREES=500
SEED=42
SKIP_SHAP=false
SKIP_INSTALL=false
FORCE_RERUN=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)        usage ;;
    --levels)         LEVELS_ARG="$2"; shift 2 ;;
    --num-trees)      NUM_TREES="$2";  shift 2 ;;
    --seed)           SEED="$2";       shift 2 ;;
    --skip-shap)      SKIP_SHAP=true;   shift ;;
    --skip-install)   SKIP_INSTALL=true; shift ;;
    --force-rerun)    FORCE_RERUN=true;  shift ;;
    --dry-run)        DRY_RUN=true;      shift ;;
    -*)               err "Unknown option: $1\n  Run: mbx_ml_classifier_run.sh --help" ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Multiple positional arguments — only one MBX_OUT_DIR expected.\n  Got extra: $1"
      fi
      ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "Missing required argument: <mbX_pro_outputs_dir>
  Run: mbx_ml_classifier_run.sh --help"

[[ -d "$MBX_OUT_DIR" ]] || err "Not a directory: $MBX_OUT_DIR
  Did you pass the path to mbX_pro_outputs_<TIMESTAMP>/?"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# Parse levels list — default = all 7
if [[ -z "$LEVELS_ARG" ]]; then
  LEVELS=("d" "p" "c" "o" "f" "g" "s")
else
  # split on comma into bash 3.2-compatible array
  _OLDIFS="$IFS"; IFS=','; LEVELS=( $LEVELS_ARG ); IFS="$_OLDIFS"
  for _l in "${LEVELS[@]}"; do
    case "$_l" in
      d|p|c|o|f|g|s) ;;
      *) err "Unknown level '$_l' in --levels.  Allowed: d,p,c,o,f,g,s" ;;
    esac
  done
fi

# Map level letters to full name + the prefix used in mbX cleaned file naming
_level_full() {
  case "$1" in
    d) echo "domain" ;;  p) echo "phylum" ;;  c) echo "class" ;;
    o) echo "order"  ;;  f) echo "family" ;;  g) echo "genus" ;;
    s) echo "species" ;;
  esac
}
_level_dirname() {
  # The mbX cleaned dir naming convention (plural, slightly inconsistent)
  case "$1" in
    d) echo "mbX_cleaned_domains_or_kingdoms" ;;
    p) echo "mbX_cleaned_phyla" ;;
    c) echo "mbX_cleaned_classes" ;;
    o) echo "mbX_cleaned_orders" ;;
    f) echo "mbX_cleaned_families" ;;
    g) echo "mbX_cleaned_genera" ;;
    s) echo "mbX_cleaned_species" ;;
  esac
}

# ── Auto-detect CPU count ─────────────────────────────────────────────────────
if command -v nproc &>/dev/null; then N_JOBS="$(nproc)"
elif command -v sysctl &>/dev/null; then N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else N_JOBS=1; fi
[[ "$N_JOBS" -gt 0 ]] || N_JOBS=1

NOW="$(date '+%Y-%m-%d %H:%M:%S')"

# ─────────────────────────────────────────────────────────────────────────────
step "1/7 — Read input paths from previous steps"
# ─────────────────────────────────────────────────────────────────────────────

TAXONOMY_INFO="${MBX_OUT_DIR}/7_taxonomy_csv/mbx_taxonomy_info.txt"
CLEANED_DIR="${MBX_OUT_DIR}/8_cleaned_files"
EZCLEAN_INFO="${CLEANED_DIR}/mbx_ezclean_info.txt"

[[ -f "$TAXONOMY_INFO" ]] || err "Missing: $TAXONOMY_INFO
  -> Run mbx_taxonomy_run.sh (step 7) first."
[[ -d "$CLEANED_DIR" ]] || err "Missing directory: $CLEANED_DIR
  -> Run mbx_ezclean_all_levels.sh (step 8) first."

METADATA_TXT="$(_read_key METADATA_TXT "$TAXONOMY_INFO")"
[[ -n "$METADATA_TXT" && -f "$METADATA_TXT" ]] \
  || err "METADATA_TXT key missing or file not found in $TAXONOMY_INFO"

info "Metadata file       : $METADATA_TXT"
info "Cleaned files dir   : $CLEANED_DIR"
if [[ -f "$EZCLEAN_INFO" ]]; then
  info "ezclean info file   : $EZCLEAN_INFO"
else
  warn "ezclean info file not found ($EZCLEAN_INFO)"
  warn "  -> Falling back to direct disk discovery of cleaned xlsx files."
  warn "  -> (This is fine; older runs of step 8 didn't write an info file.)"
fi
info "Levels requested    : ${LEVELS[*]}"
info "Random forest trees : $NUM_TREES"
info "Random seed         : $SEED"
info "CPU cores available : $N_JOBS"
info "Force re-run        : $FORCE_RERUN"
info "Dry-run             : $DRY_RUN"

# Map each requested level letter to its cleaned xlsx path.
# Strategy: prefer the info file (canonical), but fall back to directly
# searching 8_cleaned_files/ for the standard mbX_cleaned_<plural>_level-7
# subdirectory + matching xlsx (older step-8 runs didn't write the info file).
LEVEL_XLSX_PATHS=()
for _lvl in "${LEVELS[@]}"; do
  _UPLVL="$(echo "$_lvl" | tr '[:lower:]' '[:upper:]')"
  _DIRNAME="$(_level_dirname "$_lvl")"
  _XLSX=""

  if [[ -f "$EZCLEAN_INFO" ]]; then
    _XLSX="$(_read_key "LEVEL_${_UPLVL}_XLSX" "$EZCLEAN_INFO")"
    if [[ "$_XLSX" == "FAILED" || "$_XLSX" == "NOT_FOUND" ]]; then
      _XLSX=""
    fi
  fi

  # Fallback: direct disk discovery
  if [[ -z "$_XLSX" || ! -f "$_XLSX" ]]; then
    # Try the standard naming pattern: mbX_cleaned_<plural>_level-7/...level-7.xlsx
    for _src_dir in "${CLEANED_DIR}/${_DIRNAME}_level-7" "${CLEANED_DIR}/${_DIRNAME}"; do
      if [[ -d "$_src_dir" ]]; then
        _candidate="$(ls "$_src_dir"/*.xlsx 2>/dev/null | head -1)"
        if [[ -n "$_candidate" && -f "$_candidate" ]]; then
          _XLSX="$_candidate"
          break
        fi
      fi
    done
  fi

  if [[ -n "$_XLSX" && -f "$_XLSX" ]]; then
    LEVEL_XLSX_PATHS+=("$_XLSX")
  else
    warn "Level '$_lvl' (${_DIRNAME}): no cleaned xlsx found.  This level will be skipped."
    LEVEL_XLSX_PATHS+=("SKIP")
  fi
done

# ─────────────────────────────────────────────────────────────────────────────
step "2/7 — Locate Rscript (system, outside conda)"
# ─────────────────────────────────────────────────────────────────────────────
# Prefer absolute SYSTEM paths first; fall back to PATH lookup last so that
# `command -v Rscript` returning conda's R doesn't pre-empt system R.
RSCRIPT_CMD=""
for _c in \
    "/opt/homebrew/bin/Rscript" \
    "/usr/local/bin/Rscript" \
    "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
    "/usr/bin/Rscript" \
    "$(command -v Rscript 2>/dev/null || true)"; do
  [[ -n "$_c" && -x "$_c" ]] && { RSCRIPT_CMD="$_c"; break; }
done
[[ -n "$RSCRIPT_CMD" ]] || err "Rscript not found.
  -> Run mbx_ezclean_all_levels.sh (step 8) first to install R, or: brew install r"

R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
ok "Rscript : $RSCRIPT_CMD  (R $R_VERSION)"

# ─────────────────────────────────────────────────────────────────────────────
step "3/7 — Check / install required R packages"
# ─────────────────────────────────────────────────────────────────────────────
# Required packages and what they're for:
#   ranger     fast Random Forest implementation (CRAN)
#   pROC       ROC curves + AUC (multi-class via one-vs-rest)
#   openxlsx   write .xlsx outputs (no Java dependency)
#   ggplot2    plotting
#   pheatmap   confusion matrix heatmap, SHAP heatmap
#   reshape2   wide<->long for ggplot input

_TMPID="${$}_$(date +%s)"
PKG_CHECK_R="/tmp/mbx_ml_pkgchk_${_TMPID}.R"
trap 'rm -f /tmp/mbx_ml_*_${_TMPID}*.R /tmp/mbx_ml_*_${_TMPID}*.txt' EXIT

cat > "$PKG_CHECK_R" << 'RPKG'
required <- c("ranger", "pROC", "openxlsx", "ggplot2", "pheatmap", "reshape2")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) == 0) {
  cat("[OK]    All required R packages already installed.\n")
  cat(sprintf("[INFO]  ranger version: %s\n",
              tryCatch(packageVersion("ranger"), error = function(e) "unknown")))
  quit(status = 0)
}
cat(sprintf("[INFO]  Missing R packages: %s\n", paste(missing, collapse = ", ")))
quit(status = 2)   # 2 => "needs install"
RPKG

if $DRY_RUN; then
  warn "[DRY-RUN] Would check/install R packages: ranger pROC openxlsx ggplot2 pheatmap reshape2"
else
  set +e
  _R --vanilla "$PKG_CHECK_R"
  RC=$?
  set -e
  if [[ $RC -eq 0 ]]; then
    :
  elif [[ $RC -eq 2 ]]; then
    if $SKIP_INSTALL; then
      err "Required R packages are missing and --skip-install was set.
  -> Install manually:
       Rscript -e 'install.packages(c(\"ranger\",\"pROC\",\"openxlsx\",\"ggplot2\",\"pheatmap\",\"reshape2\"))'
  -> Then re-run this script."
    fi
    info "Installing missing R packages from CRAN (may take 1-3 minutes)..."
    PKG_INSTALL_R="/tmp/mbx_ml_pkginst_${_TMPID}.R"
    cat > "$PKG_INSTALL_R" << 'RPKGI'
required <- c("ranger", "pROC", "openxlsx", "ggplot2", "pheatmap", "reshape2")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org",
                   quiet = FALSE, dependencies = TRUE)
  failed <- missing[!sapply(missing, requireNamespace, quietly = TRUE)]
  if (length(failed) > 0) {
    stop(sprintf("Failed to install: %s", paste(failed, collapse = ", ")))
  }
}
cat("[OK]    All required R packages installed.\n")
RPKGI
    _R --vanilla "$PKG_INSTALL_R" \
      || err "R package installation failed.  See output above."
    rm -f "$PKG_INSTALL_R"
    ok "All required R packages installed."
  else
    err "R package check failed (exit $RC)."
  fi
fi
rm -f "$PKG_CHECK_R"

# ─────────────────────────────────────────────────────────────────────────────
step "4/7 — Detect categorical metadata columns"
# ─────────────────────────────────────────────────────────────────────────────
# Same logic as ezviz / ezstat:
#   * skip QIIME2 #q2:types row
#   * skip the sample-id column
#   * skip numeric columns
#   * skip all-unique columns (free-text / IDs)
#   * skip constant / single-value columns
#   * require at least one group with >= 2 samples (RF needs replication)

CAT_DETECT_R="/tmp/mbx_ml_detect_${_TMPID}.R"
CATS_STDOUT="/tmp/mbx_ml_cats_${_TMPID}.out"
CATS_STDERR="/tmp/mbx_ml_cats_${_TMPID}.err"

# The R script writes column names ONLY between the BEGIN/END sentinels on
# stdout.  All [INFO] and warnings go to stderr.  Bash parses only the lines
# strictly between the sentinels — this is bulletproof against R printing
# stray NULLs, package-loaded messages, or version warnings.
cat > "$CAT_DETECT_R" << RDETECT
suppressPackageStartupMessages({ })
mfile <- "${METADATA_TXT}"
ext   <- tolower(tools::file_ext(mfile))
if (ext == "csv") {
  meta <- read.csv(mfile, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, comment.char = "")
} else {
  meta <- read.delim(mfile, header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE, comment.char = "")
}
if (nrow(meta) > 0 && grepl("^#q2:types", as.character(meta[1,1]),
                            ignore.case = TRUE)) {
  meta <- meta[-1, , drop = FALSE]
}
n_samples <- nrow(meta)
col_names <- names(meta)
sid_re <- "^(sample[-_ ]?id|sampleid|id|featureid|feature[-_ ]id)\$"
sid_idx <- which(grepl(sid_re, col_names, ignore.case = TRUE))
if (length(sid_idx) == 0) sid_idx <- 1
candidates <- col_names[-sid_idx]

categorical <- character(0)
for (col in candidates) {
  vals_raw <- meta[[col]]
  vals     <- vals_raw[!is.na(vals_raw) & vals_raw != ""]
  if (length(vals) == 0) next
  num_test <- suppressWarnings(as.numeric(vals))
  if (!any(is.na(num_test))) next
  n_unique <- length(unique(trimws(as.character(vals))))
  if (n_unique <= 1)         next
  if (n_unique == n_samples) next
  group_counts <- table(trimws(as.character(vals)))
  if (max(group_counts) < 2) next
  categorical <- c(categorical, col)
}

cat(sprintf("[INFO]  Metadata: %d samples, %d candidate columns.\n",
            n_samples, length(candidates)), file = stderr())
cat(sprintf("[INFO]  Categorical columns kept (%d): %s\n",
            length(categorical), paste(categorical, collapse = ", ")),
    file = stderr())

cat("MBX_CATS_BEGIN\n")
for (cc in categorical) cat(cc, "\n", sep = "")
cat("MBX_CATS_END\n")
invisible(NULL)
RDETECT

if $DRY_RUN; then
  warn "[DRY-RUN] Would detect categorical columns from $METADATA_TXT"
  CATEGORICAL_COLS=("Treatment" "SampleType")
  info "Dry-run placeholder columns: ${CATEGORICAL_COLS[*]}"
else
  _R --vanilla "$CAT_DETECT_R" > "$CATS_STDOUT" 2> "$CATS_STDERR" \
    || { cat "$CATS_STDERR" >&2; err "Categorical column detection failed.
  -> Check your metadata file: $METADATA_TXT"; }

  grep "^\[" "$CATS_STDERR" || true

  # Parse ONLY the lines strictly between MBX_CATS_BEGIN and MBX_CATS_END
  CATEGORICAL_COLS=()
  _IN_BLOCK=false
  while IFS= read -r line; do
    if [[ "$line" == "MBX_CATS_BEGIN" ]]; then _IN_BLOCK=true;  continue; fi
    if [[ "$line" == "MBX_CATS_END"   ]]; then _IN_BLOCK=false; continue; fi
    [[ "$_IN_BLOCK" == true ]] || continue
    [[ -z "$line" ]] && continue
    [[ "$line" == "NULL" ]] && continue   # belt-and-suspenders
    CATEGORICAL_COLS+=("$line")
  done < "$CATS_STDOUT"

  if [[ ${#CATEGORICAL_COLS[@]} -eq 0 ]]; then
    err "No categorical columns found in metadata.
  -> Need at least one column with >= 2 distinct non-numeric values
     and at least one group with 2+ samples."
  fi
  ok "Categorical columns: ${CATEGORICAL_COLS[*]}"
fi
rm -f "$CAT_DETECT_R" "$CATS_STDOUT" "$CATS_STDERR"

# ─────────────────────────────────────────────────────────────────────────────
step "5/7 — Create output directories"
# ─────────────────────────────────────────────────────────────────────────────
ML_DIR="${MBX_OUT_DIR}/16_ml_biomarkers"
WORK_DIR="${ML_DIR}/working_dir_ml"

if $DRY_RUN; then
  info "[DRY-RUN] Would create: $ML_DIR/"
  info "[DRY-RUN] Would create: $WORK_DIR/"
  for col in "${CATEGORICAL_COLS[@]}"; do
    info "[DRY-RUN] Would create: 16_ml_biomarkers/$(_sanitize_dirname "$col")/"
  done
else
  mkdir -p "$WORK_DIR" \
    || err "Could not create: $WORK_DIR — check permissions."
  for col in "${CATEGORICAL_COLS[@]}"; do
    _SC="$(_sanitize_dirname "$col")"
    mkdir -p "${ML_DIR}/${_SC}" \
      || err "Could not create: ${ML_DIR}/${_SC}"
    info "  Created: 16_ml_biomarkers/${_SC}/"
  done
fi

# ─────────────────────────────────────────────────────────────────────────────
step "6/7 — Train Random Forest for each level x variable"
# ─────────────────────────────────────────────────────────────────────────────

_TOTAL_RUNS=$(( ${#LEVELS[@]} * ${#CATEGORICAL_COLS[@]} ))
info "Total RF runs planned : $_TOTAL_RUNS  (${#LEVELS[@]} levels x ${#CATEGORICAL_COLS[@]} variables)"
info "Estimated runtime     : ~5-30 s per run (i.e. ~$(( _TOTAL_RUNS / 10 + 1 )) minutes total)"

_RUN_IDX=0
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
RESULT_LINES=""

# ── Build the R analysis script template (used for every combination) ────────
# This is large; written once to a temp file and reused 21+ times via env vars.
R_SCRIPT_TEMPLATE="/tmp/mbx_ml_run_${_TMPID}.R"
cat > "$R_SCRIPT_TEMPLATE" << 'RTPL'
suppressPackageStartupMessages({
  library(openxlsx)
  library(ranger)
  library(pROC)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
})

# Bash injects the parameters via env vars (so the template is reusable)
CLEANED_XLSX <- Sys.getenv("MBX_ML_CLEANED_XLSX")
META_TXT     <- Sys.getenv("MBX_ML_META_TXT")
VAR_COL      <- Sys.getenv("MBX_ML_VAR_COL")
LEVEL_LET    <- Sys.getenv("MBX_ML_LEVEL_LET")
LEVEL_NAME   <- Sys.getenv("MBX_ML_LEVEL_NAME")
OUT_DIR      <- Sys.getenv("MBX_ML_OUT_DIR")
NUM_TREES    <- as.integer(Sys.getenv("MBX_ML_NUM_TREES", "500"))
SEED         <- as.integer(Sys.getenv("MBX_ML_SEED", "42"))
DO_SHAP      <- toupper(Sys.getenv("MBX_ML_DO_SHAP", "TRUE")) == "TRUE"
N_JOBS       <- as.integer(Sys.getenv("MBX_ML_N_JOBS", "1"))

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 1. Read cleaned table (samples in rows; col 1 = sample-id) ──────────────
df <- read.xlsx(CLEANED_XLSX, check.names = FALSE)
if (nrow(df) < 4)
  stop(sprintf("Too few samples in %s (%d) — need >= 4", CLEANED_XLSX, nrow(df)))

sid_col <- names(df)[1]
sids    <- as.character(df[[sid_col]])

# Read metadata to identify which columns in df are metadata vs taxa
ext <- tolower(tools::file_ext(META_TXT))
meta <- if (ext == "csv") {
  read.csv(META_TXT, header = TRUE, check.names = FALSE,
           stringsAsFactors = FALSE, comment.char = "")
} else {
  read.delim(META_TXT, header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, comment.char = "")
}
if (nrow(meta) > 0 && grepl("^#q2:types", as.character(meta[1,1]),
                            ignore.case = TRUE)) {
  meta <- meta[-1, , drop = FALSE]
}
meta_cols <- intersect(names(meta), names(df))   # mbX cleaned merged metadata in
if (!(VAR_COL %in% names(df))) {
  if (VAR_COL %in% names(meta)) {
    df <- merge(df, meta[, c(names(meta)[1], VAR_COL), drop = FALSE],
                by.x = sid_col, by.y = names(meta)[1], all.x = TRUE)
    meta_cols <- c(meta_cols, VAR_COL)
  } else {
    stop(sprintf("Variable '%s' not found in cleaned xlsx OR metadata.", VAR_COL))
  }
}

# Class label
y_raw <- as.character(df[[VAR_COL]])
keep  <- !(is.na(y_raw) | trimws(y_raw) == "")
if (sum(keep) < 4)
  stop(sprintf("Only %d non-missing labels for '%s' — need >= 4 to fit RF",
               sum(keep), VAR_COL))
df    <- df[keep, , drop = FALSE]
sids  <- sids[keep]
y     <- factor(y_raw[keep])

# Feature matrix = everything except sample-id and metadata cols
feat_cols <- setdiff(names(df), c(sid_col, meta_cols))
if (length(feat_cols) < 2)
  stop(sprintf("Only %d taxa columns at level '%s' — need >= 2 for RF",
               length(feat_cols), LEVEL_LET))

X <- as.data.frame(lapply(df[, feat_cols, drop = FALSE], function(v) {
  v <- suppressWarnings(as.numeric(as.character(v)))
  v[is.na(v)] <- 0
  v
}), check.names = FALSE)
rownames(X) <- sids

# Drop zero-variance taxa (no information)
v_var <- apply(X, 2, function(x) stats::var(x, na.rm = TRUE))
v_var[is.na(v_var)] <- 0
X <- X[, v_var > 0, drop = FALSE]
if (ncol(X) < 2)
  stop(sprintf("After zero-variance filter, only %d taxa remain.", ncol(X)))

n_samp <- nrow(X)
n_feat <- ncol(X)
n_class <- nlevels(y)
class_counts <- table(y)

cat(sprintf("[INFO]  N samples  : %d\n", n_samp))
cat(sprintf("[INFO]  N features : %d\n", n_feat))
cat(sprintf("[INFO]  N classes  : %d  (%s)\n", n_class,
            paste(sprintf("%s=%d", names(class_counts), as.integer(class_counts)),
                  collapse = ", ")))

# ── 2. Cross-validation strategy ─────────────────────────────────────────────
# 5-fold stratified when N >= 20; else LOOCV.  If any class has < 2 samples,
# stratification falls back to LOOCV automatically (RF can't be evaluated
# meaningfully with a singleton class otherwise).
if (n_samp >= 20 && min(class_counts) >= 2) {
  cv_method <- "5-fold stratified"
  set.seed(SEED)
  fold_ids <- integer(n_samp)
  for (lev in levels(y)) {
    ix <- which(y == lev)
    fold_ids[ix] <- ((sample(seq_along(ix)) - 1) %% 5) + 1
  }
  K <- 5
} else {
  cv_method <- "leave-one-out (LOOCV)"
  fold_ids  <- seq_len(n_samp)
  K <- n_samp
}
cat(sprintf("[INFO]  CV method  : %s\n", cv_method))

# Class weights for imbalance (inverse class frequency)
imbalance_ratio <- max(class_counts) / max(min(class_counts), 1)
class_w <- if (imbalance_ratio > 2) {
  cat(sprintf("[INFO]  Class imbalance ratio = %.1fx -> applying inverse-frequency weights\n",
              imbalance_ratio))
  setNames(as.numeric(1 / class_counts), names(class_counts))
} else {
  NULL
}

# ── 3. Cross-validation loop ─────────────────────────────────────────────────
prob_mat <- matrix(NA_real_, nrow = n_samp, ncol = n_class,
                   dimnames = list(sids, levels(y)))
pred_lab <- character(n_samp)
oob_errs <- numeric(K)

case_w_fn <- function(y_train) {
  if (is.null(class_w)) return(NULL)
  unname(class_w[as.character(y_train)])
}

for (k in seq_len(K)) {
  test_ix  <- which(fold_ids == k)
  if (length(test_ix) == 0) next
  train_ix <- setdiff(seq_len(n_samp), test_ix)
  if (length(unique(y[train_ix])) < 2) {
    # Can happen with LOOCV when held-out is a singleton class — predict mode
    pred_lab[test_ix] <- as.character(names(sort(table(y[train_ix]),
                                                 decreasing = TRUE)[1]))
    next
  }
  fit_k <- ranger(
    x             = X[train_ix, , drop = FALSE],
    y             = y[train_ix],
    num.trees     = NUM_TREES,
    probability   = TRUE,
    importance    = "none",
    case.weights  = case_w_fn(y[train_ix]),
    num.threads   = N_JOBS,
    seed          = SEED + k
  )
  pp <- predict(fit_k, data = X[test_ix, , drop = FALSE])$predictions
  # Align prediction columns to canonical level order
  prob_mat[test_ix, colnames(pp)] <- pp
  pred_lab[test_ix] <- colnames(pp)[max.col(pp, ties.method = "first")]
  oob_errs[k] <- fit_k$prediction.error
}
pred_fac <- factor(pred_lab, levels = levels(y))
mean_oob <- mean(oob_errs[oob_errs > 0 | TRUE], na.rm = TRUE)

# ── 4. Metrics ───────────────────────────────────────────────────────────────
cm <- table(actual = y, predicted = pred_fac)
acc <- sum(diag(cm)) / sum(cm)

# Per-class metrics (one-vs-rest)
per_class <- data.frame(
  class = levels(y),
  n     = as.integer(class_counts[levels(y)]),
  sensitivity = NA_real_,
  specificity = NA_real_,
  precision   = NA_real_,
  f1          = NA_real_,
  AUC         = NA_real_,
  stringsAsFactors = FALSE
)
for (i in seq_len(n_class)) {
  cls   <- levels(y)[i]
  TP <- cm[cls, cls]
  FN <- sum(cm[cls, ]) - TP
  FP <- sum(cm[, cls]) - TP
  TN <- sum(cm) - TP - FN - FP
  per_class$sensitivity[i] <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  per_class$specificity[i] <- if (TN + FP > 0) TN / (TN + FP) else NA_real_
  per_class$precision[i]   <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  prec <- per_class$precision[i]; sens <- per_class$sensitivity[i]
  per_class$f1[i] <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
                      2 * prec * sens / (prec + sens) else NA_real_
  # AUC via pROC, one-vs-rest
  y_bin <- as.integer(y == cls)
  auc_k <- tryCatch(
    suppressMessages(as.numeric(pROC::auc(pROC::roc(y_bin, prob_mat[, cls],
                                                    quiet = TRUE)))),
    error = function(e) NA_real_)
  per_class$AUC[i] <- auc_k
}

# Macro-AUC (mean of per-class AUCs ignoring NAs)
macro_auc <- mean(per_class$AUC, na.rm = TRUE)
macro_f1  <- mean(per_class$f1,  na.rm = TRUE)

cat(sprintf("[OK]    Accuracy = %.3f, macro-AUC = %.3f, macro-F1 = %.3f\n",
            acc, macro_auc, macro_f1))

# ── 5. Final model on ALL data (for importance + reuse) ──────────────────────
fit_full <- ranger(
  x            = X,
  y            = y,
  num.trees    = NUM_TREES,
  probability  = TRUE,
  importance   = "permutation",
  case.weights = case_w_fn(y),
  local.importance = DO_SHAP,
  num.threads  = N_JOBS,
  seed         = SEED
)

# Permutation importance (full model)
imp <- ranger::importance(fit_full)
imp_df <- data.frame(taxon = names(imp),
                     importance = as.numeric(imp),
                     stringsAsFactors = FALSE)
imp_df <- imp_df[order(-imp_df$importance), , drop = FALSE]

# ── 6. Write metrics xlsx ────────────────────────────────────────────────────
metrics_overall <- data.frame(
  metric = c("N_samples", "N_features", "N_classes", "CV_method",
             "CV_accuracy", "CV_macro_AUC", "CV_macro_F1",
             "Mean_fold_OOB_error", "imbalance_ratio",
             "num_trees", "seed"),
  value  = c(n_samp, n_feat, n_class, cv_method,
             sprintf("%.4f", acc), sprintf("%.4f", macro_auc),
             sprintf("%.4f", macro_f1),
             sprintf("%.4f", mean_oob), sprintf("%.2f", imbalance_ratio),
             NUM_TREES, SEED),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overall");    writeData(wb, "overall",    metrics_overall)
addWorksheet(wb, "per_class");  writeData(wb, "per_class",  per_class)
addWorksheet(wb, "confusion");  writeData(wb, "confusion",
                                          as.data.frame.matrix(cm),
                                          rowNames = TRUE)
saveWorkbook(wb, file.path(OUT_DIR, "model_metrics.xlsx"), overwrite = TRUE)

# ── 7. Confusion matrix heatmap ──────────────────────────────────────────────
cm_pct <- sweep(as.matrix(cm), 1, pmax(rowSums(cm), 1), "/") * 100
disp_lab <- matrix(sprintf("%d\n(%.0f%%)", as.integer(cm), as.numeric(cm_pct)),
                   nrow = nrow(cm))
title_str <- sprintf("Confusion Matrix -- %s by %s\nAccuracy = %.1f%% (CV = %s)",
                     LEVEL_NAME, VAR_COL, acc * 100, cv_method)
png(file.path(OUT_DIR, "confusion_matrix.png"),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(cm_pct, color = colorRampPalette(c("white", "#3b73af"))(100),
         display_numbers = disp_lab, fontsize_number = 11,
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = title_str)
dev.off()
pdf(file.path(OUT_DIR, "confusion_matrix.pdf"), width = 7, height = 6)
pheatmap(cm_pct, color = colorRampPalette(c("white", "#3b73af"))(100),
         display_numbers = disp_lab, fontsize_number = 11,
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = title_str)
dev.off()

# ── 8. ROC curves (one-vs-rest, all classes overlaid) ────────────────────────
roc_df_list <- list()
for (cls in levels(y)) {
  y_bin <- as.integer(y == cls)
  if (length(unique(y_bin)) < 2) next
  rk <- tryCatch(suppressMessages(pROC::roc(y_bin, prob_mat[, cls], quiet = TRUE)),
                 error = function(e) NULL)
  if (is.null(rk)) next
  roc_df_list[[cls]] <- data.frame(
    FPR = 1 - rk$specificities,
    TPR = rk$sensitivities,
    class = sprintf("%s (AUC = %.3f)", cls, as.numeric(pROC::auc(rk))),
    stringsAsFactors = FALSE
  )
}
if (length(roc_df_list) > 0) {
  roc_df <- do.call(rbind, roc_df_list)
  p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR, color = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_step(linewidth = 1) +
    coord_equal() +
    scale_x_continuous("False positive rate", limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous("True positive rate",  limits = c(0, 1), expand = c(0, 0)) +
    labs(title = sprintf("ROC curves -- %s by %s", LEVEL_NAME, VAR_COL),
         subtitle = sprintf("Cross-validation: %s", cv_method),
         color = "Class") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40"))
  ggsave(file.path(OUT_DIR, "roc_curves.png"), p_roc,
         width = 7, height = 6, dpi = 300)
  ggsave(file.path(OUT_DIR, "roc_curves.pdf"), p_roc,
         width = 7, height = 6)
}

# ── 9. Feature importance: full xlsx + top-20 barplot ────────────────────────
write.xlsx(imp_df, file.path(OUT_DIR, "feature_importance.xlsx"),
           overwrite = TRUE)

top_n <- min(20, nrow(imp_df))
top_df <- imp_df[seq_len(top_n), , drop = FALSE]
top_df$taxon <- factor(top_df$taxon, levels = rev(top_df$taxon))
p_imp <- ggplot(top_df, aes(x = importance, y = taxon)) +
  geom_col(fill = "#3b73af") +
  labs(title = sprintf("Top %d taxa by permutation importance", top_n),
       subtitle = sprintf("%s by %s -- accuracy = %.1f%%, macro-AUC = %.3f",
                          LEVEL_NAME, VAR_COL, acc * 100, macro_auc),
       x = "Mean decrease in accuracy (permutation)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40"),
        axis.text.y   = element_text(size = 10))
ggsave(file.path(OUT_DIR, "top20_importance.png"), p_imp,
       width = 9, height = max(5, 0.32 * top_n), dpi = 300)
ggsave(file.path(OUT_DIR, "top20_importance.pdf"), p_imp,
       width = 9, height = max(5, 0.32 * top_n))

# ── 10. SHAP-style per-sample contributions ──────────────────────────────────
# ranger returns local.importance: a (samples x features) matrix where each
# entry is the per-sample permutation importance contribution.  This is the
# Saabas-equivalent precursor to TreeSHAP and is mathematically interpretable
# as "how much this taxon's value pushed THIS sample's prediction away from
# the dataset mean prediction" (sign + magnitude).
if (DO_SHAP && !is.null(fit_full$variable.importance.local)) {
  loc_imp <- fit_full$variable.importance.local
  # Restrict to the same top-20 taxa we plotted above
  shap_top <- intersect(as.character(top_df$taxon), colnames(loc_imp))
  shap_mat <- loc_imp[, shap_top, drop = FALSE]
  rownames(shap_mat) <- sids
  # Annotate samples with their true class for visual grouping
  ann_row <- data.frame(class = y, row.names = sids)
  shap_title <- sprintf(
    "Per-sample feature contributions (Saabas/SHAP-style)\n%s by %s -- top %d taxa",
    LEVEL_NAME, VAR_COL, length(shap_top))
  if (length(shap_top) >= 2 && nrow(shap_mat) >= 2) {
    png(file.path(OUT_DIR, "shap_per_sample.png"),
        width = 11, height = max(5, 0.3 * nrow(shap_mat) + 2),
        units = "in", res = 300)
    pheatmap(t(shap_mat),
             color = colorRampPalette(c("#3b73af", "white", "#c0392b"))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             annotation_col = ann_row,
             show_rownames = TRUE, show_colnames = TRUE,
             fontsize_row = 9, fontsize_col = 8,
             main = shap_title)
    dev.off()
    pdf(file.path(OUT_DIR, "shap_per_sample.pdf"),
        width = 11, height = max(5, 0.3 * nrow(shap_mat) + 2))
    pheatmap(t(shap_mat),
             color = colorRampPalette(c("#3b73af", "white", "#c0392b"))(100),
             cluster_rows = TRUE, cluster_cols = TRUE,
             annotation_col = ann_row,
             show_rownames = TRUE, show_colnames = TRUE,
             fontsize_row = 9, fontsize_col = 8,
             main = shap_title)
    dev.off()
    write.xlsx(as.data.frame(shap_mat, check.names = FALSE),
               file.path(OUT_DIR, "shap_per_sample.xlsx"),
               rowNames = TRUE, overwrite = TRUE)
  } else {
    cat("[INFO]  SHAP heatmap skipped (need >= 2 samples and >= 2 top taxa)\n")
  }
} else if (DO_SHAP) {
  cat("[INFO]  SHAP local importance unavailable for this model.\n")
}

# ── 11. Predicted vs actual (per-sample table) ──────────────────────────────
pred_df <- data.frame(
  sample_id     = sids,
  actual_class  = as.character(y),
  predicted_class = as.character(pred_fac),
  correct       = as.character(y) == as.character(pred_fac),
  stringsAsFactors = FALSE
)
prob_mat_df <- as.data.frame(prob_mat, check.names = FALSE)
names(prob_mat_df) <- paste0("prob_", names(prob_mat_df))
pred_df <- cbind(pred_df, prob_mat_df)
write.xlsx(pred_df, file.path(OUT_DIR, "predicted_vs_actual.xlsx"),
           overwrite = TRUE)

# ── 12. Save model object ────────────────────────────────────────────────────
saveRDS(fit_full, file.path(OUT_DIR, "model.rds"))

# ── 13. One-line summary written to stdout for the bash caller to pick up ────
top_taxa_str <- paste(head(imp_df$taxon, 5), collapse = " | ")
cat(sprintf("MBX_RF_SUMMARY\tlevel=%s\tlevel_name=%s\tvariable=%s\tn=%d\tn_feat=%d\tn_class=%d\tacc=%.4f\tauc=%.4f\tf1=%.4f\toob=%.4f\tcv=%s\ttop5=%s\n",
            LEVEL_LET, LEVEL_NAME, VAR_COL, n_samp, n_feat, n_class,
            acc, macro_auc, macro_f1, mean_oob, cv_method, top_taxa_str))

cat("[OK]    All RF outputs written.\n")
RTPL

# ── Loop over (level x variable) ─────────────────────────────────────────────
SUMMARY_LINES_FILE="${WORK_DIR}/_summary_lines.tsv"
: > "$SUMMARY_LINES_FILE"

for j in "${!CATEGORICAL_COLS[@]}"; do
  COL="${CATEGORICAL_COLS[$j]}"
  COL_DIR="$(_sanitize_dirname "$COL")"
  for i in "${!LEVELS[@]}"; do
    LVL="${LEVELS[$i]}"
    LVL_NAME="$(_level_full "$LVL")"
    XLSX="${LEVEL_XLSX_PATHS[$i]}"
    _RUN_IDX=$(( _RUN_IDX + 1 ))

    info "[$_RUN_IDX/$_TOTAL_RUNS] RF -- level '${LVL}' (${LVL_NAME}) | variable '${COL}'"

    if [[ "$XLSX" == "SKIP" ]]; then
      warn "  Skipping (no cleaned xlsx for level '${LVL}')"
      SKIP_COUNT=$(( SKIP_COUNT + 1 ))
      RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} : SKIPPED (missing input)\n"
      continue
    fi

    OUT_RUN_DIR="${ML_DIR}/${COL_DIR}/RF_${LVL_NAME}_by_${COL_DIR}"
    MODEL_RDS="${OUT_RUN_DIR}/model.rds"

    if [[ -f "$MODEL_RDS" && "$FORCE_RERUN" == false ]]; then
      skipped "  Existing model.rds (use --force-rerun to recompute)"
      PASS_COUNT=$(( PASS_COUNT + 1 ))
      continue
    fi

    if $DRY_RUN; then
      warn "  [DRY-RUN] Would train RF using:"
      echo "    cleaned xlsx : $XLSX"
      echo "    output dir   : $OUT_RUN_DIR"
      echo "    num.trees    : $NUM_TREES,  seed = $SEED,  threads = $N_JOBS"
      RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} : DRY_RUN\n"
      continue
    fi

    mkdir -p "$OUT_RUN_DIR" \
      || { warn "  Could not create $OUT_RUN_DIR"; FAIL_COUNT=$(( FAIL_COUNT + 1 )); continue; }

    timer_start
    set +e
    MBX_ML_CLEANED_XLSX="$XLSX" \
    MBX_ML_META_TXT="$METADATA_TXT" \
    MBX_ML_VAR_COL="$COL" \
    MBX_ML_LEVEL_LET="$LVL" \
    MBX_ML_LEVEL_NAME="$LVL_NAME" \
    MBX_ML_OUT_DIR="$OUT_RUN_DIR" \
    MBX_ML_NUM_TREES="$NUM_TREES" \
    MBX_ML_SEED="$SEED" \
    MBX_ML_DO_SHAP="$([[ "$SKIP_SHAP" == true ]] && echo FALSE || echo TRUE)" \
    MBX_ML_N_JOBS="$N_JOBS" \
      _R --vanilla "$R_SCRIPT_TEMPLATE" 2>&1 | tee "${OUT_RUN_DIR}/_run.log"
    RC=${PIPESTATUS[0]}
    set -e
    timer_end

    if [[ $RC -eq 0 ]]; then
      ok "  ✔  ${COL} x ${LVL_NAME}"
      PASS_COUNT=$(( PASS_COUNT + 1 ))
      # Capture the MBX_RF_SUMMARY line for later aggregation
      grep '^MBX_RF_SUMMARY' "${OUT_RUN_DIR}/_run.log" >> "$SUMMARY_LINES_FILE" 2>/dev/null || true
      RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} : OK\n"
    else
      warn "  ✘  ${COL} x ${LVL_NAME} -- see ${OUT_RUN_DIR}/_run.log"
      FAIL_COUNT=$(( FAIL_COUNT + 1 ))
      RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} : FAILED\n"
    fi
  done
done

rm -f "$R_SCRIPT_TEMPLATE"

# ─────────────────────────────────────────────────────────────────────────────
step "7/7 — Build per-variable summary xlsx"
# ─────────────────────────────────────────────────────────────────────────────
# For each categorical variable, aggregate one row per level so the user can
# eyeball "which taxonomic level is most predictive for this variable".

if $DRY_RUN; then
  warn "[DRY-RUN] Would aggregate per-variable summary xlsx files."
elif [[ ! -s "$SUMMARY_LINES_FILE" ]]; then
  warn "No successful RF runs -- no Summary_RF_<variable>.xlsx files written."
else
  AGG_R="/tmp/mbx_ml_agg_${_TMPID}.R"
  cat > "$AGG_R" << RAGG
suppressPackageStartupMessages({ library(openxlsx) })
LINES_FILE <- "${SUMMARY_LINES_FILE}"
ML_DIR     <- "${ML_DIR}"
sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^[._-]+|[._-]+\$", "", x)
  x
}

lines <- readLines(LINES_FILE)
lines <- lines[nchar(lines) > 0 & startsWith(lines, "MBX_RF_SUMMARY")]
if (length(lines) == 0) {
  cat("[INFO]  No summary lines found.\n"); quit(status = 0)
}

# Each line: TAB-separated key=value pairs after the leading tag
parse_line <- function(s) {
  parts <- strsplit(s, "\t", fixed = TRUE)[[1]][-1]
  kv <- strsplit(parts, "=", fixed = TRUE)
  setNames(sapply(kv, function(x) paste(x[-1], collapse = "=")),
           sapply(kv, '[', 1))
}
rows <- lapply(lines, parse_line)
all_keys <- unique(unlist(lapply(rows, names)))
df <- do.call(rbind, lapply(rows, function(r) {
  v <- setNames(rep(NA_character_, length(all_keys)), all_keys)
  v[names(r)] <- r
  as.data.frame(as.list(v), stringsAsFactors = FALSE)
}))

# Coerce numeric columns
for (k in c("n", "n_feat", "n_class", "acc", "auc", "f1", "oob")) {
  if (k %in% names(df)) df[[k]] <- suppressWarnings(as.numeric(df[[k]]))
}

# Order taxonomic levels d -> p -> c -> o -> f -> g -> s
lvl_order <- c("d", "p", "c", "o", "f", "g", "s")
df\$level <- factor(df\$level, levels = lvl_order)

vars <- unique(df\$variable)
for (v in vars) {
  sub <- df[df\$variable == v, , drop = FALSE]
  sub <- sub[order(sub\$level), , drop = FALSE]
  out <- data.frame(
    taxonomic_level = paste0(sub\$level, " (", sub\$level_name, ")"),
    n_samples       = sub\$n,
    n_features      = sub\$n_feat,
    n_classes       = sub\$n_class,
    accuracy        = round(sub\$acc, 4),
    macro_AUC       = round(sub\$auc, 4),
    macro_F1        = round(sub\$f1,  4),
    mean_OOB        = round(sub\$oob, 4),
    cv_method       = sub\$cv,
    top5_taxa       = sub\$top5,
    stringsAsFactors = FALSE
  )
  best_ix <- which.max(out\$accuracy)
  cat(sprintf("[INFO]  Variable '%s': best level = %s (acc = %.3f, AUC = %.3f)\n",
              v, out\$taxonomic_level[best_ix],
              out\$accuracy[best_ix], out\$macro_AUC[best_ix]))
  vd <- file.path(ML_DIR, sanitize(v))
  dir.create(vd, showWarnings = FALSE, recursive = TRUE)
  write.xlsx(out, file.path(vd, sprintf("Summary_RF_%s.xlsx", sanitize(v))),
             overwrite = TRUE)
}
cat("[OK]    Per-variable summaries written.\n")
RAGG

  _R --vanilla "$AGG_R" || warn "Per-variable summary aggregation failed (non-fatal)."
  rm -f "$AGG_R"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Write provenance file
# ─────────────────────────────────────────────────────────────────────────────
INFO_TXT="${ML_DIR}/mbx_ml_biomarkers_info.txt"
cat > "$INFO_TXT" << INFO
# ============================================================================
# mbx_ml_biomarkers_info.txt
# Generated by mbx_ml_classifier_run.sh   (step 16)
# Date : $NOW
# ============================================================================

# ── Inputs used ──────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
METADATA_TXT=$METADATA_TXT
EZCLEAN_INFO=$EZCLEAN_INFO
R_VERSION=$R_VERSION
RANGER_NUM_TREES=$NUM_TREES
RANGER_SEED=$SEED
SHAP_ENABLED=$([[ "$SKIP_SHAP" == true ]] && echo FALSE || echo TRUE)
N_JOBS=$N_JOBS

# ── Output structure ─────────────────────────────────────────────────────────
ML_DIR=$ML_DIR
WORK_DIR=$WORK_DIR

# ── Per-(level x variable) outputs ───────────────────────────────────────────
# For every successful run, the following files exist inside
#   16_ml_biomarkers/<sanitized_variable>/RF_<level_name>_by_<sanitized_variable>/
#
#   model_metrics.xlsx          accuracy, AUC, F1, sens/spec per class
#   confusion_matrix.{png,pdf}  per-class confusion heatmap
#   roc_curves.{png,pdf}        one-vs-rest ROC + AUC labels
#   feature_importance.xlsx     full permutation importance table
#   top20_importance.{png,pdf}  ranked horizontal barplot
#   shap_per_sample.{png,pdf}   per-sample local importance heatmap (if SHAP enabled)
#   shap_per_sample.xlsx        raw per-sample contributions (if SHAP enabled)
#   predicted_vs_actual.xlsx    per-sample predicted class + class probabilities
#   model.rds                   ranger fitted model object (reusable in R)
#   _run.log                    full R stdout/stderr for this run

# ── Variables analysed ───────────────────────────────────────────────────────
INFO
for col in "${CATEGORICAL_COLS[@]}"; do
  echo "VARIABLE=$col" >> "$INFO_TXT"
done
cat >> "$INFO_TXT" << INFO2

# ── Levels analysed ──────────────────────────────────────────────────────────
LEVELS_REQUESTED=${LEVELS[*]}

# ── Run summary ──────────────────────────────────────────────────────────────
TOTAL_RUNS=$_TOTAL_RUNS
PASSED_RUNS=$PASS_COUNT
FAILED_RUNS=$FAIL_COUNT
SKIPPED_RUNS=$SKIP_COUNT
INFO2

ok "Provenance written -> $INFO_TXT"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete -- no R code was executed."
else
  echo ""
  echo "  ╔══════════════════════════════════════════════════════════════╗"
  echo "  ║  Random Forest biomarker analysis complete                   ║"
  echo "  ╚══════════════════════════════════════════════════════════════╝"
  echo ""
  echo "  Output directory : $ML_DIR"
  echo "  Total runs       : $_TOTAL_RUNS"
  echo "    - Succeeded    : $PASS_COUNT"
  echo "    - Failed       : $FAIL_COUNT"
  echo "    - Skipped      : $SKIP_COUNT"
  echo ""
  if [[ -n "$RESULT_LINES" ]]; then
    echo "  Per-combination outcome:"
    printf '%b' "$RESULT_LINES" | sed 's/^/    /'
  fi
  echo ""
  echo "  Inspect the per-variable summary xlsx files first:"
  for col in "${CATEGORICAL_COLS[@]}"; do
    SC="$(_sanitize_dirname "$col")"
    echo "    ${ML_DIR}/${SC}/Summary_RF_${SC}.xlsx"
  done
  echo ""
  echo "  Provenance file : $INFO_TXT"
  echo ""
fi
sep
