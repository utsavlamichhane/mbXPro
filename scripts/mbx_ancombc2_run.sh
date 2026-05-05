#!/usr/bin/env bash
# =============================================================================
#  mbx_ancombc2_run.sh   —   Step 14 of the mbX Pro pipeline
# =============================================================================
#  PURPOSE
#  -------
#  Differential abundance analysis with ANCOMBC2 (Lin & Peddada 2024,
#  Bioconductor) at multiple taxonomic levels for every categorical metadata
#  variable.  Always runs with pairwise = TRUE and global = TRUE so users get:
#      • global (reference-free) omnibus test  → which taxa differ across
#        ANY combination of groups
#      • pairwise contrasts (mdFDR-corrected)  → which specific groups differ
#  Reference level for each categorical variable = alphabetically first level
#  (R's default behaviour for factor()).
#
#  Reads from previous step's info file:
#      <mbX_pro_outputs_*>/13_beta_diversity_results/mbx_beta_diversity_info.txt
#  (also reads alpha + pre-diversity info files chained from there.)
#
#  Inputs:
#      • Filtered feature table  (from step 7 — mito/chloro removed)
#      • Taxonomy classification (from step 6)
#      • Metadata
#  ANCOMBC2 needs RAW counts, NOT rarefied — it does its own bias correction.
#
#  WHAT THIS SCRIPT PRODUCES
#  -------------------------
#  <mbX_pro_outputs_*>/
#  └── 14_differential_abundance_ANCOMBC2/
#      ├── ANCOMBC2_phylum/                        (and class, order, family,
#      │   ├── ancombc2_<Variable>/                 genus, species)
#      │   │   ├── ancombc2_primary_<Variable>.xlsx
#      │   │   ├── ancombc2_pairwise_<Variable>.xlsx
#      │   │   ├── ancombc2_global_<Variable>.xlsx
#      │   │   ├── ancombc2_structural_zeros_<Variable>.xlsx
#      │   │   ├── ancombc2_summary_<Variable>.xlsx
#      │   │   ├── volcano_pairwise_<Variable>.png
#      │   │   └── heatmap_significant_<Variable>.png
#      │   └── Summary_<level>_all_variables.xlsx
#      ├── working_dir_differential_abundance/
#      │   ├── collapsed_tables/                    table_collapsed_L<n>.qza
#      │   ├── exported_tables/L<n>/                feature-table.{biom,tsv}
#      │   ├── ancombc2_raw_objects/                <level>_<var>.rds
#      │   └── ancombc2_run_logs/                   <level>_<var>.log
#      ├── Summary_all_levels_all_variables.xlsx
#      └── mbx_ancombc2_info.txt
#
#  GATING
#  ------
#  Refuses to run if step 11 reported OVERALL_STATUS=FAIL or
#  READY_FOR_DIVERSITY=no, unless --force is passed.
#
#  Compatible with bash 3.2+ (macOS default shell).
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
  grep "^${key}=" "$file" 2>/dev/null | cut -d= -f2- | head -1
}

_sanitize() {
  printf '%s' "$1" \
    | tr ' ' '_' \
    | tr -d '()[]{}/<>|\\:*?"' \
    | sed 's/^[_.-]*//;s/[_.-]*$//'
}

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

mbx_ancombc2_run.sh — Differential abundance with ANCOMBC2 at every
                      taxonomic level for every categorical metadata variable.
                      Runs pairwise + global tests.  Reference = alphabetical
                      first level.

USAGE:
  mbx_ancombc2_run.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads <mbX_pro_outputs_dir>/13_beta_diversity_results/mbx_beta_diversity_info.txt
  Then for every (taxonomic level × categorical variable):
      1. Collapses feature table to that level (qiime taxa collapse).
      2. Builds a phyloseq object  (taxa = collapsed taxonomy strings).
      3. Runs ancombc2() with:
            fix_formula = <variable>
            group       = <variable>
            pairwise    = TRUE
            global      = TRUE
            reference   = alphabetical first level (R default)
      4. Writes 5 xlsx tables + 2 PNG plots + 1 RDS.

OPTIONS:
  --levels LIST          Comma-separated levels to run.
                         Choices: domain, phylum, class, order, family, genus, species
                         Default: phylum,class,order,family,genus,species
                         (Skip domain — typically only Bacteria/Archaea.)
  --variables LIST       Comma-separated categorical variables (must exist in metadata).
                         Default: ALL auto-detected categorical columns.
  --p-adjust METHOD      P-value adjustment for primary tests
                         Choices: holm, BH, BY, bonferroni  (default: holm)
  --alpha N              Significance threshold (default: 0.05)
  --prv-cut N            Prevalence filter (default: 0.10 = features in <10% of samples dropped)
  --lib-cut N            Library size filter (default: 1000 = samples with <1000 reads dropped)
  --permutations N       For mdFDR pairwise control (default: 100)
  --skip-deps-check      Don't check / install R packages (assume they're available)
  --skip-plots           Don't generate volcano + heatmap PNGs
  --force                Run even if step 11 OVERALL_STATUS != PASS / PASS_WITH_WARNINGS
  --force-rerun          Recompute everything (ignore existing outputs)
  --rscript PATH         Override Rscript location (default: smart-detect)
  --dry-run              Print commands; do not execute
  -h, --help             Show this help

EXAMPLES:
  mbx_ancombc2_run.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_ancombc2_run.sh /path/to/mbX_pro_outputs_20260417_121431 --levels genus,family
  mbx_ancombc2_run.sh /path/to/mbX_pro_outputs_20260417_121431 --variables Treatment

PROVENANCE:
  All paths come from 13_beta_diversity_results/mbx_beta_diversity_info.txt.
  Re-runs are safe — completed (level × variable) combinations are skipped.

DEPENDENCIES:
  System    : QIIME2 conda env  (active during invocation)
              + bash 3.2+, biom, awk, sed
  R         : auto-detects an Rscript with ANCOMBC pre-installed; otherwise
              installs ANCOMBC via BiocManager (requires Xcode CLT on macOS,
              build-essential on Linux, RTools on Windows).
  R packages: ANCOMBC, phyloseq, openxlsx, ggplot2, pheatmap, RColorBrewer,
              dplyr, tidyr, BiocManager.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
USER_LEVELS=""
USER_VARIABLES=""
P_ADJUST="holm"
ALPHA="0.05"
PRV_CUT="0.10"
LIB_CUT="1000"
B_PERMS="100"
SKIP_DEPS_CHECK=false
SKIP_PLOTS=false
FORCE_GATE=false
FORCE_RERUN=false
USER_RSCRIPT=""
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)         usage ;;
    --levels)          USER_LEVELS="${2:-}";    [[ -z "$USER_LEVELS" ]]    && err "--levels requires a value"; shift 2 ;;
    --variables)       USER_VARIABLES="${2:-}"; [[ -z "$USER_VARIABLES" ]] && err "--variables requires a value"; shift 2 ;;
    --p-adjust)        P_ADJUST="${2:-}";       [[ -z "$P_ADJUST" ]]       && err "--p-adjust requires a value"; shift 2 ;;
    --alpha)           ALPHA="${2:-}";          [[ -z "$ALPHA" ]]          && err "--alpha requires a value"; shift 2 ;;
    --prv-cut)         PRV_CUT="${2:-}";        [[ -z "$PRV_CUT" ]]        && err "--prv-cut requires a value"; shift 2 ;;
    --lib-cut)         LIB_CUT="${2:-}";        [[ -z "$LIB_CUT" ]]        && err "--lib-cut requires a value"; shift 2 ;;
    --permutations)    B_PERMS="${2:-}";        case "$B_PERMS" in ''|*[!0-9]*) err "--permutations must be a positive integer" ;; esac; shift 2 ;;
    --skip-deps-check) SKIP_DEPS_CHECK=true; shift ;;
    --skip-plots)      SKIP_PLOTS=true; shift ;;
    --force)           FORCE_GATE=true; shift ;;
    --force-rerun)     FORCE_RERUN=true; shift ;;
    --rscript)         USER_RSCRIPT="${2:-}";   [[ -z "$USER_RSCRIPT" ]]   && err "--rscript requires a value"; shift 2 ;;
    --dry-run)         DRY_RUN=true; shift ;;
    -*)                err "Unknown option: '$1'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '$1'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -d "$MBX_OUT_DIR" ]] || err "Directory does not exist: '$MBX_OUT_DIR'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

case "$(basename "$MBX_OUT_DIR")" in
  mbX_pro_outputs_*|mbx_pro_outputs_*) : ;;
  *) warn "Directory name does not match 'mbX_pro_outputs_*': $(basename "$MBX_OUT_DIR")
        Continuing anyway." ;;
esac

case "$P_ADJUST" in
  holm|BH|BY|bonferroni) : ;;
  *) err "--p-adjust must be one of: holm, BH, BY, bonferroni  (got: '$P_ADJUST')" ;;
esac

# ── CPU detection ─────────────────────────────────────────────────────────────
if   command -v nproc   &>/dev/null; then N_JOBS="$(nproc)"
elif command -v sysctl  &>/dev/null; then N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else                                       N_JOBS=1
fi

# ── Detect platform (Apple Silicon / Intel Mac / Linux / Windows) ────────────
_PLATFORM_OS="$(uname -s 2>/dev/null || echo Unknown)"
_PLATFORM_ARCH="$(uname -m 2>/dev/null || echo unknown)"
case "$_PLATFORM_OS" in
  Darwin)
    case "$_PLATFORM_ARCH" in
      arm64|aarch64) PLATFORM_LABEL="macOS Apple Silicon" ;;
      x86_64)        PLATFORM_LABEL="macOS Intel" ;;
      *)             PLATFORM_LABEL="macOS ($_PLATFORM_ARCH)" ;;
    esac ;;
  Linux)             PLATFORM_LABEL="Linux ($_PLATFORM_ARCH)" ;;
  MINGW*|MSYS*|CYGWIN*) PLATFORM_LABEL="Windows ($_PLATFORM_OS)" ;;
  *)                 PLATFORM_LABEL="$_PLATFORM_OS ($_PLATFORM_ARCH)" ;;
esac

# ── R / conda-env interaction (CRITICAL — read carefully) ────────────────────
# This script needs two DIFFERENT R installations to potentially work together:
#
#   1. System R (e.g. /usr/local/bin/Rscript)
#      Default for the rest of the pipeline. Compiled extensions = .so (Apple
#      Silicon native). May or may not have ANCOMBC2 — install via BiocManager
#      if missing.
#
#   2. Conda env's R ($CONDA_PREFIX/lib/R)
#      Used INTERNALLY by qiime — q2-composition imports rpy2, which boots an
#      embedded R. Conda's R has compiled extensions named .dylib, AND has
#      ANCOMBC pre-installed in the qiime2-amplicon-2025.4 env.  Smart-detect
#      will prefer this R for ANCOMBC2 since it avoids a 10-20 min install.
#
# When you `conda activate qiime2-amplicon-2025.4`, conda exports
# R_LIBS_USER=$CONDA_PREFIX/lib/R/library.  This poisons the system R: it tries
# to load methods.dylib from a conda lib path that doesn't have it (with .so
# instead).  Result:  "shared object 'methods.dylib' not found".
#
# We solve this with TWO wrappers used everywhere:
#   • _R     : strips every R_* env var so SYSTEM R uses its native libPaths
#   • _QIIME : sets R_HOME to CONDA's R so qiime's embedded R finds .dylib
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

_CONDA_R_HOME=""
if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib/R" ]]; then
  _CONDA_R_HOME="${CONDA_PREFIX}/lib/R"
fi

# Helper: invoke an Rscript with all R env vars stripped
_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}

# _QIIME wraps `qiime` so its rpy2-embedded R points at conda's R installation
_QIIME() {
  if [[ -n "$_CONDA_R_HOME" ]]; then
    R_HOME="$_CONDA_R_HOME" \
    R_LIBS_USER="$_CONDA_R_HOME/library/" \
      qiime "$@"
  else
    qiime "$@"
  fi
}

# ── R dependency manifest (used by both probe and installer) ──────────────────
# All packages needed at runtime, partitioned by repo so the installer can use
# the right install function (install.packages vs BiocManager::install).
REQUIRED_PKGS="ANCOMBC phyloseq openxlsx ggplot2 pheatmap microbiome mia TreeSummarizedExperiment dplyr tidyr RColorBrewer patchwork"
CRAN_ONLY_PKGS="openxlsx ggplot2 pheatmap dplyr tidyr RColorBrewer patchwork"
BIOC_PKGS="ANCOMBC phyloseq microbiome mia TreeSummarizedExperiment"

# Probe: returns space-separated list of missing packages on stdout
_probe_missing() {
  local rcmd="$1"
  local rhome="${2:-}"
  local probe="
req <- strsplit('${REQUIRED_PKGS}', ' ', fixed=TRUE)[[1]]
miss <- req[!sapply(req, requireNamespace, quietly=TRUE)]
cat(paste(miss, collapse=' '))
"
  if [[ -n "$rhome" ]]; then
    R_HOME="$rhome" R_LIBS_USER="$rhome/library/" \
      "$rcmd" --vanilla -e "$probe" 2>/dev/null
  else
    _strip_env "$rcmd" --vanilla -e "$probe" 2>/dev/null
  fi
}

# Probe: just check if Rscript runs
_probe_basic() {
  "$1" --version >/dev/null 2>&1
}

# Smart Rscript locator — ranks candidates by # missing packages.
# Tie-breaker: earlier in the list (conda R is listed first because it has
# the heavy Bioc deps like 'microbiome' pre-installed in the qiime2 env).
RSCRIPT_CMD=""
RSCRIPT_USES_CONDA_R=false   # if true, must invoke with conda's R_HOME
RSCRIPT_NEEDS_INSTALL=false  # if true, must install missing packages
MISSING_PKGS_LIST=""

# Build candidate list (parallel arrays — bash 3.2)
CAND_RCMDS=()
CAND_RHOMES=()
_add_cand() {
  local rcmd="$1" rhome="${2:-}"
  [[ -n "$rcmd" && -x "$rcmd" ]] || return 0
  if [[ ${#CAND_RCMDS[@]} -gt 0 ]]; then
    for existing in "${CAND_RCMDS[@]}"; do
      [[ "$existing" == "$rcmd" ]] && return 0   # dedupe
    done
  fi
  CAND_RCMDS+=( "$rcmd" )
  CAND_RHOMES+=( "$rhome" )
}

if [[ -n "$USER_RSCRIPT" ]]; then
  [[ -x "$USER_RSCRIPT" ]] || err "--rscript path is not executable: '$USER_RSCRIPT'"
  _add_cand "$USER_RSCRIPT" ""
else
  # Conda R first — ANCOMBC's heavy 'microbiome' dep is pre-installed in qiime2 env
  [[ -n "$_CONDA_R_HOME" ]] && _add_cand "${_CONDA_R_HOME}/bin/Rscript" "$_CONDA_R_HOME"
  _add_cand "/opt/homebrew/bin/Rscript" ""
  _add_cand "/usr/local/bin/Rscript" ""
  _add_cand "/Library/Frameworks/R.framework/Resources/bin/Rscript" ""
  _add_cand "/usr/bin/Rscript" ""
  _add_cand "$(command -v Rscript 2>/dev/null || true)" ""
fi

[[ ${#CAND_RCMDS[@]} -gt 0 ]] || err "No working Rscript found.

  Install R appropriate for your platform: $PLATFORM_LABEL
  • macOS (Apple Silicon/Intel):  brew install r
  • Linux (Debian/Ubuntu):        sudo apt install r-base
  • Linux (RHEL/Fedora):          sudo dnf install R
  • Windows:                      https://cran.r-project.org/

  Then re-run this script."

# Probe each candidate; rank by INSTALL COST not raw count.
# A missing CRAN package = 1 point  (binary install ~30 s).
# A missing Bioconductor package = 100 points  (compile from source, ~2-15 min each).
# A perfect candidate (0 missing) wins immediately.
BEST_IDX=-1
BEST_COST=99999
BEST_COUNT=0
BEST_MISSING=""
for i in "${!CAND_RCMDS[@]}"; do
  rcmd="${CAND_RCMDS[$i]}"
  rhome="${CAND_RHOMES[$i]}"
  _probe_basic "$rcmd" || continue
  miss="$(_probe_missing "$rcmd" "$rhome" || true)"
  count=0
  cost=0
  if [[ -n "$miss" ]]; then
    for _pkg in $miss; do
      count=$(( count + 1 ))
      _is_bioc=false
      for _b in $BIOC_PKGS; do
        [[ "$_pkg" == "$_b" ]] && { _is_bioc=true; break; }
      done
      if $_is_bioc; then cost=$(( cost + 100 ))
      else               cost=$(( cost + 1   )); fi
    done
  fi
  if [[ "$cost" -lt "$BEST_COST" ]]; then
    BEST_IDX=$i
    BEST_COST=$cost
    BEST_COUNT=$count
    BEST_MISSING="$miss"
  fi
  if [[ "$cost" -eq 0 ]]; then break; fi   # perfect candidate
done

[[ $BEST_IDX -ge 0 ]] || err "No working Rscript found among candidates."

RSCRIPT_CMD="${CAND_RCMDS[$BEST_IDX]}"
[[ -n "${CAND_RHOMES[$BEST_IDX]}" ]] && RSCRIPT_USES_CONDA_R=true
[[ "$BEST_COST" -gt 0 ]] && RSCRIPT_NEEDS_INSTALL=true
BEST_MISSING_COUNT="$BEST_COUNT"
MISSING_PKGS_LIST="$BEST_MISSING"

# Dispatcher: invoke the chosen Rscript with the right env
_R() {
  if $RSCRIPT_USES_CONDA_R; then
    R_HOME="$_CONDA_R_HOME" \
    R_LIBS_USER="$_CONDA_R_HOME/library/" \
      "$RSCRIPT_CMD" "$@"
  else
    _strip_env "$RSCRIPT_CMD" "$@"
  fi
}

R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"

# ── PID-based temp file IDs (avoids mktemp suffix issues on macOS bash 3.2) ───
_TMPID="${$}_$(date +%s)"
trap 'rm -f /tmp/mbx_ancombc2_${_TMPID}*.R /tmp/mbx_ancombc2_${_TMPID}*.txt' EXIT

# =============================================================================
step "1/9 — Read mbx_beta_diversity_info.txt and chained info files"
# =============================================================================
BETA_INFO="${MBX_OUT_DIR}/13_beta_diversity_results/mbx_beta_diversity_info.txt"
ALPHA_INFO="${MBX_OUT_DIR}/12_alpha_diversity_results/mbx_alpha_diversity_info.txt"
PRE_DIV_INFO="${MBX_OUT_DIR}/11_pre_diversity/mbx_pre_diversity_info.txt"
TAXONOMY_INFO="${MBX_OUT_DIR}/7_taxonomy_csv/mbx_taxonomy_info.txt"

[[ -f "$BETA_INFO" ]] || err "Beta-diversity info file not found:
    $BETA_INFO
  → Run mbx_beta_diversity_run.sh first (step 13)."
[[ -f "$PRE_DIV_INFO" ]] || err "Pre-diversity info file not found:
    $PRE_DIV_INFO
  → Run mbx_pre_diversity_parameters.sh first (step 11)."

OVERALL_STATUS="$(_read_key OVERALL_STATUS         "$PRE_DIV_INFO")"
READY="$(        _read_key READY_FOR_DIVERSITY    "$PRE_DIV_INFO")"
METADATA_TXT="$( _read_key METADATA_TXT           "$BETA_INFO")"

# Filtered (non-rarefied) feature table — comes from step 7 (taxonomy filtered).
# Beta info doesn't store it directly; pull from alpha info or taxonomy info,
# else default to known location.
FT_FILTERED_QZA=""
if [[ -f "$ALPHA_INFO" ]]; then
  FT_FILTERED_QZA="$(_read_key FEATURE_TABLE_QZA "$ALPHA_INFO")"
fi
if [[ -z "$FT_FILTERED_QZA" || ! -f "$FT_FILTERED_QZA" ]]; then
  FT_FILTERED_QZA="${MBX_OUT_DIR}/7_taxonomy_csv/feature_table_filtered.qza"
fi
[[ -f "$FT_FILTERED_QZA" ]] || err "Filtered feature table not found:
    $FT_FILTERED_QZA
  → Run mbx_taxonomy_run.sh first (step 7)."

# Taxonomy
TAXONOMY_QZA=""
if [[ -f "$TAXONOMY_INFO" ]]; then
  TAXONOMY_QZA="$(_read_key TAXONOMY_QZA "$TAXONOMY_INFO")"
fi
if [[ -z "$TAXONOMY_QZA" || ! -f "$TAXONOMY_QZA" ]]; then
  TAXONOMY_QZA="${MBX_OUT_DIR}/6_classifier_taxonomy/taxonomy.qza"
fi
[[ -f "$TAXONOMY_QZA" ]] || err "Taxonomy classification not found:
    $TAXONOMY_QZA
  → Run mbx_classifier_run.sh first (step 6)."

[[ -n "$METADATA_TXT" && -f "$METADATA_TXT" ]] || err "METADATA_TXT not usable: '$METADATA_TXT'"

# ── Gating ────────────────────────────────────────────────────────────────────
case "$OVERALL_STATUS" in
  PASS|PASS_WITH_WARNINGS) : ;;
  *)
    if $FORCE_GATE; then
      warn "Step 11 OVERALL_STATUS = '$OVERALL_STATUS' — proceeding because --force was given."
    else
      err "Step 11 OVERALL_STATUS = '$OVERALL_STATUS'.
  → Differential abundance on a non-PASS state may produce misleading results.
  → Re-run with --force to bypass." ; fi ;;
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
info "Platform             : $PLATFORM_LABEL"
info "QIIME2 env           : ${CONDA_DEFAULT_ENV:-unknown}"
if $RSCRIPT_USES_CONDA_R; then
  info "Rscript              : $RSCRIPT_CMD  (R $R_VERSION, conda env's R)"
else
  info "Rscript              : $RSCRIPT_CMD  (R $R_VERSION, system R)"
fi
if $RSCRIPT_NEEDS_INSTALL; then
  info "R deps status        : MISSING $BEST_MISSING_COUNT — will install: $MISSING_PKGS_LIST"
else
  info "R deps status        : ALL PRESENT"
fi
info "Output root          : $MBX_OUT_DIR"
info "Metadata             : $METADATA_TXT"
info "Filtered table (raw counts): $FT_FILTERED_QZA"
info "Taxonomy             : $TAXONOMY_QZA"
info "P-adjust method      : $P_ADJUST"
info "Alpha                : $ALPHA"
info "Prevalence cutoff    : $PRV_CUT"
info "Library cutoff       : $LIB_CUT"
info "Permutations (mdFDR) : $B_PERMS"
info "CPU cores            : $N_JOBS"
info "Step-11 status       : $OVERALL_STATUS  /  ready=$READY"
sep

# ── Set up output directories ─────────────────────────────────────────────────
ANCO_DIR="${MBX_OUT_DIR}/14_differential_abundance_ANCOMBC2"
WORK_DIR="${ANCO_DIR}/working_dir_differential_abundance"
COLLAPSE_DIR="${WORK_DIR}/collapsed_tables"
EXPORT_DIR="${WORK_DIR}/exported_tables"
RAW_DIR="${WORK_DIR}/ancombc2_raw_objects"
LOG_DIR="${WORK_DIR}/ancombc2_run_logs"
INFO_FILE="${ANCO_DIR}/mbx_ancombc2_info.txt"
GLOBAL_SUM_XLSX="${ANCO_DIR}/Summary_all_levels_all_variables.xlsx"

mkdir -p "$ANCO_DIR" "$WORK_DIR" "$COLLAPSE_DIR" "$EXPORT_DIR" "$RAW_DIR" "$LOG_DIR"

# ── Tax level tables (parallel arrays — bash 3.2 compatible) ─────────────────
ALL_LEVEL_NUMS=( "1"      "2"      "3"     "4"     "5"      "6"     "7"       )
ALL_LEVEL_NAMES=( "domain" "phylum" "class" "order" "family" "genus" "species" )

# Filter to user-selected levels (default: phylum,class,order,family,genus,species)
LEVEL_NUMS=()
LEVEL_NAMES=()
if [[ -z "$USER_LEVELS" ]]; then
  USER_LEVELS="phylum,class,order,family,genus,species"
fi
IFS=',' read -r -a _user_levels <<< "$USER_LEVELS"
for _u in "${_user_levels[@]}"; do
  _u_trim="$(echo "$_u" | tr -d ' ' | tr '[:upper:]' '[:lower:]')"
  _matched=false
  for i in "${!ALL_LEVEL_NAMES[@]}"; do
    if [[ "${ALL_LEVEL_NAMES[$i]}" == "$_u_trim" ]]; then
      LEVEL_NUMS+=( "${ALL_LEVEL_NUMS[$i]}" )
      LEVEL_NAMES+=( "${ALL_LEVEL_NAMES[$i]}" )
      _matched=true; break
    fi
  done
  $_matched || err "Unknown level in --levels: '$_u'
  Choices: domain, phylum, class, order, family, genus, species"
done
[[ ${#LEVEL_NUMS[@]} -eq 0 ]] && err "No valid levels in --levels='$USER_LEVELS'"
N_LEVELS=${#LEVEL_NUMS[@]}
info "Taxonomic levels     : ${LEVEL_NAMES[*]}"
sep

# =============================================================================
step "2/9 — Check / install R dependencies (ANCOMBC2 + phyloseq + openxlsx + ...)"
# =============================================================================
if $SKIP_DEPS_CHECK; then
  warn "Skipping dependency check (--skip-deps-check)."
elif ! $RSCRIPT_NEEDS_INSTALL; then
  ok "All R dependencies already available — skipping install."
else
  info "Missing R packages in $RSCRIPT_CMD : $MISSING_PKGS_LIST"
  echo ""

  # Partition missing packages into CRAN-only vs Bioconductor
  MISS_CRAN=""
  MISS_BIOC=""
  for pkg in $MISSING_PKGS_LIST; do
    case " $CRAN_ONLY_PKGS " in *" $pkg "*) MISS_CRAN="$MISS_CRAN $pkg" ;; esac
    case " $BIOC_PKGS "      in *" $pkg "*) MISS_BIOC="$MISS_BIOC $pkg" ;; esac
  done
  MISS_CRAN="$(echo "$MISS_CRAN" | sed 's/^ *//;s/ *$//')"
  MISS_BIOC="$(echo "$MISS_BIOC" | sed 's/^ *//;s/ *$//')"

  if [[ -n "$MISS_CRAN" ]]; then
    info "Will install from CRAN     : $MISS_CRAN  (fast)"
  fi
  if [[ -n "$MISS_BIOC" ]]; then
    info "Will install from Bioconductor : $MISS_BIOC  (slow on first run, ~5-15 min)"
  fi
  echo ""
  echo "  Platform: $PLATFORM_LABEL"

  # ── Platform-specific system tool checks (only matter for source compilation)
  case "$_PLATFORM_OS" in
    Darwin)
      if ! xcode-select -p >/dev/null 2>&1; then
        err "Xcode Command Line Tools required to compile R packages on macOS.
  Install with:
      xcode-select --install
  Then re-run this script."
      fi
      info "Xcode CLT : $(xcode-select -p)"
      ;;
    Linux)
      if ! command -v gcc >/dev/null 2>&1; then
        err "GCC required to compile R packages on Linux.
  Install with:
      sudo apt install build-essential          # Debian/Ubuntu
      sudo dnf groupinstall 'Development Tools' # RHEL/Fedora
  Then re-run this script."
      fi
      info "GCC : $(gcc --version 2>&1 | head -1)"
      ;;
    MINGW*|MSYS*|CYGWIN*)
      warn "Detected Windows ($_PLATFORM_OS).
  Make sure RTools matching your R version is installed:
      https://cran.r-project.org/bin/windows/Rtools/"
      ;;
  esac
  echo ""

  # Build R vectors of missing pkg names
  _to_r_vec() {
    local items="$1" out="c("
    local first=true
    for p in $items; do
      if $first; then out="${out}\"$p\""; first=false
      else            out="${out}, \"$p\""; fi
    done
    echo "${out})"
  }
  MISS_CRAN_R="$(_to_r_vec "$MISS_CRAN")"
  MISS_BIOC_R="$(_to_r_vec "$MISS_BIOC")"
  REQUIRED_R="$(_to_r_vec "$REQUIRED_PKGS")"

  R_INSTALL_SCRIPT="/tmp/mbx_ancombc2_${_TMPID}_install.R"
  cat > "$R_INSTALL_SCRIPT" << RINSTALL
options(repos = c(CRAN = "https://cloud.r-project.org"))
cat(sprintf("[INFO]  R version : %s\n", R.version.string))
cat(sprintf("[INFO]  Lib paths : %s\n", paste(.libPaths(), collapse="; ")))

cran_missing <- ${MISS_CRAN_R}
bioc_missing <- ${MISS_BIOC_R}

if (length(cran_missing) > 0) {
  cat(sprintf("[INFO]  Installing CRAN: %s\n", paste(cran_missing, collapse=", ")))
  install.packages(cran_missing, quiet = TRUE)
}

if (length(bioc_missing) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("[INFO]  Installing BiocManager...\n")
    install.packages("BiocManager", quiet = TRUE)
  }
  suppressMessages(library(BiocManager))
  cat(sprintf("[INFO]  BiocManager version: %s\n",
              as.character(BiocManager::version())))
  cat(sprintf("[INFO]  Installing Bioconductor: %s\n",
              paste(bioc_missing, collapse=", ")))
  cat("[INFO]  This downloads many packages — please be patient.\n")
  BiocManager::install(bioc_missing, update = FALSE, ask = FALSE, quiet = TRUE)
}

# Final verification — re-check FULL required list, not just what we tried to install
required <- ${REQUIRED_R}
still_missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(still_missing) > 0) {
  cat(sprintf("[ERROR] Still missing after install: %s\n",
              paste(still_missing, collapse=", ")), file = stderr())
  quit(status = 1)
}
cat("\n[OK]    All R dependencies present.\n")
RINSTALL

  if $DRY_RUN; then
    warn "[DRY-RUN] Would install R packages."
  else
    timer_start
    _R --vanilla "$R_INSTALL_SCRIPT" \
      || err "R dependency installation failed.

  Common causes on $PLATFORM_LABEL:
   • Missing system libraries (e.g. libxml2-dev, libssl-dev, gfortran)
   • Network connectivity to CRAN / Bioconductor
   • Insufficient disk space (Bioc 'microbiome' install is ~500 MB)

  Try installing manually:
      $RSCRIPT_CMD -e 'install.packages(\"BiocManager\"); BiocManager::install(c(\"ANCOMBC\",\"phyloseq\",\"microbiome\"))'

  Or re-run with --rscript pointing to an Rscript that already has ANCOMBC2."
    timer_end
  fi
  rm -f "$R_INSTALL_SCRIPT"
  ok "R dependencies installed."
fi
sep

# =============================================================================
step "3/9 — Detect categorical metadata columns"
# =============================================================================
CATS_FILE="/tmp/mbx_ancombc2_${_TMPID}_cats.txt"
R_DETECT="/tmp/mbx_ancombc2_${_TMPID}_detect.R"

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
  warn "[DRY-RUN] Would detect categorical columns."
  AUTODETECTED_COLS=( "Treatment" "SampleType" )
else
  _R --vanilla "$R_DETECT" "$METADATA_TXT" > "$CATS_FILE" 2>&1 \
    || err "Categorical column detection failed.
  → Verify metadata file: $METADATA_TXT"
  grep "^\[" "$CATS_FILE" || true
  AUTODETECTED_COLS=()
  while IFS= read -r _col; do
    [[ -n "$_col" ]] && AUTODETECTED_COLS+=("$_col")
  done < <(grep -v "^\[" "$CATS_FILE" | grep -v '^$' || true)
  [[ ${#AUTODETECTED_COLS[@]} -eq 0 ]] && err "No categorical columns detected.
  → Verify your metadata has at least one grouping column."
fi
rm -f "$R_DETECT"

# Filter to user-selected variables, if any
CATEGORICAL_COLS=()
if [[ -z "$USER_VARIABLES" ]]; then
  CATEGORICAL_COLS=( "${AUTODETECTED_COLS[@]}" )
else
  IFS=',' read -r -a _uv <<< "$USER_VARIABLES"
  for _v in "${_uv[@]}"; do
    _v_trim="$(echo "$_v" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    [[ -z "$_v_trim" ]] && continue
    _found=false
    for col in "${AUTODETECTED_COLS[@]}"; do
      [[ "$col" == "$_v_trim" ]] && { CATEGORICAL_COLS+=( "$col" ); _found=true; break; }
    done
    $_found || err "Variable '$_v_trim' not found among auto-detected categorical columns: ${AUTODETECTED_COLS[*]}
  Use --variables with a subset of these (case-sensitive)."
  done
  [[ ${#CATEGORICAL_COLS[@]} -eq 0 ]] && err "No valid variables in --variables='$USER_VARIABLES'"
fi
ok "Variables to test: ${CATEGORICAL_COLS[*]}"

SANITIZED_COLS=()
for col in "${CATEGORICAL_COLS[@]}"; do
  SANITIZED_COLS+=( "$(_sanitize "$col")" )
done

# Set up per-level output directories
i=0
while [[ $i -lt $N_LEVELS ]]; do
  LNAME="${LEVEL_NAMES[$i]}"
  mkdir -p "${ANCO_DIR}/ANCOMBC2_${LNAME}"
  i=$(( i + 1 ))
done
sep

# =============================================================================
step "4/9 — Collapse feature table by taxonomy at $N_LEVELS levels"
# =============================================================================
i=0
while [[ $i -lt $N_LEVELS ]]; do
  LNUM="${LEVEL_NUMS[$i]}"
  LNAME="${LEVEL_NAMES[$i]}"
  OUT_QZA="${COLLAPSE_DIR}/table_collapsed_L${LNUM}.qza"

  if [[ -f "$OUT_QZA" && "$FORCE_RERUN" == false ]]; then
    skipped "[$((i+1))/$N_LEVELS] L${LNUM} ($LNAME) — already collapsed"
    i=$(( i + 1 )); continue
  fi

  info "[$((i+1))/$N_LEVELS] Collapsing to L${LNUM} ($LNAME)..."
  cmd_show "qiime taxa collapse" \
    "--i-table $FT_FILTERED_QZA" \
    "--i-taxonomy $TAXONOMY_QZA" \
    "--p-level $LNUM" \
    "--o-collapsed-table $OUT_QZA"

  if ! $DRY_RUN; then
    timer_start
    _QIIME taxa collapse \
      --i-table             "$FT_FILTERED_QZA" \
      --i-taxonomy          "$TAXONOMY_QZA" \
      --p-level             "$LNUM" \
      --o-collapsed-table   "$OUT_QZA" \
      || err "qiime taxa collapse failed at L${LNUM}.
  → Make sure the taxonomy contains every feature in the filtered table.
  → If they were generated at different times, re-run step 7 (mbx_taxonomy_run.sh)."
    timer_end
  fi
  ok "  → $(basename "$OUT_QZA")"
  i=$(( i + 1 ))
done
sep

# =============================================================================
step "5/9 — Export collapsed tables to TSV (via biom convert)"
# =============================================================================
i=0
while [[ $i -lt $N_LEVELS ]]; do
  LNUM="${LEVEL_NUMS[$i]}"
  LNAME="${LEVEL_NAMES[$i]}"
  IN_QZA="${COLLAPSE_DIR}/table_collapsed_L${LNUM}.qza"
  OUT_DIR="${EXPORT_DIR}/L${LNUM}"
  BIOM_FILE="${OUT_DIR}/feature-table.biom"
  TSV_FILE="${OUT_DIR}/feature-table.tsv"

  if [[ -f "$TSV_FILE" && "$FORCE_RERUN" == false ]]; then
    skipped "[$((i+1))/$N_LEVELS] L${LNUM} TSV exists"
    i=$(( i + 1 )); continue
  fi

  info "[$((i+1))/$N_LEVELS] Exporting L${LNUM} → TSV..."
  rm -rf "$OUT_DIR"; mkdir -p "$OUT_DIR"

  if ! $DRY_RUN; then
    _QIIME tools export \
      --input-path  "$IN_QZA" \
      --output-path "$OUT_DIR" \
      >/dev/null \
      || err "Failed to export $IN_QZA"
    [[ -f "$BIOM_FILE" ]] || err "Expected biom file not found: $BIOM_FILE"

    # biom convert is part of biom-format (installed in qiime2 env)
    if command -v biom >/dev/null 2>&1; then
      biom convert -i "$BIOM_FILE" -o "$TSV_FILE" --to-tsv \
        || err "biom convert failed for L${LNUM}"
    else
      err "'biom' command not found. Activate the qiime2 conda env first."
    fi
    [[ -f "$TSV_FILE" ]] || err "TSV not produced: $TSV_FILE"
  fi
  ok "  → $(basename "$TSV_FILE")"
  i=$(( i + 1 ))
done
sep

# =============================================================================
step "6/9 — Run ANCOMBC2  ($N_LEVELS levels × ${#CATEGORICAL_COLS[@]} variables)"
# =============================================================================
R_ANCO="/tmp/mbx_ancombc2_${_TMPID}_analysis.R"

# Build comma-separated quoted lists for R
VARS_R=""
for col in "${CATEGORICAL_COLS[@]}"; do
  _esc="${col//\"/\\\"}"
  if [[ -z "$VARS_R" ]]; then VARS_R="\"$_esc\""
  else                        VARS_R="$VARS_R, \"$_esc\""; fi
done
LEVEL_NUMS_R=""
LEVEL_NAMES_R=""
i=0
while [[ $i -lt $N_LEVELS ]]; do
  if [[ -z "$LEVEL_NUMS_R" ]]; then
    LEVEL_NUMS_R="${LEVEL_NUMS[$i]}"
    LEVEL_NAMES_R="\"${LEVEL_NAMES[$i]}\""
  else
    LEVEL_NUMS_R="${LEVEL_NUMS_R}, ${LEVEL_NUMS[$i]}"
    LEVEL_NAMES_R="${LEVEL_NAMES_R}, \"${LEVEL_NAMES[$i]}\""
  fi
  i=$(( i + 1 ))
done

if $SKIP_PLOTS; then SKIP_PLOTS_R="TRUE"; else SKIP_PLOTS_R="FALSE"; fi
if $FORCE_RERUN; then FORCE_RERUN_R="TRUE"; else FORCE_RERUN_R="FALSE"; fi

cat > "$R_ANCO" << RANCO
# =============================================================================
# ANCOMBC2 differential abundance analysis  (per level × variable)
# =============================================================================
suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Inputs ──────────────────────────────────────────────────────────────────
METADATA       <- "${METADATA_TXT}"
EXPORT_DIR     <- "${EXPORT_DIR}"
ANCO_DIR       <- "${ANCO_DIR}"
RAW_DIR        <- "${RAW_DIR}"
LOG_DIR        <- "${LOG_DIR}"
GLOBAL_SUM_XLSX<- "${GLOBAL_SUM_XLSX}"
P_ADJUST       <- "${P_ADJUST}"
ALPHA          <- ${ALPHA}
PRV_CUT        <- ${PRV_CUT}
LIB_CUT        <- ${LIB_CUT}
B_PERMS        <- ${B_PERMS}
N_JOBS         <- ${N_JOBS}
SKIP_PLOTS     <- ${SKIP_PLOTS_R}
FORCE_RERUN    <- ${FORCE_RERUN_R}
VARS           <- c(${VARS_R})
LEVEL_NUMS     <- c(${LEVEL_NUMS_R})
LEVEL_NAMES    <- c(${LEVEL_NAMES_R})

set.seed(42)

# ── Helpers ─────────────────────────────────────────────────────────────────
sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^[._-]+|[._-]+\$", "", x)
  x
}

read_metadata <- function(path) {
  ext <- tolower(tools::file_ext(path))
  md  <- if (ext == "csv") {
    read.csv(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  } else {
    read.delim(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE,
               na.strings=c("","NA","N/A","na","n/a"))
  }
  q2 <- grepl("^#", md[[1]])
  if (any(q2)) md <- md[!q2, , drop=FALSE]
  names(md)[1] <- "sample-id"
  md[["sample-id"]] <- trimws(as.character(md[["sample-id"]]))
  md
}

# Read taxa-collapsed feature table from biom-converted TSV.
# TSV format (biom convert --to-tsv):
#   # Constructed from biom file
#   #OTU ID  Sample1  Sample2  ...
#   k__Bact;p__Firmicutes  100  200  ...
read_collapsed_tsv <- function(tsv_path) {
  if (!file.exists(tsv_path)) stop(sprintf("TSV not found: %s", tsv_path))
  raw <- readLines(tsv_path)
  hdr_idx <- grep("^#OTU", raw, ignore.case = FALSE)
  if (length(hdr_idx) == 0) hdr_idx <- grep("^#?[Ss]ample[ -]?[Ii][Dd]", raw)
  if (length(hdr_idx) == 0) hdr_idx <- which(!grepl("^#", raw))[1] - 1
  if (length(hdr_idx) == 0 || hdr_idx < 1)
    stop(sprintf("Could not locate header row in %s", tsv_path))
  hdr <- strsplit(sub("^#", "", raw[hdr_idx]), "\t", fixed=TRUE)[[1]]
  body_lines <- raw[-seq_len(hdr_idx)]
  body_lines <- body_lines[nchar(body_lines) > 0 & !grepl("^#", body_lines)]
  body <- do.call(rbind, lapply(strsplit(body_lines, "\t", fixed=TRUE), function(x) {
    if (length(x) < length(hdr)) x <- c(x, rep("", length(hdr) - length(x)))
    x[seq_along(hdr)]
  }))
  colnames(body) <- hdr
  body <- as.data.frame(body, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(body) <- body[[1]]
  body[[1]] <- NULL
  for (j in seq_len(ncol(body))) body[[j]] <- as.numeric(body[[j]])
  body <- as.matrix(body)
  body[is.na(body)] <- 0
  body
}

cat(sprintf("[INFO]  ANCOMBC version: %s\n", as.character(packageVersion("ANCOMBC"))))
cat(sprintf("[INFO]  Reading metadata: %s\n", METADATA))
md <- read_metadata(METADATA)

# Bookkeeping for global summary
global_summary_rows <- list()
fail_log <- character(0)

# ───────────────────────────────────────────────────────────────────────────
# Loop over (level × variable)
# ───────────────────────────────────────────────────────────────────────────
n_combos <- length(LEVEL_NAMES) * length(VARS)
combo_i <- 0L

for (li in seq_along(LEVEL_NAMES)) {
  lnum  <- LEVEL_NUMS[li]
  lname <- LEVEL_NAMES[li]
  level_dir <- file.path(ANCO_DIR, sprintf("ANCOMBC2_%s", lname))
  dir.create(level_dir, recursive = TRUE, showWarnings = FALSE)

  tsv_path <- file.path(EXPORT_DIR, sprintf("L%d", lnum), "feature-table.tsv")
  cat(sprintf("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))
  cat(sprintf("  Level %d  (%s)   →   %s/\n", lnum, lname, level_dir))
  cat(sprintf("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"))

  otu_mat <- tryCatch(read_collapsed_tsv(tsv_path),
                      error = function(e) {
                        fail_log <<- c(fail_log,
                          sprintf("L%d (%s): could not read TSV — %s",
                                  lnum, lname, conditionMessage(e)))
                        NULL
                      })
  if (is.null(otu_mat)) next

  cat(sprintf("[INFO]  L%d table: %d taxa × %d samples\n",
              lnum, nrow(otu_mat), ncol(otu_mat)))

  level_summary_rows <- list()

  for (var in VARS) {
    combo_i <- combo_i + 1L
    safe_var <- sanitize(var)
    out_dir  <- file.path(level_dir, sprintf("ancombc2_%s", safe_var))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    rds_path <- file.path(RAW_DIR, sprintf("L%d_%s.rds", lnum, safe_var))
    log_path <- file.path(LOG_DIR, sprintf("L%d_%s.log", lnum, safe_var))

    primary_xlsx  <- file.path(out_dir, sprintf("ancombc2_primary_%s.xlsx", safe_var))
    pairwise_xlsx <- file.path(out_dir, sprintf("ancombc2_pairwise_%s.xlsx", safe_var))
    global_xlsx   <- file.path(out_dir, sprintf("ancombc2_global_%s.xlsx", safe_var))
    sszero_xlsx   <- file.path(out_dir, sprintf("ancombc2_structural_zeros_%s.xlsx", safe_var))
    summary_xlsx  <- file.path(out_dir, sprintf("ancombc2_summary_%s.xlsx", safe_var))
    volcano_png   <- file.path(out_dir, sprintf("volcano_pairwise_%s.png", safe_var))
    heatmap_png   <- file.path(out_dir, sprintf("heatmap_significant_%s.png", safe_var))

    cat(sprintf("\n[%d/%d]  L%d (%s)  ×  %s\n",
                combo_i, n_combos, lnum, lname, var))

    # Idempotency: skip if all main outputs exist
    main_outputs <- c(primary_xlsx, pairwise_xlsx, global_xlsx, summary_xlsx)
    if (!FORCE_RERUN && all(file.exists(main_outputs))) {
      cat(sprintf("[SKIP]  L%d × %s — all xlsx outputs exist (use --force-rerun to recompute)\n",
                  lnum, var))
      # Still need to load summary into level_summary_rows for level summary
      if (file.exists(summary_xlsx)) {
        s <- tryCatch(read.xlsx(summary_xlsx), error = function(e) NULL)
        if (!is.null(s)) {
          s\$Level <- lname
          level_summary_rows[[length(level_summary_rows)+1]] <- s
          global_summary_rows[[length(global_summary_rows)+1]] <- s
        }
      }
      next
    }

    # ── Filter metadata to samples present in the OTU table and to non-NA var ─
    samp_in_otu <- colnames(otu_mat)
    md_var <- md[md[["sample-id"]] %in% samp_in_otu, , drop = FALSE]
    md_var <- md_var[!is.na(md_var[[var]]) &
                     nchar(trimws(as.character(md_var[[var]]))) > 0, , drop = FALSE]
    md_var[[var]] <- factor(trimws(as.character(md_var[[var]])),
                            levels = sort(unique(trimws(as.character(md_var[[var]])))))

    # Drop singleton groups (n<2 won't fit ANCOMBC2 model)
    grp_n <- table(md_var[[var]])
    bad   <- names(grp_n)[grp_n < 2]
    if (length(bad) > 0) {
      md_var <- md_var[!md_var[[var]] %in% bad, , drop = FALSE]
      md_var[[var]] <- droplevels(md_var[[var]])
    }
    if (length(unique(md_var[[var]])) < 2) {
      msg <- sprintf("  <2 valid groups after filtering — skipping %s × %s.", var, lname)
      cat(sprintf("[WARN]  %s\n", msg))
      fail_log <<- c(fail_log, sprintf("L%d × %s: <2 valid groups", lnum, var))
      next
    }
    small_groups <- names(grp_n[grp_n >= 2 & grp_n < 5])
    if (length(small_groups) > 0)
      cat(sprintf("[WARN]  L%d × %s: groups n<5 (%s) — ANCOMBC2 results may be unstable.\n",
                  lnum, var, paste(small_groups, collapse=",")))

    samp_keep <- md_var[["sample-id"]]
    otu_sub   <- otu_mat[, samp_keep, drop = FALSE]
    cat(sprintf("[INFO]  Samples kept: %d  |  Groups: %s\n",
                ncol(otu_sub),
                paste(sprintf("%s(n=%d)", names(table(md_var[[var]])),
                              as.integer(table(md_var[[var]]))),
                      collapse=", ")))

    # ── Drop zero-variance taxa (would crash ANCOMBC2's bias-correction) ─────
    # ANCOMBC2 errors hard on taxa whose variance is zero WITHIN ANY GROUP
    # (e.g., all-zero counts in one group and constant non-zero in another).
    # Its built-in prv_cut = 0.10 only filters by overall presence, so we
    # must add this stronger filter ourselves.
    # NOTE: must use stats::var explicitly because the loop variable 'var'
    # holds the metadata column name as a string and shadows R's var().
    grp_vec   <- as.character(md_var[[var]])
    grp_var_mat <- sapply(unique(grp_vec), function(g) {
      sub_cols <- which(grp_vec == g)
      if (length(sub_cols) < 2) return(rep(NA_real_, nrow(otu_sub)))
      apply(otu_sub[, sub_cols, drop = FALSE], 1, stats::var)
    })
    if (is.null(dim(grp_var_mat))) grp_var_mat <- matrix(grp_var_mat, nrow = nrow(otu_sub))
    min_within <- suppressWarnings(apply(grp_var_mat, 1, min, na.rm = TRUE))
    bad_taxa   <- which(is.na(min_within) | !is.finite(min_within) | min_within == 0)
    if (length(bad_taxa) > 0) {
      cat(sprintf("[INFO]  Dropping %d zero-within-group-variance taxa (would crash ANCOMBC2)\n",
                  length(bad_taxa)))
      otu_sub <- otu_sub[-bad_taxa, , drop = FALSE]
    }
    if (nrow(otu_sub) < 5) {
      msg <- sprintf("Only %d taxa with non-zero within-group variance — not enough for ANCOMBC2",
                     nrow(otu_sub))
      cat(sprintf("[WARN]  %s\n", msg))
      fail_log <<- c(fail_log, sprintf("L%d × %s: %s", lnum, var, msg))
      next
    }

    # Build phyloseq object (no tax_table — taxa = collapsed taxonomy strings)
    sam_df <- as.data.frame(md_var, check.names = FALSE)
    rownames(sam_df) <- sam_df[["sample-id"]]
    ps <- tryCatch(
      phyloseq::phyloseq(
        phyloseq::otu_table(otu_sub, taxa_are_rows = TRUE),
        phyloseq::sample_data(sam_df)),
      error = function(e) NULL)
    if (is.null(ps)) {
      fail_log <<- c(fail_log, sprintf("L%d × %s: phyloseq build failed", lnum, var))
      cat(sprintf("[FAIL]  Could not build phyloseq object for L%d × %s.\n",
                  lnum, var))
      next
    }

    # ── Run ANCOMBC2 ─────────────────────────────────────────────────────────
    cat(sprintf("[RUN]   ancombc2(fix_formula = '%s', group = '%s', pairwise = TRUE, global = TRUE) ...\n",
                var, var))
    t0 <- Sys.time()
    out <- tryCatch(
      withCallingHandlers(
        ANCOMBC::ancombc2(
          data            = ps,
          fix_formula     = var,
          group           = var,
          p_adj_method    = P_ADJUST,
          pseudo_sens     = TRUE,
          prv_cut         = PRV_CUT,
          lib_cut         = LIB_CUT,
          struc_zero      = TRUE,
          neg_lb          = TRUE,
          alpha           = ALPHA,
          n_cl            = N_JOBS,
          verbose         = FALSE,
          global          = TRUE,
          pairwise        = TRUE,
          dunnet          = FALSE,
          trend           = FALSE,
          mdfdr_control   = list(fwer_ctrl_method = P_ADJUST,
                                 B = B_PERMS)),
        warning = function(w) {
          message(sprintf("[ANCOMBC WARN]  %s", conditionMessage(w)))
          invokeRestart("muffleWarning")
        }),
      error = function(e) {
        msg <- sprintf("L%d × %s: ANCOMBC2 error — %s", lnum, var,
                       conditionMessage(e))
        fail_log <<- c(fail_log, msg)
        cat(sprintf("[FAIL]  %s\n", msg))
        NULL
      })
    if (is.null(out)) next
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("[OK]    ancombc2 done in %.1f s\n", elapsed))

    # ── Save raw RDS for advanced users ──────────────────────────────────────
    saveRDS(out, file = rds_path)

    # ── Tidy each result table and write xlsx ────────────────────────────────
    res_primary <- out\$res        # main effect of each coefficient
    res_global  <- out\$res_global # global F-test
    res_pair    <- out\$res_pair   # all pairwise contrasts

    # 1. Primary results
    if (!is.null(res_primary) && nrow(res_primary) > 0) {
      pr <- as.data.frame(res_primary, check.names = FALSE)
      # Reorder so taxon is first
      if ("taxon" %in% names(pr))
        pr <- pr[, c("taxon", setdiff(names(pr), "taxon")), drop = FALSE]
      wb <- createWorkbook()
      addWorksheet(wb, sprintf("primary_%s", substr(safe_var, 1, 18)))
      writeData(wb, 1, pr,
                headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
      freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 2)
      setColWidths(wb, 1, cols = 1:ncol(pr), widths = "auto")
      saveWorkbook(wb, primary_xlsx, overwrite = TRUE)
      cat(sprintf("[OK]    primary  → %s  (%d taxa × %d cols)\n",
                  basename(primary_xlsx), nrow(pr), ncol(pr)))
    } else {
      cat("[WARN]  Primary results empty.\n")
    }

    # 2. Pairwise results
    n_pair_sig <- 0
    if (!is.null(res_pair) && nrow(res_pair) > 0) {
      pw <- as.data.frame(res_pair, check.names = FALSE)
      if ("taxon" %in% names(pw))
        pw <- pw[, c("taxon", setdiff(names(pw), "taxon")), drop = FALSE]
      wb <- createWorkbook()
      addWorksheet(wb, "pairwise")
      writeData(wb, 1, pw,
                headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
      freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 2)
      setColWidths(wb, 1, cols = 1:ncol(pw), widths = "auto")
      saveWorkbook(wb, pairwise_xlsx, overwrite = TRUE)
      diff_cols <- grep("^diff_", names(pw), value = TRUE)
      n_pair_sig <- if (length(diff_cols) > 0) {
        sum(rowSums(sapply(diff_cols, function(c) isTRUE(pw[[c]]) | pw[[c]] %in% TRUE), na.rm=TRUE) > 0)
      } else 0L
      cat(sprintf("[OK]    pairwise → %s  (%d sig taxa across %d contrasts)\n",
                  basename(pairwise_xlsx), n_pair_sig, length(diff_cols)))
    } else {
      cat("[WARN]  Pairwise results empty.\n")
    }

    # 3. Global results
    n_global_sig <- 0
    if (!is.null(res_global) && nrow(res_global) > 0) {
      gl <- as.data.frame(res_global, check.names = FALSE)
      if ("taxon" %in% names(gl))
        gl <- gl[, c("taxon", setdiff(names(gl), "taxon")), drop = FALSE]
      wb <- createWorkbook()
      addWorksheet(wb, "global")
      writeData(wb, 1, gl,
                headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
      freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 2)
      setColWidths(wb, 1, cols = 1:ncol(gl), widths = "auto")
      saveWorkbook(wb, global_xlsx, overwrite = TRUE)
      n_global_sig <- if ("diff_abn" %in% names(gl))
        sum(gl\$diff_abn %in% TRUE, na.rm = TRUE) else 0L
      cat(sprintf("[OK]    global   → %s  (%d sig taxa)\n",
                  basename(global_xlsx), n_global_sig))
    } else {
      cat("[WARN]  Global results empty.\n")
    }

    # 4. Structural zeros
    if (!is.null(out\$zero_ind) && nrow(out\$zero_ind) > 0) {
      sz <- as.data.frame(out\$zero_ind, check.names = FALSE)
      if ("taxon" %in% names(sz))
        sz <- sz[, c("taxon", setdiff(names(sz), "taxon")), drop = FALSE]
      sz_cols <- setdiff(names(sz), "taxon")
      sz_any  <- if (length(sz_cols) > 0)
        rowSums(sapply(sz_cols, function(c) sz[[c]] %in% TRUE), na.rm=TRUE) > 0 else FALSE
      sz_only <- sz[sz_any, , drop = FALSE]
      if (nrow(sz_only) > 0) {
        wb <- createWorkbook()
        addWorksheet(wb, "structural_zeros")
        writeData(wb, 1, sz_only,
                  headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
        setColWidths(wb, 1, cols = 1:ncol(sz_only), widths = "auto")
        saveWorkbook(wb, sszero_xlsx, overwrite = TRUE)
        cat(sprintf("[OK]    str.zero → %s  (%d taxa absent in >=1 group)\n",
                    basename(sszero_xlsx), nrow(sz_only)))
      } else {
        cat("[INFO]  No structural zeros — all taxa observed in every group.\n")
      }
    }

    # 5. Per-(level × var) summary xlsx
    summ_df <- data.frame(
      Level                    = lname,
      Variable                 = var,
      Reference_level          = levels(md_var[[var]])[1],
      All_groups               = paste(levels(md_var[[var]]), collapse=", "),
      n_samples                = ncol(otu_sub),
      n_groups                 = length(levels(md_var[[var]])),
      n_taxa_input             = nrow(otu_sub),
      n_taxa_after_prv_filter  = if (!is.null(out\$feature_table)) nrow(out\$feature_table) else NA_integer_,
      n_significant_global     = n_global_sig,
      n_significant_pairwise   = n_pair_sig,
      p_adjust_method          = P_ADJUST,
      alpha                    = ALPHA,
      prv_cut                  = PRV_CUT,
      lib_cut                  = LIB_CUT,
      runtime_sec              = round(elapsed, 1),
      ancombc_version          = as.character(packageVersion("ANCOMBC")),
      stringsAsFactors         = FALSE)
    wb <- createWorkbook()
    addWorksheet(wb, "summary")
    writeData(wb, 1, summ_df,
              headerStyle = createStyle(textDecoration = "bold"))
    setColWidths(wb, 1, cols = 1:ncol(summ_df), widths = "auto")
    saveWorkbook(wb, summary_xlsx, overwrite = TRUE)

    level_summary_rows[[length(level_summary_rows)+1]] <- summ_df
    global_summary_rows[[length(global_summary_rows)+1]] <- summ_df

    # ── Plots ────────────────────────────────────────────────────────────────
    if (!SKIP_PLOTS) {
      # Volcano: faceted across all pairwise comparisons
      if (!is.null(res_pair) && nrow(res_pair) > 0) {
        pw <- as.data.frame(res_pair, check.names = FALSE)
        # Find paired (lfc, q) column groups by stripping the "lfc_" prefix
        lfc_cols <- grep("^lfc_", names(pw), value = TRUE)
        contrasts <- sub("^lfc_", "", lfc_cols)

        long_rows <- list()
        for (cn in contrasts) {
          lfc_col <- paste0("lfc_", cn)
          q_col   <- paste0("q_",   cn)
          if (!(q_col %in% names(pw))) next
          long_rows[[cn]] <- data.frame(
            taxon       = pw\$taxon,
            contrast    = cn,
            lfc         = suppressWarnings(as.numeric(pw[[lfc_col]])),
            q           = suppressWarnings(as.numeric(pw[[q_col]])),
            stringsAsFactors = FALSE)
        }
        if (length(long_rows) > 0) {
          long_df <- do.call(rbind, long_rows)
          long_df\$neg_log10_q <- -log10(pmax(long_df\$q, 1e-300))
          long_df\$signif <- !is.na(long_df\$q) & long_df\$q < ALPHA &
                            !is.na(long_df\$lfc) & abs(long_df\$lfc) > 1

          ncols_facet <- max(1, min(3, ceiling(sqrt(length(unique(long_df\$contrast))))))
          p <- ggplot(long_df, aes(x = lfc, y = neg_log10_q,
                                   colour = signif)) +
            geom_point(size = 1.6, alpha = 0.7) +
            scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#cc1f1a"),
                                labels = c("FALSE" = "ns",
                                           "TRUE" = sprintf("q<%.2g & |lfc|>1", ALPHA))) +
            geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
                       colour = "grey40", linewidth = 0.4) +
            geom_vline(xintercept = c(-1, 1), linetype = "dashed",
                       colour = "grey40", linewidth = 0.4) +
            facet_wrap(~ contrast, ncol = ncols_facet, scales = "free") +
            labs(title    = sprintf("ANCOMBC2 volcano  --  L%d (%s)  by  %s",
                                    lnum, lname, var),
                 subtitle = sprintf("Reference: %s   |   p-adj: %s   |   alpha: %.2g",
                                    levels(md_var[[var]])[1], P_ADJUST, ALPHA),
                 x = "log fold change", y = expression(-log[10](q)),
                 colour = NULL) +
            theme_classic(base_size = 11) +
            theme(plot.title    = element_text(face = "bold"),
                  plot.subtitle = element_text(colour = "grey30"),
                  legend.position = "bottom",
                  strip.text    = element_text(face = "bold", size = 9))
          tryCatch({
            ggsave(volcano_png, p,
                   width  = max(8, 4 * ncols_facet),
                   height = max(5, 3.5 * ceiling(length(unique(long_df\$contrast)) / ncols_facet)),
                   dpi = 300, bg = "white")
            cat(sprintf("[OK]    volcano  → %s\n", basename(volcano_png)))
          }, error = function(e)
            cat(sprintf("[WARN]  Volcano save failed: %s\n", conditionMessage(e))))
        }
      }

      # Heatmap of significant taxa (any test)
      tryCatch({
        sig_taxa <- character(0)
        if (!is.null(res_global) && "diff_abn" %in% names(res_global)) {
          sig_taxa <- c(sig_taxa,
                        as.character(res_global\$taxon[res_global\$diff_abn %in% TRUE]))
        }
        if (!is.null(res_pair)) {
          diff_cols_pw <- grep("^diff_", names(res_pair), value = TRUE)
          if (length(diff_cols_pw) > 0) {
            mask <- rowSums(sapply(diff_cols_pw,
                                   function(c) res_pair[[c]] %in% TRUE),
                            na.rm = TRUE) > 0
            sig_taxa <- c(sig_taxa, as.character(res_pair\$taxon[mask]))
          }
        }
        sig_taxa <- unique(sig_taxa)
        sig_taxa <- intersect(sig_taxa, rownames(otu_sub))
        if (length(sig_taxa) >= 2) {
          # Use bias-corrected log abundance if available, else log1p of counts
          mat <- if (!is.null(out\$bias_correct_log_table)) {
            m <- as.matrix(out\$bias_correct_log_table)
            common <- intersect(sig_taxa, rownames(m))
            m[common, , drop = FALSE]
          } else {
            log1p(otu_sub[sig_taxa, , drop = FALSE])
          }
          # Cap at top 60 taxa for readability (most variable)
          if (nrow(mat) > 60) {
            v <- apply(mat, 1, var, na.rm = TRUE)
            mat <- mat[order(-v)[1:60], , drop = FALSE]
          }
          # Sort samples by group
          samp_order <- order(md_var[[var]])
          mat <- mat[, samp_order, drop = FALSE]
          mat[is.na(mat)] <- 0

          # Truncate long taxon names for display
          short_names <- sapply(rownames(mat), function(s) {
            parts <- strsplit(s, ";", fixed = TRUE)[[1]]
            tail(parts, 2) |> paste(collapse = ";")
          })
          rownames(mat) <- short_names

          annot_col <- data.frame(grp = md_var[[var]][samp_order])
          names(annot_col) <- var
          rownames(annot_col) <- colnames(mat)

          n_lvl <- length(levels(annot_col[[var]]))
          pal <- if (n_lvl <= 8) brewer.pal(max(3, n_lvl), "Set2")[seq_len(n_lvl)]
                 else if (n_lvl <= 12) brewer.pal(n_lvl, "Set3")
                 else colorRampPalette(brewer.pal(8, "Set2"))(n_lvl)
          annot_colors <- list(setNames(pal, levels(annot_col[[var]])))
          names(annot_colors) <- var

          pheatmap(mat,
                   cluster_rows      = TRUE,
                   cluster_cols      = FALSE,
                   show_rownames     = TRUE,
                   show_colnames     = TRUE,
                   annotation_col    = annot_col,
                   annotation_colors = annot_colors,
                   color             = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
                   main              = sprintf("Significant taxa (any test) -- L%d (%s)  x  %s",
                                              lnum, lname, var),
                   fontsize_row      = 7,
                   fontsize_col      = 8,
                   filename          = heatmap_png,
                   width             = max(8, 0.25 * ncol(mat) + 5),
                   height            = max(6, 0.18 * nrow(mat) + 3))
          cat(sprintf("[OK]    heatmap  → %s\n", basename(heatmap_png)))
        } else {
          cat("[INFO]  <2 significant taxa — heatmap not generated.\n")
        }
      }, error = function(e)
        cat(sprintf("[WARN]  Heatmap failed: %s\n", conditionMessage(e))))
    }
  }   # end variable loop

  # Per-level summary
  if (length(level_summary_rows) > 0) {
    ldf <- do.call(rbind, level_summary_rows); rownames(ldf) <- NULL
    out_xlsx <- file.path(level_dir, sprintf("Summary_%s_all_variables.xlsx", lname))
    wb <- createWorkbook()
    addWorksheet(wb, sprintf("L%d_%s", lnum, substr(lname, 1, 24)))
    writeData(wb, 1, ldf,
              headerStyle = createStyle(textDecoration = "bold"))
    setColWidths(wb, 1, cols = 1:ncol(ldf), widths = "auto")
    saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    cat(sprintf("[OK]    Level summary → %s\n", out_xlsx))
  }
}   # end level loop

# ── Global summary across all levels × variables ─────────────────────────────
if (length(global_summary_rows) > 0) {
  gdf <- do.call(rbind, global_summary_rows); rownames(gdf) <- NULL
  wb <- createWorkbook()
  addWorksheet(wb, "all_combos")
  writeData(wb, 1, gdf,
            headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb, 1, firstActiveRow = 2)
  setColWidths(wb, 1, cols = 1:ncol(gdf), widths = "auto")
  saveWorkbook(wb, GLOBAL_SUM_XLSX, overwrite = TRUE)
  cat(sprintf("\n[OK]    Global summary → %s  (%d combinations)\n",
              GLOBAL_SUM_XLSX, nrow(gdf)))
}

# Failures
if (length(fail_log) > 0) {
  cat("\n[SUMMARY]  Failures during analysis:\n")
  for (f in fail_log) cat(sprintf("  - %s\n", f))
}
cat("\n[DONE]  ANCOMBC2 analysis complete.\n")
RANCO

if $DRY_RUN; then
  warn "[DRY-RUN] Would run ANCOMBC2 R analysis."
else
  _R --vanilla "$R_ANCO" \
    || err "R ANCOMBC2 analysis failed — see output above."
fi
rm -f "$R_ANCO"
sep

# =============================================================================
step "7/9 — Verify output files"
# =============================================================================
N_PRIMARY=0
N_PAIRWISE=0
N_GLOBAL=0
N_PNG=0
N_RDS=0
i=0
while [[ $i -lt $N_LEVELS ]]; do
  LNAME="${LEVEL_NAMES[$i]}"
  LDIR="${ANCO_DIR}/ANCOMBC2_${LNAME}"
  if [[ -d "$LDIR" ]]; then
    N_PRIMARY=$((  N_PRIMARY  + $(find "$LDIR" -name 'ancombc2_primary_*.xlsx'  | wc -l) ))
    N_PAIRWISE=$(( N_PAIRWISE + $(find "$LDIR" -name 'ancombc2_pairwise_*.xlsx' | wc -l) ))
    N_GLOBAL=$((   N_GLOBAL   + $(find "$LDIR" -name 'ancombc2_global_*.xlsx'   | wc -l) ))
    N_PNG=$((      N_PNG      + $(find "$LDIR" -name '*.png'                     | wc -l) ))
  fi
  i=$(( i + 1 ))
done
N_RDS=$(find "$RAW_DIR" -name '*.rds' 2>/dev/null | wc -l)

info "primary xlsx   : $N_PRIMARY"
info "pairwise xlsx  : $N_PAIRWISE"
info "global xlsx    : $N_GLOBAL"
info "PNG plots      : $N_PNG"
info "raw RDS objects: $N_RDS"
sep

# =============================================================================
step "8/9 — Write mbx_ancombc2_info.txt"
# =============================================================================
QIIME_VER="$(qiime --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo unknown)"
PY_BIN="$(command -v python || true)"
PY_VER="$($PY_BIN -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>/dev/null || echo unknown)"
ANCOMBC_VER="$(_R -e 'cat(as.character(packageVersion("ANCOMBC")))' 2>/dev/null | tail -1)"

VAR_LINES=""
for col in "${CATEGORICAL_COLS[@]}"; do
  VAR_LINES="${VAR_LINES}VARIABLE=${col}
"
done
LEVEL_LINES=""
i=0
while [[ $i -lt $N_LEVELS ]]; do
  LEVEL_LINES="${LEVEL_LINES}LEVEL=${LEVEL_NAMES[$i]}
"
  i=$(( i + 1 ))
done

cat > "$INFO_FILE" << INFO
# ============================================================================
# mbx_ancombc2_info.txt
# Generated by mbx_ancombc2_run.sh
# Date : $NOW_PRETTY
# ============================================================================
#
# Do NOT edit the key=value lines below — they are parsed programmatically.

# ── Inputs ───────────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
BETA_INFO=$BETA_INFO
ALPHA_INFO=$ALPHA_INFO
PRE_DIVERSITY_INFO=$PRE_DIV_INFO
METADATA_TXT=$METADATA_TXT
FEATURE_TABLE_FILTERED_QZA=$FT_FILTERED_QZA
TAXONOMY_QZA=$TAXONOMY_QZA

# ── Output locations ─────────────────────────────────────────────────────────
ANCOMBC2_OUT_DIR=$ANCO_DIR
WORK_DIR=$WORK_DIR
COLLAPSE_DIR=$COLLAPSE_DIR
EXPORT_DIR=$EXPORT_DIR
RAW_OBJECTS_DIR=$RAW_DIR
LOGS_DIR=$LOG_DIR
GLOBAL_SUMMARY_XLSX=$GLOBAL_SUM_XLSX

# ── Settings ─────────────────────────────────────────────────────────────────
P_ADJUST_METHOD=$P_ADJUST
ALPHA=$ALPHA
PREVALENCE_CUTOFF=$PRV_CUT
LIBRARY_CUTOFF=$LIB_CUT
PERMUTATIONS=$B_PERMS

# ── Levels analysed ──────────────────────────────────────────────────────────
N_LEVELS=$N_LEVELS
${LEVEL_LINES}

# ── Categorical variables analysed ───────────────────────────────────────────
N_VARIABLES=${#CATEGORICAL_COLS[@]}
${VAR_LINES}

# ── Output counts ────────────────────────────────────────────────────────────
N_PRIMARY_XLSX=$N_PRIMARY
N_PAIRWISE_XLSX=$N_PAIRWISE
N_GLOBAL_XLSX=$N_GLOBAL
N_PNG_PLOTS=$N_PNG
N_RDS_OBJECTS=$N_RDS

# ── Provenance ───────────────────────────────────────────────────────────────
GENERATED_AT=$NOW_PRETTY
SCRIPT_NAME=mbx_ancombc2_run.sh
PLATFORM=$PLATFORM_LABEL
QIIME2_VERSION=$QIIME_VER
PYTHON_VERSION=$PY_VER
R_VERSION=$R_VERSION
ANCOMBC_VERSION=$ANCOMBC_VER
RSCRIPT_USED=$RSCRIPT_CMD
RSCRIPT_USED_CONDA_R=$RSCRIPT_USES_CONDA_R
N_JOBS=$N_JOBS
SKIPPED_DEPS_CHECK=$SKIP_DEPS_CHECK
SKIPPED_PLOTS=$SKIP_PLOTS
INVOCATION_USER=${USER:-unknown}
INVOCATION_CWD=$(pwd)
INVOCATION_ARGV=$0 $*
INFO

ok "Info file → $INFO_FILE"
sep

# =============================================================================
step "9/9 — Done"
# =============================================================================
echo ""
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║  STEP 14 COMPLETE  —  DIFFERENTIAL ABUNDANCE (ANCOMBC2)     ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "  Per-level results : $ANCO_DIR/ANCOMBC2_<level>/"
echo "  Working files     : $WORK_DIR/"
echo "  Levels            : ${LEVEL_NAMES[*]}"
echo "  Variables         : ${CATEGORICAL_COLS[*]}"
echo ""
echo "  Each ANCOMBC2_<level>/ancombc2_<var>/ folder contains:"
echo "    • ancombc2_primary_<var>.xlsx          (model coefficients: lfc, q, diff per coef)"
echo "    • ancombc2_pairwise_<var>.xlsx         (all pairwise contrasts, mdFDR-adjusted)"
echo "    • ancombc2_global_<var>.xlsx           (omnibus test, reference-free)"
echo "    • ancombc2_structural_zeros_<var>.xlsx (taxa absent from a whole group)"
echo "    • ancombc2_summary_<var>.xlsx          (high-level counts, settings, runtime)"
if ! $SKIP_PLOTS; then
echo "    • volcano_pairwise_<var>.png           (faceted across all pairs)"
echo "    • heatmap_significant_<var>.png        (DA taxa × samples)"
fi
echo ""
echo "  Global summary  : $GLOBAL_SUM_XLSX"
echo "  Open it in Excel:"
echo "    open '$GLOBAL_SUM_XLSX'"
echo ""
echo "  Raw RDS objects (for advanced custom analysis):"
echo "    $RAW_DIR/<level>_<var>.rds"
echo ""
echo "  Provenance (machine-readable) → $INFO_FILE"
echo ""
