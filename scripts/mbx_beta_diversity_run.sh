#!/usr/bin/env bash
# =============================================================================
#  mbx_beta_diversity_run.sh   —   Step 13 of the mbX Pro pipeline
# =============================================================================
#  PURPOSE
#  -------
#  Compute beta diversity (4 distance metrics), run PERMANOVA + PERMDISP +
#  pairwise PERMANOVA + Adonis (Type III + univariate) per categorical
#  metadata column, and produce publication-ready PCoAs, distance heatmaps,
#  and UPGMA dendrograms.
#
#  All inputs are read from the previous step's machine-readable info file:
#      <mbX_pro_outputs_*>/12_alpha_diversity_results/mbx_alpha_diversity_info.txt
#  The user provides only the mbX_pro_outputs directory path.
#
#  WHAT THIS SCRIPT PRODUCES
#  -------------------------
#  <mbX_pro_outputs_*>/
#  └── 13_beta_diversity_results/
#      ├── results_by_categorical_variables/
#      │   └── <Variable>/
#      │       ├── PCoA_<Metric>_<Variable>.png                 (4 metrics)
#      │       ├── PCoA_panel_4metrics_<Variable>.png            bonus
#      │       ├── Boxplot_DistanceToCentroid_<Metric>_<Variable>.png
#      │       ├── PERMANOVA_results_<Variable>.xlsx            1 sheet, 4 rows
#      │       ├── Pairwise_PERMANOVA_<Variable>.xlsx           only if k>2
#      │       ├── PERMDISP_results_<Variable>.xlsx
#      │       └── qzv/
#      │           ├── permanova_<Metric>_<Variable>.qzv         (4)
#      │           └── permdisp_<Metric>_<Variable>.qzv          (4)
#      ├── all_samples_beta_diversity/
#      │   ├── Adonis_multivariable_PERMANOVA.xlsx              4 sheets
#      │   ├── Distance_heatmap_<Metric>.png                    (4)
#      │   ├── UPGMA_dendrogram_<Metric>.png                    (4)
#      │   └── qzv/
#      │       └── emperor_<Metric>.qzv                         (4)
#      ├── working_dir_beta_diversity/
#      │   ├── distance_matrices/<metric>_distance_matrix.qza   (4)
#      │   ├── pcoa/<metric>_pcoa.qza                            (4)
#      │   └── exported/{distance_matrices,pcoa}/<metric>/...
#      └── mbx_beta_diversity_info.txt
#
#  GATING
#  ------
#  Refuses to run if step 11 reported OVERALL_STATUS=FAIL or
#  READY_FOR_DIVERSITY=no, unless --force is passed.
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
  grep "^${key}=" "$file" 2>/dev/null | cut -d= -f2- | head -1
}

# Sanitize a string for use as a directory / filename.
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

mbx_beta_diversity_run.sh — Beta diversity:  4 metrics × N categorical vars.
                            PERMANOVA + PERMDISP + pairwise + Adonis + PCoA
                            + distance heatmap + UPGMA dendrogram.

USAGE:
  mbx_beta_diversity_run.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads <mbX_pro_outputs_dir>/12_alpha_diversity_results/mbx_alpha_diversity_info.txt
  to discover:
      • Rarefied feature table (from step 12)
      • Rooted phylogenetic tree
      • Metadata file path
  Then computes:
      1. 4 distance matrices: Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac
      2. PCoA + emperor (interactive) per metric
      3. PERMANOVA + PERMDISP QZVs per (variable × metric)
      4. R-tabular: PERMANOVA, pairwise PERMANOVA (BH-adj), PERMDISP, distance-to-centroid
      5. PCoA PNG with 95% ellipses per (variable × metric)
      6. All-samples: multivariable Adonis (Type III + univariate),
                      Bray-Curtis distance heatmap (annotated by all categorical vars),
                      UPGMA dendrogram (one per metric).

OPTIONS:
  --metrics LIST         Comma-separated metric subset.  Default: all 4.
                         Choices: braycurtis,jaccard,weighted_unifrac,unweighted_unifrac
  --permutations N       PERMANOVA / PERMDISP permutations (default: 999)
  --skip-qzvs            Skip QIIME2 QZV generation (only R-tabular + PNGs)
  --skip-stats           Skip per-variable R analysis (only QZVs + heatmap/dendrogram)
  --skip-heatmap         Skip the all-samples heatmap + dendrogram
  --force                Run even if step 11 OVERALL_STATUS != PASS / PASS_WITH_WARNINGS
  --force-rerun          Recompute everything (ignore existing outputs)
  --dry-run              Print commands; do not execute
  -h, --help             Show this help

EXAMPLES:
  mbx_beta_diversity_run.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_beta_diversity_run.sh /path/to/mbX_pro_outputs_20260417_121431 --permutations 9999
  mbx_beta_diversity_run.sh /path/to/mbX_pro_outputs_20260417_121431 --metrics braycurtis,weighted_unifrac

PROVENANCE:
  All paths are read from 12_alpha_diversity_results/mbx_alpha_diversity_info.txt.
  Re-runs are safe — completed steps are skipped automatically.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
USER_METRICS=""
PERMS=999
SKIP_QZVS=false
SKIP_STATS=false
SKIP_HEATMAP=false
FORCE_GATE=false
FORCE_RERUN=false
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)        usage ;;
    --metrics)        USER_METRICS="${2:-}"; [[ -z "$USER_METRICS" ]] && err "--metrics requires a value"; shift 2 ;;
    --permutations)   PERMS="${2:-}"; case "$PERMS" in ''|*[!0-9]*) err "--permutations must be a positive integer";; esac; shift 2 ;;
    --skip-qzvs)      SKIP_QZVS=true; shift ;;
    --skip-stats)     SKIP_STATS=true; shift ;;
    --skip-heatmap)   SKIP_HEATMAP=true; shift ;;
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

case "$(basename "$MBX_OUT_DIR")" in
  mbX_pro_outputs_*|mbx_pro_outputs_*) : ;;
  *) warn "Directory name does not match 'mbX_pro_outputs_*': $(basename "$MBX_OUT_DIR")
        Continuing anyway." ;;
esac

# ── CPU detection ─────────────────────────────────────────────────────────────
if   command -v nproc   &>/dev/null; then N_JOBS="$(nproc)"
elif command -v sysctl  &>/dev/null; then N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else                                       N_JOBS=1
fi

# ── R / conda-env interaction (CRITICAL — read carefully) ────────────────────
#
# This script needs two DIFFERENT R installations to work simultaneously:
#
#   1. System R (e.g. /usr/local/bin/Rscript)
#      Used directly by us for tabular stats + plotting.  Has packages
#      installed by mbx_ezclean_all_levels.sh  (vegan, openxlsx, etc.).
#      Its compiled extensions use .so (Apple Silicon native convention).
#
#   2. Conda env's R ($CONDA_PREFIX/lib/R)
#      Used INTERNALLY by qiime — q2-composition imports rpy2, which boots
#      an embedded R.  Conda's R has compiled extensions named .dylib.
#
# When you `conda activate qiime2-amplicon-2025.4`, conda exports:
#     R_LIBS_USER=$CONDA_PREFIX/lib/R/library/
# This poisons the system R: it tries to load compiled extensions
# (methods.so, Rcpp.so, ...) from a path that contains .dylib variants.
# Result:  "shared object 'methods.dylib' not found"  or  "Rcpp.so not found".
#
# We solve this with TWO wrappers used everywhere in the script:
#   • _R     : strips every R_* env var so SYSTEM R uses its native libPaths
#   • _QIIME : sets R_HOME to CONDA's R so qiime's embedded R uses .dylib paths
#
# We also clear R_LIBS_USER from the parent shell's env so child processes
# (e.g. nested python -> rpy2) start with a clean slate.
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

# Detect the conda env's R home (used by _QIIME wrapper below).
_CONDA_R_HOME=""
if [[ -n "${CONDA_PREFIX:-}" && -d "${CONDA_PREFIX}/lib/R" ]]; then
  _CONDA_R_HOME="${CONDA_PREFIX}/lib/R"
fi

# ── Locate Rscript with working packages, strip conda R env vars on every call ─
# (Same fix as alpha-diversity:  conda activates set R_LIBS_USER which poisons
#  the system R's library path search → Rcpp.so / openxlsx etc. fail to load.)
_R() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$RSCRIPT_CMD" "$@"
}

# _QIIME wraps `qiime` so its rpy2-embedded R points at conda's R installation
# (which has matching .dylib compiled extensions).  Without this, q2-composition
# fails to boot rpy2 -> "shared object 'methods.dylib' not found".
_QIIME() {
  if [[ -n "$_CONDA_R_HOME" ]]; then
    R_HOME="$_CONDA_R_HOME" \
    R_LIBS_USER="$_CONDA_R_HOME/library/" \
      qiime "$@"
  else
    qiime "$@"
  fi
}
_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}
RSCRIPT_CMD=""
RPROBE='if (!requireNamespace("vegan", quietly=TRUE) ||
            !requireNamespace("openxlsx", quietly=TRUE)) quit(status=1); cat("ok")'
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
  → Then run mbx_ezclean_all_levels.sh once to install required R packages."
R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"

# ── PID-based temp file IDs (avoids mktemp suffix issues on macOS bash 3.2) ───
_TMPID="${$}_$(date +%s)"
trap 'rm -f /tmp/mbx_beta_${_TMPID}*.R /tmp/mbx_beta_${_TMPID}*.txt' EXIT

# =============================================================================
step "1/9 — Read mbx_alpha_diversity_info.txt"
# =============================================================================
ALPHA_INFO="${MBX_OUT_DIR}/12_alpha_diversity_results/mbx_alpha_diversity_info.txt"
PRE_DIV_INFO="${MBX_OUT_DIR}/11_pre_diversity/mbx_pre_diversity_info.txt"

[[ -f "$ALPHA_INFO" ]] || err "Alpha-diversity info file not found:
    $ALPHA_INFO
  → Run mbx_alpha_diversity_run.sh first (step 12)."
[[ -f "$PRE_DIV_INFO" ]] || err "Pre-diversity info file not found:
    $PRE_DIV_INFO
  → Run mbx_pre_diversity_parameters.sh first (step 11)."

OVERALL_STATUS="$(_read_key OVERALL_STATUS         "$PRE_DIV_INFO")"
READY="$(        _read_key READY_FOR_DIVERSITY    "$PRE_DIV_INFO")"
METADATA_TXT="$( _read_key METADATA_TXT           "$ALPHA_INFO")"
RAREFIED_QZA="$( _read_key RAREFIED_TABLE_QZA     "$ALPHA_INFO")"
TREE_QZA="$(     _read_key ROOTED_TREE_QZA        "$ALPHA_INFO")"
REC_DEPTH="$(    _read_key RECOMMENDED_DEPTH      "$ALPHA_INFO")"

[[ -n "$METADATA_TXT" && -f "$METADATA_TXT" ]] || err "METADATA_TXT not usable: '$METADATA_TXT'"
[[ -n "$RAREFIED_QZA" && -f "$RAREFIED_QZA" ]] || err "RAREFIED_TABLE_QZA not usable: '$RAREFIED_QZA'
  → Re-run mbx_alpha_diversity_run.sh."
[[ -n "$TREE_QZA"     && -f "$TREE_QZA"     ]] || err "ROOTED_TREE_QZA not usable: '$TREE_QZA'"

# ── Gating ────────────────────────────────────────────────────────────────────
case "$OVERALL_STATUS" in
  PASS|PASS_WITH_WARNINGS) : ;;
  *)
    if $FORCE_GATE; then
      warn "Step 11 OVERALL_STATUS = '$OVERALL_STATUS' — proceeding because --force was given."
    else
      err "Step 11 OVERALL_STATUS = '$OVERALL_STATUS'.
  → Diversity analysis on a non-PASS state may produce misleading results.
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
info "QIIME2 env           : ${CONDA_DEFAULT_ENV:-unknown}"
info "Rscript              : $RSCRIPT_CMD  (R $R_VERSION)"
info "Output root          : $MBX_OUT_DIR"
info "Metadata             : $METADATA_TXT"
info "Rarefied table       : $RAREFIED_QZA"
info "Rooted tree          : $TREE_QZA"
info "Sampling depth (q11) : $REC_DEPTH"
info "CPU cores            : $N_JOBS"
info "Permutations         : $PERMS"
info "Step-11 status       : $OVERALL_STATUS  /  ready=$READY"
sep

# ── Distance metric tables ───────────────────────────────────────────────────
# parallel arrays: machine name (QIIME2/dir), pretty label, safe filename token
ALL_QIIME=(   "braycurtis" "jaccard" "weighted_unifrac" "unweighted_unifrac" )
ALL_LABELS=(  "Bray-Curtis" "Jaccard" "Weighted UniFrac" "Unweighted UniFrac" )
ALL_SAFE=(    "BrayCurtis"  "Jaccard" "WeightedUniFrac"  "UnweightedUniFrac"  )
ALL_DIRTOK=(  "bray_curtis" "jaccard" "weighted_unifrac" "unweighted_unifrac" )

# Filter to user-selected metrics if --metrics was given
QIIME_METRICS=()
LABELS=()
SAFE_NAMES=()
DIR_TOKENS=()
if [[ -z "$USER_METRICS" ]]; then
  QIIME_METRICS=( "${ALL_QIIME[@]}" )
  LABELS=(        "${ALL_LABELS[@]}" )
  SAFE_NAMES=(    "${ALL_SAFE[@]}" )
  DIR_TOKENS=(    "${ALL_DIRTOK[@]}" )
else
  IFS=',' read -r -a _user <<< "$USER_METRICS"
  for _u in "${_user[@]}"; do
    _u_trimmed="$(echo "$_u" | tr -d ' ' | tr '[:upper:]' '[:lower:]')"
    _found=false
    for i in "${!ALL_QIIME[@]}"; do
      if [[ "${ALL_QIIME[$i]}" == "$_u_trimmed" ]]; then
        QIIME_METRICS+=( "${ALL_QIIME[$i]}" )
        LABELS+=(        "${ALL_LABELS[$i]}" )
        SAFE_NAMES+=(    "${ALL_SAFE[$i]}" )
        DIR_TOKENS+=(    "${ALL_DIRTOK[$i]}" )
        _found=true
        break
      fi
    done
    $_found || err "Unknown metric in --metrics: '$_u'
  Choices: braycurtis, jaccard, weighted_unifrac, unweighted_unifrac"
  done
  [[ ${#QIIME_METRICS[@]} -eq 0 ]] && err "No valid metrics in --metrics='$USER_METRICS'"
fi
N_METRICS=${#QIIME_METRICS[@]}
info "Distance metrics : ${QIIME_METRICS[*]}"
sep

# ── Set up output directories ─────────────────────────────────────────────────
BETA_DIR="${MBX_OUT_DIR}/13_beta_diversity_results"
PER_VAR_DIR="${BETA_DIR}/results_by_categorical_variables"
ALL_SAMP_DIR="${BETA_DIR}/all_samples_beta_diversity"
ALL_SAMP_QZV="${ALL_SAMP_DIR}/qzv"
WORK_DIR="${BETA_DIR}/working_dir_beta_diversity"
DM_DIR="${WORK_DIR}/distance_matrices"
PCOA_DIR="${WORK_DIR}/pcoa"
EXPORT_DIR="${WORK_DIR}/exported"
EXPORT_DM="${EXPORT_DIR}/distance_matrices"
EXPORT_PCOA="${EXPORT_DIR}/pcoa"
INFO_FILE="${BETA_DIR}/mbx_beta_diversity_info.txt"

mkdir -p "$BETA_DIR" "$PER_VAR_DIR" "$ALL_SAMP_DIR" "$ALL_SAMP_QZV" \
         "$WORK_DIR" "$DM_DIR" "$PCOA_DIR" "$EXPORT_DM" "$EXPORT_PCOA"

# =============================================================================
step "2/9 — Compute ${N_METRICS} distance matrices (QIIME2)"
# =============================================================================
i=0
while [[ $i -lt $N_METRICS ]]; do
  METRIC="${QIIME_METRICS[$i]}"
  LABEL="${LABELS[$i]}"
  TOKEN="${DIR_TOKENS[$i]}"
  OUT_QZA="${DM_DIR}/${TOKEN}_distance_matrix.qza"

  if [[ -f "$OUT_QZA" && "$FORCE_RERUN" == false ]]; then
    skipped "[$((i+1))/$N_METRICS] $METRIC — distance matrix already exists"
    i=$(( i + 1 )); continue
  fi

  info "[$((i+1))/$N_METRICS] Computing $METRIC ($LABEL)..."

  if [[ "$METRIC" == "weighted_unifrac" || "$METRIC" == "unweighted_unifrac" ]]; then
    cmd_show "qiime diversity beta-phylogenetic" \
      "--i-table $RAREFIED_QZA" \
      "--i-phylogeny $TREE_QZA" \
      "--p-metric $METRIC" \
      "--p-threads $N_JOBS" \
      "--o-distance-matrix $OUT_QZA"

    if ! $DRY_RUN; then
      timer_start
      _QIIME diversity beta-phylogenetic \
        --i-table          "$RAREFIED_QZA" \
        --i-phylogeny      "$TREE_QZA" \
        --p-metric         "$METRIC" \
        --p-threads        "$N_JOBS" \
        --o-distance-matrix "$OUT_QZA" \
        || err "qiime diversity beta-phylogenetic ($METRIC) failed.
  → Make sure the rooted tree contains every feature in the rarefied table.
  → If feature IDs have drifted between steps, re-run step 11 to regenerate the tree."
      timer_end
    fi
  else
    cmd_show "qiime diversity beta" \
      "--i-table $RAREFIED_QZA" \
      "--p-metric $METRIC" \
      "--p-n-jobs $N_JOBS" \
      "--o-distance-matrix $OUT_QZA"

    if ! $DRY_RUN; then
      timer_start
      _QIIME diversity beta \
        --i-table          "$RAREFIED_QZA" \
        --p-metric         "$METRIC" \
        --p-n-jobs         "$N_JOBS" \
        --o-distance-matrix "$OUT_QZA" \
        || err "qiime diversity beta ($METRIC) failed."
      timer_end
    fi
  fi
  ok "  → $(basename "$OUT_QZA")"
  i=$(( i + 1 ))
done
sep

# =============================================================================
step "3/9 — PCoA (QIIME2) + Emperor visualizations"
# =============================================================================
i=0
while [[ $i -lt $N_METRICS ]]; do
  METRIC="${QIIME_METRICS[$i]}"
  LABEL="${LABELS[$i]}"
  TOKEN="${DIR_TOKENS[$i]}"
  SAFE="${SAFE_NAMES[$i]}"
  DM_QZA="${DM_DIR}/${TOKEN}_distance_matrix.qza"
  PCOA_QZA="${PCOA_DIR}/${TOKEN}_pcoa.qza"
  EMP_QZV="${ALL_SAMP_QZV}/emperor_${SAFE}.qzv"

  # PCoA
  if [[ -f "$PCOA_QZA" && "$FORCE_RERUN" == false ]]; then
    skipped "[$((i+1))/$N_METRICS] PCoA exists for $METRIC"
  else
    info "[$((i+1))/$N_METRICS] PCoA for $METRIC..."
    cmd_show "qiime diversity pcoa" \
      "--i-distance-matrix $DM_QZA" \
      "--o-pcoa $PCOA_QZA"
    if ! $DRY_RUN; then
      _QIIME diversity pcoa \
        --i-distance-matrix "$DM_QZA" \
        --o-pcoa            "$PCOA_QZA" \
        || err "qiime diversity pcoa ($METRIC) failed."
    fi
    ok "  PCoA → $(basename "$PCOA_QZA")"
  fi

  # Emperor
  if $SKIP_QZVS; then
    skipped "[--skip-qzvs] emperor for $METRIC"
  elif [[ -f "$EMP_QZV" && "$FORCE_RERUN" == false ]]; then
    skipped "[$((i+1))/$N_METRICS] Emperor exists for $METRIC"
  else
    info "[$((i+1))/$N_METRICS] Emperor plot for $METRIC..."
    cmd_show "qiime emperor plot" \
      "--i-pcoa $PCOA_QZA" \
      "--m-metadata-file $METADATA_TXT" \
      "--o-visualization $EMP_QZV"
    if ! $DRY_RUN; then
      _QIIME emperor plot \
        --i-pcoa            "$PCOA_QZA" \
        --m-metadata-file   "$METADATA_TXT" \
        --o-visualization   "$EMP_QZV" \
        || err "qiime emperor plot ($METRIC) failed."
    fi
    ok "  Emperor → $(basename "$EMP_QZV")"
  fi
  i=$(( i + 1 ))
done
sep

# =============================================================================
step "4/9 — Detect categorical metadata columns"
# =============================================================================
CATS_FILE="/tmp/mbx_beta_${_TMPID}_cats.txt"
R_DETECT="/tmp/mbx_beta_${_TMPID}_detect.R"

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
  warn "[DRY-RUN] Would detect categorical columns from $METADATA_TXT"
  CATEGORICAL_COLS=( "Treatment" "SampleType" )
else
  _R --vanilla "$R_DETECT" "$METADATA_TXT" > "$CATS_FILE" 2>&1 \
    || err "Categorical column detection failed.
  → Verify metadata file: $METADATA_TXT"
  grep "^\[" "$CATS_FILE" || true
  CATEGORICAL_COLS=()
  while IFS= read -r _col; do
    [[ -n "$_col" ]] && CATEGORICAL_COLS+=("$_col")
  done < <(grep -v "^\[" "$CATS_FILE" | grep -v '^$' || true)
  [[ ${#CATEGORICAL_COLS[@]} -eq 0 ]] && err "No categorical columns detected.
  → Verify your metadata has at least one grouping column."
fi
rm -f "$R_DETECT"

ok "Categorical columns: ${CATEGORICAL_COLS[*]}"

SANITIZED_COLS=()
for col in "${CATEGORICAL_COLS[@]}"; do
  SANITIZED_COLS+=( "$(_sanitize "$col")" )
  mkdir -p "${PER_VAR_DIR}/$(_sanitize "$col")/qzv"
done
sep

# =============================================================================
step "5/9 — Per-variable QIIME2 QZVs (PERMANOVA + PERMDISP)"
# =============================================================================
if $SKIP_QZVS; then
  warn "Skipping per-variable QIIME2 QZV generation (--skip-qzvs)."
else
  N_VARS=${#CATEGORICAL_COLS[@]}
  TOTAL_QZV=$(( N_VARS * N_METRICS * 2 ))
  RUN_NUMBER=0

  for v_idx in "${!CATEGORICAL_COLS[@]}"; do
    COL="${CATEGORICAL_COLS[$v_idx]}"
    SAFE_VAR="${SANITIZED_COLS[$v_idx]}"
    QZV_DIR="${PER_VAR_DIR}/${SAFE_VAR}/qzv"

    echo ""
    echo "  ── Variable: $COL  →  $QZV_DIR/"

    i=0
    while [[ $i -lt $N_METRICS ]]; do
      METRIC="${QIIME_METRICS[$i]}"
      SAFE_M="${SAFE_NAMES[$i]}"
      TOKEN="${DIR_TOKENS[$i]}"
      DM_QZA="${DM_DIR}/${TOKEN}_distance_matrix.qza"
      PERMA_QZV="${QZV_DIR}/permanova_${SAFE_M}_${SAFE_VAR}.qzv"
      PERMD_QZV="${QZV_DIR}/permdisp_${SAFE_M}_${SAFE_VAR}.qzv"

      # PERMANOVA
      RUN_NUMBER=$(( RUN_NUMBER + 1 ))
      if [[ -f "$PERMA_QZV" && "$FORCE_RERUN" == false ]]; then
        skipped "[$RUN_NUMBER/$TOTAL_QZV] permanova $METRIC × $COL — already exists"
      else
        info "[$RUN_NUMBER/$TOTAL_QZV] PERMANOVA  $METRIC × $COL ..."
        if ! $DRY_RUN; then
          _QIIME diversity beta-group-significance \
            --i-distance-matrix  "$DM_QZA" \
            --m-metadata-file    "$METADATA_TXT" \
            --m-metadata-column  "$COL" \
            --p-method           permanova \
            --p-pairwise         \
            --p-permutations     "$PERMS" \
            --o-visualization    "$PERMA_QZV" \
            >/dev/null 2>&1 \
            || warn "  ✘  permanova $METRIC × $COL failed (probably <2 valid groups). Skipping."
        fi
        [[ -f "$PERMA_QZV" ]] && ok "    permanova → $(basename "$PERMA_QZV")"
      fi

      # PERMDISP
      RUN_NUMBER=$(( RUN_NUMBER + 1 ))
      if [[ -f "$PERMD_QZV" && "$FORCE_RERUN" == false ]]; then
        skipped "[$RUN_NUMBER/$TOTAL_QZV] permdisp $METRIC × $COL — already exists"
      else
        info "[$RUN_NUMBER/$TOTAL_QZV] PERMDISP   $METRIC × $COL ..."
        if ! $DRY_RUN; then
          _QIIME diversity beta-group-significance \
            --i-distance-matrix  "$DM_QZA" \
            --m-metadata-file    "$METADATA_TXT" \
            --m-metadata-column  "$COL" \
            --p-method           permdisp \
            --p-permutations     "$PERMS" \
            --o-visualization    "$PERMD_QZV" \
            >/dev/null 2>&1 \
            || warn "  ✘  permdisp $METRIC × $COL failed. Skipping."
        fi
        [[ -f "$PERMD_QZV" ]] && ok "    permdisp  → $(basename "$PERMD_QZV")"
      fi
      i=$(( i + 1 ))
    done
  done
fi
sep

# =============================================================================
step "6/9 — Export distance matrices + PCoA artifacts to TSV"
# =============================================================================
i=0
while [[ $i -lt $N_METRICS ]]; do
  METRIC="${QIIME_METRICS[$i]}"
  TOKEN="${DIR_TOKENS[$i]}"

  DM_QZA="${DM_DIR}/${TOKEN}_distance_matrix.qza"
  DM_OUT="${EXPORT_DM}/${TOKEN}"
  DM_TSV="${DM_OUT}/distance-matrix.tsv"

  if [[ -f "$DM_TSV" && "$FORCE_RERUN" == false ]]; then
    skipped "DM TSV exists: $METRIC"
  else
    rm -rf "$DM_OUT"; mkdir -p "$DM_OUT"
    if ! $DRY_RUN; then
      _QIIME tools export --input-path "$DM_QZA" --output-path "$DM_OUT" >/dev/null \
        || err "Failed to export $DM_QZA"
      [[ -f "$DM_TSV" ]] || err "Expected DM TSV not found: $DM_TSV"
    fi
  fi

  PCOA_QZA="${PCOA_DIR}/${TOKEN}_pcoa.qza"
  PCOA_OUT="${EXPORT_PCOA}/${TOKEN}"
  PCOA_TXT="${PCOA_OUT}/ordination.txt"
  if [[ -f "$PCOA_TXT" && "$FORCE_RERUN" == false ]]; then
    skipped "PCoA ordination exists: $METRIC"
  else
    rm -rf "$PCOA_OUT"; mkdir -p "$PCOA_OUT"
    if ! $DRY_RUN; then
      _QIIME tools export --input-path "$PCOA_QZA" --output-path "$PCOA_OUT" >/dev/null \
        || err "Failed to export $PCOA_QZA"
    fi
  fi
  i=$(( i + 1 ))
done
ok "Distance matrices + PCoA exports under: $EXPORT_DIR"
sep

# =============================================================================
step "7/9 — R analysis (per-variable + all-samples)"
# =============================================================================
if $SKIP_STATS && $SKIP_HEATMAP; then
  warn "Skipping all R analysis (--skip-stats and --skip-heatmap)."
else
  R_BETA="/tmp/mbx_beta_${_TMPID}_analysis.R"

  # ── Build R-side variable lists ─────────────────────────────────────────────
  VARS_R=""
  for col in "${CATEGORICAL_COLS[@]}"; do
    _esc="${col//\"/\\\"}"
    if [[ -z "$VARS_R" ]]; then VARS_R="\"$_esc\""
    else                        VARS_R="$VARS_R, \"$_esc\""; fi
  done
  METRICS_R=""
  LABELS_R=""
  SAFE_R=""
  TOKENS_R=""
  i=0
  while [[ $i -lt $N_METRICS ]]; do
    if [[ -z "$METRICS_R" ]]; then
      METRICS_R="\"${QIIME_METRICS[$i]}\""
      LABELS_R="\"${LABELS[$i]}\""
      SAFE_R="\"${SAFE_NAMES[$i]}\""
      TOKENS_R="\"${DIR_TOKENS[$i]}\""
    else
      METRICS_R="$METRICS_R, \"${QIIME_METRICS[$i]}\""
      LABELS_R="$LABELS_R, \"${LABELS[$i]}\""
      SAFE_R="$SAFE_R, \"${SAFE_NAMES[$i]}\""
      TOKENS_R="$TOKENS_R, \"${DIR_TOKENS[$i]}\""
    fi
    i=$(( i + 1 ))
  done

  if $SKIP_STATS;   then SKIP_STATS_R="TRUE";   else SKIP_STATS_R="FALSE";   fi
  if $SKIP_HEATMAP; then SKIP_HEATMAP_R="TRUE"; else SKIP_HEATMAP_R="FALSE"; fi

  cat > "$R_BETA" << RBETA
# =============================================================================
# Beta-diversity analysis  (per-variable + all-samples)
# =============================================================================
suppressPackageStartupMessages({
  required_pkgs <- c("vegan","ape","ggplot2","openxlsx","dplyr","pheatmap",
                     "RColorBrewer","patchwork")
  to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly=TRUE)]
  if (length(to_install) > 0) {
    cat(sprintf("[INFO]  Installing R package(s): %s\n", paste(to_install, collapse=", ")))
    install.packages(to_install, repos="https://cloud.r-project.org", quiet=TRUE)
  }
  library(vegan); library(ape); library(ggplot2); library(openxlsx)
  library(dplyr); library(pheatmap); library(RColorBrewer); library(patchwork)
})

# ── Inputs ──────────────────────────────────────────────────────────────────
METADATA       <- "${METADATA_TXT}"
EXPORT_DIR     <- "${EXPORT_DIR}"
PER_VAR_DIR    <- "${PER_VAR_DIR}"
ALL_SAMP_DIR   <- "${ALL_SAMP_DIR}"
PERMS          <- ${PERMS}
SKIP_STATS     <- ${SKIP_STATS_R}
SKIP_HEATMAP   <- ${SKIP_HEATMAP_R}
VARS           <- c(${VARS_R})
METRICS        <- c(${METRICS_R})
LABELS         <- c(${LABELS_R})
SAFE_M         <- c(${SAFE_R})
TOKENS         <- c(${TOKENS_R})
names(LABELS)  <- METRICS
names(SAFE_M)  <- METRICS
names(TOKENS)  <- METRICS

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
load_dm <- function(token) {
  tsv <- file.path(EXPORT_DIR, "distance_matrices", token, "distance-matrix.tsv")
  if (!file.exists(tsv)) stop(sprintf("Missing distance-matrix TSV: %s", tsv))
  m <- read.table(tsv, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  rownames(m) <- trimws(rownames(m))
  colnames(m) <- trimws(colnames(m))
  as.dist(as.matrix(m))
}

cat(sprintf("[INFO]  Reading metadata: %s\n", METADATA))
md <- read_metadata(METADATA)

# Restrict metadata to samples actually present in the (rarefied) distance matrix.
dm0  <- load_dm(TOKENS[METRICS[1]])
keep_ids <- attr(dm0, "Labels")
md_full  <- md
md       <- md[md[["sample-id"]] %in% keep_ids, , drop=FALSE]
rownames(md) <- md[["sample-id"]]
md <- md[match(keep_ids, md[["sample-id"]]), , drop=FALSE]
rownames(md) <- md[["sample-id"]]

dropped <- setdiff(md_full[["sample-id"]], keep_ids)
if (length(dropped) > 0)
  cat(sprintf("[INFO]  %d sample(s) in metadata not in distance matrix (rarefaction dropped): %s\n",
              length(dropped), paste(head(dropped, 5), collapse=", ")))
cat(sprintf("[INFO]  Beta-diversity samples: %d  |  variables: %d  |  metrics: %d\n",
            nrow(md), length(VARS), length(METRICS)))

# Ensure VARS exist in metadata
VARS <- VARS[VARS %in% names(md)]
if (length(VARS) == 0) stop("No requested variables found in metadata after sample filtering.")

# ── Pre-load every distance matrix once ─────────────────────────────────────
DMS <- list()
for (m in METRICS) {
  DMS[[m]] <- load_dm(TOKENS[m])
  cat(sprintf("[INFO]  Loaded %-22s (%d samples, %d pairwise distances)\n",
              LABELS[m], attr(DMS[[m]], "Size"), length(DMS[[m]])))
}

# Reasonable colour palette per group
qual_colours <- function(n) {
  if (n <= 8) RColorBrewer::brewer.pal(max(n,3), "Set2")[seq_len(n)]
  else if (n <= 12) RColorBrewer::brewer.pal(n, "Set3")
  else colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
}

# ───────────────────────────────────────────────────────────────────────────
# PER-VARIABLE ANALYSIS
# ───────────────────────────────────────────────────────────────────────────
permanova_rows  <- list()  # global, key = variable
pairwise_rows   <- list()
permdisp_rows   <- list()
fail_log        <- character(0)

if (!SKIP_STATS) {
  for (var in VARS) {
    safe_var <- sanitize(var)
    var_dir  <- file.path(PER_VAR_DIR, safe_var)
    dir.create(var_dir, recursive=TRUE, showWarnings=FALSE)
    cat(sprintf("\n  ── %s  →  %s/\n", var, var_dir))

    # Per-variable working frame
    df <- md[!is.na(md[[var]]) & nchar(trimws(as.character(md[[var]]))) > 0, , drop=FALSE]
    df[[var]] <- factor(trimws(as.character(df[[var]])))
    # drop singletons
    grp_n <- table(df[[var]])
    bad   <- names(grp_n)[grp_n < 2]
    if (length(bad) > 0) {
      df <- df[!df[[var]] %in% bad, , drop=FALSE]
      df[[var]] <- droplevels(df[[var]])
    }
    if (length(unique(df[[var]])) < 2) {
      cat(sprintf("[WARN]  %s: <2 valid groups after filtering — skipping.\n", var))
      next
    }
    keep_ids_var <- df[["sample-id"]]
    k     <- length(unique(df[[var]]))
    n_var <- nrow(df)
    cat(sprintf("[INFO]  %s: n=%d, k=%d groups, dropped=%s\n",
                var, n_var, k,
                if (length(bad)) paste(bad, collapse=",") else "none"))

    # Per-variable PERMANOVA / PERMDISP / pairwise across metrics
    perm_var      <- list()
    pairwise_var  <- list()
    permdisp_var  <- list()
    panel_plots   <- list()
    pal_var       <- setNames(qual_colours(k), levels(df[[var]]))

    for (metric in METRICS) {
      label  <- LABELS[metric];  safe_m <- SAFE_M[metric]
      dm_full <- DMS[[metric]]

      # Subset distance matrix to the samples retained for this variable
      idx     <- match(keep_ids_var, attr(dm_full, "Labels"))
      idx     <- idx[!is.na(idx)]
      m_full  <- as.matrix(dm_full)
      dm_sub  <- as.dist(m_full[idx, idx])
      df_sub  <- df[match(rownames(m_full)[idx], df[["sample-id"]]), , drop=FALSE]

      # ── PERMANOVA (omnibus) ──────────────────────────────────────────────
      adn <- tryCatch(
        adonis2(as.formula(sprintf("dm_sub ~ \`%s\`", var)),
                data = df_sub, permutations = PERMS, by = "terms"),
        error = function(e) { fail_log <<- c(fail_log,
          sprintf("PERMANOVA %s × %s: %s", var, metric, conditionMessage(e))); NULL })
      if (is.null(adn)) next

      F_stat  <- adn\$F[1]
      R2      <- adn\$R2[1]
      p_val   <- adn\$\`Pr(>F)\`[1]
      df_b    <- adn\$Df[1]
      df_w    <- adn\$Df[2]
      perm_var[[metric]] <- data.frame(
        Metric            = label,
        Variable          = var,
        n_samples         = n_var,
        n_groups          = k,
        Pseudo_F          = round(F_stat, 4),
        R_squared         = round(R2, 4),
        Variance_explained_pct = round(R2 * 100, 2),
        df_between        = df_b,
        df_within         = df_w,
        Permutations      = PERMS,
        p_value           = signif(p_val, 6),
        Significant_at_05 = !is.na(p_val) && p_val < 0.05,
        stringsAsFactors  = FALSE
      )
      cat(sprintf("[OK]    %-18s × %-14s  F=%.3f  R^2=%.4f (%.2f%%)  p=%s\n",
                  var, label, F_stat, R2, R2*100,
                  format.pval(p_val, digits=3, eps=1e-6)))

      # ── Pairwise PERMANOVA (only if k > 2) ───────────────────────────────
      if (k > 2) {
        groups   <- levels(df_sub[[var]])
        pairs    <- combn(groups, 2, simplify=FALSE)
        rows     <- list()
        for (pp in pairs) {
          mask <- df_sub[[var]] %in% pp
          d2   <- df_sub[mask, , drop=FALSE]
          d2[[var]] <- droplevels(factor(d2[[var]]))
          ids2 <- d2[["sample-id"]]
          ix   <- match(ids2, attr(dm_sub, "Labels"))
          ix   <- ix[!is.na(ix)]
          m2   <- as.matrix(dm_sub)
          dm2  <- as.dist(m2[ix, ix])
          adn2 <- tryCatch(
            adonis2(as.formula(sprintf("dm2 ~ \`%s\`", var)),
                    data=d2, permutations=PERMS, by="terms"),
            error = function(e) NULL)
          if (is.null(adn2)) next
          rows[[length(rows)+1]] <- data.frame(
            Metric          = label,
            Variable        = var,
            comparison      = paste(pp, collapse=" vs "),
            n_pair          = nrow(d2),
            Pseudo_F        = round(adn2\$F[1], 4),
            R_squared       = round(adn2\$R2[1], 4),
            Variance_explained_pct = round(adn2\$R2[1] * 100, 2),
            p_unadjusted    = signif(adn2\$\`Pr(>F)\`[1], 6),
            stringsAsFactors = FALSE
          )
        }
        if (length(rows) > 0) {
          pw <- do.call(rbind, rows)
          pw\$p_BH_adjusted    <- signif(p.adjust(pw\$p_unadjusted, method="BH"), 6)
          pw\$Significant_BH05 <- !is.na(pw\$p_BH_adjusted) & pw\$p_BH_adjusted < 0.05
          pairwise_var[[metric]] <- pw
        }
      }

      # ── PERMDISP ─────────────────────────────────────────────────────────
      bd <- tryCatch(betadisper(dm_sub, df_sub[[var]], type="centroid"),
                     error = function(e) { fail_log <<- c(fail_log,
                       sprintf("PERMDISP %s × %s: %s", var, metric,
                               conditionMessage(e))); NULL })
      if (!is.null(bd)) {
        pd <- tryCatch(permutest(bd, permutations=PERMS),
                       error = function(e) NULL)
        if (!is.null(pd)) {
          permdisp_var[[metric]] <- data.frame(
            Metric            = label,
            Variable          = var,
            n_samples         = n_var,
            n_groups          = k,
            F_stat            = round(pd\$tab[1, "F"], 4),
            df_between        = pd\$tab[1, "Df"],
            df_within         = pd\$tab[2, "Df"],
            Permutations      = PERMS,
            p_value           = signif(pd\$tab[1, "Pr(>F)"], 6),
            Significant_at_05 = !is.na(pd\$tab[1, "Pr(>F)"]) && pd\$tab[1, "Pr(>F)"] < 0.05,
            note              = "PERMDISP p>0.05 means dispersions are HOMOGENEOUS (good — PERMANOVA is valid)",
            stringsAsFactors  = FALSE)

          # Distance-to-centroid boxplot
          dcent <- data.frame(group = bd\$group, distance = bd\$distances)
          box_p <- ggplot(dcent, aes(x = group, y = distance, fill = group)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.85, colour = "#1f2d3d") +
            geom_jitter(width = 0.15, height = 0, size = 1.6,
                        alpha = 0.7, colour = "#1f2d3d") +
            scale_fill_manual(values = pal_var, guide = "none") +
            labs(title = sprintf("Distance-to-centroid (%s) by %s", label, var),
                 subtitle = sprintf("PERMDISP F = %.3f, p = %s, perm = %d",
                                    pd\$tab[1, "F"],
                                    format.pval(pd\$tab[1, "Pr(>F)"], digits=3, eps=1e-6),
                                    PERMS),
                 x = var, y = sprintf("Distance to %s centroid", label)) +
            theme_classic(base_size = 13) +
            theme(plot.title = element_text(face="bold"),
                  axis.text.x = element_text(angle = 25, hjust = 1))
          ggsave(file.path(var_dir,
                  sprintf("Boxplot_DistanceToCentroid_%s_%s.png", safe_m, safe_var)),
                 box_p, width = 7, height = 5, dpi = 300, bg = "white")
        }
      }

      # ── PCoA plot (R cmdscale; same eigenvalues QIIME2 produces) ─────────
      pco <- tryCatch(cmdscale(dm_sub, k = 3, eig = TRUE, add = TRUE),
                      error = function(e) NULL)
      if (!is.null(pco)) {
        eig_pos <- pco\$eig[pco\$eig > 0]
        pct12   <- round(100 * pco\$eig[1:2] / sum(eig_pos), 1)
        pts <- data.frame(
          SampleID = rownames(pco\$points),
          PC1 = pco\$points[, 1],
          PC2 = pco\$points[, 2],
          stringsAsFactors = FALSE)
        pts <- merge(pts, df_sub, by.x="SampleID", by.y="sample-id", all.x=TRUE)

        # Build R^2 / p annotation line
        ann <- sprintf("PERMANOVA: F = %.2f, R^2 = %.3f (%.1f%%), p = %s, perm = %d",
                       F_stat, R2, R2*100,
                       format.pval(p_val, digits=3, eps=1e-6), PERMS)

        p <- ggplot(pts, aes(x = PC1, y = PC2, colour = .data[[var]])) +
          geom_point(size = 3, alpha = 0.85) +
          stat_ellipse(aes(group = .data[[var]]), type = "norm", level = 0.95,
                       linewidth = 0.6) +
          scale_colour_manual(values = pal_var) +
          labs(title    = sprintf("%s PCoA -- coloured by %s", label, var),
               subtitle = ann,
               x = sprintf("PC1 (%.1f%%)", pct12[1]),
               y = sprintf("PC2 (%.1f%%)", pct12[2]),
               colour = var) +
          theme_classic(base_size = 13) +
          theme(plot.title    = element_text(face = "bold"),
                plot.subtitle = element_text(colour = "grey30", size = 10),
                legend.position = "right")

        ggsave(file.path(var_dir,
                sprintf("PCoA_%s_%s.png", safe_m, safe_var)),
               p, width = 7.5, height = 5.5, dpi = 300, bg = "white")
        panel_plots[[metric]] <- p
      }
    }   # end metric loop

    # ── Write per-variable xlsx (PERMANOVA, Pairwise, PERMDISP) ────────────
    if (length(perm_var) > 0) {
      perm_df <- do.call(rbind, perm_var)
      write.xlsx(perm_df,
        file.path(var_dir, sprintf("PERMANOVA_results_%s.xlsx", safe_var)),
        overwrite = TRUE,
        headerStyle = createStyle(textDecoration = "bold"))
      permanova_rows[[var]] <- perm_df
    }
    if (length(pairwise_var) > 0) {
      pw_df <- do.call(rbind, pairwise_var)
      write.xlsx(pw_df,
        file.path(var_dir, sprintf("Pairwise_PERMANOVA_%s.xlsx", safe_var)),
        overwrite = TRUE,
        headerStyle = createStyle(textDecoration = "bold"))
      pairwise_rows[[var]] <- pw_df
    }
    if (length(permdisp_var) > 0) {
      pd_df <- do.call(rbind, permdisp_var)
      write.xlsx(pd_df,
        file.path(var_dir, sprintf("PERMDISP_results_%s.xlsx", safe_var)),
        overwrite = TRUE,
        headerStyle = createStyle(textDecoration = "bold"))
      permdisp_rows[[var]] <- pd_df
    }

    # ── 4-metric PCoA panel (one PNG per variable) ─────────────────────────
    if (length(panel_plots) >= 2) {
      panel <- Reduce(\`+\`, panel_plots) +
        plot_layout(ncol = 2, guides = "collect") +
        plot_annotation(title = sprintf("PCoA panel by %s (%d metrics)",
                                        var, length(panel_plots)),
                        theme = theme(plot.title = element_text(face="bold", size=14)))
      ggsave(file.path(var_dir,
              sprintf("PCoA_panel_4metrics_%s.png", safe_var)),
             panel, width = 14, height = 11, dpi = 300, bg = "white")
    }
  }   # end variable loop

  # ── Global summary across all variables × metrics ────────────────────────
  if (length(permanova_rows) > 0) {
    glob_perm <- do.call(rbind, permanova_rows); rownames(glob_perm) <- NULL
    write.xlsx(glob_perm,
      file.path(PER_VAR_DIR, "Summary_PERMANOVA_all_variables_all_metrics.xlsx"),
      overwrite = TRUE,
      headerStyle = createStyle(textDecoration = "bold"))
    cat(sprintf("\n[OK]    Global PERMANOVA summary  → %s\n",
                file.path(PER_VAR_DIR, "Summary_PERMANOVA_all_variables_all_metrics.xlsx")))
  }
}

# ───────────────────────────────────────────────────────────────────────────
# ALL-SAMPLES ANALYSIS:  Adonis (Type III + univariate), heatmap, dendrogram
# ───────────────────────────────────────────────────────────────────────────
if (!SKIP_HEATMAP) {
  dir.create(ALL_SAMP_DIR, recursive=TRUE, showWarnings=FALSE)

  # Filter metadata to only samples present in the distance matrix and only
  # to columns we are actually testing (drops sample-id and any non-categorical).
  md_adonis <- md[, c("sample-id", VARS), drop=FALSE]
  for (v in VARS) md_adonis[[v]] <- factor(trimws(as.character(md_adonis[[v]])))
  rownames(md_adonis) <- md_adonis[["sample-id"]]

  # Drop samples with NA in any tested variable (adonis can't handle them)
  complete_idx <- complete.cases(md_adonis[, VARS, drop=FALSE])
  if (any(!complete_idx)) {
    cat(sprintf("[INFO]  Dropping %d sample(s) with missing values in tested variables.\n",
                sum(!complete_idx)))
    md_adonis <- md_adonis[complete_idx, , drop=FALSE]
  }

  if (nrow(md_adonis) < 3) {
    cat("[WARN]  <3 samples remain for Adonis after NA filtering — skipping.\n")
  } else {
    # ── Multivariable Adonis (Type III + univariate) ─────────────────────
    adonis_xlsx <- file.path(ALL_SAMP_DIR, "Adonis_multivariable_PERMANOVA.xlsx")
    wb <- createWorkbook()
    rhs_full <- paste(sprintf("\`%s\`", VARS), collapse = " + ")

    for (metric in METRICS) {
      label <- LABELS[metric]
      dm_full <- DMS[[metric]]
      ids_md  <- md_adonis[["sample-id"]]
      ids_dm  <- attr(dm_full, "Labels")
      common  <- intersect(ids_md, ids_dm)
      m_mat   <- as.matrix(dm_full)
      dm_use  <- as.dist(m_mat[common, common])
      md_use  <- md_adonis[match(common, md_adonis[["sample-id"]]), , drop=FALSE]

      # Type III (margin) — each variable adjusted for the others
      adn3 <- tryCatch(
        adonis2(as.formula(sprintf("dm_use ~ %s", rhs_full)),
                data=md_use, permutations=PERMS, by="margin"),
        error = function(e) {
          cat(sprintf("[WARN]  Type III adonis failed for %s: %s\n", label,
                      conditionMessage(e))); NULL })
      # Type I sequential (just to also report cumulative R^2 path)
      adn1 <- tryCatch(
        adonis2(as.formula(sprintf("dm_use ~ %s", rhs_full)),
                data=md_use, permutations=PERMS, by="terms"),
        error = function(e) NULL)
      # Univariate (each var on its own)
      uni_rows <- list()
      for (v in VARS) {
        au <- tryCatch(
          adonis2(as.formula(sprintf("dm_use ~ \`%s\`", v)),
                  data=md_use, permutations=PERMS, by="terms"),
          error = function(e) NULL)
        if (!is.null(au)) {
          uni_rows[[v]] <- data.frame(
            Variable = v,
            Test     = "Univariate (each variable alone)",
            n_samples = nrow(md_use),
            df       = au\$Df[1],
            Pseudo_F = round(au\$F[1], 4),
            R_squared = round(au\$R2[1], 4),
            Variance_explained_pct = round(au\$R2[1] * 100, 2),
            p_value  = signif(au\$\`Pr(>F)\`[1], 6),
            Significant_at_05 = !is.na(au\$\`Pr(>F)\`[1]) && au\$\`Pr(>F)\`[1] < 0.05,
            stringsAsFactors = FALSE)
        }
      }

      # Combine into a single sheet per metric
      out_blocks <- list()

      if (!is.null(adn3)) {
        a3 <- as.data.frame(adn3)
        a3\$Predictor <- rownames(a3)
        a3 <- a3[!a3\$Predictor %in% c("Total"), , drop=FALSE]
        a3 <- a3[, c("Predictor","Df","SumOfSqs","R2","F","Pr(>F)")]
        names(a3) <- c("Predictor","Df","SumOfSqs","R_squared","Pseudo_F","p_value")
        a3\$Variance_explained_pct <- round(a3\$R_squared * 100, 2)
        a3\$Test <- "Type III (margin) — each var ADJUSTED for the others"
        a3\$Significant_at_05 <- !is.na(a3\$p_value) & a3\$p_value < 0.05
        a3 <- a3[, c("Test","Predictor","Df","SumOfSqs","Pseudo_F",
                     "R_squared","Variance_explained_pct","p_value","Significant_at_05")]
        out_blocks\$type3 <- a3
      }
      if (!is.null(adn1)) {
        a1 <- as.data.frame(adn1)
        a1\$Predictor <- rownames(a1)
        a1 <- a1[!a1\$Predictor %in% c("Total"), , drop=FALSE]
        a1 <- a1[, c("Predictor","Df","SumOfSqs","R2","F","Pr(>F)")]
        names(a1) <- c("Predictor","Df","SumOfSqs","R_squared","Pseudo_F","p_value")
        a1\$Variance_explained_pct <- round(a1\$R_squared * 100, 2)
        a1\$Test <- "Type I (sequential) — order matters"
        a1\$Significant_at_05 <- !is.na(a1\$p_value) & a1\$p_value < 0.05
        a1 <- a1[, c("Test","Predictor","Df","SumOfSqs","Pseudo_F",
                     "R_squared","Variance_explained_pct","p_value","Significant_at_05")]
        out_blocks\$type1 <- a1
      }
      if (length(uni_rows) > 0) {
        u <- do.call(rbind, uni_rows); rownames(u) <- NULL
        names(u)[names(u) == "Variable"] <- "Predictor"
        u\$Test <- "Univariate (each variable alone)"
        u\$SumOfSqs <- NA_real_
        u <- u[, c("Test","Predictor","df","SumOfSqs","Pseudo_F",
                   "R_squared","Variance_explained_pct","p_value","Significant_at_05")]
        names(u)[3] <- "Df"
        out_blocks\$uni <- u
      }

      sheet_df <- do.call(rbind, out_blocks)
      sheet_name <- substr(label, 1, 31)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, sheet_df, headerStyle = createStyle(textDecoration="bold"))
      freezePane(wb, sheet_name, firstActiveRow = 2)
      setColWidths(wb, sheet_name, cols = 1:ncol(sheet_df), widths = "auto")

      cat(sprintf("[OK]    Adonis (%s):  %d rows  →  sheet '%s'\n",
                  label, nrow(sheet_df), sheet_name))
    }
    saveWorkbook(wb, adonis_xlsx, overwrite = TRUE)
    cat(sprintf("[OK]    Multivariable Adonis xlsx → %s\n", adonis_xlsx))

    # ── Distance heatmaps + UPGMA dendrograms (per metric) ───────────────
    annot_df <- md_adonis[, VARS, drop=FALSE]
    rownames(annot_df) <- md_adonis[["sample-id"]]

    # Build named colour list for pheatmap annotations
    annot_colors <- list()
    for (v in VARS) {
      lv <- levels(annot_df[[v]])
      annot_colors[[v]] <- setNames(qual_colours(length(lv)), lv)
    }

    for (metric in METRICS) {
      label <- LABELS[metric];  safe_m <- SAFE_M[metric]
      dm_full <- DMS[[metric]]
      ids_md  <- md_adonis[["sample-id"]]
      common  <- intersect(ids_md, attr(dm_full, "Labels"))
      m_mat   <- as.matrix(dm_full)[common, common]
      annot_use <- annot_df[common, , drop=FALSE]

      # Heatmap — clustered with UPGMA on both axes
      heat_path <- file.path(ALL_SAMP_DIR, sprintf("Distance_heatmap_%s.png", safe_m))
      tryCatch({
        pheatmap(m_mat,
                 clustering_method     = "average",
                 clustering_distance_rows = as.dist(m_mat),
                 clustering_distance_cols = as.dist(m_mat),
                 annotation_row        = annot_use,
                 annotation_col        = annot_use,
                 annotation_colors     = annot_colors,
                 color                 = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
                 main                  = sprintf("%s -- sample-by-sample distances (UPGMA)", label),
                 fontsize              = 9,
                 filename              = heat_path,
                 width                 = max(8, 0.35 * length(common) + 4),
                 height                = max(7, 0.35 * length(common) + 3))
        cat(sprintf("[OK]    Heatmap → %s\n", basename(heat_path)))
      }, error = function(e)
        cat(sprintf("[WARN]  Heatmap failed for %s: %s\n", label, conditionMessage(e))))

      # UPGMA dendrogram (clean, separate from heatmap)
      dend_path <- file.path(ALL_SAMP_DIR, sprintf("UPGMA_dendrogram_%s.png", safe_m))
      tryCatch({
        hc      <- hclust(as.dist(m_mat), method = "average")
        phylo   <- ape::as.phylo(hc)

        # Color tip labels by the variable that explains the most variance
        # for THIS metric (compute univariate Adonis with light permutations).
        primary_var <- VARS[1]
        best_r2 <- -1
        for (v in VARS) {
          au <- tryCatch(
            adonis2(as.formula(sprintf("as.dist(m_mat) ~ \`%s\`", v)),
                    data = annot_use, permutations = 99, by = "terms"),
            error = function(e) NULL)
          if (!is.null(au) && !is.na(au\$R2[1]) && au\$R2[1] > best_r2) {
            best_r2 <- au\$R2[1]; primary_var <- v
          }
        }

        tip_groups <- annot_use[phylo\$tip.label, primary_var]
        tip_pal    <- annot_colors[[primary_var]]
        tip_cols   <- as.character(tip_pal[as.character(tip_groups)])
        tip_cols[is.na(tip_cols)] <- "#777777"

        png(dend_path, width = 9 * 300,
            height = max(5, 0.22 * length(common) + 2) * 300, res = 300, bg = "white")
        par(mar = c(4, 1, 3, max(8, 0.4 * max(nchar(common)) + 2)))
        plot(phylo, type = "phylogram",
             tip.color = tip_cols, edge.width = 1.2,
             label.offset = 0.005,
             main = sprintf("UPGMA dendrogram -- %s   (tips coloured by %s)",
                            label, primary_var),
             cex = 0.7)
        legend("topleft",
               legend = names(tip_pal),
               fill   = unname(tip_pal),
               cex    = 0.7, bty = "n",
               title  = primary_var)
        dev.off()
        cat(sprintf("[OK]    UPGMA dendrogram → %s\n", basename(dend_path)))
      }, error = function(e)
        cat(sprintf("[WARN]  Dendrogram failed for %s: %s\n", label, conditionMessage(e))))
    }
  }
}

if (length(fail_log) > 0) {
  cat("\n[SUMMARY] Failures during analysis:\n")
  for (f in fail_log) cat(sprintf("  - %s\n", f))
}
cat("\n[DONE]  R analysis complete.\n")
RBETA

  if $DRY_RUN; then
    warn "[DRY-RUN] Would run R analysis."
  else
    _R --vanilla "$R_BETA" \
      || err "R beta-diversity analysis failed — see output above."
  fi
  rm -f "$R_BETA"
fi
sep

# =============================================================================
step "8/9 — Reorganise per-variable QZVs (cosmetic — they live in qzv/)"
# =============================================================================
ok "QZV files: <var>/qzv/permanova_<metric>_<var>.qzv  +  permdisp_<metric>_<var>.qzv"
ok "Emperors : $ALL_SAMP_QZV/emperor_<metric>.qzv"
sep

# =============================================================================
step "9/9 — Write mbx_beta_diversity_info.txt"
# =============================================================================
QIIME_VER="$(qiime --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo unknown)"
PY_BIN="$(command -v python || true)"
PY_VER="$($PY_BIN -c 'import sys;print(".".join(map(str,sys.version_info[:3])))' 2>/dev/null || echo unknown)"

# Build per-line entries for variables and metrics
VAR_LINES=""
for col in "${CATEGORICAL_COLS[@]}"; do
  VAR_LINES="${VAR_LINES}VARIABLE=${col}
"
done
METRIC_LINES=""
for m in "${QIIME_METRICS[@]}"; do
  METRIC_LINES="${METRIC_LINES}METRIC=${m}
"
done

cat > "$INFO_FILE" << INFO
# ============================================================================
# mbx_beta_diversity_info.txt
# Generated by mbx_beta_diversity_run.sh
# Date : $NOW_PRETTY
# ============================================================================
#
# Do NOT edit the key=value lines below — they are parsed programmatically.

# ── Inputs ───────────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
ALPHA_INFO=$ALPHA_INFO
PRE_DIVERSITY_INFO=$PRE_DIV_INFO
METADATA_TXT=$METADATA_TXT
RAREFIED_TABLE_QZA=$RAREFIED_QZA
ROOTED_TREE_QZA=$TREE_QZA
SAMPLING_DEPTH=$REC_DEPTH
PERMUTATIONS=$PERMS

# ── Output locations ─────────────────────────────────────────────────────────
BETA_OUT_DIR=$BETA_DIR
PER_VAR_DIR=$PER_VAR_DIR
ALL_SAMP_DIR=$ALL_SAMP_DIR
WORK_DIR=$WORK_DIR
DM_QZA_DIR=$DM_DIR
PCOA_QZA_DIR=$PCOA_DIR
EXPORT_DIR=$EXPORT_DIR
EMPEROR_QZV_DIR=$ALL_SAMP_QZV
ADONIS_XLSX=${ALL_SAMP_DIR}/Adonis_multivariable_PERMANOVA.xlsx
GLOBAL_PERMANOVA_XLSX=${PER_VAR_DIR}/Summary_PERMANOVA_all_variables_all_metrics.xlsx

# ── Distance metrics analysed ────────────────────────────────────────────────
N_METRICS=${N_METRICS}
${METRIC_LINES}

# ── Categorical variables analysed ───────────────────────────────────────────
N_VARIABLES=${#CATEGORICAL_COLS[@]}
${VAR_LINES}

# ── Provenance ───────────────────────────────────────────────────────────────
GENERATED_AT=$NOW_PRETTY
SCRIPT_NAME=mbx_beta_diversity_run.sh
QIIME2_VERSION=$QIIME_VER
PYTHON_VERSION=$PY_VER
R_VERSION=$R_VERSION
N_JOBS=$N_JOBS
SKIPPED_QZVS=$SKIP_QZVS
SKIPPED_STATS=$SKIP_STATS
SKIPPED_HEATMAP=$SKIP_HEATMAP
INVOCATION_USER=${USER:-unknown}
INVOCATION_CWD=$(pwd)
INVOCATION_ARGV=$0 $*
INFO

ok "Info file → $INFO_FILE"
sep

# ── Final banner ──────────────────────────────────────────────────────────────
echo ""
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║  STEP 13 COMPLETE  —  BETA DIVERSITY                        ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "  Per-variable results : $PER_VAR_DIR/<Variable>/"
echo "  All-samples results  : $ALL_SAMP_DIR/"
echo "  Working files        : $WORK_DIR/"
echo "  Variables            : ${CATEGORICAL_COLS[*]}"
echo "  Metrics              : ${QIIME_METRICS[*]}"
echo ""
echo "  Each <Variable>/ folder contains:"
echo "    • PCoA_<Metric>_<Variable>.png            (one per metric, with 95% ellipses)"
echo "    • PCoA_panel_4metrics_<Variable>.png      (publication-style 4-metric panel)"
echo "    • Boxplot_DistanceToCentroid_*.png        (visual PERMDISP)"
echo "    • PERMANOVA_results_<Variable>.xlsx       (Pseudo-F, R^2, p — per metric)"
echo "    • Pairwise_PERMANOVA_<Variable>.xlsx      (only if k>2)"
echo "    • PERMDISP_results_<Variable>.xlsx        (homogeneity-of-dispersion test)"
echo "    • qzv/permanova_<Metric>_<Variable>.qzv   (interactive QIIME2 view)"
echo "    • qzv/permdisp_<Metric>_<Variable>.qzv"
echo ""
echo "  $(basename "$ALL_SAMP_DIR")/ contains:"
echo "    • Adonis_multivariable_PERMANOVA.xlsx     (Type III + univariate)"
echo "    • Distance_heatmap_<Metric>.png           (UPGMA-clustered, annotated)"
echo "    • UPGMA_dendrogram_<Metric>.png           (clean tree, tips coloured by primary var)"
echo "    • qzv/emperor_<Metric>.qzv                (interactive 3D PCoA)"
echo ""
echo "  Open the global summary:"
echo "    open '$PER_VAR_DIR/Summary_PERMANOVA_all_variables_all_metrics.xlsx'"
echo "    open '$ALL_SAMP_DIR/Adonis_multivariable_PERMANOVA.xlsx'"
echo ""
echo "  Provenance (machine-readable) → $INFO_FILE"
echo ""
