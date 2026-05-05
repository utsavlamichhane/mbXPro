#!/usr/bin/env bash
# =============================================================================
#  mbx_network_run.sh   (step 17)
#  Co-occurrence / correlation networks for microbial taxa
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  PURPOSE
#    For each (taxonomic_level x categorical_variable) combination, build:
#      1) a GLOBAL network using all samples (the "ecological backbone")
#      2) one PER-GROUP network for each level of the categorical variable
#         (only when N_samples_in_group >= --min-group-n)
#    Plus a side-by-side group_comparison.xlsx + multi-panel plot.
#
#  ALGORITHM (compositional-aware, Gloor et al. 2017)
#    1. Prevalence filter   (drop taxa present in < --prevalence-threshold)
#    2. Pseudocount + CLR transform                  (centered log-ratio)
#    3. Spearman correlation on CLR values           (rank-based, robust)
#    4. Benjamini-Hochberg FDR adjustment of pairwise p-values
#    5. Filter edges: |rho| >= --rho-threshold AND q <= --q-threshold
#    6. Build undirected graph; compute node-level metrics
#    7. Module detection via Louvain (max-modularity)
#    8. Hub taxa: top 10 % by combined degree + betweenness + hub_score
#
#  WHY this approach?
#    - SparCC is the textbook "compositional-aware" method but ships only as
#      Python; SpiecEasi is excellent but pulls in heavy Bioconductor deps.
#    - CLR-transform + Spearman is the recommended lightweight alternative
#      (Gloor et al. 2017 "Microbiome datasets are compositional") that
#      properly handles the simplex constraint without iterative bootstrapping
#      and runs entirely on CRAN packages.
#
#  OUTPUTS PER NETWORK
#    edges.tsv, nodes.tsv, network.graphml,
#    network_plot.{png,pdf}, network_summary.xlsx,
#    hub_taxa.xlsx, modules.tsv
#
#  GATING (reads from previous steps' info files)
#    7_taxonomy_csv/mbx_taxonomy_info.txt   -> METADATA_TXT
#    8_cleaned_files/mbx_ezclean_info.txt   -> 7 cleaned xlsx paths
#      (with disk-discovery fallback if info file is missing)
#
#  OUTPUT STRUCTURE
#    17_co_occurrence_networks/
#      working_dir_networks/
#      global_networks/
#         network_<level>/                        (1 per level requested)
#      <variable>/
#         network_<level>_by_<variable>/
#            per_group/<group>/...                (1 per group)
#            group_comparison.xlsx
#            multi_panel_plot.{png,pdf}
#      mbx_networks_info.txt
#
#  RUNTIME
#    ~5-15 s per network on N=20.  Default scope (g/f/s x N_vars x ~3 groups
#    + global) ~ 30-60 networks => 2-10 minutes total.
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

_sanitize_dirname() {
  printf '%s' "$1" \
    | tr ' ' '_' \
    | tr -d '()[]{}/<>|\\:*?"' \
    | sed 's/^[_.-]*//;s/[_.-]*$//'
}

# ── R env wrapper — strip conda's R_LIBS_USER pollution ──────────────────────
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}
_R() { _strip_env "$RSCRIPT_CMD" "$@"; }

# ── Level mapping helpers ────────────────────────────────────────────────────
_level_full() {
  case "$1" in
    d) echo "domain" ;;  p) echo "phylum" ;;  c) echo "class" ;;
    o) echo "order"  ;;  f) echo "family" ;;  g) echo "genus" ;;
    s) echo "species" ;;
  esac
}
_level_dirname() {
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

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_network_run.sh  --  Co-occurrence / correlation networks (step 17)

USAGE:
  mbx_network_run.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Builds CLR-transformed Spearman correlation networks of taxa, with
  Benjamini-Hochberg FDR correction.  For every requested taxonomic level:
    - one GLOBAL network using all samples
    - one PER-GROUP network for each level of every categorical variable
      (skipped automatically when group has fewer than --min-group-n samples)
  Plus side-by-side multi-panel plot and group_comparison.xlsx per variable.

OPTIONS:
  --levels <list>            Comma-separated subset of {d,p,c,o,f,g,s}
                             (default: g,f,s -- only the levels with enough
                             nodes for meaningful networks)
  --rho-threshold <val>      Min |Spearman rho| to keep an edge   (default: 0.6)
  --q-threshold <val>        Max BH-adjusted q-value for an edge  (default: 0.05)
  --prevalence-threshold <p> Min fraction of samples a taxon must
                             be present in to be retained         (default: 0.30)
  --min-group-n <N>          Min samples per group to compute a
                             per-group network                    (default: 5)
  --skip-per-group           Compute global networks only (faster)
  --skip-install             Don't auto-install missing R packages
  --force-rerun              Recompute even if outputs already exist
  --dry-run                  Print planned R calls without running them
  -h, --help                 Show this help and exit

OUTPUT STRUCTURE
  17_co_occurrence_networks/
    global_networks/
      network_<level>/  edges.tsv, nodes.tsv, network.graphml,
                        network_plot.{png,pdf}, network_summary.xlsx,
                        hub_taxa.xlsx, modules.tsv
    <variable>/
      network_<level>_by_<variable>/
        per_group/<group>/...   (same files as above, per group)
        group_comparison.xlsx
        multi_panel_plot.{png,pdf}
    mbx_networks_info.txt

EXAMPLES
  mbx_network_run.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_network_run.sh /path/to/mbX_pro_outputs_20260417_121431 --levels g,f
  mbx_network_run.sh /path/to/mbX_pro_outputs_20260417_121431 \
    --rho-threshold 0.5 --q-threshold 0.10 --skip-per-group
EOF
  exit 0
}

# ── Argument parsing ──────────────────────────────────────────────────────────
MBX_OUT_DIR=""
LEVELS_ARG=""
RHO_THRESHOLD=0.6
Q_THRESHOLD=0.05
PREV_THRESHOLD=0.30
MIN_GROUP_N=5
SKIP_PER_GROUP=false
SKIP_INSTALL=false
FORCE_RERUN=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)              usage ;;
    --levels)               LEVELS_ARG="$2";        shift 2 ;;
    --rho-threshold)        RHO_THRESHOLD="$2";     shift 2 ;;
    --q-threshold)          Q_THRESHOLD="$2";       shift 2 ;;
    --prevalence-threshold) PREV_THRESHOLD="$2";    shift 2 ;;
    --min-group-n)          MIN_GROUP_N="$2";       shift 2 ;;
    --skip-per-group)       SKIP_PER_GROUP=true;     shift ;;
    --skip-install)         SKIP_INSTALL=true;       shift ;;
    --force-rerun)          FORCE_RERUN=true;        shift ;;
    --dry-run)              DRY_RUN=true;            shift ;;
    -*)                     err "Unknown option: $1\n  Run: mbx_network_run.sh --help" ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Multiple positional arguments — only one MBX_OUT_DIR expected.\n  Got extra: $1"
      fi
      ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "Missing required argument: <mbX_pro_outputs_dir>
  Run: mbx_network_run.sh --help"
[[ -d "$MBX_OUT_DIR" ]] || err "Not a directory: $MBX_OUT_DIR
  Did you pass the path to mbX_pro_outputs_<TIMESTAMP>/?"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# Parse levels list — default = g, f, s (most ecologically meaningful for nets)
if [[ -z "$LEVELS_ARG" ]]; then
  LEVELS=("g" "f" "s")
else
  _OLDIFS="$IFS"; IFS=','; LEVELS=( $LEVELS_ARG ); IFS="$_OLDIFS"
  for _l in "${LEVELS[@]}"; do
    case "$_l" in
      d|p|c|o|f|g|s) ;;
      *) err "Unknown level '$_l' in --levels.  Allowed: d,p,c,o,f,g,s" ;;
    esac
  done
fi

# Auto-detect CPU count
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

info "Metadata file              : $METADATA_TXT"
info "Cleaned files dir          : $CLEANED_DIR"
[[ -f "$EZCLEAN_INFO" ]] && info "ezclean info file          : $EZCLEAN_INFO"
[[ -f "$EZCLEAN_INFO" ]] || warn "ezclean info file missing  -> will use disk discovery"
info "Levels requested           : ${LEVELS[*]}"
info "rho threshold              : $RHO_THRESHOLD"
info "q-value threshold (BH)     : $Q_THRESHOLD"
info "Prevalence threshold       : $PREV_THRESHOLD"
info "Min samples per group      : $MIN_GROUP_N"
info "Skip per-group networks    : $SKIP_PER_GROUP"
info "Force re-run               : $FORCE_RERUN"
info "Dry-run                    : $DRY_RUN"
info "CPU cores available        : $N_JOBS"

# Map each requested level letter to its cleaned xlsx path (with fallback)
LEVEL_XLSX_PATHS=()
for _lvl in "${LEVELS[@]}"; do
  _UPLVL="$(echo "$_lvl" | tr '[:lower:]' '[:upper:]')"
  _DIRNAME="$(_level_dirname "$_lvl")"
  _XLSX=""

  if [[ -f "$EZCLEAN_INFO" ]]; then
    _XLSX="$(_read_key "LEVEL_${_UPLVL}_XLSX" "$EZCLEAN_INFO")"
    [[ "$_XLSX" == "FAILED" || "$_XLSX" == "NOT_FOUND" ]] && _XLSX=""
  fi

  if [[ -z "$_XLSX" || ! -f "$_XLSX" ]]; then
    for _src_dir in "${CLEANED_DIR}/${_DIRNAME}_level-7" "${CLEANED_DIR}/${_DIRNAME}"; do
      if [[ -d "$_src_dir" ]]; then
        _candidate="$(ls "$_src_dir"/*.xlsx 2>/dev/null | head -1)"
        if [[ -n "$_candidate" && -f "$_candidate" ]]; then
          _XLSX="$_candidate"; break
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
# Required:
#   igraph        graph operations + Louvain modularity + base plots
#   psych         corr.test() returns p-value matrix alongside correlation matrix
#   openxlsx      xlsx I/O
#   ggplot2       comparison plots + multi-panel
#   RColorBrewer  module / sign palettes

_TMPID="${$}_$(date +%s)"
PKG_CHECK_R="/tmp/mbx_net_pkgchk_${_TMPID}.R"
trap 'rm -f /tmp/mbx_net_*_${_TMPID}*.R /tmp/mbx_net_*_${_TMPID}*.txt /tmp/mbx_net_*_${_TMPID}*.out /tmp/mbx_net_*_${_TMPID}*.err' EXIT

cat > "$PKG_CHECK_R" << 'RPKG'
required <- c("igraph", "psych", "openxlsx", "ggplot2", "RColorBrewer")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) == 0) {
  cat("[OK]    All required R packages already installed.\n")
  cat(sprintf("[INFO]  igraph version: %s\n",
              tryCatch(packageVersion("igraph"), error = function(e) "unknown")))
  quit(status = 0)
}
cat(sprintf("[INFO]  Missing R packages: %s\n", paste(missing, collapse = ", ")))
quit(status = 2)
RPKG

if $DRY_RUN; then
  warn "[DRY-RUN] Would check/install R packages: igraph psych openxlsx ggplot2 RColorBrewer"
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
       Rscript -e 'install.packages(c(\"igraph\",\"psych\",\"openxlsx\",\"ggplot2\",\"RColorBrewer\"))'
  -> Then re-run this script."
    fi
    info "Installing missing R packages from CRAN (may take 1-3 minutes)..."
    PKG_INSTALL_R="/tmp/mbx_net_pkginst_${_TMPID}.R"
    cat > "$PKG_INSTALL_R" << 'RPKGI'
required <- c("igraph", "psych", "openxlsx", "ggplot2", "RColorBrewer")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org",
                   quiet = FALSE, dependencies = TRUE)
  failed <- missing[!sapply(missing, requireNamespace, quietly = TRUE)]
  if (length(failed) > 0) stop(sprintf("Failed to install: %s",
                                       paste(failed, collapse = ", ")))
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
# Same sentinel-based detection used in step 16: prints column names ONLY
# between MBX_CATS_BEGIN / MBX_CATS_END markers on stdout.
CAT_DETECT_R="/tmp/mbx_net_detect_${_TMPID}.R"
CATS_STDOUT="/tmp/mbx_net_cats_${_TMPID}.out"
CATS_STDERR="/tmp/mbx_net_cats_${_TMPID}.err"

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

  CATEGORICAL_COLS=()
  _IN_BLOCK=false
  while IFS= read -r line; do
    if [[ "$line" == "MBX_CATS_BEGIN" ]]; then _IN_BLOCK=true;  continue; fi
    if [[ "$line" == "MBX_CATS_END"   ]]; then _IN_BLOCK=false; continue; fi
    [[ "$_IN_BLOCK" == true ]] || continue
    [[ -z "$line" ]] && continue
    [[ "$line" == "NULL" ]] && continue
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
NET_DIR="${MBX_OUT_DIR}/17_co_occurrence_networks"
WORK_DIR="${NET_DIR}/working_dir_networks"
GLOBAL_DIR="${NET_DIR}/global_networks"

if $DRY_RUN; then
  info "[DRY-RUN] Would create: $NET_DIR/"
  info "[DRY-RUN] Would create: $WORK_DIR/, $GLOBAL_DIR/"
  if ! $SKIP_PER_GROUP; then
    for col in "${CATEGORICAL_COLS[@]}"; do
      info "[DRY-RUN] Would create: 17_co_occurrence_networks/$(_sanitize_dirname "$col")/"
    done
  fi
else
  mkdir -p "$WORK_DIR" "$GLOBAL_DIR" \
    || err "Could not create output directories under $NET_DIR — check permissions."
  if ! $SKIP_PER_GROUP; then
    for col in "${CATEGORICAL_COLS[@]}"; do
      _SC="$(_sanitize_dirname "$col")"
      mkdir -p "${NET_DIR}/${_SC}" \
        || err "Could not create: ${NET_DIR}/${_SC}"
      info "  Created: 17_co_occurrence_networks/${_SC}/"
    done
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "6/7 — Build co-occurrence networks"
# ─────────────────────────────────────────────────────────────────────────────

# Counters
_GLOBAL_PASS=0; _GLOBAL_FAIL=0; _GLOBAL_SKIP=0
_GROUP_PASS=0;  _GROUP_FAIL=0;  _GROUP_SKIP=0
RESULT_LINES=""

# ── R analysis script template (used for every network) ─────────────────────
# Bash injects per-call parameters via env vars so we write the template once
# and reuse it for every (level x scope) combination.  This template handles
# BOTH "global" mode (all samples, no group filter) and "per_group" mode
# (only samples in one group of one variable) via env vars.
R_SCRIPT_TEMPLATE="/tmp/mbx_net_run_${_TMPID}.R"
cat > "$R_SCRIPT_TEMPLATE" << 'RTPL'
suppressPackageStartupMessages({
  library(openxlsx)
  library(igraph)
  library(psych)
  library(ggplot2)
  library(RColorBrewer)
})

CLEANED_XLSX <- Sys.getenv("MBX_NET_CLEANED_XLSX")
META_TXT     <- Sys.getenv("MBX_NET_META_TXT")
LEVEL_LET    <- Sys.getenv("MBX_NET_LEVEL_LET")
LEVEL_NAME   <- Sys.getenv("MBX_NET_LEVEL_NAME")
SCOPE        <- Sys.getenv("MBX_NET_SCOPE")            # "global" or "per_group"
VAR_COL      <- Sys.getenv("MBX_NET_VAR_COL", "")      # only used in per_group
GROUP        <- Sys.getenv("MBX_NET_GROUP", "")        # only used in per_group
OUT_DIR      <- Sys.getenv("MBX_NET_OUT_DIR")
RHO_THR      <- as.numeric(Sys.getenv("MBX_NET_RHO_THR", "0.6"))
Q_THR        <- as.numeric(Sys.getenv("MBX_NET_Q_THR",   "0.05"))
PREV_THR     <- as.numeric(Sys.getenv("MBX_NET_PREV_THR","0.30"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
context_str <- if (SCOPE == "global") {
  sprintf("Global (all samples) -- %s", LEVEL_NAME)
} else {
  sprintf("%s = %s -- %s", VAR_COL, GROUP, LEVEL_NAME)
}

# ── 1. Read cleaned table ───────────────────────────────────────────────────
df <- read.xlsx(CLEANED_XLSX, check.names = FALSE)
sid_col <- names(df)[1]
sids    <- as.character(df[[sid_col]])

# Identify metadata vs taxa columns by intersecting with the metadata file
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
meta_cols <- intersect(names(meta), names(df))

# Filter to requested scope
if (SCOPE == "per_group") {
  if (!(VAR_COL %in% names(df))) {
    if (VAR_COL %in% names(meta)) {
      df <- merge(df, meta[, c(names(meta)[1], VAR_COL), drop = FALSE],
                  by.x = sid_col, by.y = names(meta)[1], all.x = TRUE)
      meta_cols <- c(meta_cols, VAR_COL)
    } else {
      stop(sprintf("Variable '%s' not found in cleaned xlsx OR metadata.", VAR_COL))
    }
  }
  keep <- !is.na(df[[VAR_COL]]) &
          trimws(as.character(df[[VAR_COL]])) == GROUP
  if (sum(keep) < 4) {
    stop(sprintf("Group '%s' has only %d samples (need >= 4 for network).",
                 GROUP, sum(keep)))
  }
  df   <- df[keep, , drop = FALSE]
  sids <- sids[keep]
}

feat_cols <- setdiff(names(df), c(sid_col, meta_cols))
if (length(feat_cols) < 3)
  stop(sprintf("Only %d taxa columns at level '%s' -- need >= 3 for a network.",
               length(feat_cols), LEVEL_LET))

X <- as.data.frame(lapply(df[, feat_cols, drop = FALSE], function(v) {
  v <- suppressWarnings(as.numeric(as.character(v)))
  v[is.na(v)] <- 0
  v
}), check.names = FALSE)
rownames(X) <- sids
n_samp <- nrow(X)

# ── 2. Prevalence filter ────────────────────────────────────────────────────
prev <- colSums(X > 0) / n_samp
keep_taxa <- prev >= PREV_THR
X <- X[, keep_taxa, drop = FALSE]
prev <- prev[keep_taxa]
if (ncol(X) < 3) {
  cat(sprintf("[WARN]  After prevalence filter (>= %.0f%%), only %d taxa remain.  Skipping.\n",
              100 * PREV_THR, ncol(X)))
  # Write a marker so bash knows we deliberately skipped
  writeLines(sprintf("Skipped: too few taxa after prevalence filter (%d)",
                     ncol(X)), file.path(OUT_DIR, "_skipped.txt"))
  cat(sprintf("MBX_NET_SUMMARY\tcontext=%s\tstatus=SKIPPED_FEW_TAXA\tn_samples=%d\tn_taxa=%d\n",
              context_str, n_samp, ncol(X)))
  quit(status = 0)
}

mean_abun <- colMeans(X)

# ── 3. CLR transform (Aitchison) ────────────────────────────────────────────
# Pseudocount = 0.5 of the smallest non-zero value (or 0.5 absolute, whichever
# is smaller).  Then x_clr = log(x_i) - mean_log(x_row).
nonzero_min <- min(X[X > 0], na.rm = TRUE)
pseudo <- min(0.5, nonzero_min * 0.5)
Xp <- X + pseudo
log_X <- log(as.matrix(Xp))
clr_X <- sweep(log_X, 1, rowMeans(log_X), "-")

# ── 4. Spearman + p-values + BH adjustment ──────────────────────────────────
# psych::corr.test gives both r and p (raw) matrices.
# Suppress its 'sample size in pairs' message; keep ci=FALSE for speed.
cor_obj <- suppressWarnings(
  psych::corr.test(clr_X, use = "pairwise", method = "spearman",
                   adjust = "none", ci = FALSE)
)
R_mat <- cor_obj$r
P_mat <- cor_obj$p   # raw p-values (lower triangle = adjusted; we re-do this)

# Get all unique pair indices and BH-adjust ourselves
n_taxa <- ncol(R_mat)
ix <- which(upper.tri(R_mat), arr.ind = TRUE)
edge_df <- data.frame(
  taxon_1 = colnames(R_mat)[ix[, 1]],
  taxon_2 = colnames(R_mat)[ix[, 2]],
  rho     = R_mat[upper.tri(R_mat)],
  p_raw   = P_mat[upper.tri(P_mat)],
  stringsAsFactors = FALSE
)
edge_df$q_BH <- p.adjust(edge_df$p_raw, method = "BH")
edge_df$sign <- ifelse(edge_df$rho > 0, "positive", "negative")

# Significant edges
sig_df <- edge_df[abs(edge_df$rho) >= RHO_THR & edge_df$q_BH <= Q_THR, , drop = FALSE]
sig_df <- sig_df[order(-abs(sig_df$rho)), , drop = FALSE]

cat(sprintf("[INFO]  %s\n", context_str))
cat(sprintf("[INFO]  N samples = %d, N taxa after prev filter = %d, candidate pairs = %d\n",
            n_samp, n_taxa, nrow(edge_df)))
cat(sprintf("[INFO]  Significant edges  : %d  (|rho| >= %.2f, q <= %.3f)\n",
            nrow(sig_df), RHO_THR, Q_THR))

# Write full edge table even if no significant edges (for transparency)
write.table(sig_df, file.path(OUT_DIR, "edges.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 5. Build graph + node metrics ───────────────────────────────────────────
if (nrow(sig_df) == 0) {
  cat("[WARN]  No significant edges -- network is empty.\n")
  writeLines("No significant edges.", file.path(OUT_DIR, "_no_edges.txt"))
  cat(sprintf("MBX_NET_SUMMARY\tcontext=%s\tstatus=NO_EDGES\tn_samples=%d\tn_taxa=%d\tn_edges=0\n",
              context_str, n_samp, n_taxa))
  quit(status = 0)
}

g <- graph_from_data_frame(
  d = sig_df[, c("taxon_1", "taxon_2", "rho", "q_BH", "sign")],
  directed = FALSE,
  vertices = data.frame(name = colnames(R_mat), stringsAsFactors = FALSE)
)
# Drop isolated nodes (those without any significant edge)
isolated <- which(degree(g) == 0)
if (length(isolated) > 0) g <- delete_vertices(g, isolated)
n_nodes <- vcount(g)
n_edges <- ecount(g)
if (n_nodes < 2 || n_edges < 1) {
  cat("[WARN]  Network too small after pruning isolated nodes.\n")
  cat(sprintf("MBX_NET_SUMMARY\tcontext=%s\tstatus=TOO_SMALL\tn_samples=%d\tn_taxa=%d\tn_edges=%d\n",
              context_str, n_samp, n_taxa, n_edges))
  quit(status = 0)
}

# Node-level metrics
deg     <- degree(g)
betw    <- betweenness(g, normalized = TRUE)
clos    <- tryCatch(closeness(g, normalized = TRUE), error = function(e) rep(NA, n_nodes))
hub     <- tryCatch(hub_score(g)$vector, error = function(e) rep(NA, n_nodes))

# Module detection -- Louvain (max modularity)
modules_obj <- cluster_louvain(g, weights = abs(E(g)$rho))
mod_membership <- membership(modules_obj)
modularity_val <- modularity(modules_obj)

node_df <- data.frame(
  taxon          = V(g)$name,
  degree         = deg,
  betweenness    = round(betw, 5),
  closeness      = round(clos, 5),
  hub_score      = round(hub,  5),
  mean_abundance = round(mean_abun[V(g)$name], 5),
  prevalence     = round(prev[V(g)$name],      3),
  module         = as.integer(mod_membership),
  stringsAsFactors = FALSE
)
node_df <- node_df[order(-node_df$degree, -node_df$betweenness), , drop = FALSE]

write.table(node_df, file.path(OUT_DIR, "nodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Hub taxa: top 10 % by combined rank of degree + betweenness + hub_score
combo_score <- rank(-deg) + rank(-betw) + rank(-hub)
hub_cutoff  <- max(2, ceiling(0.1 * n_nodes))
hub_idx     <- order(combo_score)[seq_len(min(20, hub_cutoff, n_nodes))]
hub_df      <- node_df[node_df$taxon %in% V(g)$name[hub_idx], , drop = FALSE]
hub_df      <- hub_df[order(-hub_df$degree), , drop = FALSE]
write.xlsx(hub_df, file.path(OUT_DIR, "hub_taxa.xlsx"), overwrite = TRUE)

# Modules export
modules_tsv <- data.frame(
  taxon  = V(g)$name,
  module = as.integer(mod_membership),
  stringsAsFactors = FALSE
)
modules_tsv <- modules_tsv[order(modules_tsv$module, modules_tsv$taxon), , drop = FALSE]
write.table(modules_tsv, file.path(OUT_DIR, "modules.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Network-level summary
n_pos <- sum(E(g)$sign == "positive")
n_neg <- sum(E(g)$sign == "negative")
sum_df <- data.frame(
  metric = c("context", "n_samples", "n_taxa_after_prev_filter",
             "n_nodes", "n_edges", "n_positive_edges", "n_negative_edges",
             "density", "mean_degree", "modularity", "n_modules",
             "global_transitivity", "rho_threshold", "q_threshold",
             "prevalence_threshold"),
  value  = c(context_str, as.character(n_samp), as.character(n_taxa),
             as.character(n_nodes), as.character(n_edges),
             as.character(n_pos), as.character(n_neg),
             sprintf("%.4f", edge_density(g)),
             sprintf("%.3f", mean(deg)),
             sprintf("%.4f", modularity_val),
             as.character(length(unique(mod_membership))),
             sprintf("%.4f", transitivity(g, type = "global")),
             sprintf("%.2f", RHO_THR), sprintf("%.3f", Q_THR),
             sprintf("%.2f", PREV_THR)),
  stringsAsFactors = FALSE
)
write.xlsx(sum_df, file.path(OUT_DIR, "network_summary.xlsx"),
           overwrite = TRUE)

# ── 6. GraphML export (Cytoscape / Gephi) ───────────────────────────────────
# Attach attributes that downstream visualisers care about
V(g)$mean_abundance <- mean_abun[V(g)$name]
V(g)$prevalence     <- prev[V(g)$name]
V(g)$module         <- as.integer(mod_membership)
V(g)$degree         <- deg
V(g)$is_hub         <- V(g)$name %in% hub_df$taxon
write_graph(g, file.path(OUT_DIR, "network.graphml"), format = "graphml")

# ── 7. Plot ─────────────────────────────────────────────────────────────────
# Layout: Fruchterman-Reingold seeded for reproducibility.
set.seed(42)
lay <- layout_with_fr(g, weights = abs(E(g)$rho))

# Edge style by sign (red = negative, blue = positive)
ecol <- ifelse(E(g)$sign == "positive", "#3b73af", "#c0392b")
ewid <- 0.3 + 2.5 * (abs(E(g)$rho) - RHO_THR) / max(1e-6, 1 - RHO_THR)

# Node style: size by mean abundance (sqrt-scaled), color by module
# Cap node size between 3 and 14
abu <- mean_abun[V(g)$name]
nsize <- 3 + 11 * sqrt((abu - min(abu)) / max(1e-6, diff(range(abu))))
nsize[is.na(nsize)] <- 4

n_mod <- length(unique(mod_membership))
mod_pal <- if (n_mod <= 12) {
  RColorBrewer::brewer.pal(max(3, min(12, n_mod)), "Set3")[seq_len(n_mod)]
} else {
  rep(RColorBrewer::brewer.pal(12, "Set3"), length.out = n_mod)
}
ncol_v <- mod_pal[as.integer(mod_membership)]

# Label only hub taxa (truncate long names)
labels <- ifelse(V(g)$name %in% hub_df$taxon,
                 substr(V(g)$name, 1, 22), NA)

plot_title <- sprintf("%s\n%d nodes, %d edges (%d+ / %d-)\nmodularity = %.3f, density = %.3f",
                      context_str, n_nodes, n_edges, n_pos, n_neg,
                      modularity_val, edge_density(g))

png(file.path(OUT_DIR, "network_plot.png"),
    width = 11, height = 9, units = "in", res = 300)
par(mar = c(1, 1, 4, 1))
plot(g, layout = lay,
     vertex.color = ncol_v, vertex.size = nsize,
     vertex.label = labels, vertex.label.cex = 0.7,
     vertex.label.color = "black", vertex.label.family = "sans",
     vertex.frame.color = "grey30",
     edge.color = ecol, edge.width = ewid,
     main = plot_title)
legend("bottomleft", legend = c("positive (rho > 0)", "negative (rho < 0)"),
       col = c("#3b73af", "#c0392b"), lty = 1, lwd = 2, bty = "n",
       cex = 0.85)
dev.off()

pdf(file.path(OUT_DIR, "network_plot.pdf"), width = 11, height = 9)
par(mar = c(1, 1, 4, 1))
plot(g, layout = lay,
     vertex.color = ncol_v, vertex.size = nsize,
     vertex.label = labels, vertex.label.cex = 0.7,
     vertex.label.color = "black", vertex.label.family = "sans",
     vertex.frame.color = "grey30",
     edge.color = ecol, edge.width = ewid,
     main = plot_title)
legend("bottomleft", legend = c("positive (rho > 0)", "negative (rho < 0)"),
       col = c("#3b73af", "#c0392b"), lty = 1, lwd = 2, bty = "n",
       cex = 0.85)
dev.off()

cat(sprintf("[OK]    Network written: %d nodes, %d edges, modularity = %.3f\n",
            n_nodes, n_edges, modularity_val))
cat(sprintf("MBX_NET_SUMMARY\tcontext=%s\tstatus=OK\tn_samples=%d\tn_taxa=%d\tn_nodes=%d\tn_edges=%d\tn_pos=%d\tn_neg=%d\tdensity=%.4f\tmodularity=%.4f\tn_modules=%d\thubs=%s\n",
            context_str, n_samp, n_taxa, n_nodes, n_edges, n_pos, n_neg,
            edge_density(g), modularity_val, length(unique(mod_membership)),
            paste(head(hub_df$taxon, 5), collapse = "|")))
RTPL

# Helper: invoke the R template with the right env vars
_run_one_network() {
  local cleaned_xlsx="$1" level_let="$2" level_name="$3"
  local scope="$4" var_col="$5" group="$6" out_dir="$7"
  local marker="${out_dir}/network_summary.xlsx"
  local marker2="${out_dir}/_no_edges.txt"
  local marker3="${out_dir}/_skipped.txt"

  if [[ ( -f "$marker" || -f "$marker2" || -f "$marker3" ) && "$FORCE_RERUN" == false ]]; then
    skipped "    Existing output (use --force-rerun to recompute)"
    return 2   # 2 = skipped
  fi

  if $DRY_RUN; then
    echo "    [DRY-RUN] cleaned_xlsx=$cleaned_xlsx scope=$scope var=$var_col group=$group out=$out_dir"
    return 3   # 3 = dry-run
  fi

  mkdir -p "$out_dir" || { warn "    Could not create $out_dir"; return 1; }
  set +e
  MBX_NET_CLEANED_XLSX="$cleaned_xlsx" \
  MBX_NET_META_TXT="$METADATA_TXT" \
  MBX_NET_LEVEL_LET="$level_let" \
  MBX_NET_LEVEL_NAME="$level_name" \
  MBX_NET_SCOPE="$scope" \
  MBX_NET_VAR_COL="$var_col" \
  MBX_NET_GROUP="$group" \
  MBX_NET_OUT_DIR="$out_dir" \
  MBX_NET_RHO_THR="$RHO_THRESHOLD" \
  MBX_NET_Q_THR="$Q_THRESHOLD" \
  MBX_NET_PREV_THR="$PREV_THRESHOLD" \
    _R --vanilla "$R_SCRIPT_TEMPLATE" 2>&1 | tee "${out_dir}/_run.log"
  local rc=${PIPESTATUS[0]}
  set -e
  return $rc
}

# Aggregate summary lines for later cross-network reporting
SUMMARY_LINES_FILE="${WORK_DIR}/_summary_lines.tsv"
[[ -d "$WORK_DIR" ]] && : > "$SUMMARY_LINES_FILE" || true

# ── 6a. Global networks (one per level, all samples) ─────────────────────────
sep
info "── Global networks (all samples) ──────────────────────────────────────"
for i in "${!LEVELS[@]}"; do
  LVL="${LEVELS[$i]}"
  LVL_NAME="$(_level_full "$LVL")"
  XLSX="${LEVEL_XLSX_PATHS[$i]}"

  info "  Global :: level '${LVL}' (${LVL_NAME})"
  if [[ "$XLSX" == "SKIP" ]]; then
    warn "    Skipping (no cleaned xlsx for this level)"
    _GLOBAL_SKIP=$(( _GLOBAL_SKIP + 1 ))
    RESULT_LINES="${RESULT_LINES}- global x ${LVL_NAME} : SKIPPED (missing input)\n"
    continue
  fi

  OUT_RUN_DIR="${GLOBAL_DIR}/network_${LVL_NAME}"
  timer_start
  set +e
  _run_one_network "$XLSX" "$LVL" "$LVL_NAME" "global" "" "" "$OUT_RUN_DIR"
  RC=$?
  set -e
  timer_end

  case $RC in
    0) ok "    ✔  global x ${LVL_NAME}"
       _GLOBAL_PASS=$(( _GLOBAL_PASS + 1 ))
       grep '^MBX_NET_SUMMARY' "${OUT_RUN_DIR}/_run.log" >> "$SUMMARY_LINES_FILE" 2>/dev/null || true
       RESULT_LINES="${RESULT_LINES}- global x ${LVL_NAME} : OK\n" ;;
    2) _GLOBAL_PASS=$(( _GLOBAL_PASS + 1 ))   # already done; count as success
       RESULT_LINES="${RESULT_LINES}- global x ${LVL_NAME} : ALREADY_DONE\n" ;;
    3) RESULT_LINES="${RESULT_LINES}- global x ${LVL_NAME} : DRY_RUN\n" ;;
    *) warn "    ✘  global x ${LVL_NAME}  -- see ${OUT_RUN_DIR}/_run.log"
       _GLOBAL_FAIL=$(( _GLOBAL_FAIL + 1 ))
       RESULT_LINES="${RESULT_LINES}- global x ${LVL_NAME} : FAILED\n" ;;
  esac
done

# ── 6b. Per-group networks (skipped if --skip-per-group) ────────────────────
if $SKIP_PER_GROUP; then
  sep
  info "Skipping per-group networks (--skip-per-group set)."
else
  for j in "${!CATEGORICAL_COLS[@]}"; do
    COL="${CATEGORICAL_COLS[$j]}"
    COL_DIR="$(_sanitize_dirname "$COL")"
    sep
    info "── Per-group networks for variable '$COL' ──"

    # Discover groups for this variable, with their N counts
    GROUPS_TSV="/tmp/mbx_net_groups_${_TMPID}_${j}.tsv"
    GROUPS_DETECT_R="/tmp/mbx_net_groups_${_TMPID}_${j}.R"
    cat > "$GROUPS_DETECT_R" << RGRP
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
                            ignore.case = TRUE)) meta <- meta[-1, , drop = FALSE]
v <- trimws(as.character(meta[["${COL}"]]))
v <- v[!is.na(v) & v != ""]
tbl <- sort(table(v), decreasing = TRUE)
for (lab in names(tbl)) cat(lab, "\t", as.integer(tbl[lab]), "\n", sep = "")
RGRP

    # NB: GROUPS (plural) is a built-in read-only bash array containing the
    # user's UNIX group IDs.  Assigning to it is silently absorbed and array
    # operations behave unpredictably.  Use GRP_LIST instead.
    if $DRY_RUN; then
      info "  [DRY-RUN] Would discover groups for '$COL'"
      GRP_LIST=("group_a" "group_b")
      GRP_NS=(7 6)
    else
      _R --vanilla "$GROUPS_DETECT_R" > "$GROUPS_TSV" 2>/dev/null \
        || { warn "  Group enumeration for '$COL' failed -- skipping this variable."; continue; }
      GRP_LIST=()
      GRP_NS=()
      while IFS=$'\t' read -r g n; do
        if [[ -z "$g" ]]; then continue; fi
        GRP_LIST+=("$g")
        GRP_NS+=("$n")
      done < "$GROUPS_TSV"
      if [[ ${#GRP_LIST[@]} -eq 0 ]]; then
        warn "  No groups detected for '$COL' -- skipping."
        rm -f "$GROUPS_DETECT_R" "$GROUPS_TSV"
        continue
      fi
      _NICE_LIST=""
      for k in "${!GRP_LIST[@]}"; do
        _NICE_LIST="${_NICE_LIST}${GRP_LIST[$k]}=${GRP_NS[$k]} "
      done
      info "  Groups detected (${#GRP_LIST[@]}): $_NICE_LIST"
    fi
    rm -f "$GROUPS_DETECT_R" "$GROUPS_TSV"

    for i in "${!LEVELS[@]}"; do
      LVL="${LEVELS[$i]}"
      LVL_NAME="$(_level_full "$LVL")"
      XLSX="${LEVEL_XLSX_PATHS[$i]}"

      if [[ "$XLSX" == "SKIP" ]]; then
        _GROUP_SKIP=$(( _GROUP_SKIP + 1 ))
        RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} : ALL_GROUPS_SKIPPED (missing input)\n"
        continue
      fi

      LEVEL_DIR="${NET_DIR}/${COL_DIR}/network_${LVL_NAME}_by_${COL_DIR}"
      PER_GROUP_DIR="${LEVEL_DIR}/per_group"
      $DRY_RUN || mkdir -p "$PER_GROUP_DIR"

      info "  ${COL} :: level '${LVL}' (${LVL_NAME})"
      for k in "${!GRP_LIST[@]}"; do
        GRP="${GRP_LIST[$k]}"
        N_GRP="${GRP_NS[$k]:-0}"

        if [[ "$N_GRP" -lt "$MIN_GROUP_N" ]]; then
          warn "    Group '${GRP}' has only ${N_GRP} samples (< --min-group-n=${MIN_GROUP_N}) -- skipping."
          _GROUP_SKIP=$(( _GROUP_SKIP + 1 ))
          RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} x ${GRP} : SKIPPED (n=${N_GRP} < ${MIN_GROUP_N})\n"
          continue
        fi

        GRP_DIR="${PER_GROUP_DIR}/$(_sanitize_dirname "$GRP")"
        info "    group '${GRP}' (n=${N_GRP})"
        timer_start
        set +e
        _run_one_network "$XLSX" "$LVL" "$LVL_NAME" "per_group" "$COL" "$GRP" "$GRP_DIR"
        RC=$?
        set -e
        timer_end

        case $RC in
          0) ok "      ✔  ${COL} x ${LVL_NAME} x ${GRP}"
             _GROUP_PASS=$(( _GROUP_PASS + 1 ))
             grep '^MBX_NET_SUMMARY' "${GRP_DIR}/_run.log" >> "$SUMMARY_LINES_FILE" 2>/dev/null || true
             RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} x ${GRP} : OK\n" ;;
          2) _GROUP_PASS=$(( _GROUP_PASS + 1 ))
             RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} x ${GRP} : ALREADY_DONE\n" ;;
          3) RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} x ${GRP} : DRY_RUN\n" ;;
          *) warn "      ✘  ${COL} x ${LVL_NAME} x ${GRP}  -- see ${GRP_DIR}/_run.log"
             _GROUP_FAIL=$(( _GROUP_FAIL + 1 ))
             RESULT_LINES="${RESULT_LINES}- ${COL} x ${LVL_NAME} x ${GRP} : FAILED\n" ;;
        esac
      done
    done
  done
fi

rm -f "$R_SCRIPT_TEMPLATE"

# ─────────────────────────────────────────────────────────────────────────────
step "7/7 — Build per-variable group_comparison.xlsx + multi-panel plots"
# ─────────────────────────────────────────────────────────────────────────────

if $DRY_RUN; then
  warn "[DRY-RUN] Would aggregate group_comparison.xlsx + multi-panel plots."
elif $SKIP_PER_GROUP; then
  info "Skipping (--skip-per-group set)."
else
  # Rebuild the summary lines file from EVERY _run.log on disk so that
  # idempotent re-runs (where most networks are skipped) still see the
  # complete picture rather than only freshly-computed networks.
  : > "$SUMMARY_LINES_FILE"
  while IFS= read -r _logf; do
    grep '^MBX_NET_SUMMARY' "$_logf" >> "$SUMMARY_LINES_FILE" 2>/dev/null || true
  done < <(find "$NET_DIR" -mindepth 2 -name "_run.log" 2>/dev/null)

  if [[ ! -s "$SUMMARY_LINES_FILE" ]]; then
    warn "No network _run.log files found -- nothing to compare."
  else
    info "Aggregating from $(wc -l < "$SUMMARY_LINES_FILE" | tr -d ' ') summary lines."
  fi

if [[ -s "$SUMMARY_LINES_FILE" ]]; then
  AGG_R="/tmp/mbx_net_agg_${_TMPID}.R"
  cat > "$AGG_R" << RAGG
suppressPackageStartupMessages({
  library(openxlsx)
  library(igraph)
  library(ggplot2)
  library(RColorBrewer)
})
NET_DIR    <- "${NET_DIR}"
LINES_FILE <- "${SUMMARY_LINES_FILE}"
sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^[._-]+|[._-]+\$", "", x)
  x
}

lines <- readLines(LINES_FILE)
lines <- lines[startsWith(lines, "MBX_NET_SUMMARY")]
if (length(lines) == 0) { cat("[INFO]  No summary lines.\n"); quit(status = 0) }

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
for (k in c("n_samples","n_taxa","n_nodes","n_edges","n_pos","n_neg",
            "density","modularity","n_modules")) {
  if (k %in% names(df)) df[[k]] <- suppressWarnings(as.numeric(df[[k]]))
}

# Parse the context column to extract variable / group / level
df\$variable <- NA_character_; df\$group <- NA_character_; df\$level <- NA_character_
for (i in seq_len(nrow(df))) {
  ctx <- df\$context[i]
  if (grepl("^Global", ctx)) {
    df\$variable[i] <- "GLOBAL"
    df\$group[i]    <- "all_samples"
    df\$level[i]    <- sub("^Global \\\\(all samples\\\\) -- ", "", ctx)
  } else {
    parts <- strsplit(ctx, " -- ", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      df\$level[i] <- parts[2]
      vg <- strsplit(parts[1], " = ", fixed = TRUE)[[1]]
      if (length(vg) == 2) {
        df\$variable[i] <- vg[1]
        df\$group[i]    <- vg[2]
      }
    }
  }
}

# For each (variable, level) group: write group_comparison.xlsx + multi-panel
# plot.  Keep NO_EDGES rows -- a "no significant co-occurrence" finding is
# itself a meaningful biological result and should appear in the comparison
# (with 0 nodes / 0 edges) so the user can see WHY some groups have no
# network rather than seeing them silently disappear.
per_var <- df[df\$variable != "GLOBAL", , drop = FALSE]
per_var <- per_var[per_var\$status %in% c("OK", "NO_EDGES", "TOO_SMALL",
                                          "SKIPPED_FEW_TAXA"), , drop = FALSE]
if (nrow(per_var) == 0) {
  cat("[INFO]  No per-group networks to aggregate.\n"); quit(status = 0)
}
# Replace NA numeric fields (NO_EDGES rows) with 0 for clean xlsx display
for (k in c("n_nodes", "n_edges", "n_pos", "n_neg",
            "density", "modularity", "n_modules")) {
  if (k %in% names(per_var))
    per_var[[k]][is.na(per_var[[k]])] <- 0
}
if (!("hubs" %in% names(per_var))) per_var\$hubs <- NA_character_
per_var\$hubs[is.na(per_var\$hubs)] <- ""

vars <- unique(per_var\$variable)
for (v in vars) {
  for (lvl in unique(per_var\$level[per_var\$variable == v])) {
    sub <- per_var[per_var\$variable == v & per_var\$level == lvl, , drop = FALSE]
    if (nrow(sub) < 1) next

    sv <- sanitize(v); slvl <- sanitize(lvl)
    level_dir <- file.path(NET_DIR, sv, sprintf("network_%s_by_%s", slvl, sv))
    if (!dir.exists(level_dir)) next

    # ---- group_comparison.xlsx -----------------------------------------------
    cmp <- data.frame(
      group         = sub\$group,
      status        = sub\$status,
      n_samples     = as.integer(sub\$n_samples),
      n_taxa_filt   = as.integer(sub\$n_taxa),
      n_nodes       = as.integer(sub\$n_nodes),
      n_edges       = as.integer(sub\$n_edges),
      n_positive    = as.integer(sub\$n_pos),
      n_negative    = as.integer(sub\$n_neg),
      density       = round(sub\$density,    4),
      modularity    = round(sub\$modularity, 4),
      n_modules     = as.integer(sub\$n_modules),
      hub_taxa_top5 = sub\$hubs,
      stringsAsFactors = FALSE
    )

    # Hub-overlap matrix (Jaccard between groups)
    hub_lists <- strsplit(sub\$hubs, "|", fixed = TRUE)
    names(hub_lists) <- sub\$group
    G <- length(hub_lists)
    if (G >= 2) {
      jacc <- matrix(NA_real_, G, G, dimnames = list(sub\$group, sub\$group))
      for (a in seq_len(G)) for (b in seq_len(G)) {
        A <- hub_lists[[a]]; B <- hub_lists[[b]]
        u <- union(A, B); jacc[a, b] <- if (length(u) > 0) length(intersect(A, B)) / length(u) else NA_real_
      }
      jacc_df <- as.data.frame(round(jacc, 3), check.names = FALSE)
    } else {
      jacc_df <- data.frame()
    }

    wb <- createWorkbook()
    addWorksheet(wb, "metrics_per_group");  writeData(wb, "metrics_per_group", cmp)
    if (nrow(jacc_df) > 0) {
      addWorksheet(wb, "hub_overlap_jaccard")
      writeData(wb, "hub_overlap_jaccard", jacc_df, rowNames = TRUE)
    }
    saveWorkbook(wb, file.path(level_dir, "group_comparison.xlsx"), overwrite = TRUE)

    # ---- multi-panel plot ----------------------------------------------------
    # Bar chart: density + modularity per group (side-by-side)
    long_df <- data.frame(
      group   = rep(cmp\$group, 2),
      metric  = factor(rep(c("density", "modularity"), each = nrow(cmp)),
                       levels = c("density", "modularity")),
      value   = c(cmp\$density, cmp\$modularity),
      stringsAsFactors = FALSE
    )
    p_bar <- ggplot(long_df, aes(x = group, y = value, fill = metric)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = sprintf("%.3f", value)),
                position = position_dodge(width = 0.8), vjust = -0.4,
                size = 3.2) +
      scale_fill_manual(values = c(density = "#3b73af", modularity = "#7f8c8d")) +
      labs(title    = sprintf("Network metrics per group  --  %s by %s",
                              lvl, v),
           subtitle = "Higher density = more connections.  Higher modularity = more compartmentalised.",
           x = NULL, y = NULL, fill = NULL) +
      theme_bw(base_size = 12) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(color = "grey40"),
            axis.text.x   = element_text(angle = 25, hjust = 1),
            legend.position = "top")
    ggsave(file.path(level_dir, "network_metrics_comparison.png"),
           p_bar, width = 9, height = 5.5, dpi = 300)
    ggsave(file.path(level_dir, "network_metrics_comparison.pdf"),
           p_bar, width = 9, height = 5.5)

    # ---- multi-panel network plot (one igraph per group) -------------------
    grps  <- cmp\$group
    n_pan <- length(grps)
    n_col <- min(3, n_pan)
    n_row <- ceiling(n_pan / n_col)
    panel_w <- 5.5; panel_h <- 5
    png(file.path(level_dir, "multi_panel_plot.png"),
        width = panel_w * n_col, height = panel_h * n_row,
        units = "in", res = 300)
    par(mfrow = c(n_row, n_col), mar = c(1, 1, 3, 1))
    for (g in grps) {
      gd <- file.path(level_dir, "per_group", sanitize(g))
      gml <- file.path(gd, "network.graphml")
      if (file.exists(gml)) {
        gi <- read_graph(gml, format = "graphml")
        if (vcount(gi) >= 2 && ecount(gi) >= 1) {
          set.seed(42)
          lay <- layout_with_fr(gi, weights = abs(E(gi)\$rho))
          ecol <- ifelse(E(gi)\$sign == "positive", "#3b73af", "#c0392b")
          plot(gi, layout = lay,
               vertex.label = NA, vertex.size = 4,
               vertex.color = "#7f8c8d", vertex.frame.color = "grey30",
               edge.color = ecol, edge.width = 0.6,
               main = sprintf("%s  (n=%d, edges=%d, mod=%.2f)",
                              g,
                              cmp\$n_samples[cmp\$group == g],
                              cmp\$n_edges[cmp\$group == g],
                              cmp\$modularity[cmp\$group == g]))
        } else {
          plot.new(); title(main = sprintf("%s\n(empty)", g))
        }
      } else {
        plot.new(); title(main = sprintf("%s\n(no graph)", g))
      }
    }
    mtext(sprintf("%s by %s", lvl, v), outer = TRUE, line = -1.2, font = 2)
    dev.off()
    pdf(file.path(level_dir, "multi_panel_plot.pdf"),
        width = panel_w * n_col, height = panel_h * n_row)
    par(mfrow = c(n_row, n_col), mar = c(1, 1, 3, 1))
    for (g in grps) {
      gd <- file.path(level_dir, "per_group", sanitize(g))
      gml <- file.path(gd, "network.graphml")
      if (file.exists(gml)) {
        gi <- read_graph(gml, format = "graphml")
        if (vcount(gi) >= 2 && ecount(gi) >= 1) {
          set.seed(42)
          lay <- layout_with_fr(gi, weights = abs(E(gi)\$rho))
          ecol <- ifelse(E(gi)\$sign == "positive", "#3b73af", "#c0392b")
          plot(gi, layout = lay,
               vertex.label = NA, vertex.size = 4,
               vertex.color = "#7f8c8d", vertex.frame.color = "grey30",
               edge.color = ecol, edge.width = 0.6,
               main = sprintf("%s  (n=%d, edges=%d, mod=%.2f)",
                              g, cmp\$n_samples[cmp\$group == g],
                              cmp\$n_edges[cmp\$group == g],
                              cmp\$modularity[cmp\$group == g]))
        } else { plot.new(); title(main = sprintf("%s\n(empty)", g)) }
      } else { plot.new(); title(main = sprintf("%s\n(no graph)", g)) }
    }
    mtext(sprintf("%s by %s", lvl, v), outer = TRUE, line = -1.2, font = 2)
    dev.off()

    cat(sprintf("[OK]    %s by %s: comparison files written.\n", lvl, v))
  }
}
RAGG

  _R --vanilla "$AGG_R" || warn "Per-variable comparison aggregation failed (non-fatal)."
  rm -f "$AGG_R"
  fi   # end of: if [[ -s "$SUMMARY_LINES_FILE" ]]
fi     # end of: if/elif/else for skip-per-group / dry-run

# ─────────────────────────────────────────────────────────────────────────────
# Write provenance file + final summary
# ─────────────────────────────────────────────────────────────────────────────
INFO_TXT="${NET_DIR}/mbx_networks_info.txt"

if $DRY_RUN; then
  warn "[DRY-RUN] Would write provenance file: $INFO_TXT"
else
cat > "$INFO_TXT" << INFO
# ============================================================================
# mbx_networks_info.txt
# Generated by mbx_network_run.sh   (step 17)
# Date : $NOW
# ============================================================================

# ── Inputs ────────────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
METADATA_TXT=$METADATA_TXT
EZCLEAN_INFO=$EZCLEAN_INFO
R_VERSION=$R_VERSION

# ── Algorithm + thresholds ────────────────────────────────────────────────────
METHOD=CLR_transform_+_Spearman_+_BH_FDR
RHO_THRESHOLD=$RHO_THRESHOLD
Q_THRESHOLD=$Q_THRESHOLD
PREVALENCE_THRESHOLD=$PREV_THRESHOLD
MIN_GROUP_N=$MIN_GROUP_N
MODULE_DETECTION=Louvain

# ── Output directory ──────────────────────────────────────────────────────────
NET_DIR=$NET_DIR
GLOBAL_DIR=$GLOBAL_DIR
WORK_DIR=$WORK_DIR

# ── Levels analysed ───────────────────────────────────────────────────────────
LEVELS_REQUESTED=${LEVELS[*]}

# ── Variables analysed ────────────────────────────────────────────────────────
INFO
for col in "${CATEGORICAL_COLS[@]}"; do
  echo "VARIABLE=$col" >> "$INFO_TXT"
done
cat >> "$INFO_TXT" << INFO2

# ── Run summary ───────────────────────────────────────────────────────────────
GLOBAL_PASS=$_GLOBAL_PASS
GLOBAL_FAIL=$_GLOBAL_FAIL
GLOBAL_SKIP=$_GLOBAL_SKIP
PERGROUP_PASS=$_GROUP_PASS
PERGROUP_FAIL=$_GROUP_FAIL
PERGROUP_SKIP=$_GROUP_SKIP

# ── Output files (per network) ────────────────────────────────────────────────
# edges.tsv               (taxon_1, taxon_2, rho, p_raw, q_BH, sign)
# nodes.tsv               (taxon, degree, betweenness, closeness, hub_score,
#                          mean_abundance, prevalence, module)
# network.graphml         (Cytoscape / Gephi import)
# network_plot.{png,pdf}  (FR layout, edges by sign, nodes by module)
# network_summary.xlsx    (n_nodes, n_edges, density, modularity, etc.)
# hub_taxa.xlsx           (top hub taxa ranked)
# modules.tsv             (taxon -> module_id)

# ── Per-variable comparison files ─────────────────────────────────────────────
# group_comparison.xlsx          (metrics_per_group + hub_overlap_jaccard)
# multi_panel_plot.{png,pdf}     (one igraph panel per group)
# network_metrics_comparison.{png,pdf}  (density + modularity bar chart)
INFO2

  ok "Provenance written -> $INFO_TXT"
fi

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete -- no R code was executed."
else
  echo ""
  echo "  ╔══════════════════════════════════════════════════════════════╗"
  echo "  ║  Co-occurrence network analysis complete                     ║"
  echo "  ╚══════════════════════════════════════════════════════════════╝"
  echo ""
  echo "  Output directory : $NET_DIR"
  echo ""
  echo "  Global networks  :  ${_GLOBAL_PASS} passed, ${_GLOBAL_FAIL} failed, ${_GLOBAL_SKIP} skipped"
  if ! $SKIP_PER_GROUP; then
    echo "  Per-group nets   :  ${_GROUP_PASS} passed, ${_GROUP_FAIL} failed, ${_GROUP_SKIP} skipped"
  fi
  echo ""
  if [[ -n "$RESULT_LINES" ]]; then
    echo "  Per-network outcome:"
    printf '%b' "$RESULT_LINES" | sed 's/^/    /'
  fi
  echo ""
  echo "  Inspect first:"
  echo "    ${GLOBAL_DIR}/network_<level>/network_plot.png   (the ecological backbone)"
  if ! $SKIP_PER_GROUP; then
    for col in "${CATEGORICAL_COLS[@]}"; do
      SC="$(_sanitize_dirname "$col")"
      echo "    ${NET_DIR}/${SC}/network_<level>_by_${SC}/group_comparison.xlsx"
      echo "    ${NET_DIR}/${SC}/network_<level>_by_${SC}/multi_panel_plot.png"
    done
  fi
  echo ""
  echo "  Provenance file : $INFO_TXT"
  echo ""
fi
sep
