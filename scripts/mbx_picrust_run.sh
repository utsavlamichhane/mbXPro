#!/usr/bin/env bash
# =============================================================================
#  mbx_picrust_run.sh   —   Step 15 of the mbX Pro pipeline
# =============================================================================
#  PURPOSE
#  -------
#  Predict metagenomic content from 16S rRNA marker gene sequences using
#  PICRUSt2 (Douglas et al. 2020).  Generates KO, EC, COG, and MetaCyc pathway
#  abundance tables, performs differential abundance analysis (Kruskal-Wallis
#  + Dunn) per categorical variable, generates publication-ready plots
#  (stacked bars, lfc heatmaps, functional PCoA, genera × pathway correlation),
#  and writes a self-contained HTML report.
#
#  ╔═══════════════════════════════════════════════════════════════════════╗
#  ║  CRITICAL: This script does NOT modify your QIIME2 conda environment.  ║
#  ║  It uses a SEPARATE conda env named "mbx_picrust2" (auto-created on    ║
#  ║  first run) and invokes picrust2 commands via `conda run` (no activate).║
#  ╚═══════════════════════════════════════════════════════════════════════╝
#
#  Reads from previous steps:
#      • 4_dada2_outputs/representative_sequences.qza       (ASV rep seqs)
#      • 7_taxonomy_csv/feature_table_filtered.qza          (filtered counts)
#      • 13_beta_diversity_results/mbx_beta_diversity_info.txt (metadata path)
#      • 14_differential_abundance_ANCOMBC2/ (optional, for genus×pathway corr)
#
#  WHAT THIS SCRIPT PRODUCES
#  -------------------------
#  <mbX_pro_outputs_*>/
#  └── 15_picrust2/
#      ├── all_picrust2_outputs/
#      │   ├── KO_metagenome/           pred_metagenome_unstrat[_descrip].tsv.gz
#      │   ├── EC_metagenome/           pred_metagenome_unstrat[_descrip].tsv.gz
#      │   ├── COG_metagenome/          pred_metagenome_unstrat[_descrip].tsv.gz
#      │   ├── pathways_metacyc/        path_abun_unstrat[_descrip].tsv.gz
#      │   ├── nsti/                    marker NSTI tables, weighted NSTI,
#      │   │                            nsti_filtering_summary.txt
#      │   └── place_seqs/              placed_seqs.tre, place_seqs.log
#      ├── picrust2_<Variable>/         (one per categorical variable)
#      │   ├── DA_KO_KW_<Var>.xlsx                (+ pairwise xlsx)
#      │   ├── DA_EC_KW_<Var>.xlsx                (+ pairwise xlsx)
#      │   ├── DA_COG_KW_<Var>.xlsx               (+ pairwise xlsx)
#      │   ├── DA_metacyc_pathways_KW_<Var>.xlsx  (+ pairwise xlsx)
#      │   ├── stacked_bar_top20_metacyc_<Var>.{png,pdf}
#      │   ├── heatmap_DA_pathways_<Var>.{png,pdf}    (lfc + q-stars)
#      │   ├── PCoA_BrayCurtis_FUNCTIONAL_<Var>.{png,pdf}
#      │   └── correlation_DAgenera_x_DApathways_<Var>.{tsv,png,pdf}
#      ├── working_dir_picrust2/
#      │   ├── inputs/                  rep_seqs.fasta, feature_table.{biom,tsv}
#      │   ├── intermediates/           HSP outputs, pathway intermediates
#      │   ├── R_DA_results/            RDS objects per (var × database)
#      │   └── picrust2_run_logs/
#      ├── picrust2_report.html         (self-contained, embedded plots)
#      ├── Summary_picrust2_NSTI.xlsx   (per-sample NSTI + sample flags)
#      └── mbx_picrust2_info.txt
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

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_picrust_run.sh — Predict and analyse functional profiles from 16S data
                     using PICRUSt2 (Douglas et al. 2020).  Generates KO, EC,
                     COG, MetaCyc pathway abundance tables; differential
                     abundance per categorical variable; publication-ready
                     plots; and a self-contained HTML report.

USAGE:
  mbx_picrust_run.sh <mbX_pro_outputs_dir> [OPTIONS]

OPTIONS:
  --nsti N               NSTI threshold (default: 2.0).
                         ASVs with NSTI > N are dropped from metagenome
                         predictions.  Default 2.0 is appropriate for
                         standard gut/soil 16S data.  Rumen, termite gut,
                         and other unusual environments may need adjustment.
                         (Future: this will be read from <mbX_outputs>/config.yaml.)

  --variables LIST       Comma-separated categorical variables (default: all
                         auto-detected).
  --picrust-env NAME     Conda env name for PICRUSt2 (default: mbx_picrust2).
                         Searched first in standard conda/mamba env locations;
                         if missing AND no other env named 'picrust2' is found,
                         a new env is auto-created.

  --picrust-env-path P   Explicit absolute path to an existing PICRUSt2 env's
                         root (the dir containing bin/place_seqs.py).
                         Overrides --picrust-env.

  --skip-deps-check      Skip conda env creation/check (assume already set up)
  --skip-plots           Don't generate volcano/heatmap/PCoA PNGs (tables only)
  --skip-correlation     Skip the genus × pathway Spearman correlation
                         (default: try if step 14 ANCOMBC2 outputs exist;
                         skip silently otherwise)
  --skip-html            Don't generate HTML report
  --force                Run even if step 11 OVERALL_STATUS != PASS
  --force-rerun          Recompute everything (ignore existing outputs)
  --threads N            Override CPU thread count (default: auto-detect)
  --dry-run              Print commands; do not execute
  -h, --help             Show this help

EXAMPLES:
  mbx_picrust_run.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_picrust_run.sh /path/to/mbX_pro_outputs_20260417_121431 --nsti 1.5
  mbx_picrust_run.sh /path/to/mbX_pro_outputs_20260417_121431 --variables Treatment

PROVENANCE:
  All paths come from previous step info files.
  Re-runs are safe — completed steps are skipped.

DEPENDENCIES:
  System  : QIIME2 conda env active
            + conda  (creates a separate 'mbx_picrust2' env on first run;
                      ~10-15 min one-time install of picrust2 + its deps)
  R       : system R (or conda env's R) with:
            openxlsx, ggplot2, dplyr, tidyr, dunn.test, vegan, pheatmap,
            RColorBrewer, patchwork
            (auto-installs missing CRAN packages.)

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
NSTI_THRESHOLD="2.0"
USER_VARIABLES=""
PICRUST_ENV="mbx_picrust2"
PICRUST_ENV_PATH=""
SKIP_DEPS_CHECK=false
SKIP_PLOTS=false
SKIP_CORRELATION_USER=""    # auto-decide unless user sets --skip-correlation
SKIP_HTML=false
FORCE_GATE=false
FORCE_RERUN=false
USER_THREADS=""
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)         usage ;;
    --nsti)            NSTI_THRESHOLD="${2:-}"; [[ -z "$NSTI_THRESHOLD" ]] && err "--nsti requires a value"; shift 2 ;;
    --variables)       USER_VARIABLES="${2:-}"; [[ -z "$USER_VARIABLES" ]] && err "--variables requires a value"; shift 2 ;;
    --picrust-env)     PICRUST_ENV="${2:-}"; [[ -z "$PICRUST_ENV" ]] && err "--picrust-env requires a value"; shift 2 ;;
    --picrust-env-path) PICRUST_ENV_PATH="${2:-}"; [[ -z "$PICRUST_ENV_PATH" ]] && err "--picrust-env-path requires a value"; shift 2 ;;
    --skip-deps-check) SKIP_DEPS_CHECK=true; shift ;;
    --skip-plots)      SKIP_PLOTS=true; shift ;;
    --skip-correlation) SKIP_CORRELATION_USER=true; shift ;;
    --skip-html)       SKIP_HTML=true; shift ;;
    --force)           FORCE_GATE=true; shift ;;
    --force-rerun)     FORCE_RERUN=true; shift ;;
    --threads)         USER_THREADS="${2:-}"; case "$USER_THREADS" in ''|*[!0-9]*) err "--threads must be a positive integer" ;; esac; shift 2 ;;
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

# Validate NSTI threshold (must be positive number)
if ! awk -v n="$NSTI_THRESHOLD" 'BEGIN{exit !(n+0 == n && n+0 > 0)}'; then
  err "--nsti must be a positive number (got: '$NSTI_THRESHOLD')"
fi

# ── CPU detection ─────────────────────────────────────────────────────────────
if [[ -n "$USER_THREADS" ]]; then
  N_JOBS="$USER_THREADS"
elif command -v nproc   &>/dev/null; then N_JOBS="$(nproc)"
elif command -v sysctl  &>/dev/null; then N_JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)"
else                                       N_JOBS=1
fi

# ── Detect platform (Apple Silicon / Intel / Linux / Windows) ────────────────
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

# ── R env wrapper (strip conda's R_LIBS_USER from system R) ──────────────────
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}

# Locate Rscript (this script needs CRAN-only packages, so simple is fine)
RSCRIPT_CMD=""
for _c in \
    "/opt/homebrew/bin/Rscript" \
    "/usr/local/bin/Rscript" \
    "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
    "/usr/bin/Rscript" \
    "$(command -v Rscript 2>/dev/null || true)"; do
  [[ -n "$_c" && -x "$_c" ]] && { RSCRIPT_CMD="$_c"; break; }
done
[[ -n "$RSCRIPT_CMD" ]] || err "No Rscript found.

  Install R for $PLATFORM_LABEL:
  • macOS:    brew install r
  • Linux:    sudo apt install r-base   (Debian/Ubuntu)
  • Windows:  https://cran.r-project.org/

  Then re-run."

_R() { _strip_env "$RSCRIPT_CMD" "$@"; }
R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"

# ── PICRUSt2 wrapper ─────────────────────────────────────────────────────────
# Set later, after env detection.  We invoke picrust2 commands by prepending
# the env's bin/ to PATH for that single command — no `conda activate`,
# no `conda run` (its --no-capture-output is awkward and mamba's run rejects
# the same flags).  picrust2's shebangs use absolute env-python paths, so
# PATH-prefixing is bullet-proof.
PICRUST_BIN_DIR=""
_PICRUST() {
  PATH="${PICRUST_BIN_DIR}:${PATH}" "$@"
}

# ── PID-based temp file IDs ──────────────────────────────────────────────────
_TMPID="${$}_$(date +%s)"
trap 'rm -f /tmp/mbx_picrust_${_TMPID}*.R /tmp/mbx_picrust_${_TMPID}*.txt /tmp/mbx_picrust_${_TMPID}*.html' EXIT

# =============================================================================
step "1/12 — Read previous step info files"
# =============================================================================
BETA_INFO="${MBX_OUT_DIR}/13_beta_diversity_results/mbx_beta_diversity_info.txt"
PRE_DIV_INFO="${MBX_OUT_DIR}/11_pre_diversity/mbx_pre_diversity_info.txt"
ANCO_INFO="${MBX_OUT_DIR}/14_differential_abundance_ANCOMBC2/mbx_ancombc2_info.txt"

[[ -f "$BETA_INFO" ]]    || err "Beta-diversity info file not found:
    $BETA_INFO
  → Run mbx_beta_diversity_run.sh first (step 13)."
[[ -f "$PRE_DIV_INFO" ]] || err "Pre-diversity info file not found:
    $PRE_DIV_INFO
  → Run mbx_pre_diversity_parameters.sh first (step 11)."

OVERALL_STATUS="$(_read_key OVERALL_STATUS      "$PRE_DIV_INFO")"
READY="$(        _read_key READY_FOR_DIVERSITY "$PRE_DIV_INFO")"
METADATA_TXT="$( _read_key METADATA_TXT        "$BETA_INFO")"

# Rep seqs from step 4 (DADA2)
REP_SEQS_QZA="${MBX_OUT_DIR}/4_dada2_outputs/representative_sequences.qza"
[[ -f "$REP_SEQS_QZA" ]] || err "Representative sequences not found:
    $REP_SEQS_QZA
  → Run mbx_dada2_run.sh first (step 4)."

# Filtered feature table from step 7 (taxonomy filtered, mito/chloro removed)
FT_FILTERED_QZA="${MBX_OUT_DIR}/7_taxonomy_csv/feature_table_filtered.qza"
[[ -f "$FT_FILTERED_QZA" ]] || err "Filtered feature table not found:
    $FT_FILTERED_QZA
  → Run mbx_taxonomy_run.sh first (step 7)."

[[ -n "$METADATA_TXT" && -f "$METADATA_TXT" ]] || err "METADATA_TXT not usable: '$METADATA_TXT'"

# Optional: ANCOMBC2 outputs for genus × pathway correlation
HAS_ANCOMBC=false
ANCOMBC_DIR=""
if [[ -f "$ANCO_INFO" ]]; then
  ANCOMBC_DIR="$(_read_key ANCOMBC2_OUT_DIR "$ANCO_INFO")"
  if [[ -n "$ANCOMBC_DIR" && -d "$ANCOMBC_DIR" ]]; then
    HAS_ANCOMBC=true
  fi
fi

# Decide whether to attempt correlation
if [[ -n "$SKIP_CORRELATION_USER" ]]; then
  DO_CORRELATION=false
elif $HAS_ANCOMBC; then
  DO_CORRELATION=true
else
  DO_CORRELATION=false
fi

# ── Gating ────────────────────────────────────────────────────────────────────
case "$OVERALL_STATUS" in
  PASS|PASS_WITH_WARNINGS) : ;;
  *)
    if $FORCE_GATE; then
      warn "Step 11 OVERALL_STATUS = '$OVERALL_STATUS' — proceeding because --force was given."
    else
      err "Step 11 OVERALL_STATUS = '$OVERALL_STATUS'.
  → Functional prediction on a non-PASS state may produce misleading results.
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
info "QIIME2 env (current) : ${CONDA_DEFAULT_ENV:-unknown}"
info "PICRUSt2 env (target): $PICRUST_ENV"
info "Rscript              : $RSCRIPT_CMD  (R $R_VERSION)"
info "Output root          : $MBX_OUT_DIR"
info "Metadata             : $METADATA_TXT"
info "Rep seqs (input)     : $REP_SEQS_QZA"
info "Filtered table       : $FT_FILTERED_QZA"
info "NSTI threshold       : $NSTI_THRESHOLD"
info "CPU threads          : $N_JOBS"
info "ANCOMBC2 results     : $($HAS_ANCOMBC && echo "found at $ANCOMBC_DIR" || echo "not found — will skip genus × pathway correlation")"
info "Correlation analysis : $($DO_CORRELATION && echo enabled || echo disabled)"
info "Step-11 status       : $OVERALL_STATUS  /  ready=$READY"
sep

# ── Set up output directories ─────────────────────────────────────────────────
PICRUST_DIR="${MBX_OUT_DIR}/15_picrust2"
ALL_DIR="${PICRUST_DIR}/all_picrust2_outputs"
WORK_DIR="${PICRUST_DIR}/working_dir_picrust2"
INPUTS_DIR="${WORK_DIR}/inputs"
INTER_DIR="${WORK_DIR}/intermediates"
RDS_DIR="${WORK_DIR}/R_DA_results"
LOG_DIR="${WORK_DIR}/picrust2_run_logs"
KO_OUT_DIR="${ALL_DIR}/KO_metagenome"
EC_OUT_DIR="${ALL_DIR}/EC_metagenome"
COG_OUT_DIR="${ALL_DIR}/COG_metagenome"
PATH_OUT_DIR="${ALL_DIR}/pathways_metacyc"
PLACE_DIR="${ALL_DIR}/place_seqs"
NSTI_DIR="${ALL_DIR}/nsti"
INFO_FILE="${PICRUST_DIR}/mbx_picrust2_info.txt"
NSTI_XLSX="${PICRUST_DIR}/Summary_picrust2_NSTI.xlsx"
HTML_REPORT="${PICRUST_DIR}/picrust2_report.html"

mkdir -p "$PICRUST_DIR" "$ALL_DIR" "$WORK_DIR" "$INPUTS_DIR" "$INTER_DIR" \
         "$RDS_DIR" "$LOG_DIR" "$KO_OUT_DIR" "$EC_OUT_DIR" "$COG_OUT_DIR" \
         "$PATH_OUT_DIR" "$PLACE_DIR" "$NSTI_DIR"

# =============================================================================
step "2/12 — Detect / install PICRUSt2 conda environment"
# =============================================================================
#
# Strategy:
#   1. If --picrust-env-path was given, validate it directly.
#   2. Otherwise search every conda/mamba install on this machine for an env
#      whose bin/ contains place_seqs.py.  Try in order:
#         a. env named exactly $PICRUST_ENV      (default 'mbx_picrust2')
#         b. env named exactly 'picrust2'         (common pre-existing env)
#         c. any env whose name ends in 'picrust2' / 'picrust'
#   3. If none found, install via mamba (fastest, handles Apple Silicon
#      compiler post-link issues better than conda) — falling back to conda.
#
# We never `conda activate` — we run picrust2 commands by prepending the
# env's bin/ to PATH for the duration of one call (see _PICRUST wrapper).
# =============================================================================

# All conda/mamba env-search roots on this machine (order matters: prefer
# user's primary conda first)
_CANDIDATE_ENV_ROOTS=()
for _root in \
    "$HOME/miniconda3/envs" \
    "$HOME/miniforge3/envs" \
    "$HOME/anaconda3/envs" \
    "$HOME/.local/share/mamba/envs" \
    "/opt/miniconda3/envs" \
    "/opt/anaconda3/envs" \
    "/opt/homebrew/Caskroom/miniconda/base/envs"; do
  [[ -d "$_root" ]] && _CANDIDATE_ENV_ROOTS+=( "$_root" )
done
if command -v conda >/dev/null 2>&1; then
  while IFS= read -r _l; do
    case "$_l" in '' | \#* ) continue ;; esac
    # `conda env list` columns: name [optional *] path
    _p="$(echo "$_l" | awk '{print $NF}')"
    [[ -d "$_p" ]] && [[ -d "$_p/.." ]] && _CANDIDATE_ENV_ROOTS+=( "$(cd "$_p/.." && pwd)" )
  done < <(conda env list 2>/dev/null)
fi
# Dedup
_CANDIDATE_ENV_ROOTS=($(printf '%s\n' "${_CANDIDATE_ENV_ROOTS[@]}" | awk '!seen[$0]++'))

# Probe a candidate env path for working picrust2
_validate_picrust_env() {
  local env_root="$1"
  [[ -d "$env_root" ]] || return 1
  [[ -x "$env_root/bin/place_seqs.py" ]] || return 1
  [[ -x "$env_root/bin/hsp.py" ]]        || return 1
  [[ -x "$env_root/bin/metagenome_pipeline.py" ]] || return 1
  [[ -x "$env_root/bin/pathway_pipeline.py" ]]    || return 1
  return 0
}

# Search for usable env
_FOUND_ENV_PATH=""
_FOUND_ENV_NAME=""

if [[ -n "$PICRUST_ENV_PATH" ]]; then
  if _validate_picrust_env "$PICRUST_ENV_PATH"; then
    _FOUND_ENV_PATH="$PICRUST_ENV_PATH"
    _FOUND_ENV_NAME="(custom: $(basename "$PICRUST_ENV_PATH"))"
  else
    err "--picrust-env-path '$PICRUST_ENV_PATH' is not a valid PICRUSt2 env.
  Expected to find:
      $PICRUST_ENV_PATH/bin/place_seqs.py
      $PICRUST_ENV_PATH/bin/hsp.py
      $PICRUST_ENV_PATH/bin/metagenome_pipeline.py
      $PICRUST_ENV_PATH/bin/pathway_pipeline.py"
  fi
elif ! $SKIP_DEPS_CHECK; then
  # Search names in priority order
  for _name in "$PICRUST_ENV" "picrust2" "mbx_picrust2"; do
    for _root in "${_CANDIDATE_ENV_ROOTS[@]}"; do
      _p="$_root/$_name"
      if _validate_picrust_env "$_p"; then
        _FOUND_ENV_PATH="$_p"
        _FOUND_ENV_NAME="$_name"
        break 2
      fi
    done
  done

  # Last resort: search every env across every root for one with place_seqs.py
  if [[ -z "$_FOUND_ENV_PATH" ]]; then
    for _root in "${_CANDIDATE_ENV_ROOTS[@]}"; do
      for _envdir in "$_root"/*; do
        if _validate_picrust_env "$_envdir"; then
          _FOUND_ENV_PATH="$_envdir"
          _FOUND_ENV_NAME="$(basename "$_envdir") (auto-discovered)"
          break 2
        fi
      done
    done
  fi
fi

if [[ -n "$_FOUND_ENV_PATH" ]]; then
  PICRUST_BIN_DIR="$_FOUND_ENV_PATH/bin"
  ok "Found existing PICRUSt2 env: $_FOUND_ENV_NAME"
  info "  Path: $_FOUND_ENV_PATH"
  info "  → Will use this env (no install needed)."
elif $SKIP_DEPS_CHECK; then
  warn "Skipping conda env check (--skip-deps-check)."
  err "PICRUSt2 env still not found.  Either install it or pass --picrust-env-path."
else
  # ── No existing env found — install one ──────────────────────────────────
  info "No existing PICRUSt2 env found.  Will install one now."
  echo ""
  echo "  This will take ~10-15 minutes and downloads ~2 GB of files."
  echo "  PICRUSt2 + dependencies will be installed in a NEW env named '$PICRUST_ENV'"
  echo "  WITHOUT touching your current QIIME2 environment."
  echo ""

  if ! command -v conda >/dev/null 2>&1; then
    err "'conda' command not found.

  PICRUSt2 requires conda for isolated installation.
  Install Miniconda for $PLATFORM_LABEL:
  • macOS:   https://docs.conda.io/projects/miniconda/en/latest/
  • Linux:   https://docs.conda.io/projects/miniconda/en/latest/
  • Windows: https://docs.conda.io/projects/miniconda/en/latest/

  Then re-run this script."
  fi
  info "conda  : $(conda --version 2>&1 | head -1)"
  if command -v mamba >/dev/null 2>&1; then
    info "mamba  : $(mamba --version 2>&1 | head -1)  ← will use (faster solver)"
    INSTALLER="mamba"
  else
    info "mamba  : not installed"
    INSTALLER="conda"
  fi

  # Apple Silicon needs osx-64 emulation (no native arm64 builds for gappa/epa-ng)
  USE_ROSETTA=false
  if [[ "$_PLATFORM_OS" == "Darwin" && "$_PLATFORM_ARCH" == "arm64" ]]; then
    USE_ROSETTA=true
    warn "Apple Silicon detected — will install osx-64 packages via Rosetta 2"
    warn "  Reason: bioconda lacks native osx-arm64 builds for gappa/epa-ng."
    if ! /usr/bin/arch -x86_64 /usr/bin/true 2>/dev/null; then
      warn "Rosetta 2 not yet installed.  If install fails:
      softwareupdate --install-rosetta --agree-to-license"
    fi
  fi
  case "$_PLATFORM_OS" in
    Darwin)
      xcode-select -p >/dev/null 2>&1 || warn "Xcode CLT not detected — install if needed: xcode-select --install" ;;
    Linux)
      command -v gcc >/dev/null 2>&1 || warn "GCC not detected — install build-essential if compile errors arise" ;;
    MINGW*|MSYS*|CYGWIN*)
      warn "Windows: bioconda binaries may need WSL2.  If install fails, use WSL2." ;;
  esac

  PICRUST_VERSION_PIN="picrust2=2.6.3"

  if $DRY_RUN; then
    warn "[DRY-RUN] Would install via $INSTALLER (Rosetta=$USE_ROSETTA): $PICRUST_VERSION_PIN"
  else
    timer_start
    if $USE_ROSETTA; then
      CONDA_SUBDIR=osx-64 "$INSTALLER" create -n "$PICRUST_ENV" \
        -c bioconda -c conda-forge "$PICRUST_VERSION_PIN" -y \
        || err "Failed to create PICRUSt2 env via $INSTALLER (Rosetta 2 mode).

  If conda failed but mamba is available, install mamba and retry:
      conda install -n base -c conda-forge mamba

  Or install manually:
      CONDA_SUBDIR=osx-64 mamba create -n $PICRUST_ENV -c bioconda -c conda-forge $PICRUST_VERSION_PIN -y

  Then re-run this script."
    else
      "$INSTALLER" create -n "$PICRUST_ENV" \
        -c bioconda -c conda-forge "$PICRUST_VERSION_PIN" -y \
        || err "Failed to create PICRUSt2 env via $INSTALLER."
    fi
    timer_end
  fi

  # Re-discover the env we just created (mamba may put it in ~/.local/share/mamba/envs)
  for _root in "${_CANDIDATE_ENV_ROOTS[@]}" "$HOME/.local/share/mamba/envs"; do
    _p="$_root/$PICRUST_ENV"
    if _validate_picrust_env "$_p"; then
      PICRUST_BIN_DIR="$_p/bin"
      _FOUND_ENV_PATH="$_p"
      _FOUND_ENV_NAME="$PICRUST_ENV (newly installed)"
      break
    fi
  done
  if [[ -z "$PICRUST_BIN_DIR" ]]; then
    err "Newly installed env not found in any expected location.
  Searched: ${_CANDIDATE_ENV_ROOTS[*]}
  Run 'conda env list' or 'mamba env list' to locate it,
  then re-run with --picrust-env-path /full/path/to/env."
  fi
  ok "Created env at: $_FOUND_ENV_PATH"
fi

# Print PICRUSt2 version
if ! $DRY_RUN && [[ -n "$PICRUST_BIN_DIR" ]]; then
  PICRUST_VERSION_STR="$(_PICRUST place_seqs.py -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 || echo 'unknown')"
  info "PICRUSt2 version: $PICRUST_VERSION_STR"
fi
sep

# =============================================================================
step "3/12 — Check / install R packages"
# =============================================================================
R_PROBE="/tmp/mbx_picrust_${_TMPID}_probe.R"
cat > "$R_PROBE" << 'RPROBE'
required <- c("openxlsx","ggplot2","dplyr","tidyr","dunn.test",
              "vegan","pheatmap","RColorBrewer","patchwork")
miss <- required[!sapply(required, requireNamespace, quietly=TRUE)]
cat(paste(miss, collapse=" "))
RPROBE

MISSING_R_PKGS="$(_R --vanilla "$R_PROBE" 2>/dev/null || true)"
rm -f "$R_PROBE"

if $SKIP_DEPS_CHECK; then
  warn "Skipping R package check (--skip-deps-check)."
elif [[ -z "$MISSING_R_PKGS" ]]; then
  ok "All R packages already available."
else
  info "Installing missing R packages: $MISSING_R_PKGS"
  R_INSTALL="/tmp/mbx_picrust_${_TMPID}_install.R"
  # Build R vector of missing pkgs
  MISS_R="c("
  _first=true
  for _p in $MISSING_R_PKGS; do
    if $_first; then MISS_R="${MISS_R}\"$_p\""; _first=false
    else             MISS_R="${MISS_R}, \"$_p\""; fi
  done
  MISS_R="${MISS_R})"

  cat > "$R_INSTALL" << RINSTALL
options(repos = c(CRAN = "https://cloud.r-project.org"))
to_install <- ${MISS_R}
cat(sprintf("[INFO]  Installing CRAN: %s\n", paste(to_install, collapse=", ")))
install.packages(to_install, quiet = TRUE)
miss <- to_install[!sapply(to_install, requireNamespace, quietly=TRUE)]
if (length(miss) > 0) {
  cat(sprintf("[ERROR] Still missing: %s\n", paste(miss, collapse=", ")), file=stderr())
  quit(status = 1)
}
cat("[OK]    R packages installed.\n")
RINSTALL

  if $DRY_RUN; then
    warn "[DRY-RUN] Would install: $MISSING_R_PKGS"
  else
    timer_start
    _R --vanilla "$R_INSTALL" \
      || err "R package install failed.
  Try manually:
      $RSCRIPT_CMD -e 'install.packages(c($(echo $MISSING_R_PKGS | sed \"s/ /, /g\")))'"
    timer_end
  fi
  rm -f "$R_INSTALL"
fi
sep

# =============================================================================
step "4/12 — Detect categorical metadata columns"
# =============================================================================
CATS_FILE="/tmp/mbx_picrust_${_TMPID}_cats.txt"
R_DETECT="/tmp/mbx_picrust_${_TMPID}_detect.R"

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
if (any(q2_rows)) md <- md[!q2_rows, , drop = FALSE]
n_samples <- nrow(md)
cat(sprintf("[INFO]  Metadata: %d samples, %d columns.\n",
            n_samples, ncol(md)), file = stderr())

col_names <- names(md)[-1]
categorical <- character(0)

for (col in col_names) {
  vals_raw <- md[[col]]
  vals     <- vals_raw[!is.na(vals_raw) & nchar(trimws(as.character(vals_raw))) > 0]
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

if (length(categorical) == 0) {
  cat("[ERROR] No categorical columns found.\n", file = stderr())
  quit(status = 1)
}
cat(sprintf("[INFO]  Categorical columns: %d → %s\n",
            length(categorical), paste(categorical, collapse=", ")), file=stderr())
cat(paste(categorical, collapse="\n"), "\n", sep="")
RDETECT

if $DRY_RUN; then
  warn "[DRY-RUN] Would detect categorical columns."
  AUTODETECTED_COLS=( "Treatment" )
else
  _R --vanilla "$R_DETECT" "$METADATA_TXT" > "$CATS_FILE" 2>&1 \
    || err "Categorical column detection failed."
  grep "^\[" "$CATS_FILE" || true
  AUTODETECTED_COLS=()
  while IFS= read -r _col; do
    [[ -n "$_col" ]] && AUTODETECTED_COLS+=("$_col")
  done < <(grep -v "^\[" "$CATS_FILE" | grep -v '^$' || true)
  [[ ${#AUTODETECTED_COLS[@]} -eq 0 ]] && err "No categorical columns detected."
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
    $_found || err "Variable '$_v_trim' not in auto-detected list: ${AUTODETECTED_COLS[*]}"
  done
fi
ok "Variables to test: ${CATEGORICAL_COLS[*]}"

# Set up per-variable output dirs
for col in "${CATEGORICAL_COLS[@]}"; do
  _safe="$(_sanitize "$col")"
  mkdir -p "${PICRUST_DIR}/picrust2_${_safe}"
done
sep

# =============================================================================
step "5/12 — Extract rep sequences and feature table from QZAs"
# =============================================================================
REP_FASTA="${INPUTS_DIR}/representative_sequences.fasta"
FT_BIOM="${INPUTS_DIR}/feature_table.biom"
FT_TSV="${INPUTS_DIR}/feature_table.tsv"

# Rep sequences
if [[ -f "$REP_FASTA" && "$FORCE_RERUN" == false ]]; then
  skipped "rep seqs FASTA already extracted: $REP_FASTA"
else
  info "Extracting rep sequences..."
  _TMP_EXT="${INPUTS_DIR}/_tmp_repseqs"
  rm -rf "$_TMP_EXT"; mkdir -p "$_TMP_EXT"
  if ! $DRY_RUN; then
    qiime tools export --input-path "$REP_SEQS_QZA" --output-path "$_TMP_EXT" >/dev/null \
      || err "qiime tools export failed for rep seqs."
    [[ -f "$_TMP_EXT/dna-sequences.fasta" ]] || err "Expected dna-sequences.fasta not produced."
    mv "$_TMP_EXT/dna-sequences.fasta" "$REP_FASTA"
    rm -rf "$_TMP_EXT"
    N_ASVS=$(grep -c '^>' "$REP_FASTA" || echo 0)
    ok "Rep seqs: $N_ASVS ASVs → $(basename "$REP_FASTA")"
  fi
fi

# Feature table
if [[ -f "$FT_BIOM" && -f "$FT_TSV" && "$FORCE_RERUN" == false ]]; then
  skipped "feature table already extracted"
else
  info "Extracting feature table..."
  _TMP_EXT="${INPUTS_DIR}/_tmp_ft"
  rm -rf "$_TMP_EXT"; mkdir -p "$_TMP_EXT"
  if ! $DRY_RUN; then
    qiime tools export --input-path "$FT_FILTERED_QZA" --output-path "$_TMP_EXT" >/dev/null \
      || err "qiime tools export failed for feature table."
    [[ -f "$_TMP_EXT/feature-table.biom" ]] || err "Expected feature-table.biom not produced."
    mv "$_TMP_EXT/feature-table.biom" "$FT_BIOM"
    biom convert -i "$FT_BIOM" -o "$FT_TSV" --to-tsv \
      || err "biom convert failed."
    rm -rf "$_TMP_EXT"
    N_SAMPLES=$(awk -F'\t' 'NR==2 {print NF-1}' "$FT_TSV")
    N_FEATURES=$(($(wc -l < "$FT_TSV") - 2))
    ok "Feature table: $N_FEATURES features × $N_SAMPLES samples"
  fi
fi
sep

# =============================================================================
step "6/12 — Place sequences on reference tree (PICRUSt2 place_seqs.py)"
# =============================================================================
PLACED_TRE="${PLACE_DIR}/placed_seqs.tre"

if [[ -f "$PLACED_TRE" && "$FORCE_RERUN" == false ]]; then
  skipped "placed tree already exists: $(basename "$PLACED_TRE")"
else
  info "Running place_seqs.py (this can take 10-30 minutes for large datasets)..."
  if $DRY_RUN; then
    warn "[DRY-RUN] place_seqs.py"
  else
    timer_start
    _PICRUST place_seqs.py \
      -s "$REP_FASTA" \
      -o "$PLACED_TRE" \
      -p "$N_JOBS" \
      --intermediate "${INTER_DIR}/place_seqs_intermediate" \
      2>&1 | tee "${LOG_DIR}/place_seqs.log" \
      || err "place_seqs.py failed.  See: ${LOG_DIR}/place_seqs.log"
    timer_end
  fi
  ok "→ $(basename "$PLACED_TRE")"
fi
sep

# =============================================================================
step "7/12 — Hidden state prediction (HSP) for 16S, KO, EC, COG"
# =============================================================================
NSTI_TSV="${NSTI_DIR}/marker_predicted_and_nsti.tsv.gz"
KO_PRED="${INTER_DIR}/KO_predicted.tsv.gz"
EC_PRED="${INTER_DIR}/EC_predicted.tsv.gz"
COG_PRED="${INTER_DIR}/COG_predicted.tsv.gz"

# 16S marker prediction (with NSTI)
if [[ -f "$NSTI_TSV" && "$FORCE_RERUN" == false ]]; then
  skipped "16S NSTI predictions already exist"
else
  info "[1/4] hsp.py for 16S marker (computes NSTI per ASV)..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST hsp.py \
      -i 16S \
      -t "$PLACED_TRE" \
      -o "$NSTI_TSV" \
      -p "$N_JOBS" \
      -n \
      2>&1 | tee "${LOG_DIR}/hsp_16S.log" \
      || err "hsp.py for 16S failed.  See: ${LOG_DIR}/hsp_16S.log"
    timer_end
  fi
fi

# KO prediction
if [[ -f "$KO_PRED" && "$FORCE_RERUN" == false ]]; then
  skipped "KO predictions already exist"
else
  info "[2/4] hsp.py for KO (KEGG Orthology)..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST hsp.py -i KO -t "$PLACED_TRE" -o "$KO_PRED" -p "$N_JOBS" \
      2>&1 | tee "${LOG_DIR}/hsp_KO.log" \
      || err "hsp.py for KO failed."
    timer_end
  fi
fi

# EC prediction
if [[ -f "$EC_PRED" && "$FORCE_RERUN" == false ]]; then
  skipped "EC predictions already exist"
else
  info "[3/4] hsp.py for EC (Enzyme Commission)..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST hsp.py -i EC -t "$PLACED_TRE" -o "$EC_PRED" -p "$N_JOBS" \
      2>&1 | tee "${LOG_DIR}/hsp_EC.log" \
      || err "hsp.py for EC failed."
    timer_end
  fi
fi

# COG prediction.
#
# IMPORTANT: PICRUSt2 v2.6.x has a known issue with COG support:
#   • 'COG' appears in argparse choices for hsp.py
#   • but it's missing from default_tables_bac in default.py
#   • the bundled cog.txt.gz uses old IMG reference IDs incompatible with the
#     current default reference tree
#
# We auto-attempt: (1) -i COG, (2) --observed_trait_table cog.txt.gz, (3) skip.
# If the user has a newer picrust2 (v2.7+) where COG is fully restored, the
# first attempt will succeed.  Otherwise we skip cleanly.
if [[ -f "$COG_PRED" && "$FORCE_RERUN" == false ]]; then
  skipped "COG predictions already exist"
elif [[ -f "${COG_PRED}.skipped" && "$FORCE_RERUN" == false ]]; then
  skipped "COG was previously marked unavailable (delete $(basename "$COG_PRED").skipped to retry)"
else
  info "[4/4] hsp.py for COG (Clusters of Orthologous Groups)..."
  if ! $DRY_RUN; then
    timer_start
    # Attempt 1: standard -i COG (works on picrust2 v2.7+ if/when restored)
    if _PICRUST hsp.py -i COG -t "$PLACED_TRE" -o "$COG_PRED" -p "$N_JOBS" \
         2>&1 | tee "${LOG_DIR}/hsp_COG.log" \
       && [[ -s "$COG_PRED" ]]; then
      ok "  COG predictions produced via -i COG"
    else
      # Attempt 2: locate bundled cog.txt.gz, supply via --observed_trait_table
      _PIC_MOD_DIR="$(_PICRUST python -c "import picrust2, os, sys; sys.stdout.write(os.path.dirname(picrust2.__file__))" 2>/dev/null || true)"
      _COG_TT=""
      if [[ -n "$_PIC_MOD_DIR" && -f "${_PIC_MOD_DIR}/default_files/prokaryotic/cog.txt.gz" ]]; then
        _COG_TT="${_PIC_MOD_DIR}/default_files/prokaryotic/cog.txt.gz"
      fi
      if [[ -n "$_COG_TT" ]] && \
         _PICRUST hsp.py --observed_trait_table "$_COG_TT" \
                 -t "$PLACED_TRE" -o "$COG_PRED" -p "$N_JOBS" \
                 >> "${LOG_DIR}/hsp_COG.log" 2>&1 \
         && [[ -s "$COG_PRED" ]]; then
        ok "  COG predictions produced via custom trait table"
      else
        warn "COG database not available in this PICRUSt2 install.

  Cause:    PICRUSt2 v2.6.x has incomplete COG support.  The bundled
            cog.txt.gz uses reference IDs that don't match the placement tree.
  Effect:   COG analysis will be skipped.  KO, EC, and MetaCyc results are
            complete and unaffected.
  To enable COG: upgrade PICRUSt2 to v2.7+ (when released), or install a
            version that has matching cog.txt.gz reference IDs."
        # Mark as skipped (zero-byte file with .skipped extension)
        rm -f "$COG_PRED"
        : > "${COG_PRED}.skipped"
      fi
    fi
    timer_end
  fi
fi
sep

# =============================================================================
step "8/12 — Predict per-sample metagenomes (with NSTI filter = $NSTI_THRESHOLD)"
# =============================================================================
KO_UNSTRAT="${KO_OUT_DIR}/pred_metagenome_unstrat.tsv.gz"
EC_UNSTRAT="${EC_OUT_DIR}/pred_metagenome_unstrat.tsv.gz"
COG_UNSTRAT="${COG_OUT_DIR}/pred_metagenome_unstrat.tsv.gz"

# KO metagenome
if [[ -f "$KO_UNSTRAT" && "$FORCE_RERUN" == false ]]; then
  skipped "KO metagenome already predicted"
else
  info "[1/3] metagenome_pipeline.py for KO..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST metagenome_pipeline.py \
      -i "$FT_BIOM" \
      -m "$NSTI_TSV" \
      -f "$KO_PRED" \
      -o "$KO_OUT_DIR" \
      --max_nsti "$NSTI_THRESHOLD" \
      2>&1 | tee "${LOG_DIR}/metagenome_KO.log" \
      || err "metagenome_pipeline for KO failed."
    timer_end
  fi
fi

# EC metagenome
if [[ -f "$EC_UNSTRAT" && "$FORCE_RERUN" == false ]]; then
  skipped "EC metagenome already predicted"
else
  info "[2/3] metagenome_pipeline.py for EC..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST metagenome_pipeline.py \
      -i "$FT_BIOM" \
      -m "$NSTI_TSV" \
      -f "$EC_PRED" \
      -o "$EC_OUT_DIR" \
      --max_nsti "$NSTI_THRESHOLD" \
      2>&1 | tee "${LOG_DIR}/metagenome_EC.log" \
      || err "metagenome_pipeline for EC failed."
    timer_end
  fi
fi

# COG metagenome (skip if COG HSP wasn't possible)
if [[ -f "$COG_UNSTRAT" && "$FORCE_RERUN" == false ]]; then
  skipped "COG metagenome already predicted"
elif [[ -f "${COG_PRED}.skipped" ]]; then
  warn "COG HSP was skipped (database missing) — skipping COG metagenome too."
elif [[ ! -f "$COG_PRED" ]]; then
  warn "COG predictions not produced — skipping COG metagenome."
else
  info "[3/3] metagenome_pipeline.py for COG..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST metagenome_pipeline.py \
      -i "$FT_BIOM" \
      -m "$NSTI_TSV" \
      -f "$COG_PRED" \
      -o "$COG_OUT_DIR" \
      --max_nsti "$NSTI_THRESHOLD" \
      2>&1 | tee "${LOG_DIR}/metagenome_COG.log" \
      || err "metagenome_pipeline for COG failed."
    timer_end
  fi
fi

# Save weighted NSTI files into the dedicated NSTI folder for easy access
for _src in "$KO_OUT_DIR/weighted_nsti.tsv.gz" "$EC_OUT_DIR/weighted_nsti.tsv.gz" "$COG_OUT_DIR/weighted_nsti.tsv.gz"; do
  [[ -f "$_src" ]] || continue
  cp -n "$_src" "$NSTI_DIR/$(basename "$(dirname "$_src")")_weighted_nsti.tsv.gz" 2>/dev/null || true
done
sep

# =============================================================================
step "9/12 — MetaCyc pathway abundance (pathway_pipeline.py)"
# =============================================================================
PATH_UNSTRAT="${PATH_OUT_DIR}/path_abun_unstrat.tsv.gz"

if [[ -f "$PATH_UNSTRAT" && "$FORCE_RERUN" == false ]]; then
  skipped "MetaCyc pathway abundance already computed"
else
  info "Running pathway_pipeline.py (uses EC predictions to infer MetaCyc pathways)..."
  if ! $DRY_RUN; then
    timer_start
    _PICRUST pathway_pipeline.py \
      -i "$EC_UNSTRAT" \
      -o "$PATH_OUT_DIR" \
      --intermediate "${INTER_DIR}/pathway_intermediate" \
      -p "$N_JOBS" \
      2>&1 | tee "${LOG_DIR}/pathway_pipeline.log" \
      || err "pathway_pipeline.py failed."
    timer_end
  fi
fi
sep

# =============================================================================
step "10/12 — Add functional descriptions (KO names, EC names, MetaCyc names)"
# =============================================================================
_add_descriptions() {
  local input="$1" mtype="$2" output="$3" log="$4"
  if [[ -f "$output" && "$FORCE_RERUN" == false ]]; then
    skipped "  $(basename "$output") already exists"
    return 0
  fi
  if ! $DRY_RUN; then
    _PICRUST add_descriptions.py \
      -i "$input" -m "$mtype" -o "$output" \
      2>&1 | tee "$log" \
      || err "add_descriptions.py failed for $mtype"
  fi
  ok "  → $(basename "$output")"
}

info "[1/4] KO descriptions..."
_add_descriptions "$KO_UNSTRAT"   "KO"      "${KO_OUT_DIR}/pred_metagenome_unstrat_descrip.tsv.gz" "${LOG_DIR}/desc_KO.log"
info "[2/4] EC descriptions..."
_add_descriptions "$EC_UNSTRAT"   "EC"      "${EC_OUT_DIR}/pred_metagenome_unstrat_descrip.tsv.gz" "${LOG_DIR}/desc_EC.log"
if [[ -f "$COG_UNSTRAT" ]]; then
  info "[3/4] COG descriptions..."
  _add_descriptions "$COG_UNSTRAT"  "COG"     "${COG_OUT_DIR}/pred_metagenome_unstrat_descrip.tsv.gz" "${LOG_DIR}/desc_COG.log"
else
  warn "[3/4] Skipping COG descriptions (COG predictions not available)."
fi
info "[4/4] MetaCyc descriptions..."
_add_descriptions "$PATH_UNSTRAT" "METACYC" "${PATH_OUT_DIR}/path_abun_unstrat_descrip.tsv.gz" "${LOG_DIR}/desc_METACYC.log"
sep

# =============================================================================
step "11/12 — NSTI analysis + DA + plots + correlation (R)"
# =============================================================================
R_ANALYSIS="/tmp/mbx_picrust_${_TMPID}_analysis.R"

# Build R-side variable list
VARS_R=""
for col in "${CATEGORICAL_COLS[@]}"; do
  _esc="${col//\"/\\\"}"
  if [[ -z "$VARS_R" ]]; then VARS_R="\"$_esc\""
  else                        VARS_R="$VARS_R, \"$_esc\""; fi
done

if $SKIP_PLOTS;       then SKIP_PLOTS_R="TRUE";        else SKIP_PLOTS_R="FALSE";        fi
if $DO_CORRELATION;   then DO_CORR_R="TRUE";           else DO_CORR_R="FALSE";           fi
if $FORCE_RERUN;      then FORCE_RERUN_R="TRUE";       else FORCE_RERUN_R="FALSE";       fi
if $HAS_ANCOMBC;      then HAS_ANCOMBC_R="TRUE";       else HAS_ANCOMBC_R="FALSE";       fi

# ANCOMBC dir for correlation; safe to pass even if not used
ANCOMBC_DIR_R="${ANCOMBC_DIR:-}"

cat > "$R_ANALYSIS" << RANALYSIS
# =============================================================================
# PICRUSt2 downstream analysis  (NSTI + DA + plots + correlation)
# =============================================================================
suppressPackageStartupMessages({
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(dunn.test)
  library(vegan)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
})

# ── Inputs (injected by bash) ───────────────────────────────────────────────
METADATA      <- "${METADATA_TXT}"
PICRUST_DIR   <- "${PICRUST_DIR}"
KO_FILE       <- "${KO_OUT_DIR}/pred_metagenome_unstrat_descrip.tsv.gz"
EC_FILE       <- "${EC_OUT_DIR}/pred_metagenome_unstrat_descrip.tsv.gz"
COG_FILE      <- "${COG_OUT_DIR}/pred_metagenome_unstrat_descrip.tsv.gz"
PATH_FILE     <- "${PATH_OUT_DIR}/path_abun_unstrat_descrip.tsv.gz"
NSTI_FILE     <- "${NSTI_TSV}"
FT_TSV        <- "${FT_TSV}"
WEIGHTED_KO   <- "${KO_OUT_DIR}/weighted_nsti.tsv.gz"
WEIGHTED_EC   <- "${EC_OUT_DIR}/weighted_nsti.tsv.gz"
NSTI_DIR      <- "${NSTI_DIR}"
NSTI_XLSX     <- "${NSTI_XLSX}"
RDS_DIR       <- "${RDS_DIR}"
NSTI_FILT_LOG <- "${NSTI_DIR}/nsti_filtering_summary.txt"
NSTI_THRESH   <- ${NSTI_THRESHOLD}

VARS         <- c(${VARS_R})
SKIP_PLOTS   <- ${SKIP_PLOTS_R}
DO_CORR      <- ${DO_CORR_R}
HAS_ANCOMBC  <- ${HAS_ANCOMBC_R}
ANCOMBC_DIR  <- "${ANCOMBC_DIR_R}"
FORCE_RERUN  <- ${FORCE_RERUN_R}

PICRUST_SUBTITLE <- "Predicted functional profiles - inferred from 16S rRNA data"

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

# Read PICRUSt2 unstrat_descrip TSV — first col = function ID, second = description,
# remaining cols = samples.
read_picrust_table <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  df <- read.delim(con, header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE, comment.char = "")
  if (ncol(df) < 3) stop(sprintf("Expected >=3 columns in %s", path))
  desc_col <- if ("description" %in% tolower(names(df))) {
    which(tolower(names(df)) == "description")[1]
  } else 2L
  ids   <- as.character(df[[1]])
  descs <- as.character(df[[desc_col]])
  vals  <- df[, -c(1, desc_col), drop = FALSE]
  for (j in seq_len(ncol(vals))) vals[[j]] <- as.numeric(vals[[j]])
  abun  <- as.matrix(vals)
  rownames(abun) <- ids
  list(ids = ids, descs = descs, abun = abun)
}

# Save a ggplot as both PNG and PDF
save_ggplot <- function(p, base_path, w = 9, h = 6.5) {
  tryCatch({
    ggsave(paste0(base_path, ".png"), p, width = w, height = h,
           dpi = 300, bg = "white")
    ggsave(paste0(base_path, ".pdf"), p, width = w, height = h,
           device = cairo_pdf, bg = "white")
  }, error = function(e) {
    # fallback: try regular pdf device if cairo unavailable
    ggsave(paste0(base_path, ".png"), p, width = w, height = h,
           dpi = 300, bg = "white")
    ggsave(paste0(base_path, ".pdf"), p, width = w, height = h,
           bg = "white")
  })
}

# Save pheatmap as both PNG and PDF (pheatmap takes filename argument)
save_pheatmap <- function(mat, filename_base, ..., w = 10, h = 7) {
  png_file <- paste0(filename_base, ".png")
  pdf_file <- paste0(filename_base, ".pdf")
  pheatmap(mat, ..., filename = png_file, width = w, height = h)
  pheatmap(mat, ..., filename = pdf_file, width = w, height = h)
}

# ───────────────────────────────────────────────────────────────────────────
# Part A — NSTI analysis
# ───────────────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Part A: NSTI filtering analysis\n")
cat("════════════════════════════════════════════════════════════════════\n")

con <- gzfile(NSTI_FILE, "rt")
nsti_df <- read.delim(con, header = TRUE, check.names = FALSE,
                      stringsAsFactors = FALSE, comment.char = "")
close(con)
# Expected columns: sequence, 16S_rRNA_Count, metadata_NSTI
seq_col  <- 1L
nsti_col <- which(tolower(names(nsti_df)) == "metadata_nsti")[1]
if (is.na(nsti_col))
  nsti_col <- which(grepl("nsti", names(nsti_df), ignore.case = TRUE))[1]
if (is.na(nsti_col))
  stop(sprintf("Could not find NSTI column in %s; columns: %s",
               NSTI_FILE, paste(names(nsti_df), collapse=", ")))

names(nsti_df)[seq_col]  <- "sequence"
names(nsti_df)[nsti_col] <- "metadata_NSTI"
nsti_df\$metadata_NSTI <- as.numeric(nsti_df\$metadata_NSTI)

# Read feature table (skip biom convert comment line)
ft_lines <- readLines(FT_TSV)
ft_hdr_idx <- grep("^#OTU", ft_lines, ignore.case = FALSE)
if (length(ft_hdr_idx) == 0) ft_hdr_idx <- which(!grepl("^#", ft_lines))[1] - 1
hdr <- strsplit(sub("^#", "", ft_lines[ft_hdr_idx]), "\t", fixed = TRUE)[[1]]
body_lines <- ft_lines[-seq_len(ft_hdr_idx)]
body_lines <- body_lines[nchar(body_lines) > 0 & !grepl("^#", body_lines)]
body <- do.call(rbind, lapply(strsplit(body_lines, "\t", fixed = TRUE), function(x) {
  if (length(x) < length(hdr)) x <- c(x, rep("", length(hdr) - length(x)))
  x[seq_along(hdr)]
}))
colnames(body) <- hdr
body <- as.data.frame(body, stringsAsFactors = FALSE, check.names = FALSE)
rownames(body) <- body[[1]]
body[[1]] <- NULL
for (j in seq_len(ncol(body))) body[[j]] <- as.numeric(body[[j]])
ft_mat <- as.matrix(body)
ft_mat[is.na(ft_mat)] <- 0

# Filter NSTI
total_asvs       <- nrow(nsti_df)
mean_nsti_before <- mean(nsti_df\$metadata_NSTI, na.rm = TRUE)
max_nsti_before  <- max(nsti_df\$metadata_NSTI, na.rm = TRUE)
median_nsti_before <- median(nsti_df\$metadata_NSTI, na.rm = TRUE)

dropped <- nsti_df[nsti_df\$metadata_NSTI > NSTI_THRESH & !is.na(nsti_df\$metadata_NSTI), ]
kept    <- nsti_df[nsti_df\$metadata_NSTI <= NSTI_THRESH & !is.na(nsti_df\$metadata_NSTI), ]
n_dropped <- nrow(dropped)
n_kept    <- nrow(kept)
mean_nsti_after <- if (n_kept > 0) mean(kept\$metadata_NSTI) else NA
max_nsti_after  <- if (n_kept > 0) max(kept\$metadata_NSTI) else NA

# % reads dropped overall and per-sample
samples_in_ft <- colnames(ft_mat)
asvs_dropped  <- as.character(dropped\$sequence)
asvs_in_ft_dropped <- intersect(asvs_dropped, rownames(ft_mat))

reads_dropped_per_sample <- if (length(asvs_in_ft_dropped) > 0) {
  colSums(ft_mat[asvs_in_ft_dropped, , drop = FALSE])
} else setNames(rep(0L, length(samples_in_ft)), samples_in_ft)

total_reads_per_sample <- colSums(ft_mat)
pct_lost_per_sample <- 100 * reads_dropped_per_sample / pmax(total_reads_per_sample, 1)

total_reads        <- sum(total_reads_per_sample)
reads_dropped_total <- sum(reads_dropped_per_sample)
pct_reads_dropped  <- 100 * reads_dropped_total / max(total_reads, 1)

# Top 5 most-affected samples
affected_idx <- order(-pct_lost_per_sample)
top_affected <- head(affected_idx, 5)

# Per-sample weighted NSTI (from PICRUSt2 EC weighted_nsti file — equivalent across functional categories)
weighted_per_sample <- NULL
if (file.exists(WEIGHTED_EC)) {
  con <- gzfile(WEIGHTED_EC, "rt")
  wn  <- read.delim(con, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, comment.char = "")
  close(con)
  # Look for the NSTI column (commonly "weighted_NSTI")
  vcol <- which(grepl("nsti", names(wn), ignore.case = TRUE))
  scol <- 1L
  if (length(vcol) > 0) {
    weighted_per_sample <- setNames(as.numeric(wn[[vcol[1]]]),
                                    as.character(wn[[scol]]))
  }
}

# Build per-sample summary table
samples <- samples_in_ft
sample_df <- data.frame(
  sample_id              = samples,
  total_reads            = total_reads_per_sample[samples],
  reads_after_NSTI       = total_reads_per_sample[samples] - reads_dropped_per_sample[samples],
  reads_dropped          = reads_dropped_per_sample[samples],
  pct_reads_dropped      = round(pct_lost_per_sample[samples], 2),
  weighted_mean_NSTI     = if (!is.null(weighted_per_sample))
                              round(weighted_per_sample[samples], 4) else NA_real_,
  flag_high_NSTI         = if (!is.null(weighted_per_sample))
                              ifelse(weighted_per_sample[samples] > 1.0, "YES", ""),
  flag_50pct_reads_lost  = ifelse(pct_lost_per_sample[samples] > 50, "YES", ""),
  stringsAsFactors       = FALSE)
rownames(sample_df) <- NULL

# Write per-sample NSTI summary xlsx
wb <- createWorkbook()
addWorksheet(wb, "per_sample_NSTI")
high_style <- createStyle(fgFill = "#ffe0b3", textDecoration = "bold")
critical_style <- createStyle(fgFill = "#ff9999", textDecoration = "bold")
writeData(wb, 1, sample_df,
          headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
# Highlight high-NSTI rows
high_rows <- which(sample_df\$flag_high_NSTI == "YES")
crit_rows <- which(sample_df\$flag_50pct_reads_lost == "YES")
if (length(high_rows) > 0)
  addStyle(wb, 1, high_style, rows = high_rows + 1, cols = 1:ncol(sample_df), gridExpand = TRUE)
if (length(crit_rows) > 0)
  addStyle(wb, 1, critical_style, rows = crit_rows + 1, cols = 1:ncol(sample_df), gridExpand = TRUE)
freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 2)
setColWidths(wb, 1, cols = 1:ncol(sample_df), widths = "auto")
saveWorkbook(wb, NSTI_XLSX, overwrite = TRUE)

# Write nsti_filtering_summary.txt
sink(NSTI_FILT_LOG)
cat("# ============================================================================\n")
cat("# NSTI filtering summary  —  PICRUSt2 step 15\n")
cat(sprintf("# Generated : %s\n", format(Sys.time())))
cat(sprintf("# Threshold : %.3f  (ASVs with NSTI > threshold are dropped)\n", NSTI_THRESH))
cat("# ============================================================================\n\n")

cat(sprintf("Total ASVs                       : %d\n", total_asvs))
cat(sprintf("ASVs dropped (NSTI > %.2f)       : %d  (%.2f%%)\n",
            NSTI_THRESH, n_dropped, 100 * n_dropped / total_asvs))
cat(sprintf("ASVs kept                        : %d  (%.2f%%)\n",
            n_kept, 100 * n_kept / total_asvs))
cat(sprintf("Total reads (all samples)        : %d\n", total_reads))
cat(sprintf("Reads dropped overall            : %d  (%.2f%%)\n",
            reads_dropped_total, pct_reads_dropped))
cat("\n")
cat(sprintf("Mean NSTI before filter          : %.4f\n", mean_nsti_before))
cat(sprintf("Median NSTI before filter        : %.4f\n", median_nsti_before))
cat(sprintf("Max  NSTI before filter          : %.4f\n", max_nsti_before))
cat(sprintf("Mean NSTI after  filter          : %.4f\n", mean_nsti_after))
cat(sprintf("Max  NSTI after  filter          : %.4f\n", max_nsti_after))
cat("\n")
cat("Top 5 most-affected samples (by % reads dropped):\n")
for (i in top_affected) {
  cat(sprintf("  %-30s  %.2f%% reads dropped\n",
              samples[i], pct_lost_per_sample[i]))
}
cat("\n")

if (!is.null(weighted_per_sample)) {
  high_samples <- names(weighted_per_sample)[weighted_per_sample > 1.0]
  if (length(high_samples) > 0) {
    cat("⚠ HIGH-NSTI SAMPLES (weighted mean NSTI > 1.0):\n")
    cat("   Functional predictions for these samples are UNRELIABLE.\n")
    cat("   Interpret results with caution — consider excluding them.\n")
    for (s in high_samples)
      cat(sprintf("   - %-30s  weighted NSTI = %.4f\n",
                  s, weighted_per_sample[s]))
    cat("\n")
  } else {
    cat("✓ No samples with weighted mean NSTI > 1.0 — all predictions reliable.\n\n")
  }
}

samples_50pct <- samples[pct_lost_per_sample > 50]
if (length(samples_50pct) > 0) {
  cat("⚠ SAMPLES LOSING >50% OF READS AFTER NSTI FILTERING:\n")
  cat("   Strongly consider excluding these from downstream DA analysis.\n")
  for (s in samples_50pct)
    cat(sprintf("   - %-30s  %.2f%% reads dropped\n", s, pct_lost_per_sample[s]))
  cat("\n")
}

cat(sprintf("Per-sample table written to: %s\n", basename(NSTI_XLSX)))
sink()

# Echo summary to console
cat(readLines(NSTI_FILT_LOG), sep = "\n")
cat("\n")

# ───────────────────────────────────────────────────────────────────────────
# Part B — Read all 4 functional databases, prep DA loop
# ───────────────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Part B: Loading functional abundance tables\n")
cat("════════════════════════════════════════════════════════════════════\n")

databases <- list(
  KO        = list(file = KO_FILE,   label = "KO",      mtype = "KO"),
  EC        = list(file = EC_FILE,   label = "EC",      mtype = "EC"),
  COG       = list(file = COG_FILE,  label = "COG",     mtype = "COG"),
  metacyc   = list(file = PATH_FILE, label = "MetaCyc", mtype = "METACYC"))

# Drop any database whose file is missing (e.g. COG on older picrust2 versions)
keep_dbs <- character(0)
for (db_key in names(databases)) {
  db <- databases[[db_key]]
  if (!file.exists(db\$file)) {
    cat(sprintf("[WARN]  %s file not found (%s) — DROPPING from analysis.\n",
                db\$label, basename(db\$file)))
    next
  }
  cat(sprintf("[INFO]  Loading %s table from %s\n",
              db\$label, basename(db\$file)))
  pt <- read_picrust_table(db\$file)
  databases[[db_key]]\$ids   <- pt\$ids
  databases[[db_key]]\$descs <- pt\$descs
  databases[[db_key]]\$abun  <- pt\$abun
  cat(sprintf("        -> %d features x %d samples\n",
              nrow(pt\$abun), ncol(pt\$abun)))
  keep_dbs <- c(keep_dbs, db_key)
}
databases <- databases[keep_dbs]
if (!"metacyc" %in% names(databases))
  stop("[FATAL] MetaCyc pathway table is required but missing.")

# Read metadata
md <- read_metadata(METADATA)

# Bookkeeping for HTML report
html_sections <- list()

# ───────────────────────────────────────────────────────────────────────────
# Part C — DA + plots per (variable × database)
# ───────────────────────────────────────────────────────────────────────────
cat("\n════════════════════════════════════════════════════════════════════\n")
cat("  Part C: Differential abundance per variable × database\n")
cat("════════════════════════════════════════════════════════════════════\n")

# Globally store DA pathways per variable for later genus × pathway correlation
da_path_per_var <- list()

for (var in VARS) {
  safe_var <- sanitize(var)
  var_dir  <- file.path(PICRUST_DIR, sprintf("picrust2_%s", safe_var))
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n━━━ Variable: %s  →  %s/  ━━━\n", var, basename(var_dir)))

  # Filter metadata to samples in feature tables; drop NA in this variable
  abun_samples <- colnames(databases\$KO\$abun)
  md_var <- md[md[["sample-id"]] %in% abun_samples, , drop = FALSE]
  md_var <- md_var[!is.na(md_var[[var]]) &
                   nchar(trimws(as.character(md_var[[var]]))) > 0, , drop = FALSE]
  md_var[[var]] <- factor(trimws(as.character(md_var[[var]])),
                          levels = sort(unique(trimws(as.character(md_var[[var]])))))
  grp_n <- table(md_var[[var]])
  bad_groups <- names(grp_n)[grp_n < 2]
  if (length(bad_groups) > 0) {
    md_var <- md_var[!md_var[[var]] %in% bad_groups, , drop = FALSE]
    md_var[[var]] <- droplevels(md_var[[var]])
  }
  if (length(unique(md_var[[var]])) < 2) {
    cat(sprintf("[WARN]  <2 valid groups for %s — skipping.\n", var))
    next
  }
  groups <- levels(md_var[[var]])
  k <- length(groups)
  samp_keep <- md_var[["sample-id"]]
  cat(sprintf("[INFO]  %d samples, %d groups (%s)\n",
              length(samp_keep), k, paste(groups, collapse = ", ")))

  # Per-database DA
  da_results <- list()
  for (db_key in names(databases)) {
    db <- databases[[db_key]]
    abun <- db\$abun[, samp_keep, drop = FALSE]
    abun[is.na(abun)] <- 0
    # Drop features with zero total abundance
    zero_features <- which(rowSums(abun) == 0)
    if (length(zero_features) > 0) abun <- abun[-zero_features, , drop = FALSE]
    # Drop features with zero variance within ANY group
    grp_vec <- as.character(md_var[[var]])
    g_var <- sapply(unique(grp_vec), function(g) {
      idx <- which(grp_vec == g)
      if (length(idx) < 2) return(rep(NA_real_, nrow(abun)))
      apply(abun[, idx, drop = FALSE], 1, stats::var)
    })
    if (is.null(dim(g_var))) g_var <- matrix(g_var, nrow = nrow(abun))
    min_within <- suppressWarnings(apply(g_var, 1, min, na.rm = TRUE))
    bad_features <- which(is.na(min_within) | !is.finite(min_within) | min_within == 0)
    # For DA we keep these (as they're informative for KW), only filter zero-everywhere

    cat(sprintf("[RUN]   %s  : %d features × %d samples — running KW + Dunn ...\n",
                db\$label, nrow(abun), ncol(abun)))

    # KW + post-hoc per feature
    n_features <- nrow(abun)
    pvals <- numeric(n_features)
    eps2  <- numeric(n_features)
    grp_means <- matrix(NA_real_, nrow = n_features, ncol = k,
                        dimnames = list(rownames(abun), groups))
    for (i in seq_len(n_features)) {
      vals <- as.numeric(abun[i, ])
      kw <- tryCatch(kruskal.test(vals ~ md_var[[var]]),
                     error = function(e) NULL)
      pvals[i] <- if (is.null(kw)) NA_real_ else kw\$p.value
      # epsilon-squared (effect size)
      n_total <- length(vals)
      eps2[i] <- if (is.null(kw)) NA_real_
                 else (kw\$statistic - k + 1) / (n_total - k)
      for (gi in seq_along(groups)) {
        grp_means[i, gi] <- mean(vals[grp_vec == groups[gi]], na.rm = TRUE)
      }
    }
    qvals <- p.adjust(pvals, method = "BH")

    # LFC (log2 of group mean / overall mean)
    overall_mean <- pmax(rowMeans(abun, na.rm = TRUE), 1e-10)
    log2_fc <- log2(pmax(grp_means, 1e-10) / overall_mean)
    colnames(log2_fc) <- paste0("log2FC_", colnames(log2_fc))

    # Pairwise (k≥3 use Dunn; k=2 use Wilcoxon)
    if (k >= 3) {
      # Build pairwise data
      pair_names <- combn(groups, 2, paste, collapse = "__vs__")
      n_pairs <- length(pair_names)
      pair_p <- matrix(NA_real_, nrow = n_features, ncol = n_pairs,
                       dimnames = list(rownames(abun), pair_names))
      pair_lfc <- matrix(NA_real_, nrow = n_features, ncol = n_pairs,
                         dimnames = list(rownames(abun), pair_names))
      for (i in seq_len(n_features)) {
        vals <- as.numeric(abun[i, ])
        if (length(unique(vals)) < 2) next
        dt <- tryCatch(suppressWarnings(suppressMessages(
                 dunn.test::dunn.test(vals, md_var[[var]], method = "bh", kw = FALSE))),
               error = function(e) NULL)
        if (is.null(dt)) next
        # dt\$comparisons format: "groupA - groupB" — match to pair_names "groupA__vs__groupB"
        for (cp in seq_along(dt\$comparisons)) {
          parts <- strsplit(dt\$comparisons[cp], " - ", fixed = TRUE)[[1]]
          if (length(parts) != 2) next
          a <- parts[1]; b <- parts[2]
          # Find matching column (ordered alphabetically by combn)
          key <- if (a < b) paste0(a, "__vs__", b) else paste0(b, "__vs__", a)
          if (!key %in% pair_names) next
          pair_p[i, key]   <- dt\$P.adjusted[cp]
          pair_lfc[i, key] <- log2(pmax(grp_means[i, a], 1e-10) /
                                   pmax(grp_means[i, b], 1e-10))
        }
      }
    } else {
      # k = 2: Wilcoxon
      pair_names <- paste0(groups[1], "__vs__", groups[2])
      pair_p <- matrix(NA_real_, nrow = n_features, ncol = 1,
                       dimnames = list(rownames(abun), pair_names))
      pair_lfc <- matrix(NA_real_, nrow = n_features, ncol = 1,
                         dimnames = list(rownames(abun), pair_names))
      for (i in seq_len(n_features)) {
        vals <- as.numeric(abun[i, ])
        a_vals <- vals[grp_vec == groups[1]]
        b_vals <- vals[grp_vec == groups[2]]
        wt <- tryCatch(suppressWarnings(wilcox.test(a_vals, b_vals)),
                       error = function(e) NULL)
        if (!is.null(wt)) {
          pair_p[i, 1]   <- wt\$p.value
          pair_lfc[i, 1] <- log2(pmax(mean(a_vals, na.rm=TRUE), 1e-10) /
                                 pmax(mean(b_vals, na.rm=TRUE), 1e-10))
        }
      }
      # BH-adjust
      pair_p[, 1] <- p.adjust(pair_p[, 1], method = "BH")
    }

    # Build KW result table
    kw_df <- data.frame(
      feature_id  = rownames(abun),
      description = db\$descs[match(rownames(abun), db\$ids)],
      KW_pvalue   = pvals,
      KW_qvalue   = qvals,
      epsilon_sq  = eps2,
      stringsAsFactors = FALSE,
      check.names      = FALSE)
    kw_df <- cbind(kw_df, as.data.frame(grp_means, check.names = FALSE),
                          as.data.frame(log2_fc,   check.names = FALSE))
    kw_df\$signif_q05 <- ifelse(!is.na(kw_df\$KW_qvalue) & kw_df\$KW_qvalue < 0.05,
                                "YES", "")
    kw_df <- kw_df[order(kw_df\$KW_qvalue), ]

    # Build pairwise result table
    pw_df <- data.frame(
      feature_id  = rownames(abun),
      description = db\$descs[match(rownames(abun), db\$ids)],
      stringsAsFactors = FALSE, check.names = FALSE)
    pw_df <- cbind(pw_df,
                   as.data.frame(pair_lfc, check.names = FALSE) |>
                     setNames(paste0("log2FC_", colnames(pair_lfc))),
                   as.data.frame(pair_p, check.names = FALSE) |>
                     setNames(paste0("q_", colnames(pair_p))))
    pw_df <- pw_df[order(rowMeans(pair_p, na.rm = TRUE)), ]

    # Save xlsx
    kw_xlsx <- file.path(var_dir, sprintf("DA_%s_KW_%s.xlsx", db\$mtype, safe_var))
    pw_xlsx <- file.path(var_dir, sprintf("DA_%s_pairwise_%s.xlsx", db\$mtype, safe_var))

    wb <- createWorkbook()
    addWorksheet(wb, sprintf("KW_%s", substr(db\$mtype, 1, 18)))
    sig_rows <- which(kw_df\$signif_q05 == "YES")
    writeData(wb, 1, kw_df,
              headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
    if (length(sig_rows) > 0)
      addStyle(wb, 1,
               createStyle(fgFill = "#ffe0b3"),
               rows = sig_rows + 1, cols = 1:ncol(kw_df), gridExpand = TRUE)
    freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 3)
    setColWidths(wb, 1, cols = 1:ncol(kw_df), widths = "auto")
    saveWorkbook(wb, kw_xlsx, overwrite = TRUE)

    wb <- createWorkbook()
    addWorksheet(wb, "pairwise")
    writeData(wb, 1, pw_df,
              headerStyle = createStyle(textDecoration = "bold", border = "Bottom"))
    freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 3)
    setColWidths(wb, 1, cols = 1:ncol(pw_df), widths = "auto")
    saveWorkbook(wb, pw_xlsx, overwrite = TRUE)

    n_sig <- sum(kw_df\$signif_q05 == "YES")
    cat(sprintf("[OK]    %s : %d significant features  →  %s\n",
                db\$label, n_sig, basename(kw_xlsx)))

    da_results[[db_key]] <- list(
      kw_df    = kw_df,
      pw_df    = pw_df,
      pair_lfc = pair_lfc,
      pair_p   = pair_p,
      abun     = abun,
      grp_means = grp_means,
      n_sig    = n_sig)

    # Save RDS
    saveRDS(da_results[[db_key]],
            file.path(RDS_DIR, sprintf("DA_%s_%s.rds", db\$mtype, safe_var)))
  }

  # Store DA pathways for later genus × pathway correlation
  if (!is.null(da_results\$metacyc)) {
    da_path_per_var[[safe_var]] <- list(
      var = var,
      sig_pathways = da_results\$metacyc\$kw_df\$feature_id[
        da_results\$metacyc\$kw_df\$signif_q05 == "YES"],
      kw_df = da_results\$metacyc\$kw_df,
      abun  = da_results\$metacyc\$abun,
      md_var = md_var)
  }

  # ── Plots ────────────────────────────────────────────────────────────────
  if (!SKIP_PLOTS) {
    # ── 1. Stacked bar plot of top 20 MetaCyc pathways ────────────────────
    tryCatch({
      mc <- databases\$metacyc
      mc_abun <- mc\$abun[, samp_keep, drop = FALSE]
      mc_abun[is.na(mc_abun)] <- 0
      # Compute mean abundance per group, take top 20 overall
      grp_mean_mat <- sapply(groups, function(g)
        rowMeans(mc_abun[, grp_vec == g, drop = FALSE], na.rm = TRUE))
      overall <- rowMeans(grp_mean_mat, na.rm = TRUE)
      top20_idx <- order(-overall)[seq_len(min(20, length(overall)))]
      top_paths <- rownames(mc_abun)[top20_idx]
      top_descs <- mc\$descs[match(top_paths, mc\$ids)]
      top_descs[is.na(top_descs) | top_descs == ""] <- top_paths[is.na(top_descs) | top_descs == ""]

      # Long-format normalized to 100% per group
      grp_top <- grp_mean_mat[top20_idx, , drop = FALSE]
      grp_top_norm <- t(t(grp_top) / colSums(grp_top)) * 100
      df_long <- as.data.frame(grp_top_norm) |>
        tibble::rownames_to_column(var = "feature_id") -> df_norm_wide
      # tibble may not be loaded — use base
      df_long <- data.frame(
        feature_id = rep(top_paths, ncol(grp_top_norm)),
        Group      = rep(colnames(grp_top_norm), each = nrow(grp_top_norm)),
        pct        = as.vector(grp_top_norm),
        stringsAsFactors = FALSE)
      # Truncate description for display
      desc_map <- setNames(
        ifelse(nchar(top_descs) > 50,
               paste0(substr(top_descs, 1, 47), "..."),
               top_descs),
        top_paths)
      df_long\$Pathway <- desc_map[df_long\$feature_id]
      df_long\$Pathway <- factor(df_long\$Pathway, levels = unique(desc_map[top_paths]))

      n_paths <- length(top_paths)
      pal <- if (n_paths <= 12) brewer.pal(max(3, n_paths), "Set3")[seq_len(n_paths)]
             else colorRampPalette(brewer.pal(12, "Set3"))(n_paths)

      p <- ggplot(df_long, aes(x = Group, y = pct, fill = Pathway)) +
        geom_bar(stat = "identity", colour = "white", linewidth = 0.15) +
        scale_fill_manual(values = pal) +
        labs(title    = sprintf("Top 20 MetaCyc pathways by %s", var),
             subtitle = PICRUST_SUBTITLE,
             x = NULL, y = "% of mean abundance",
             fill = "Pathway") +
        theme_classic(base_size = 11) +
        theme(plot.title    = element_text(face = "bold"),
              plot.subtitle = element_text(colour = "grey30", size = 9),
              axis.text.x   = element_text(angle = 30, hjust = 1),
              legend.text   = element_text(size = 7),
              legend.key.size = unit(0.4, "cm"),
              legend.position = "right") +
        guides(fill = guide_legend(ncol = 1, reverse = FALSE))

      base_path <- file.path(var_dir, sprintf("stacked_bar_top20_metacyc_%s", safe_var))
      save_ggplot(p, base_path,
                  w = 10, h = max(6, 0.5 + 0.6 * length(groups)))
      cat(sprintf("[OK]    stacked bar  → %s.{png,pdf}\n", basename(base_path)))
    }, error = function(e)
      cat(sprintf("[WARN]  Stacked bar failed: %s\n", conditionMessage(e))))

    # ── 2. DA pathways heatmap (lfc + q-stars) ─────────────────────────────
    tryCatch({
      mc_da <- da_results\$metacyc
      sig_features <- mc_da\$kw_df\$feature_id[mc_da\$kw_df\$signif_q05 == "YES"]
      if (length(sig_features) >= 2) {
        # Limit to top 40 by significance for readability
        sig_features <- head(sig_features, 40)
        lfc_mat <- mc_da\$pair_lfc[sig_features, , drop = FALSE]
        q_mat   <- mc_da\$pair_p[sig_features,   , drop = FALSE]
        lfc_mat[is.na(lfc_mat) | !is.finite(lfc_mat)] <- 0

        # Star annotation
        star_mat <- matrix("", nrow = nrow(q_mat), ncol = ncol(q_mat),
                           dimnames = dimnames(q_mat))
        star_mat[!is.na(q_mat) & q_mat < 0.001] <- "***"
        star_mat[!is.na(q_mat) & q_mat < 0.01  & q_mat >= 0.001] <- "**"
        star_mat[!is.na(q_mat) & q_mat < 0.05  & q_mat >= 0.01]  <- "*"

        # Truncate row names (descriptions)
        descs <- databases\$metacyc\$descs[match(sig_features, databases\$metacyc\$ids)]
        descs[is.na(descs) | descs == ""] <- sig_features[is.na(descs) | descs == ""]
        descs_trunc <- ifelse(nchar(descs) > 60,
                              paste0(substr(descs, 1, 57), "..."),
                              descs)
        rownames(lfc_mat) <- descs_trunc

        # Color: blue → white → red, centred at 0
        max_abs <- max(abs(lfc_mat), na.rm = TRUE)
        max_abs <- max(max_abs, 1)
        breaks  <- seq(-max_abs, max_abs, length.out = 101)
        colours <- colorRampPalette(c("#3b73af", "white", "#cc1f1a"))(100)

        base_path <- file.path(var_dir, sprintf("heatmap_DA_pathways_%s", safe_var))
        png_file  <- paste0(base_path, ".png")
        pdf_file  <- paste0(base_path, ".pdf")
        # Adjust size for # rows × # cols
        w <- max(8, 1.5 + 1.0 * ncol(lfc_mat))
        h <- max(6, 0.3 * nrow(lfc_mat) + 2)

        for (out_file in c(png_file, pdf_file)) {
          pheatmap(lfc_mat,
                   color             = colours,
                   breaks            = breaks,
                   cluster_rows      = TRUE,
                   cluster_cols      = FALSE,
                   show_rownames     = TRUE,
                   show_colnames     = TRUE,
                   display_numbers   = star_mat,
                   number_color      = "black",
                   fontsize_number   = 10,
                   fontsize_row      = 7,
                   fontsize_col      = 8,
                   main              = sprintf("DA MetaCyc pathways  --  %s\n%s",
                                              var, PICRUST_SUBTITLE),
                   filename          = out_file,
                   width             = w,
                   height            = h)
        }
        cat(sprintf("[OK]    DA heatmap   → %s.{png,pdf}  (%d sig pathways)\n",
                    basename(base_path), length(sig_features)))
      } else {
        cat(sprintf("[INFO]  <2 significant pathways for %s — heatmap skipped.\n", var))
      }
    }, error = function(e)
      cat(sprintf("[WARN]  DA heatmap failed: %s\n", conditionMessage(e))))

    # ── 3. Functional beta-diversity PCoA (Bray-Curtis on pathways) ───────
    tryCatch({
      mc_abun <- databases\$metacyc\$abun[, samp_keep, drop = FALSE]
      mc_abun[is.na(mc_abun)] <- 0
      # Drop zero-sum features (no info)
      mc_abun <- mc_abun[rowSums(mc_abun) > 0, , drop = FALSE]
      if (ncol(mc_abun) >= 4) {
        bc <- vegan::vegdist(t(mc_abun), method = "bray")
        pcoa <- stats::cmdscale(bc, k = 2, eig = TRUE)
        pts <- as.data.frame(pcoa\$points)
        names(pts) <- c("PC1", "PC2")
        pts\$sample_id <- rownames(pts)
        pts\$Group <- md_var[[var]][match(pts\$sample_id, md_var[["sample-id"]])]
        eig_pos <- pcoa\$eig[pcoa\$eig > 0]
        var_exp <- 100 * eig_pos / sum(eig_pos)
        pc1_lab <- sprintf("PC1 (%.1f%%)", var_exp[1])
        pc2_lab <- sprintf("PC2 (%.1f%%)", ifelse(length(var_exp) >= 2, var_exp[2], 0))

        n_lvl <- length(levels(pts\$Group))
        pal <- if (n_lvl <= 8) brewer.pal(max(3, n_lvl), "Set2")[seq_len(n_lvl)]
               else if (n_lvl <= 12) brewer.pal(n_lvl, "Set3")
               else colorRampPalette(brewer.pal(8, "Set2"))(n_lvl)

        p <- ggplot(pts, aes(x = PC1, y = PC2, colour = Group, fill = Group)) +
          stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.12,
                       linetype = 2, colour = NA) +
          geom_point(size = 3, alpha = 0.85) +
          scale_colour_manual(values = pal) +
          scale_fill_manual(values   = pal) +
          labs(title    = sprintf("Functional beta diversity (Bray-Curtis, MetaCyc)  --  by %s", var),
               subtitle = PICRUST_SUBTITLE,
               x = pc1_lab, y = pc2_lab,
               colour = var, fill = var) +
          theme_classic(base_size = 11) +
          theme(plot.title    = element_text(face = "bold", size = 12),
                plot.subtitle = element_text(colour = "grey30", size = 9),
                legend.position = "right",
                panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

        base_path <- file.path(var_dir, sprintf("PCoA_BrayCurtis_FUNCTIONAL_%s", safe_var))
        save_ggplot(p, base_path, w = 8, h = 6)
        cat(sprintf("[OK]    PCoA         → %s.{png,pdf}\n", basename(base_path)))
      } else {
        cat("[INFO]  <4 samples — PCoA skipped.\n")
      }
    }, error = function(e)
      cat(sprintf("[WARN]  PCoA failed: %s\n", conditionMessage(e))))
  }   # end SKIP_PLOTS
}     # end variable loop

# ───────────────────────────────────────────────────────────────────────────
# Part D — DA genera × DA pathways Spearman correlation
# ───────────────────────────────────────────────────────────────────────────
if (DO_CORR && HAS_ANCOMBC) {
  cat("\n════════════════════════════════════════════════════════════════════\n")
  cat("  Part D: DA genera × DA pathways Spearman correlation\n")
  cat("════════════════════════════════════════════════════════════════════\n")

  for (safe_var in names(da_path_per_var)) {
    info_var <- da_path_per_var[[safe_var]]
    var      <- info_var\$var
    var_dir  <- file.path(PICRUST_DIR, sprintf("picrust2_%s", safe_var))
    sig_paths <- info_var\$sig_pathways

    if (length(sig_paths) < 2) {
      cat(sprintf("[INFO]  %s: <2 DA pathways — correlation skipped.\n", var))
      next
    }

    # Load step-14 ANCOMBC2 genus DA results
    anco_genus_pri <- file.path(ANCOMBC_DIR, "ANCOMBC2_genus",
                                sprintf("ancombc2_%s", safe_var),
                                sprintf("ancombc2_primary_%s.xlsx", safe_var))
    anco_genus_pw <- file.path(ANCOMBC_DIR, "ANCOMBC2_genus",
                               sprintf("ancombc2_%s", safe_var),
                               sprintf("ancombc2_pairwise_%s.xlsx", safe_var))
    anco_genus_gl <- file.path(ANCOMBC_DIR, "ANCOMBC2_genus",
                               sprintf("ancombc2_%s", safe_var),
                               sprintf("ancombc2_global_%s.xlsx", safe_var))

    sig_genera <- character(0)
    if (file.exists(anco_genus_gl)) {
      gl <- tryCatch(read.xlsx(anco_genus_gl), error = function(e) NULL)
      if (!is.null(gl) && "diff_abn" %in% names(gl))
        sig_genera <- c(sig_genera,
                        as.character(gl\$taxon[gl\$diff_abn %in% TRUE]))
    }
    if (file.exists(anco_genus_pw)) {
      pw <- tryCatch(read.xlsx(anco_genus_pw), error = function(e) NULL)
      if (!is.null(pw)) {
        diff_cols <- grep("^diff_", names(pw), value = TRUE)
        if (length(diff_cols) > 0) {
          mask <- rowSums(sapply(diff_cols, function(c) pw[[c]] %in% TRUE),
                          na.rm = TRUE) > 0
          sig_genera <- c(sig_genera, as.character(pw\$taxon[mask]))
        }
      }
    }
    sig_genera <- unique(sig_genera[!is.na(sig_genera) & nchar(sig_genera) > 0])
    if (length(sig_genera) < 2) {
      cat(sprintf("[INFO]  %s: <2 DA genera — correlation skipped.\n", var))
      next
    }

    # Load genus-collapsed feature table from step 14 working dir
    genus_tsv <- file.path(ANCOMBC_DIR, "working_dir_differential_abundance",
                           "exported_tables", "L6", "feature-table.tsv")
    if (!file.exists(genus_tsv)) {
      cat(sprintf("[INFO]  %s: genus TSV not found at %s — correlation skipped.\n",
                  var, genus_tsv))
      next
    }
    g_lines <- readLines(genus_tsv)
    g_hdr_idx <- grep("^#OTU", g_lines, ignore.case = FALSE)
    if (length(g_hdr_idx) == 0) g_hdr_idx <- which(!grepl("^#", g_lines))[1] - 1
    g_hdr <- strsplit(sub("^#", "", g_lines[g_hdr_idx]), "\t", fixed=TRUE)[[1]]
    g_body <- g_lines[-seq_len(g_hdr_idx)]
    g_body <- g_body[nchar(g_body) > 0 & !grepl("^#", g_body)]
    g_mat  <- do.call(rbind, lapply(strsplit(g_body, "\t", fixed=TRUE), function(x) {
      if (length(x) < length(g_hdr)) x <- c(x, rep("", length(g_hdr) - length(x)))
      x[seq_along(g_hdr)]
    }))
    colnames(g_mat) <- g_hdr
    g_df <- as.data.frame(g_mat, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(g_df) <- g_df[[1]]; g_df[[1]] <- NULL
    for (j in seq_len(ncol(g_df))) g_df[[j]] <- as.numeric(g_df[[j]])
    g_mat <- as.matrix(g_df); g_mat[is.na(g_mat)] <- 0

    # Restrict to common samples
    md_var <- info_var\$md_var
    samp_keep <- intersect(md_var[["sample-id"]],
                  intersect(colnames(g_mat), colnames(info_var\$abun)))
    if (length(samp_keep) < 4) {
      cat(sprintf("[INFO]  %s: <4 common samples — correlation skipped.\n", var))
      next
    }
    sig_genera_in   <- intersect(sig_genera, rownames(g_mat))
    sig_paths_in    <- intersect(sig_paths,  rownames(info_var\$abun))
    if (length(sig_genera_in) < 2 || length(sig_paths_in) < 2) {
      cat(sprintf("[INFO]  %s: <2 matched genera/pathways — correlation skipped.\n", var))
      next
    }

    g_sub <- g_mat[sig_genera_in, samp_keep, drop = FALSE]
    p_sub <- info_var\$abun[sig_paths_in, samp_keep, drop = FALSE]

    # Spearman corr (genera rows × pathways cols)
    cor_mat <- suppressWarnings(
      cor(t(g_sub), t(p_sub), method = "spearman"))

    # p-values via cor.test
    p_mat <- matrix(NA_real_, nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                    dimnames = dimnames(cor_mat))
    for (i in seq_len(nrow(cor_mat))) {
      for (j in seq_len(ncol(cor_mat))) {
        ct <- tryCatch(suppressWarnings(
                cor.test(g_sub[i, ], p_sub[j, ], method = "spearman")),
              error = function(e) NULL)
        if (!is.null(ct)) p_mat[i, j] <- ct\$p.value
      }
    }
    q_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                    nrow = nrow(p_mat), dimnames = dimnames(p_mat))

    # Save TSV with genera in rows, pathways in cols
    desc_map <- setNames(databases\$metacyc\$descs, databases\$metacyc\$ids)
    cor_short <- cor_mat
    colnames(cor_short) <- ifelse(!is.na(desc_map[colnames(cor_mat)]) &
                                  desc_map[colnames(cor_mat)] != "",
                                  desc_map[colnames(cor_mat)],
                                  colnames(cor_mat))
    cor_tsv  <- file.path(var_dir,
                          sprintf("correlation_DAgenera_x_DApathways_%s.tsv", safe_var))
    write.table(cbind(genus = rownames(cor_short), cor_short),
                cor_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("[OK]    correlation TSV → %s\n", basename(cor_tsv)))

    # Heatmap with q-stars
    if (!SKIP_PLOTS) {
      tryCatch({
        # Truncate names for display
        g_short <- ifelse(nchar(rownames(cor_mat)) > 50,
                          paste0(substr(rownames(cor_mat), 1, 47), "..."),
                          rownames(cor_mat))
        p_short <- ifelse(nchar(colnames(cor_short)) > 50,
                          paste0(substr(colnames(cor_short), 1, 47), "..."),
                          colnames(cor_short))
        rownames(cor_mat) <- g_short
        colnames(cor_mat) <- p_short
        rownames(p_mat)   <- g_short
        colnames(p_mat)   <- p_short
        rownames(q_mat)   <- g_short
        colnames(q_mat)   <- p_short

        star_mat <- matrix("", nrow = nrow(q_mat), ncol = ncol(q_mat),
                           dimnames = dimnames(q_mat))
        star_mat[!is.na(q_mat) & q_mat < 0.001] <- "***"
        star_mat[!is.na(q_mat) & q_mat < 0.01  & q_mat >= 0.001] <- "**"
        star_mat[!is.na(q_mat) & q_mat < 0.05  & q_mat >= 0.01]  <- "*"

        max_abs <- max(abs(cor_mat), na.rm = TRUE)
        breaks  <- seq(-max_abs, max_abs, length.out = 101)
        colours <- colorRampPalette(c("#3b73af", "white", "#cc1f1a"))(100)

        base_path <- file.path(var_dir,
                               sprintf("correlation_DAgenera_x_DApathways_%s", safe_var))
        w <- max(8, 0.4 * ncol(cor_mat) + 4)
        h <- max(6, 0.3 * nrow(cor_mat) + 2)
        for (out_file in c(paste0(base_path, ".png"), paste0(base_path, ".pdf"))) {
          pheatmap(cor_mat,
                   color             = colours,
                   breaks            = breaks,
                   cluster_rows      = TRUE,
                   cluster_cols      = TRUE,
                   show_rownames     = TRUE,
                   show_colnames     = TRUE,
                   display_numbers   = star_mat,
                   number_color      = "black",
                   fontsize_number   = 9,
                   fontsize_row      = 7,
                   fontsize_col      = 7,
                   main              = sprintf("Spearman: DA genera x DA MetaCyc pathways  --  %s\n%s",
                                              var, PICRUST_SUBTITLE),
                   filename          = out_file,
                   width             = w,
                   height            = h)
        }
        cat(sprintf("[OK]    correlation heatmap → %s.{png,pdf}\n", basename(base_path)))
      }, error = function(e)
        cat(sprintf("[WARN]  Correlation heatmap failed: %s\n", conditionMessage(e))))
    }
  }
} else {
  cat("\n[INFO]  Skipping DA genera × DA pathway correlation (step 14 outputs not found or --skip-correlation).\n")
}

cat("\n[DONE]  Analysis complete.\n")
RANALYSIS

if $DRY_RUN; then
  warn "[DRY-RUN] Would run R analysis."
else
  _R --vanilla "$R_ANALYSIS" \
    || err "R analysis failed — see output above."
fi
rm -f "$R_ANALYSIS"
sep

# =============================================================================
step "12/12 — Generate self-contained HTML report + write info file"
# =============================================================================
if $SKIP_HTML; then
  warn "Skipping HTML report (--skip-html)."
else
  R_HTML="/tmp/mbx_picrust_${_TMPID}_html.R"

  # COG status flag for the report
  if [[ -f "${COG_PRED}.skipped" ]]; then COG_SKIPPED_R="TRUE"; else COG_SKIPPED_R="FALSE"; fi
  PICRUST2_VER_FOR_HTML="$(_PICRUST place_seqs.py -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 | awk '{print $2}' || echo 'unknown')"

  cat > "$R_HTML" << RHTML
suppressPackageStartupMessages({
  library(openxlsx)
})

PICRUST_DIR <- "${PICRUST_DIR}"
HTML_REPORT <- "${HTML_REPORT}"
NSTI_XLSX   <- "${NSTI_XLSX}"
NSTI_LOG    <- "${NSTI_DIR}/nsti_filtering_summary.txt"
NSTI_THRESH <- ${NSTI_THRESHOLD}
VARS        <- c(${VARS_R})
COG_SKIPPED <- ${COG_SKIPPED_R}
PICRUST2_VER <- "${PICRUST2_VER_FOR_HTML}"

# Helper: file paths relative to the HTML file's location (PICRUST_DIR)
rel <- function(abs_path) {
  if (!file.exists(abs_path)) return(NA_character_)
  sub(paste0(PICRUST_DIR, "/"), "", abs_path, fixed = TRUE)
}

img_or_msg <- function(path, alt) {
  rp <- rel(path)
  if (is.na(rp)) sprintf('<p class="missing">[%s not generated]</p>', alt)
  else sprintf('<img src="%s" alt="%s" loading="lazy" />', rp, alt)
}

# Read NSTI summary
nsti_text <- if (file.exists(NSTI_LOG)) {
  paste(readLines(NSTI_LOG), collapse = "\n")
} else "(NSTI summary not available)"

# Read per-sample NSTI table for embedded display
sample_html <- if (file.exists(NSTI_XLSX)) {
  df <- read.xlsx(NSTI_XLSX)
  rows <- apply(df, 1, function(r) {
    cls <- ""
    if ("flag_50pct_reads_lost" %in% names(df) && r["flag_50pct_reads_lost"] == "YES")
      cls <- "critical"
    else if ("flag_high_NSTI" %in% names(df) && r["flag_high_NSTI"] == "YES")
      cls <- "warning"
    cells <- sapply(r, function(v) sprintf("<td>%s</td>", as.character(v)))
    sprintf('<tr class="%s">%s</tr>', cls, paste(cells, collapse = ""))
  })
  hdr <- paste(sprintf("<th>%s</th>", names(df)), collapse = "")
  sprintf('<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>',
          hdr, paste(rows, collapse = "\n"))
} else "<p>(per-sample NSTI table not available)</p>"

# Per-variable HTML sections
sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^[._-]+|[._-]+\$", "", x)
  x
}

var_html <- list()
for (v in VARS) {
  sv <- sanitize(v)
  vd <- file.path(PICRUST_DIR, sprintf("picrust2_%s", sv))
  bar     <- file.path(vd, sprintf("stacked_bar_top20_metacyc_%s.png", sv))
  heat    <- file.path(vd, sprintf("heatmap_DA_pathways_%s.png", sv))
  pcoa    <- file.path(vd, sprintf("PCoA_BrayCurtis_FUNCTIONAL_%s.png", sv))
  corr    <- file.path(vd, sprintf("correlation_DAgenera_x_DApathways_%s.png", sv))
  da_xlsx <- file.path(vd, sprintf("DA_metacyc_pathways_KW_%s.xlsx", sv))
  # Note: the actual file mtype is "METACYC" (uppercase) per how add_descriptions outputs;
  # try alternate name if first not present
  alt_xlsx <- file.path(vd, sprintf("DA_METACYC_KW_%s.xlsx", sv))
  da_files_html <- ""
  for (xf in list.files(vd, pattern = "^DA_.*[.]xlsx$", full.names = TRUE)) {
    da_files_html <- paste0(da_files_html,
      sprintf('<li><a href="%s">%s</a></li>', rel(xf), basename(xf)))
  }

  var_html[[v]] <- sprintf('
<section class="variable-block">
  <h2>Variable: %s</h2>
  <h3>Top 20 MetaCyc pathways (stacked bar)</h3>
  %s
  <h3>Differentially abundant pathways (heatmap, lfc + q-stars)</h3>
  %s
  <h3>Functional beta diversity PCoA (Bray-Curtis on MetaCyc)</h3>
  %s
  <h3>DA genera x DA pathways (Spearman correlation)</h3>
  %s
  <h3>Differential abundance result tables</h3>
  <ul>%s</ul>
</section>',
    htmltools::htmlEscape(v),
    img_or_msg(bar,  "Stacked bar"),
    img_or_msg(heat, "DA heatmap"),
    img_or_msg(pcoa, "Functional PCoA"),
    img_or_msg(corr, "Genera x pathway correlation"),
    da_files_html)
}

# (Avoid hard dep on htmltools — provide simple escape)
htmltools <- list(htmlEscape = function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
})
var_html_join <- paste(unlist(var_html), collapse = "\n")

now <- format(Sys.time())

# Build "Known limitations" section for the HTML report.
# COG database is the most common gotcha — explicitly call it out when skipped.
limitations_html <- if (COG_SKIPPED) {
  sprintf('<div class="limitation-box">
  <h2>Known limitation: COG database not analysed in this run</h2>
  <p>
    Of the four functional databases this pipeline normally produces
    (<strong>KO, EC, COG, MetaCyc pathways</strong>), the
    <strong>COG (Clusters of Orthologous Groups)</strong> database
    was <strong>skipped</strong> because the installed PICRUSt2 build
    (<code>v%s</code>) has incomplete COG support:
  </p>
  <ul>
    <li><code>COG</code> appears in the <code>hsp.py -i</code> argparse choices, but it is missing from the <code>default_tables_bac</code> dictionary in PICRUSt2&apos;s source code.</li>
    <li>The bundled <code>cog.txt.gz</code> trait table uses <em>old IMG reference IDs</em> that do not match the current default placement tree, so the fallback (<code>--observed_trait_table</code>) also fails.</li>
  </ul>
  <p>
    <strong>What this means for you:</strong> KO, EC, and MetaCyc results in this report are
    complete and unaffected.  COG-level functional categorisation is unavailable for this run.
    If COG is critical to your analysis, upgrade PICRUSt2 (e.g. when v2.7+ restores COG)
    and rerun this step — the script will pick it up automatically.
  </p>
</div>', PICRUST2_VER)
} else { "" }

html <- sprintf('<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>PICRUSt2 functional prediction report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, sans-serif;
         max-width: 1100px; margin: 2em auto; padding: 0 1.5em;
         color: #222; line-height: 1.55; }
  h1 { border-bottom: 3px solid #3b73af; padding-bottom: 0.3em; }
  h2 { border-bottom: 1px solid #ccc; padding-bottom: 0.2em; margin-top: 2em; }
  h3 { color: #3b73af; margin-top: 1.5em; }
  pre { background: #f5f5f5; padding: 1em; overflow-x: auto;
        font-size: 0.85em; border-left: 4px solid #3b73af; }
  table { border-collapse: collapse; width: 100%%; margin: 1em 0;
          font-size: 0.85em; }
  th, td { border: 1px solid #ddd; padding: 0.4em 0.6em; text-align: left; }
  th { background: #f0f4f8; font-weight: 600; }
  tr.warning td { background: #fff3cd; }
  tr.critical td { background: #f8d7da; font-weight: 600; }
  img { max-width: 100%%; height: auto; margin: 1em 0;
        border: 1px solid #ddd; padding: 4px; background: white;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
  .subtitle { color: #666; font-style: italic; }
  .badge { display: inline-block; padding: 0.2em 0.6em;
           background: #3b73af; color: white; border-radius: 4px;
           font-size: 0.8em; margin-right: 0.5em; }
  .missing { color: #999; font-style: italic; }
  ul { line-height: 1.8; }
  ul a { text-decoration: none; color: #3b73af; }
  ul a:hover { text-decoration: underline; }
  .info-box { background: #e7f1f8; border-left: 4px solid #3b73af;
              padding: 1em; margin: 1em 0; border-radius: 4px; }
  .limitation-box { background: #fff3cd; border-left: 4px solid #d39e00;
                    padding: 1em 1.2em; margin: 1.5em 0; border-radius: 4px; }
  .limitation-box h2 { border: none; margin-top: 0; color: #856404;
                       font-size: 1.1em; padding-bottom: 0.2em; }
  .limitation-box code { background: #f8f1d6; padding: 0.1em 0.3em;
                         border-radius: 3px; font-size: 0.9em; }
  footer { margin-top: 3em; padding-top: 1em; border-top: 1px solid #ccc;
           color: #666; font-size: 0.85em; }
</style>
</head>
<body>

<h1>PICRUSt2 functional prediction report</h1>
<p class="subtitle">Predicted functional profiles - inferred from 16S rRNA data</p>
<p>
  <span class="badge">Step 15</span>
  <span class="badge">NSTI threshold: %.2f</span>
  <span class="badge">Generated: %s</span>
</p>

<div class="info-box">
  <strong>About this report.</strong>
  All abundance values are <em>predicted</em> from 16S rRNA gene placement
  on a reference tree.  These predictions assume that the closest reference
  genome on the tree is biochemically representative of your sample.  ASVs
  with high NSTI (Nearest Sequenced Taxon Index) values are unreliable and
  have been filtered out at the threshold shown above.
</div>

%s

<h2>NSTI filtering summary</h2>
<pre>%s</pre>

<h2>Per-sample NSTI scores</h2>
<p>
  Samples with weighted mean NSTI &gt; 1.0 are flagged in
  <span style="background:#fff3cd;padding:0 4px;">amber</span>:
  functional predictions for those samples are unreliable.  Samples that
  lost &gt;50%% of reads to NSTI filtering are flagged in
  <span style="background:#f8d7da;padding:0 4px;font-weight:600;">red</span>
  and should be considered for exclusion from downstream DA analysis.
</p>
%s

%s

<footer>
  Generated by mbx_picrust_run.sh on %s.
  PICRUSt2 reference: Douglas et al. 2020. Nature Biotechnology.
  Differential abundance: Kruskal-Wallis + Dunn (k&ge;3) or Wilcoxon (k=2),
  Benjamini-Hochberg adjusted.
</footer>

</body>
</html>',
  NSTI_THRESH, now,
  limitations_html,
  htmltools\$htmlEscape(nsti_text),
  sample_html,
  var_html_join,
  now)

writeLines(html, HTML_REPORT)
cat(sprintf("[OK]    HTML report → %s\n", HTML_REPORT))
RHTML

  if $DRY_RUN; then
    warn "[DRY-RUN] Would render HTML report."
  else
    _R --vanilla "$R_HTML" \
      || warn "HTML report rendering failed (results still complete)."
  fi
  rm -f "$R_HTML"
fi

# ── Write info file ───────────────────────────────────────────────────────────
QIIME_VER="$(qiime --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || echo unknown)"
PICRUST2_VER="$(_PICRUST place_seqs.py -h 2>&1 | grep -oE 'version [0-9.]+' | head -1 | awk '{print $2}' || echo unknown)"

VAR_LINES=""
for col in "${CATEGORICAL_COLS[@]}"; do
  VAR_LINES="${VAR_LINES}VARIABLE=${col}
"
done

cat > "$INFO_FILE" << INFO
# ============================================================================
# mbx_picrust2_info.txt
# Generated by mbx_picrust_run.sh
# Date : $NOW_PRETTY
# ============================================================================
#
# Do NOT edit the key=value lines below — they are parsed programmatically.

# ── Inputs ───────────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
BETA_INFO=$BETA_INFO
PRE_DIVERSITY_INFO=$PRE_DIV_INFO
ANCOMBC_INFO=$ANCO_INFO
METADATA_TXT=$METADATA_TXT
REP_SEQS_QZA=$REP_SEQS_QZA
FEATURE_TABLE_FILTERED_QZA=$FT_FILTERED_QZA

# ── Output locations ─────────────────────────────────────────────────────────
PICRUST2_OUT_DIR=$PICRUST_DIR
ALL_OUTPUTS_DIR=$ALL_DIR
WORK_DIR=$WORK_DIR
KO_OUT_DIR=$KO_OUT_DIR
EC_OUT_DIR=$EC_OUT_DIR
COG_OUT_DIR=$COG_OUT_DIR
PATHWAYS_OUT_DIR=$PATH_OUT_DIR
NSTI_DIR=$NSTI_DIR
RDS_DIR=$RDS_DIR
NSTI_SUMMARY_TXT=${NSTI_DIR}/nsti_filtering_summary.txt
NSTI_SAMPLE_XLSX=$NSTI_XLSX
HTML_REPORT=$HTML_REPORT

# ── Settings ─────────────────────────────────────────────────────────────────
NSTI_THRESHOLD=$NSTI_THRESHOLD
DA_METHOD=Kruskal-Wallis_BH_with_Dunn_or_Wilcoxon
PICRUST_CONDA_ENV=$PICRUST_ENV
CORRELATION_PERFORMED=$DO_CORRELATION

# ── Categorical variables analysed ───────────────────────────────────────────
N_VARIABLES=${#CATEGORICAL_COLS[@]}
${VAR_LINES}

# ── Provenance ───────────────────────────────────────────────────────────────
GENERATED_AT=$NOW_PRETTY
SCRIPT_NAME=mbx_picrust_run.sh
PLATFORM=$PLATFORM_LABEL
QIIME2_VERSION=$QIIME_VER
PICRUST2_VERSION=$PICRUST2_VER
R_VERSION=$R_VERSION
N_JOBS=$N_JOBS
SKIPPED_DEPS_CHECK=$SKIP_DEPS_CHECK
SKIPPED_PLOTS=$SKIP_PLOTS
SKIPPED_HTML=$SKIP_HTML
INVOCATION_USER=${USER:-unknown}
INVOCATION_CWD=$(pwd)
INVOCATION_ARGV=$0 $*
INFO

ok "Info file → $INFO_FILE"

# Cleanup empty COG dir if COG was skipped (cosmetic)
if [[ -f "${COG_PRED}.skipped" ]] && [[ -d "$COG_OUT_DIR" ]] \
   && [[ -z "$(ls -A "$COG_OUT_DIR" 2>/dev/null)" ]]; then
  rmdir "$COG_OUT_DIR" 2>/dev/null || true
fi

sep

# =============================================================================
step "Done"
# =============================================================================
echo ""
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║  STEP 15 COMPLETE  —  PICRUSt2 FUNCTIONAL PREDICTIONS       ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo ""
  if [[ -f "${COG_PRED}.skipped" ]]; then
    echo "  PICRUSt2 outputs   : $ALL_DIR/{KO,EC,pathways_metacyc}_metagenome/"
    echo "                       (COG skipped — see warning above; not supported in this PICRUSt2 build)"
  else
    echo "  PICRUSt2 outputs   : $ALL_DIR/{KO,EC,COG,pathways_metacyc}_metagenome/"
  fi
  echo "  NSTI summary       : $NSTI_DIR/nsti_filtering_summary.txt"
echo "  Per-sample NSTI    : $NSTI_XLSX"
echo "  Per-variable dirs  : $PICRUST_DIR/picrust2_<Variable>/"
echo "  Variables analysed : ${CATEGORICAL_COLS[*]}"
echo "  HTML report        : $HTML_REPORT"
echo "  Open it:"
echo "    open '$HTML_REPORT'"
echo ""
echo "  Each picrust2_<Var>/ folder contains:"
echo "    • DA_KO_KW_<Var>.xlsx + DA_KO_pairwise_<Var>.xlsx"
echo "    • DA_EC_KW_<Var>.xlsx + DA_EC_pairwise_<Var>.xlsx"
if [[ ! -f "${COG_PRED}.skipped" ]]; then
echo "    • DA_COG_KW_<Var>.xlsx + DA_COG_pairwise_<Var>.xlsx"
fi
echo "    • DA_METACYC_KW_<Var>.xlsx + DA_METACYC_pairwise_<Var>.xlsx"
if ! $SKIP_PLOTS; then
echo "    • stacked_bar_top20_metacyc_<Var>.{png,pdf}"
echo "    • heatmap_DA_pathways_<Var>.{png,pdf}     (lfc + q-stars)"
echo "    • PCoA_BrayCurtis_FUNCTIONAL_<Var>.{png,pdf}"
if $DO_CORRELATION; then
echo "    • correlation_DAgenera_x_DApathways_<Var>.{tsv,png,pdf}"
fi
fi
echo ""
echo "  Provenance file    : $INFO_FILE"
echo ""
