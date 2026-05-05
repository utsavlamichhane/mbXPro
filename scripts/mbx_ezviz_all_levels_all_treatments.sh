#!/usr/bin/env bash
# =============================================================================
#  mbx_ezviz_all_levels_all_treatments.sh
#  Run mbX::ezviz() for all 7 taxonomic levels × all categorical metadata columns
#
#  Compatible with bash 3.2+ (macOS default shell)
#  Requires: Rscript + mbX (both installed by mbx_ezclean_all_levels.sh)
#
#  STEPS:
#    1  Parse 7_taxonomy_csv/mbx_taxonomy_info.txt → level-7.csv + metadata
#    2  Detect categorical metadata columns (R-based, QIIME2-aware)
#    3  Create 9_visualization_entire/<variable>/ for each categorical column
#    4  For each variable × each level:
#         setwd(9_viz/<variable>/)
#         plot <- ezviz(level-7.csv, metadata, level, variable, threshold=0.5)
#         ggsave("ezviz_<level>_<variable>.png")
#    5  Write summary report mbx_ezviz_info.txt
#
#  CATEGORICAL DETECTION RULES:
#    ✓ Include : columns where values are non-numeric strings
#    ✗ Exclude : first column (sample-id)
#    ✗ Exclude : #q2:types directive row (QIIME2 metadata)
#    ✗ Exclude : all-unique columns (free text / IDs) 
#    ✗ Exclude : single-value columns (no variation to plot)
#    ✗ Exclude : columns with < 2 samples per group (not plottable)
#
#  OUTPUT STRUCTURE:
#    mbX_pro_outputs_<timestamp>/
#    └── 9_visualization_entire/
#        ├── <Variable1>/
#        │   ├── ezviz_domain_<Variable1>.png
#        │   ├── ezviz_phylum_<Variable1>.png
#        │   ├── ezviz_class_<Variable1>.png
#        │   ├── ezviz_order_<Variable1>.png
#        │   ├── ezviz_family_<Variable1>.png
#        │   ├── ezviz_genus_<Variable1>.png
#        │   └── ezviz_species_<Variable1>.png
#        ├── <Variable2>/
#        │   └── ...
#        └── mbx_ezviz_info.txt
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

# Sanitize a string for safe use as a directory name:
# replace spaces/special chars with underscores, strip leading/trailing
_sanitize_dirname() {
  printf '%s' "$1" \
    | tr ' ' '_' \
    | tr -d '()[]{}/<>|\\:*?"' \
    | sed 's/^[_.-]*//;s/[_.-]*$//'
}

# ── R env wrapper — strip conda's R_LIBS_USER pollution ──────────────────────
# When run from inside a conda env (e.g. qiime2-amplicon), conda exports
# R_LIBS_USER / R_HOME pointing at the conda env's R library tree.  System
# Rscript (Homebrew/CRAN) then tries to load packages from those paths and
# fails with errors like:
#     Error: shared object 'Rcpp.so' not found
#     Error in library.dynam: ... 'methods.dylib' not found
# Unset at script start AND strip again at every system Rscript invocation.
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}

# Defined after RSCRIPT_CMD is found
_R() { _strip_env "$RSCRIPT_CMD" "$@"; }

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_ezviz_all_levels_all_treatments.sh — Run ezviz() for all levels × all categorical variables

USAGE:
  mbx_ezviz_all_levels_all_treatments.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads level-7.csv and metadata.txt from:
    <mbX_pro_outputs_dir>/7_taxonomy_csv/mbx_taxonomy_info.txt

  Auto-detects all categorical columns in the metadata, then runs ezviz()
  for every combination of taxonomic level (d, p, c, o, f, g, s) and
  categorical variable. Each plot is saved as a 300 DPI PNG.

  level-7.csv is used for all 7 levels — it contains the full taxonomy string
  that ezviz() parses to extract each level's taxa.

RUNTIME ESTIMATE:
  Each ezviz() call re-runs ezclean() internally (~30 sec–2 min per call).
  Total = 7 levels × N categorical variables.
  Example: 4 categorical variables = 28 runs ≈ 15–60 minutes.

OPTIONS:
  --threshold <val>    Minimum % abundance threshold (default: 0.5)
                       Taxa below this are grouped into "Other"
  --top-taxa <N>       Use top N taxa instead of threshold (mutually exclusive)
  --width <inches>     Plot width in inches (default: 12)
  --height <inches>    Plot height in inches (default: 7)
  --dpi <N>            Plot resolution (default: 300)
  --dry-run            Print R commands without executing
  -h, --help           Show this help and exit

EXAMPLES:
  mbx_ezviz_all_levels_all_treatments.sh /path/to/mbX_pro_outputs_20250422_143022

  # Use top_taxa instead of threshold
  mbx_ezviz_all_levels_all_treatments.sh /path/to/mbX_pro_outputs_20250422_143022 \
    --top-taxa 20

  # Higher resolution plots
  mbx_ezviz_all_levels_all_treatments.sh /path/to/mbX_pro_outputs_20250422_143022 \
    --dpi 600 --width 14 --height 8

CATEGORICAL COLUMN DETECTION:
  A metadata column is treated as categorical if:
    • Values are non-numeric strings
    • Has ≥ 2 distinct values (not constant)
    • Is not all-unique (not a free-text / ID column)
    • Has ≥ 2 samples per group (at least one group has 2+ members)
  QIIME2 #q2:types rows are automatically removed before detection.

COMMON ERRORS:
  "The selected metadata is either not in the metadata or not a categorical value"
    → The column was detected as categorical but ezviz disagrees.
      Check the column values for mixed numeric/string content.
  "No categorical columns found"
    → All metadata columns are numeric or all-unique.
      Check your metadata file — at least one grouping column is needed.
  "object 'plot_obj' not found / ezviz returned non-plot"
    → ezviz() returned an error string instead of a plot. Check the
      level-7.csv format and metadata compatibility.
  ezviz fails for species (s) but succeeds for others
    → Normal — short reads often don't classify to species level.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
THRESHOLD="0.5"
TOP_TAXA=""
PLOT_WIDTH="12"
PLOT_HEIGHT="7"
PLOT_DPI="300"
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)      usage ;;
    --dry-run)      DRY_RUN=true; shift ;;
    --threshold)    THRESHOLD="${2:?Missing value for --threshold}"; shift 2 ;;
    --top-taxa)     TOP_TAXA="${2:?Missing value for --top-taxa}"; shift 2 ;;
    --width)        PLOT_WIDTH="${2:?Missing value for --width}"; shift 2 ;;
    --height)       PLOT_HEIGHT="${2:?Missing value for --height}"; shift 2 ;;
    --dpi)          PLOT_DPI="${2:?Missing value for --dpi}"; shift 2 ;;
    -*)  err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '${1}'"; fi ;;
  esac
done

# Validate mutually exclusive options
if [[ -n "$TOP_TAXA" && "$THRESHOLD" != "0.5" ]]; then
  err "--threshold and --top-taxa are mutually exclusive. Use one or the other."
fi

[[ -z "$MBX_OUT_DIR" ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -d "$MBX_OUT_DIR" ]] || err "Directory does not exist: '${MBX_OUT_DIR}'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# ── Locate Rscript ─────────────────────────────────────────────────────────────
# Prefer absolute SYSTEM paths first; fall back to PATH lookup last.
# Rationale: when run inside a conda env, `command -v Rscript` may return
# conda's R, which has its own (incompatible) library tree. The _strip_env
# wrapper still helps but choosing system R first is more reliable.
RSCRIPT_CMD=""
for _c in "/opt/homebrew/bin/Rscript" \
          "/usr/local/bin/Rscript" \
          "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
          "/usr/bin/Rscript" \
          "$(command -v Rscript 2>/dev/null || true)"; do
  [[ -n "$_c" && -x "$_c" ]] && { RSCRIPT_CMD="$_c"; break; }
done
[[ -n "$RSCRIPT_CMD" ]] || err "Rscript not found.
  → Run mbx_ezclean_all_levels.sh first to install R, or: brew install r"

# ─────────────────────────────────────────────────────────────────────────────
step "1/5 — Parse mbx_taxonomy_info.txt"
# ─────────────────────────────────────────────────────────────────────────────

TAXONOMY_INFO="${MBX_OUT_DIR}/7_taxonomy_csv/mbx_taxonomy_info.txt"
[[ -f "$TAXONOMY_INFO" ]] || err "mbx_taxonomy_info.txt not found:
  $TAXONOMY_INFO
  → Run mbx_taxonomy_run.sh first."

LEVEL7_CSV="$(_read_key  "LEVEL_7_CSV"   "$TAXONOMY_INFO")"
METADATA_TXT="$(_read_key "METADATA_TXT" "$TAXONOMY_INFO")"

[[ -z "$LEVEL7_CSV"   ]] && err "LEVEL_7_CSV not found in mbx_taxonomy_info.txt"
[[ -z "$METADATA_TXT" ]] && err "METADATA_TXT not found in mbx_taxonomy_info.txt"
[[ -f "$LEVEL7_CSV"   ]] || err "level-7.csv not found: $LEVEL7_CSV"
[[ -f "$METADATA_TXT" ]] || err "Metadata file not found: $METADATA_TXT"

info "level-7.csv  : $LEVEL7_CSV"
info "metadata     : $METADATA_TXT"
[[ -n "$TOP_TAXA"   ]] && info "Mode         : top_taxa = $TOP_TAXA" \
                       || info "Mode         : threshold = $THRESHOLD"
info "Plot size    : ${PLOT_WIDTH}in × ${PLOT_HEIGHT}in  @ ${PLOT_DPI} DPI"
$DRY_RUN && warn "DRY-RUN — R code will be printed but NOT executed."
sep

# ─────────────────────────────────────────────────────────────────────────────
step "2/5 — Detect categorical metadata columns"
# ─────────────────────────────────────────────────────────────────────────────

# Use PID + timestamp for unique temp filenames — avoids mktemp suffix issues on macOS
_TMPID="${$}_$(date +%s)"
CATS_FILE="/tmp/mbx_cats_${_TMPID}.txt"
trap 'rm -f "/tmp/mbx_cats_${_TMPID}.txt" "/tmp/mbx_detect_${_TMPID}.R" /tmp/mbx_ezviz_${_TMPID}_*.R' EXIT

CAT_DETECT_R="/tmp/mbx_detect_${_TMPID}.R"
cat > "$CAT_DETECT_R" << RDETECT
# Categorical column detection
# Rules:
#   Include: non-numeric string columns with 2..n-1 unique values
#            and at least one group with >= 2 members
#   Exclude: sample-id (col 1), #q2:types row, all-unique, single-value

metadata_path <- "${METADATA_TXT}"

# Read metadata — handle both tab-separated and CSV
ext <- tolower(tools::file_ext(metadata_path))
if (ext == "csv") {
  meta <- read.csv(metadata_path, header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
} else {
  meta <- read.delim(metadata_path, header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE,
                     na.strings = c("", "NA", "N/A", "na", "n/a"))
}

# Remove QIIME2 #q2:types directive row if present (row 2 in raw file)
q2_rows <- grepl("^#", meta[[1]])
if (any(q2_rows)) {
  meta <- meta[!q2_rows, , drop = FALSE]
  cat("[INFO]  Removed QIIME2 #q2:types row(s).\n", file = stderr())
}

n_samples <- nrow(meta)
cat(sprintf("[INFO]  Metadata: %d samples, %d columns.\n",
            n_samples, ncol(meta)), file = stderr())

# Skip first column (sample-id)
col_names <- names(meta)[-1]

categorical <- character(0)
skipped_numeric  <- character(0)
skipped_alluniq  <- character(0)
skipped_noval    <- character(0)
skipped_mingroupsize <- character(0)

for (col in col_names) {
  vals_raw <- meta[[col]]
  vals     <- vals_raw[!is.na(vals_raw) & nchar(trimws(as.character(vals_raw))) > 0]

  if (length(vals) == 0) {
    skipped_noval <- c(skipped_noval, col)
    next
  }

  # Check if all non-NA values are numeric
  num_test <- suppressWarnings(as.numeric(vals))
  if (!any(is.na(num_test))) {
    skipped_numeric <- c(skipped_numeric, col)
    next
  }

  n_unique <- length(unique(trimws(as.character(vals))))

  # Skip single-value columns (no variation)
  if (n_unique <= 1) {
    skipped_noval <- c(skipped_noval, col)
    next
  }

  # Skip all-unique columns (likely IDs or free text)
  if (n_unique == n_samples) {
    skipped_alluniq <- c(skipped_alluniq, col)
    next
  }

  # Skip if any group has fewer than 2 samples (ezviz can't plot singletons)
  group_counts <- table(trimws(as.character(vals)))
  if (max(group_counts) < 2) {
    skipped_mingroupsize <- c(skipped_mingroupsize, col)
    next
  }

  categorical <- c(categorical, col)
}

# Report
if (length(skipped_numeric)       > 0) cat(sprintf("[INFO]  Skipped numeric columns      : %s\n", paste(skipped_numeric, collapse=", ")), file=stderr())
if (length(skipped_alluniq)       > 0) cat(sprintf("[INFO]  Skipped all-unique columns   : %s\n", paste(skipped_alluniq, collapse=", ")), file=stderr())
if (length(skipped_noval)         > 0) cat(sprintf("[INFO]  Skipped empty/constant cols  : %s\n", paste(skipped_noval, collapse=", ")), file=stderr())
if (length(skipped_mingroupsize)  > 0) cat(sprintf("[INFO]  Skipped singleton-group cols : %s\n", paste(skipped_mingroupsize, collapse=", ")), file=stderr())

if (length(categorical) == 0) {
  cat("[ERROR] No categorical columns found in metadata.\n", file=stderr())
  quit(status=1)
}

cat(sprintf("[INFO]  Categorical columns found    : %d → %s\n",
            length(categorical), paste(categorical, collapse=", ")), file=stderr())

# Write one column name per line to stdout (read by bash)
cat(paste(categorical, collapse="\n"), "\n", sep="")
RDETECT

if $DRY_RUN; then
  warn "[DRY-RUN] Would detect categorical columns from: $METADATA_TXT"
  CATEGORICAL_COLS=("Treatment" "SampleType")   # placeholder for dry run
  info "Dry-run placeholder columns: ${CATEGORICAL_COLS[*]}"
else
  _R --vanilla "$CAT_DETECT_R" > "$CATS_FILE" 2>&1 \
    || err "Categorical column detection failed.
  → Check your metadata file: $METADATA_TXT
  → Run manually: $RSCRIPT_CMD $CAT_DETECT_R"

  # Print the detection output (contains [INFO] lines from R stderr)
  grep "^\[" "$CATS_FILE" || true

  # Read the column names (non-bracket lines = actual column names)
  # bash 3.2-safe replacement for mapfile
  CATEGORICAL_COLS=()
  while IFS= read -r _col; do
    [[ -n "$_col" ]] && CATEGORICAL_COLS+=("$_col")
  done < <(grep -v "^\[" "$CATS_FILE" | grep -v '^$' || true)

  [[ ${#CATEGORICAL_COLS[@]} -eq 0 ]] && err "No categorical columns detected.
  → Verify metadata has at least one grouping column with string values.
  → Common issue: all columns are numeric or QIIME2 #q2:types was not removed."
fi
rm -f "$CAT_DETECT_R"

ok "Categorical columns detected: ${CATEGORICAL_COLS[*]}"
sep

# ─────────────────────────────────────────────────────────────────────────────
step "3/5 — Create directory structure"
# ─────────────────────────────────────────────────────────────────────────────

VIZ_DIR="${MBX_OUT_DIR}/9_visualization_entire"
mkdir -p "$VIZ_DIR"

# Pre-create all variable subdirectories; also build sanitized-name map
declare -a SANITIZED_NAMES=()
for col in "${CATEGORICAL_COLS[@]}"; do
  SAFE="$(_sanitize_dirname "$col")"
  SANITIZED_NAMES+=("$SAFE")
  mkdir -p "${VIZ_DIR}/${SAFE}"
  info "  Created: 9_visualization_entire/${SAFE}/"
done

# Compute total runs for progress tracking
N_VARS="${#CATEGORICAL_COLS[@]}"
N_LEVELS=7
TOTAL_RUNS=$(( N_VARS * N_LEVELS ))
info "Total ezviz() calls planned: ${TOTAL_RUNS}  (${N_LEVELS} levels × ${N_VARS} variables)"
info "Estimated runtime: ~$(( TOTAL_RUNS / 2 ))–$(( TOTAL_RUNS * 2 )) minutes (depends on data size)"
sep

# ─────────────────────────────────────────────────────────────────────────────
step "4/5 — Run ezviz() for all variables × all levels"
# ─────────────────────────────────────────────────────────────────────────────

LEVELS=("d" "p" "c" "o" "f" "g" "s")
LEVEL_NAMES=("domain" "phylum" "class" "order" "family" "genus" "species")

# Build the ezviz param string (threshold OR top_taxa — mutually exclusive)
if [[ -n "$TOP_TAXA" ]]; then
  R_EZVIZ_PARAM="top_taxa = ${TOP_TAXA}"
else
  R_EZVIZ_PARAM="threshold = ${THRESHOLD}"
fi

RUN_NUMBER=0
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
FAIL_LOG=""

for v_idx in "${!CATEGORICAL_COLS[@]}"; do
  COL="${CATEGORICAL_COLS[$v_idx]}"
  SAFE="${SANITIZED_NAMES[$v_idx]}"
  VAR_DIR="${VIZ_DIR}/${SAFE}"

  echo ""
  echo "  ════════════════════════════════════════════════════════════════"
  echo "  Variable: ${COL}  →  ${VAR_DIR}"
  echo "  ════════════════════════════════════════════════════════════════"

  for l_idx in "${!LEVELS[@]}"; do
    LVL="${LEVELS[$l_idx]}"
    LVL_NAME="${LEVEL_NAMES[$l_idx]}"
    RUN_NUMBER=$(( RUN_NUMBER + 1 ))
    OUTPNG="${VAR_DIR}/ezviz_${LVL_NAME}_${SAFE}.png"

    # Skip if output already exists (idempotent re-runs)
    if [[ -f "$OUTPNG" ]]; then
      skipped "[$RUN_NUMBER/$TOTAL_RUNS] level '${LVL}' (${LVL_NAME}) for '${COL}'"
      SKIP_COUNT=$(( SKIP_COUNT + 1 ))
      continue
    fi

    info "[$RUN_NUMBER/$TOTAL_RUNS] ezviz() — level '${LVL}' (${LVL_NAME}) | variable '${COL}'"

    if $DRY_RUN; then
      echo "  [DRY-RUN] setwd(\"${VAR_DIR}\")"
      echo "  [DRY-RUN] ezviz(\"$LEVEL7_CSV\", \"$METADATA_TXT\", \"${LVL}\","
      echo "                   \"${COL}\", ${R_EZVIZ_PARAM})"
      echo "  [DRY-RUN] ggsave(\"ezviz_${LVL_NAME}_${SAFE}.png\", dpi=${PLOT_DPI})"
      PASS_COUNT=$(( PASS_COUNT + 1 ))
      continue
    fi

    # Write per-run R script
    R_SCRIPT="/tmp/mbx_ezviz_${_TMPID}_${RUN_NUMBER}.R"

    cat > "$R_SCRIPT" << RSCRIPT
suppressPackageStartupMessages({
  library(mbX)
  library(ggplot2)
})

# CRITICAL: set working directory to the variable's subdirectory
# so ezclean() intermediate files land here, not in the launch directory
setwd("${VAR_DIR}")

level7_csv   <- "${LEVEL7_CSV}"
metadata_txt <- "${METADATA_TXT}"
col_name     <- "${COL}"
lvl          <- "${LVL}"
lvl_name     <- "${LVL_NAME}"
out_png      <- "ezviz_${LVL_NAME}_${SAFE}.png"

cat(sprintf("[INFO]  cwd         : %s\n", getwd()))
cat(sprintf("[INFO]  level       : %s (%s)\n", lvl, lvl_name))
cat(sprintf("[INFO]  variable    : %s\n", col_name))

result <- tryCatch({

  plot_obj <- ezviz(
    microbiome_data   = level7_csv,
    metadata          = metadata_txt,
    level             = lvl,
    selected_metadata = col_name,
    ${R_EZVIZ_PARAM}
  )

  # ezviz() can return an error string instead of a plot on failure
  if (is.character(plot_obj)) {
    stop(sprintf("ezviz() returned an error message: %s", plot_obj))
  }

  if (!inherits(plot_obj, "gg")) {
    stop(sprintf("ezviz() did not return a ggplot object (got: %s)", class(plot_obj)[1]))
  }

  ggplot2::ggsave(
    filename = out_png,
    plot     = plot_obj,
    dpi      = ${PLOT_DPI},
    width    = ${PLOT_WIDTH},
    height   = ${PLOT_HEIGHT},
    units    = "in"
  )

  cat(sprintf("[OK]    Saved: %s\n", file.path(getwd(), out_png)))
  0L   # success exit code

}, error = function(e) {
  cat(sprintf("[ERROR] ezviz level '%s' for '%s' failed:\n        %s\n",
              lvl, col_name, conditionMessage(e)))
  1L   # failure exit code
})

quit(status = as.integer(result))
RSCRIPT

    timer_start
    if _R --vanilla "$R_SCRIPT" 2>&1; then
      timer_end
      ok "  ✔  level '${LVL}' (${LVL_NAME}) for '${COL}'"
      PASS_COUNT=$(( PASS_COUNT + 1 ))
    else
      timer_end
      warn "  ✘  level '${LVL}' (${LVL_NAME}) for '${COL}' — see output above"
      FAIL_COUNT=$(( FAIL_COUNT + 1 ))
      FAIL_LOG="${FAIL_LOG}  - ${COL} × ${LVL_NAME}\n"
    fi
    rm -f "$R_SCRIPT"

  done   # levels loop
done     # variables loop

# ─────────────────────────────────────────────────────────────────────────────
step "5/5 — Write mbx_ezviz_info.txt + summary"
# ─────────────────────────────────────────────────────────────────────────────

NOW="$(date '+%Y-%m-%d %H:%M:%S')"
INFO_TXT="${VIZ_DIR}/mbx_ezviz_info.txt"

{
  echo "# ============================================================================"
  echo "# mbx_ezviz_info.txt"
  echo "# Generated by mbx_ezviz_all_levels_all_treatments.sh"
  echo "# Date        : $NOW"
  echo "# ============================================================================"
  echo ""
  echo "# ── Inputs ──────────────────────────────────────────────────────────────────"
  echo "MBX_OUTPUTS_DIR=$MBX_OUT_DIR"
  echo "LEVEL_7_CSV=$LEVEL7_CSV"
  echo "METADATA_TXT=$METADATA_TXT"
  echo "EZVIZ_PARAM=${R_EZVIZ_PARAM}"
  echo "PLOT_DPI=${PLOT_DPI}"
  echo "PLOT_WIDTH=${PLOT_WIDTH}"
  echo "PLOT_HEIGHT=${PLOT_HEIGHT}"
  echo ""
  echo "# ── Categorical variables detected ──────────────────────────────────────────"
  for col in "${CATEGORICAL_COLS[@]}"; do
    echo "VARIABLE=${col}"
  done
  echo ""
  echo "# ── Run results ─────────────────────────────────────────────────────────────"
  echo "TOTAL_RUNS=$TOTAL_RUNS"
  echo "PASSED=$PASS_COUNT"
  echo "FAILED=$FAIL_COUNT"
  echo "SKIPPED=$SKIP_COUNT"
  echo ""
  if [[ -n "$FAIL_LOG" ]]; then
    echo "# ── Failed combinations ─────────────────────────────────────────────────────"
    printf '%b' "$FAIL_LOG"
    echo ""
  fi
  echo "# ── Output PNG files ────────────────────────────────────────────────────────"
  for v_idx in "${!CATEGORICAL_COLS[@]}"; do
    SAFE="${SANITIZED_NAMES[$v_idx]}"
    for l_idx in "${!LEVELS[@]}"; do
      LVL_NAME="${LEVEL_NAMES[$l_idx]}"
      OUTPNG="${VIZ_DIR}/${SAFE}/ezviz_${LVL_NAME}_${SAFE}.png"
      if [[ -f "$OUTPNG" ]]; then
        echo "PNG=${OUTPNG}"
      fi
    done
  done
  echo ""
  echo "# ── Using outputs with ezstat() ─────────────────────────────────────────────"
  echo "# Run ezstat() on the same inputs to get statistical comparisons:"
  echo "# library(mbX)"
  for col in "${CATEGORICAL_COLS[@]}"; do
    SAFE="$(_sanitize_dirname "$col")"
    echo "# ezstat(\"$LEVEL7_CSV\", \"$METADATA_TXT\", \"g\", \"${col}\")"
  done
} > "$INFO_TXT"

ok "Info file written → $INFO_TXT"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete — no R code was executed."
else
  ok "ezviz complete!"
  ok "  Passed  : $PASS_COUNT / $TOTAL_RUNS"
  [[ "$SKIP_COUNT" -gt 0 ]] && ok "  Skipped : $SKIP_COUNT (already existed)"
  [[ "$FAIL_COUNT" -gt 0 ]] && warn "  Failed  : $FAIL_COUNT — check output above"
fi
sep
echo ""
echo "  Output structure:"
echo "  $MBX_OUT_DIR/"
echo "  └── 9_visualization_entire/"
for v_idx in "${!CATEGORICAL_COLS[@]}"; do
  COL="${CATEGORICAL_COLS[$v_idx]}"
  SAFE="${SANITIZED_NAMES[$v_idx]}"
  printf '      ├── %s/\n' "$SAFE"
  printf '      │   ├── ezviz_domain_%s.png\n' "$SAFE"
  printf '      │   ├── ezviz_phylum_%s.png\n' "$SAFE"
  printf '      │   ├── ezviz_genus_%s.png\n'  "$SAFE"
  printf '      │   └── ezviz_species_%s.png\n' "$SAFE"
done
echo "      └── mbx_ezviz_info.txt"
echo ""
if [[ "$FAIL_COUNT" -gt 0 ]]; then
  echo "  ── Failed combinations (re-run to retry) ────────────────────────"
  printf '%b' "$FAIL_LOG"
  echo "  Species (s) failures are normal for low-depth sequencing."
  echo "  Re-running is safe — completed plots are skipped automatically."
  echo ""
fi
