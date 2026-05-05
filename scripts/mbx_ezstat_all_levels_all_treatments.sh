#!/usr/bin/env bash
# =============================================================================
#  mbx_ezstat_all_levels_all_treatments.sh
#  Run mbX::ezstat() for all 7 taxonomic levels × all categorical metadata columns
#
#  Compatible with bash 3.2+ (macOS default shell)
#  Requires: Rscript + mbX (installed by mbx_ezclean_all_levels.sh)
#
#  STEPS:
#    1  Parse 7_taxonomy_csv/mbx_taxonomy_info.txt → level-7.csv + metadata
#    2  Detect categorical metadata columns (R-based, QIIME2-aware)
#    3  Create 10_stats/<variable>/ for each categorical column
#    4  For each variable × each level:
#         setwd(10_stats/<variable>/)
#         ezstat(level-7.csv, metadata, level, variable)
#         ezstat handles all file writing internally
#    5  Write summary report mbx_ezstat_info.txt
#
#  WHAT ezstat() PRODUCES (per variable × level):
#    ezstat_KW_<level>_by_<variable>.xlsx       Kruskal-Wallis results
#    ezstat_pairwise_<level>_by_<variable>.xlsx Dunn test pairwise comparisons
#    ezstat_CLD_Summary_<level>_by_<variable>.xlsx Compact Letter Display
#    Boxplots_<variable>/                        PNG per significant taxon (p≤0.05)
#
#  OUTPUT STRUCTURE:
#    mbX_pro_outputs_<timestamp>/
#    └── 10_stats/
#        ├── <Variable1>/
#        │   ├── ezstat_KW_genera_by_<Variable1>.xlsx
#        │   ├── ezstat_pairwise_genera_by_<Variable1>.xlsx
#        │   ├── ezstat_CLD_Summary_genera_by_<Variable1>.xlsx
#        │   └── Boxplots_<Variable1>/
#        │       ├── Taxon1.png
#        │       └── Taxon2.png
#        ├── <Variable2>/
#        │   └── ...
#        └── mbx_ezstat_info.txt
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

mbx_ezstat_all_levels_all_treatments.sh — Run ezstat() for all levels × all categorical variables

USAGE:
  mbx_ezstat_all_levels_all_treatments.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads level-7.csv and metadata.txt from:
    <mbX_pro_outputs_dir>/7_taxonomy_csv/mbx_taxonomy_info.txt

  Auto-detects all categorical columns in the metadata, then runs ezstat()
  for every combination of taxonomic level (d, p, c, o, f, g, s) and
  categorical variable.

  ezstat() performs per-taxon:
    • Kruskal-Wallis test
    • Post-hoc Dunn test with BH (FDR) correction
    • Compact Letter Display (CLD) summary
    • Boxplots for significant taxa (p ≤ 0.05) at 900 DPI

  All file writing is handled internally by ezstat(). No manual ggsave needed.

RUNTIME ESTIMATE:
  Each ezstat() call re-runs ezclean() internally.
  Total = 7 levels × N categorical variables.
  Example: 3 categorical variables = 21 runs ≈ 20–60 minutes.

OPTIONS:
  --dry-run      Print R commands without executing
  -h, --help     Show this help and exit

WHAT ezstat() OUTPUTS (per variable × level, inside 10_stats/<variable>/):
  ezstat_KW_<level>_by_<variable>.xlsx          Kruskal-Wallis p-values per taxon
  ezstat_pairwise_<level>_by_<variable>.xlsx    Dunn test comparisons + FDR p-values
  ezstat_CLD_Summary_<level>_by_<variable>.xlsx CLD letters + median + mean per group
  Boxplots_<variable>/                          PNG boxplots for significant taxa only

CATEGORICAL COLUMN DETECTION (same rules as mbx_ezviz_all_levels_all_treatments.sh):
  ✓ Non-numeric string columns
  ✓ 2..n-1 unique values (not constant, not all-unique)
  ✓ At least one group with ≥ 2 members (needed for Kruskal-Wallis)
  ✗ Sample-id column, #q2:types row, numeric-only, all-unique, constant columns

COMMON ERRORS:
  "No significant taxa found" for a level
    → Normal — Kruskal-Wallis found no p ≤ 0.05 taxa. Boxplots folder may be empty.
  "ezstat failed for level s (species)"
    → Normal for sparse species-level annotations. Other levels continue.
  "argument is of length zero" or similar R error
    → Usually means a group has too few samples for pairwise comparison.
      Check that each group in the variable has ≥ 2 samples.

EXAMPLES:
  mbx_ezstat_all_levels_all_treatments.sh /path/to/mbX_pro_outputs_20250422_143022
  mbx_ezstat_all_levels_all_treatments.sh /path/to/mbX_pro_outputs_20250422_143022 --dry-run

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)  usage ;;
    --dry-run)  DRY_RUN=true; shift ;;
    -*)  err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '${1}'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -d "$MBX_OUT_DIR" ]] || err "Directory does not exist: '${MBX_OUT_DIR}'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# ── Locate Rscript (system-wide, outside conda) ───────────────────────────────
# Prefer absolute SYSTEM paths first; fall back to PATH lookup last so that
# `command -v Rscript` returning conda's R doesn't pre-empt system R.  The
# _strip_env wrapper still defends against conda env-var pollution.
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
  → Run mbx_ezclean_all_levels.sh first to install R, or: brew install r"

R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
info "Rscript : $RSCRIPT_CMD  (R $R_VERSION)"

# ── PID-based temp file IDs (avoids mktemp suffix issues on macOS bash 3.2) ───
_TMPID="${$}_$(date +%s)"
CATS_FILE="/tmp/mbx_stats_cats_${_TMPID}.txt"
CAT_DETECT_R="/tmp/mbx_stats_detect_${_TMPID}.R"
trap 'rm -f "$CATS_FILE" "$CAT_DETECT_R" /tmp/mbx_ezstat_${_TMPID}_*.R' EXIT

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
$DRY_RUN && warn "DRY-RUN — R code will be printed but NOT executed."
sep

# ─────────────────────────────────────────────────────────────────────────────
step "2/5 — Detect categorical metadata columns"
# ─────────────────────────────────────────────────────────────────────────────

cat > "$CAT_DETECT_R" << RDETECT
# Detect categorical columns — same logic as mbx_ezviz_all_levels_all_treatments.sh
metadata_path <- "${METADATA_TXT}"

ext <- tolower(tools::file_ext(metadata_path))
if (ext == "csv") {
  meta <- read.csv(metadata_path, header = TRUE,
                   check.names = FALSE, stringsAsFactors = FALSE)
} else {
  meta <- read.delim(metadata_path, header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE,
                     na.strings = c("", "NA", "N/A", "na", "n/a"))
}

# Remove QIIME2 #q2:types directive row if present
q2_rows <- grepl("^#", meta[[1]])
if (any(q2_rows)) {
  meta <- meta[!q2_rows, , drop = FALSE]
  cat("[INFO]  Removed QIIME2 #q2:types row(s).\n", file = stderr())
}

n_samples <- nrow(meta)
cat(sprintf("[INFO]  Metadata: %d samples, %d columns.\n",
            n_samples, ncol(meta)), file = stderr())

col_names <- names(meta)[-1]   # skip sample-id column

categorical          <- character(0)
skipped_numeric      <- character(0)
skipped_alluniq      <- character(0)
skipped_noval        <- character(0)
skipped_mingroupsize <- character(0)

for (col in col_names) {
  vals_raw <- meta[[col]]
  vals     <- vals_raw[!is.na(vals_raw) & nchar(trimws(as.character(vals_raw))) > 0]

  if (length(vals) == 0)                              { skipped_noval    <- c(skipped_noval, col);   next }
  num_test <- suppressWarnings(as.numeric(vals))
  if (!any(is.na(num_test)))                          { skipped_numeric  <- c(skipped_numeric, col); next }
  n_unique <- length(unique(trimws(as.character(vals))))
  if (n_unique <= 1)                                  { skipped_noval    <- c(skipped_noval, col);   next }
  if (n_unique == n_samples)                          { skipped_alluniq  <- c(skipped_alluniq, col); next }

  # Kruskal-Wallis needs at least 1 group with >= 2 samples
  group_counts <- table(trimws(as.character(vals)))
  if (max(group_counts) < 2)                         { skipped_mingroupsize <- c(skipped_mingroupsize, col); next }

  categorical <- c(categorical, col)
}

if (length(skipped_numeric)      > 0) cat(sprintf("[INFO]  Skipped numeric       : %s\n", paste(skipped_numeric, collapse=", ")),      file=stderr())
if (length(skipped_alluniq)      > 0) cat(sprintf("[INFO]  Skipped all-unique    : %s\n", paste(skipped_alluniq, collapse=", ")),      file=stderr())
if (length(skipped_noval)        > 0) cat(sprintf("[INFO]  Skipped empty/const   : %s\n", paste(skipped_noval, collapse=", ")),        file=stderr())
if (length(skipped_mingroupsize) > 0) cat(sprintf("[INFO]  Skipped singleton grp : %s\n", paste(skipped_mingroupsize, collapse=", ")), file=stderr())

if (length(categorical) == 0) {
  cat("[ERROR] No categorical columns found.\n", file = stderr())
  quit(status = 1)
}

cat(sprintf("[INFO]  Categorical columns : %d → %s\n",
            length(categorical), paste(categorical, collapse=", ")), file=stderr())

cat(paste(categorical, collapse="\n"), "\n", sep="")
RDETECT

if $DRY_RUN; then
  warn "[DRY-RUN] Would detect categorical columns from: $METADATA_TXT"
  CATEGORICAL_COLS=("Treatment" "SampleType")
  info "Dry-run placeholder columns: ${CATEGORICAL_COLS[*]}"
else
  _R --vanilla "$CAT_DETECT_R" > "$CATS_FILE" 2>&1 \
    || err "Categorical column detection failed.
  → Check your metadata: $METADATA_TXT"

  grep "^\[" "$CATS_FILE" || true

  # bash 3.2-safe array fill (no mapfile)
  CATEGORICAL_COLS=()
  while IFS= read -r _col; do
    [[ -n "$_col" ]] && CATEGORICAL_COLS+=("$_col")
  done < <(grep -v "^\[" "$CATS_FILE" | grep -v '^$' || true)

  [[ ${#CATEGORICAL_COLS[@]} -eq 0 ]] && err "No categorical columns detected.
  → Verify metadata has at least one grouping column with string values."
fi
rm -f "$CAT_DETECT_R"

ok "Categorical columns: ${CATEGORICAL_COLS[*]}"
sep

# ─────────────────────────────────────────────────────────────────────────────
step "3/5 — Create 10_stats/ directory structure"
# ─────────────────────────────────────────────────────────────────────────────

STATS_DIR="${MBX_OUT_DIR}/10_stats"
mkdir -p "$STATS_DIR"

declare -a SANITIZED_NAMES=()
for col in "${CATEGORICAL_COLS[@]}"; do
  SAFE="$(_sanitize_dirname "$col")"
  SANITIZED_NAMES+=("$SAFE")
  mkdir -p "${STATS_DIR}/${SAFE}"
  info "  Created: 10_stats/${SAFE}/"
done

N_VARS="${#CATEGORICAL_COLS[@]}"
N_LEVELS=7
TOTAL_RUNS=$(( N_VARS * N_LEVELS ))
info "Total ezstat() calls planned : ${TOTAL_RUNS}  (${N_LEVELS} levels × ${N_VARS} variables)"
info "Estimated runtime            : ~$(( TOTAL_RUNS / 3 ))–$(( TOTAL_RUNS * 2 )) minutes"
sep

# ─────────────────────────────────────────────────────────────────────────────
step "4/5 — Run ezstat() for all variables × all levels"
# ─────────────────────────────────────────────────────────────────────────────

LEVELS=("d" "p" "c" "o" "f" "g" "s")
LEVEL_NAMES=("domain" "phylum" "class" "order" "family" "genus" "species")

# Map level letter → canonical name used in ezstat output filenames
# (mirrors the mbX internal naming convention from ezclean/ezstat)
level_canonical() {
  case "$1" in
    d) echo "domains_or_kingdom" ;;
    p) echo "phylum" ;;
    c) echo "classes" ;;
    o) echo "orders" ;;
    f) echo "families" ;;
    g) echo "genera" ;;
    s) echo "species" ;;
  esac
}

RUN_NUMBER=0
PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0
FAIL_LOG=""

for v_idx in "${!CATEGORICAL_COLS[@]}"; do
  COL="${CATEGORICAL_COLS[$v_idx]}"
  SAFE="${SANITIZED_NAMES[$v_idx]}"
  VAR_DIR="${STATS_DIR}/${SAFE}"

  echo ""
  echo "  ════════════════════════════════════════════════════════════════"
  echo "  Variable: ${COL}  →  ${VAR_DIR}"
  echo "  ════════════════════════════════════════════════════════════════"

  for l_idx in "${!LEVELS[@]}"; do
    LVL="${LEVELS[$l_idx]}"
    LVL_NAME="${LEVEL_NAMES[$l_idx]}"
    LVL_CANONICAL="$(level_canonical "$LVL")"
    RUN_NUMBER=$(( RUN_NUMBER + 1 ))

    # Skip check: look for the KW output xlsx (first file ezstat writes)
    # Pattern: ezstat_KW_<canonical>_by_<variable>.xlsx
    KW_EXPECTED="${VAR_DIR}/ezstat_KW_${LVL_CANONICAL}_by_${SAFE}.xlsx"
    if [[ -f "$KW_EXPECTED" ]]; then
      skipped "[$RUN_NUMBER/$TOTAL_RUNS] level '${LVL}' (${LVL_NAME}) for '${COL}' — KW results already exist"
      SKIP_COUNT=$(( SKIP_COUNT + 1 ))
      continue
    fi

    info "[$RUN_NUMBER/$TOTAL_RUNS] ezstat() — level '${LVL}' (${LVL_NAME}) | variable '${COL}'"

    if $DRY_RUN; then
      echo "  [DRY-RUN] setwd(\"${VAR_DIR}\")"
      echo "  [DRY-RUN] library(mbX)"
      echo "  [DRY-RUN] ezstat("
      echo "  [DRY-RUN]   microbiome_data = \"$LEVEL7_CSV\","
      echo "  [DRY-RUN]   metadata        = \"$METADATA_TXT\","
      echo "  [DRY-RUN]   level           = \"${LVL}\","
      echo "  [DRY-RUN]   selected_metadata = \"${COL}\""
      echo "  [DRY-RUN] )"
      PASS_COUNT=$(( PASS_COUNT + 1 ))
      continue
    fi

    R_SCRIPT="/tmp/mbx_ezstat_${_TMPID}_${RUN_NUMBER}.R"

    cat > "$R_SCRIPT" << RSCRIPT
suppressPackageStartupMessages(library(mbX))

# CRITICAL: set working directory to the variable's subdirectory.
# ezstat() writes KW xlsx, pairwise xlsx, CLD xlsx, and Boxplots_<var>/
# all relative to the current working directory.
setwd("${VAR_DIR}")

cat(sprintf("[INFO]  cwd      : %s\n", getwd()))
cat(sprintf("[INFO]  level    : ${LVL} (${LVL_NAME})\n"))
cat(sprintf("[INFO]  variable : ${COL}\n"))

result <- tryCatch({

  out <- ezstat(
    microbiome_data   = "${LEVEL7_CSV}",
    metadata          = "${METADATA_TXT}",
    level             = "${LVL}",
    selected_metadata = "${COL}"
  )

  # ezstat() can return an error string on validation failure
  if (is.character(out)) {
    stop(sprintf("ezstat() returned an error message: %s", out))
  }

  cat(sprintf("[OK]    ezstat level '${LVL}' for '${COL}' complete.\n"))
  0L

}, error = function(e) {
  cat(sprintf("[ERROR] ezstat level '${LVL}' for '${COL}' failed:\n        %s\n",
              conditionMessage(e)))
  1L
}, warning = function(w) {
  # Catch and display warnings but continue execution
  cat(sprintf("[WARN]  ezstat warning (level '${LVL}'): %s\n", conditionMessage(w)))
  invokeRestart("muffleWarning")
})

quit(status = as.integer(result))
RSCRIPT

    timer_start
    if _R --vanilla "$R_SCRIPT" 2>&1; then
      timer_end
      ok "  ✔  level '${LVL}' (${LVL_NAME}) for '${COL}'"

      # Report what was created
      for _f in \
          "${VAR_DIR}/ezstat_KW_${LVL_CANONICAL}_by_${SAFE}.xlsx" \
          "${VAR_DIR}/ezstat_pairwise_${LVL_CANONICAL}_by_${SAFE}.xlsx" \
          "${VAR_DIR}/ezstat_CLD_Summary_${LVL_CANONICAL}_by_${SAFE}.xlsx"; do
        [[ -f "$_f" ]] && info "    ↳ $(basename "$_f")"
      done
      _BOX_DIR="${VAR_DIR}/Boxplots_${SAFE}"
      if [[ -d "$_BOX_DIR" ]]; then
        _N_BOX="$(find "$_BOX_DIR" -name "*.png" | wc -l | tr -d ' ')"
        info "    ↳ Boxplots_${SAFE}/ (${_N_BOX} significant taxa)"
      else
        info "    ↳ No significant taxa (p ≤ 0.05) at this level — no boxplots generated"
      fi

      PASS_COUNT=$(( PASS_COUNT + 1 ))
    else
      timer_end
      warn "  ✘  level '${LVL}' (${LVL_NAME}) for '${COL}' — see output above"
      FAIL_COUNT=$(( FAIL_COUNT + 1 ))
      FAIL_LOG="${FAIL_LOG}  - ${COL} × ${LVL_NAME} (${LVL})\n"
    fi
    rm -f "$R_SCRIPT"

  done   # levels loop
done     # variables loop

# ─────────────────────────────────────────────────────────────────────────────
step "5/5 — Write mbx_ezstat_info.txt + final summary"
# ─────────────────────────────────────────────────────────────────────────────

NOW="$(date '+%Y-%m-%d %H:%M:%S')"
INFO_TXT="${STATS_DIR}/mbx_ezstat_info.txt"

{
  echo "# ============================================================================"
  echo "# mbx_ezstat_info.txt"
  echo "# Generated by mbx_ezstat_all_levels_all_treatments.sh"
  echo "# Date        : $NOW"
  echo "# ============================================================================"
  echo ""
  echo "# ── Inputs ──────────────────────────────────────────────────────────────────"
  echo "MBX_OUTPUTS_DIR=$MBX_OUT_DIR"
  echo "LEVEL_7_CSV=$LEVEL7_CSV"
  echo "METADATA_TXT=$METADATA_TXT"
  echo ""
  echo "# ── Categorical variables ───────────────────────────────────────────────────"
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
  echo "# ── Output file reference ───────────────────────────────────────────────────"
  echo "# Per run, ezstat() creates inside 10_stats/<variable>/:"
  echo "#   ezstat_KW_<level>_by_<variable>.xlsx          Kruskal-Wallis p-values"
  echo "#   ezstat_pairwise_<level>_by_<variable>.xlsx    Dunn FDR pairwise results"
  echo "#   ezstat_CLD_Summary_<level>_by_<variable>.xlsx CLD letters + median + mean"
  echo "#   Boxplots_<variable>/                          PNGs for significant taxa"
  echo ""
  echo "# ── Using outputs with ezviz() for follow-up visualization ─────────────────"
  echo "# To visualize only the significant taxa found by ezstat:"
  echo "# library(mbX)"
  for col in "${CATEGORICAL_COLS[@]}"; do
    echo "# ezviz(\"$LEVEL7_CSV\", \"$METADATA_TXT\", \"g\", \"${col}\", threshold = 0.5)"
  done
} > "$INFO_TXT"

ok "Info file → $INFO_TXT"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete — no R code was executed."
else
  ok "ezstat complete!"
  ok "  Passed  : $PASS_COUNT / $TOTAL_RUNS"
  [[ "$SKIP_COUNT" -gt 0 ]] && ok "  Skipped : $SKIP_COUNT (already existed)"
  [[ "$FAIL_COUNT" -gt 0 ]] && warn "  Failed  : $FAIL_COUNT — check output above"
fi
sep
echo ""
echo "  Output structure:"
echo "  $MBX_OUT_DIR/"
echo "  └── 10_stats/"
for v_idx in "${!CATEGORICAL_COLS[@]}"; do
  COL="${CATEGORICAL_COLS[$v_idx]}"
  SAFE="${SANITIZED_NAMES[$v_idx]}"
  CANON="$(level_canonical "g")"
  printf '      ├── %s/\n' "$SAFE"
  printf '      │   ├── ezstat_KW_%s_by_%s.xlsx\n'          "$CANON" "$SAFE"
  printf '      │   ├── ezstat_pairwise_%s_by_%s.xlsx\n'    "$CANON" "$SAFE"
  printf '      │   ├── ezstat_CLD_Summary_%s_by_%s.xlsx\n' "$CANON" "$SAFE"
  printf '      │   └── Boxplots_%s/  (PNGs per significant taxon)\n' "$SAFE"
done
echo "      └── mbx_ezstat_info.txt"
echo ""
if [[ "$FAIL_COUNT" -gt 0 ]]; then
  echo "  ── Failed combinations ──────────────────────────────────────────"
  printf '%b' "$FAIL_LOG"
  echo ""
  echo "  Re-running is safe — completed levels are skipped automatically."
  echo "  Species (s) failures are normal for low-depth classification."
  echo ""
fi
