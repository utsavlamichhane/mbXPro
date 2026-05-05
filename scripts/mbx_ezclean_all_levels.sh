#!/usr/bin/env bash
# =============================================================================
#  mbx_ezclean_all_levels.sh
#  Install mbX (CRAN) and run ezclean() for all 7 taxonomic levels
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  STEPS:
#    1  Parse 7_taxonomy_csv/mbx_taxonomy_info.txt → level-7.csv + metadata
#    2  Check R installation → install via Homebrew if missing
#    3  Install mbX (CRAN v0.2.0) + dependencies if missing
#    4  Create 8_cleaned_files/
#    5  Run ezclean() for all 7 levels from inside 8_cleaned_files/
#         (working directory must be 8_cleaned_files/ — ezclean writes
#          intermediate files relative to cwd)
#    6  Verify outputs + write mbx_ezclean_info.txt
#
#  NOTE ON WORKING DIRECTORY:
#    ezclean() writes ~12 intermediate .xlsx files relative to cwd, then
#    cleans them up. It also creates one subdirectory per level:
#      8_cleaned_files/mbX_cleaned_genera_level-7/
#      └── mbX_cleaned_genera_level-7.xlsx
#
#  OUTPUT STRUCTURE:
#    mbX_pro_outputs_<timestamp>/
#    └── 8_cleaned_files/
#        ├── mbX_cleaned_domains_or_kingdom_level-7/
#        │   └── mbX_cleaned_domains_or_kingdom_level-7.xlsx
#        ├── mbX_cleaned_phylum_level-7/
#        ├── mbX_cleaned_classes_level-7/
#        ├── mbX_cleaned_orders_level-7/
#        ├── mbX_cleaned_families_level-7/
#        ├── mbX_cleaned_genera_level-7/
#        ├── mbX_cleaned_species_level-7/
#        └── mbx_ezclean_info.txt
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
sep()  { echo "────────────────────────────────────────────────────────────────"; }
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

# ── R env wrapper — strip conda's R_LIBS_USER pollution ──────────────────────
# When the user runs this script from inside a conda env (e.g. qiime2-amplicon),
# conda exports R_LIBS_USER / R_HOME / etc. that point to the conda env's R
# library tree.  System Rscript (Homebrew/CRAN) then tries to load packages
# from those paths and fails with errors like:
#     Error: shared object 'Rcpp.so' not found
#     Error in library.dynam: ... 'methods.dylib' not found
# We unset those vars at script start AND strip them again at every invocation
# of system Rscript.  This makes the script bullet-proof regardless of which
# conda env is active when it's run.
unset R_LIBS R_LIBS_USER R_LIBS_SITE \
      R_PROFILE R_PROFILE_USER R_ENVIRON R_ENVIRON_USER 2>/dev/null || true

_strip_env() {
  env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
      -u R_PROFILE -u R_PROFILE_USER -u R_ENVIRON -u R_ENVIRON_USER \
      -u R_PAPERSIZE -u R_INCLUDE_DIR -u R_DOC_DIR -u R_SHARE_DIR \
      "$@"
}

# Defined after RSCRIPT_CMD is found (see step 2)
_R() { _strip_env "$RSCRIPT_CMD" "$@"; }

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_ezclean_all_levels.sh — Run mbX::ezclean() for all 7 taxonomic levels

USAGE:
  mbx_ezclean_all_levels.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Reads level-7.csv and metadata.txt paths from:
    <mbX_pro_outputs_dir>/7_taxonomy_csv/mbx_taxonomy_info.txt

  Installs R (via Homebrew) and mbX (from CRAN) if not already present,
  then runs ezclean() for all 7 levels from inside 8_cleaned_files/:

    ezclean("level-7.csv", "metadata.txt", "d")  → Domain / Kingdom
    ezclean("level-7.csv", "metadata.txt", "p")  → Phylum
    ezclean("level-7.csv", "metadata.txt", "c")  → Class
    ezclean("level-7.csv", "metadata.txt", "o")  → Order
    ezclean("level-7.csv", "metadata.txt", "f")  → Family
    ezclean("level-7.csv", "metadata.txt", "g")  → Genus
    ezclean("level-7.csv", "metadata.txt", "s")  → Species

  level-7.csv contains all 7 taxonomic levels and is used for all runs.
  Each ezclean() call extracts the requested level's taxa from level-7.csv.

WHY level-7.csv FOR ALL LEVELS:
  The QIIME2 taxa barplot export at level-7 (species) contains all higher
  levels embedded in the taxonomy string (d__; p__; c__; o__; f__; g__; s__).
  ezclean() parses this string to extract whichever level you request.

OPTIONS:
  --dry-run      Print R commands without executing them
  --skip-install Skip R and mbX installation checks
  -h, --help     Show this help and exit

EXAMPLES:
  mbx_ezclean_all_levels.sh /path/to/mbX_pro_outputs_20250422_143022
  mbx_ezclean_all_levels.sh /path/to/mbX_pro_outputs_20250422_143022 --dry-run

OUTPUT (inside 8_cleaned_files/):
  mbX_cleaned_domains_or_kingdom_level-7/
  └── mbX_cleaned_domains_or_kingdom_level-7.xlsx
  mbX_cleaned_phylum_level-7/
  └── mbX_cleaned_phylum_level-7.xlsx
  ... (one subdirectory per level)
  mbX_cleaned_genera_level-7/
  └── mbX_cleaned_genera_level-7.xlsx     ← most commonly used downstream
  mbX_cleaned_species_level-7/
  └── mbX_cleaned_species_level-7.xlsx
  mbx_ezclean_info.txt

USING THE OUTPUT WITH mbX:
  The cleaned .xlsx files are the direct inputs for ezviz() and ezstat():
    ezviz("mbX_cleaned_genera_level-7/.../file.xlsx", "metadata.txt", "g", "Group")
    ezstat("mbX_cleaned_genera_level-7/.../file.xlsx", "metadata.txt", "g", "Group")

COMMON ERRORS:
  "brew: command not found"
    → Install Homebrew: /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  "there is no package called mbX"
    → The script will auto-install from CRAN. If it fails, run manually:
      Rscript -e 'install.packages("mbX")'
  "Please check the first header of the metadata"
    → Your metadata first column must be named: sample-id, id, sampleid, etc.
  "The file is not csv, xls, or xlsx format"
    → level-7.csv must have a .csv extension (check the file exists).
  "ezclean failed for level s (species)"
    → Species-level annotations are often sparse. The script continues with
      other levels and reports which ones succeeded.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
DRY_RUN=false
SKIP_INSTALL=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)       usage ;;
    --dry-run)       DRY_RUN=true; shift ;;
    --skip-install)  SKIP_INSTALL=true; shift ;;
    -*)  err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Unexpected extra argument: '${1}'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -d "$MBX_OUT_DIR" ]] || err "Directory does not exist: '${MBX_OUT_DIR}'"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

# ─────────────────────────────────────────────────────────────────────────────
step "1/6 — Read paths from mbx_taxonomy_info.txt"
# ─────────────────────────────────────────────────────────────────────────────

TAXONOMY_INFO="${MBX_OUT_DIR}/7_taxonomy_csv/mbx_taxonomy_info.txt"
[[ -f "$TAXONOMY_INFO" ]] || err "mbx_taxonomy_info.txt not found:
  $TAXONOMY_INFO
  → Run mbx_taxonomy_run.sh first."

LEVEL7_CSV="$(_read_key  "LEVEL_7_CSV"   "$TAXONOMY_INFO")"
METADATA_TXT="$(_read_key "METADATA_TXT" "$TAXONOMY_INFO")"

[[ -z "$LEVEL7_CSV"   ]] && err "LEVEL_7_CSV not found in mbx_taxonomy_info.txt"
[[ -z "$METADATA_TXT" ]] && err "METADATA_TXT not found in mbx_taxonomy_info.txt"
[[ -f "$LEVEL7_CSV"   ]] || err "level-7.csv not found: $LEVEL7_CSV
  → Run mbx_taxonomy_run.sh first."
[[ -f "$METADATA_TXT" ]] || err "Metadata file not found: $METADATA_TXT
  → Check that the path in mbx_taxonomy_info.txt is still correct."

info "level-7.csv  : $LEVEL7_CSV"
info "metadata.txt : $METADATA_TXT"
sep

# ─────────────────────────────────────────────────────────────────────────────
step "2/6 — Check R installation"
# ─────────────────────────────────────────────────────────────────────────────
# R is installed OUTSIDE the QIIME2 conda env (system-wide via Homebrew).
# Installing R into the QIIME2 conda env risks breaking Python dependencies.

RSCRIPT_CMD=""

# Search order: absolute SYSTEM paths first (Homebrew → /usr/local → /usr/bin),
# then PATH lookup as last resort.  This avoids picking up conda's R when
# the script is run from inside a conda env, which would re-introduce the
# library-poisoning issue we strip env vars to avoid.
for _candidate in \
    "/opt/homebrew/bin/Rscript" \
    "/usr/local/bin/Rscript" \
    "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
    "/usr/bin/Rscript" \
    "$(command -v Rscript 2>/dev/null || true)"; do
  if [[ -n "$_candidate" && -x "$_candidate" ]]; then
    RSCRIPT_CMD="$_candidate"
    break
  fi
done

if [[ -n "$RSCRIPT_CMD" ]]; then
  R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
  ok "R found: $RSCRIPT_CMD  (version $R_VERSION)"
elif $SKIP_INSTALL; then
  err "R not found and --skip-install was set.
  → Install R manually: brew install r
  → Then re-run this script."
else
  warn "R not found. Attempting to install via Homebrew..."

  # Check for Homebrew
  if ! command -v brew &>/dev/null; then
    err "Homebrew not found — cannot auto-install R.
  Install Homebrew first:
    /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"
  Then install R:
    brew install r
  Then re-run this script."
  fi

  if $DRY_RUN; then
    warn "[DRY-RUN] Would run: brew install r"
  else
    info "Running: brew install r  (this may take a few minutes)..."
    brew install r \
      || err "brew install r failed.
  → Try manually: brew install r
  → Or download R from https://cran.r-project.org/bin/macosx/"
    # Re-locate after install (absolute paths first; see note above)
    for _candidate in \
        "/opt/homebrew/bin/Rscript" \
        "/usr/local/bin/Rscript" \
        "/Library/Frameworks/R.framework/Resources/bin/Rscript" \
        "$(command -v Rscript 2>/dev/null || true)"; do
      if [[ -n "$_candidate" && -x "$_candidate" ]]; then
        RSCRIPT_CMD="$_candidate"
        break
      fi
    done
    [[ -n "$RSCRIPT_CMD" ]] || err "R was installed but Rscript still not found.
  → Open a new terminal and re-run this script."
    R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    ok "R installed: $RSCRIPT_CMD  (version $R_VERSION)"
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "3/6 — Install mbX (CRAN v0.2.0) + dependencies"
# ─────────────────────────────────────────────────────────────────────────────

if $SKIP_INSTALL; then
  skipped "mbX installation check (--skip-install)"
elif $DRY_RUN; then
  warn "[DRY-RUN] Would check/install mbX and dependencies via CRAN."
else
  info "Checking mbX installation..."

  _R - << 'RINSTALL'
# Check and install mbX + all required dependencies
required_pkgs <- c(
  "mbX",          # main package (CRAN v0.2.0)
  "openxlsx",     # read/write xlsx
  "readxl",       # read .xls / .xlsx
  "dplyr",        # data manipulation
  "tidyr",        # data reshaping
  "ggplot2",      # visualizations (ezviz / ezstat)
  "tools"         # file utilities (base, always available)
)

# Identify which packages need installing
to_install <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]

if (length(to_install) == 0) {
  cat("[OK]    All required packages already installed.\n")
} else {
  cat(sprintf("[INFO]  Installing %d package(s): %s\n",
              length(to_install), paste(to_install, collapse = ", ")))
  install.packages(
    to_install,
    repos   = "https://cloud.r-project.org",
    quiet   = FALSE,
    dependencies = TRUE
  )
  # Verify installation succeeded
  failed <- to_install[!sapply(to_install, requireNamespace, quietly = TRUE)]
  if (length(failed) > 0) {
    stop(sprintf(
      "Failed to install: %s\nTry manually: install.packages(c('%s'))",
      paste(failed, collapse = ", "),
      paste(failed, collapse = "', '")
    ))
  }
  cat("[OK]    All packages installed successfully.\n")
}

# Print installed version of mbX
mbx_ver <- tryCatch(
  packageVersion("mbX"),
  error = function(e) "unknown"
)
cat(sprintf("[INFO]  mbX version: %s\n", mbx_ver))
RINSTALL

  ok "mbX and dependencies ready."
fi

# ─────────────────────────────────────────────────────────────────────────────
step "4/6 — Create 8_cleaned_files/"
# ─────────────────────────────────────────────────────────────────────────────

CLEANED_DIR="${MBX_OUT_DIR}/8_cleaned_files"
mkdir -p "$CLEANED_DIR" \
  || err "Could not create: $CLEANED_DIR — check permissions."
ok "Output directory: $CLEANED_DIR"

# ─────────────────────────────────────────────────────────────────────────────
step "5/6 — Run ezclean() for all 7 taxonomic levels"
# ─────────────────────────────────────────────────────────────────────────────
# CRITICAL: R's working directory is set to 8_cleaned_files/ before running.
# ezclean() writes ~12 intermediate .xlsx files relative to cwd and deletes
# them on completion. The final cleaned file is placed in a subdirectory:
#   8_cleaned_files/mbX_cleaned_genera_level-7/mbX_cleaned_genera_level-7.xlsx
#
# Each level is run in a separate Rscript call so a failure at one level
# (e.g. sparse species annotations) does not abort the remaining levels.

LEVELS=("d" "p" "c" "o" "f" "g" "s")
LEVEL_NAMES=("Domain_Kingdom" "Phylum" "Class" "Order" "Family" "Genus" "Species")
RESULTS=()   # track pass/fail per level

for i in "${!LEVELS[@]}"; do
  LVL="${LEVELS[$i]}"
  LVL_NAME="${LEVEL_NAMES[$i]}"

  info "Running ezclean() — level '${LVL}' (${LVL_NAME})..."
  timer_start

  if $DRY_RUN; then
    echo ""
    echo "  [DRY-RUN] Would run inside $CLEANED_DIR:"
    echo "  Rscript -e '"
    echo "    library(mbX)"
    echo "    setwd(\"$CLEANED_DIR\")"
    echo "    ezclean("
    echo "      microbiome_data = \"$LEVEL7_CSV\","
    echo "      metadata        = \"$METADATA_TXT\","
    echo "      level           = \"${LVL}\""
    echo "    )'"
    echo ""
    RESULTS+=("${LVL}:DRY_RUN")
    continue
  fi

  # Write a per-level R script to a temp file (avoids quoting nightmares)
  R_SCRIPT="$(mktemp /tmp/mbx_ezclean_XXXXXX.R)"
  trap 'rm -f "$R_SCRIPT"' EXIT

  cat > "$R_SCRIPT" << RSCRIPT
library(mbX)

# Set working directory to 8_cleaned_files/ so ezclean() writes
# intermediate and final files here — NOT in the launch directory.
setwd("${CLEANED_DIR}")

cat(sprintf("[INFO]  Working directory: %s\n", getwd()))
cat(sprintf("[INFO]  level-7.csv : ${LEVEL7_CSV}\n"))
cat(sprintf("[INFO]  metadata    : ${METADATA_TXT}\n"))
cat(sprintf("[INFO]  level       : ${LVL} (${LVL_NAME})\n\n"))

result <- tryCatch(
  {
    out <- ezclean(
      microbiome_data = "${LEVEL7_CSV}",
      metadata        = "${METADATA_TXT}",
      level           = "${LVL}"
    )
    cat(sprintf("[OK]    ezclean level '${LVL}' complete.\n"))
    cat(sprintf("[OK]    Output file: %s\n", out))
    out
  },
  error = function(e) {
    cat(sprintf("[ERROR] ezclean level '${LVL}' failed: %s\n", conditionMessage(e)))
    NULL
  },
  warning = function(w) {
    cat(sprintf("[WARN]  ezclean level '${LVL}' warning: %s\n", conditionMessage(w))  )
    # Invoke restart to continue execution after warning
    invokeRestart("muffleWarning")
  }
)

if (is.null(result)) {
  quit(status = 1)
} else {
  quit(status = 0)
}
RSCRIPT

  # Run the R script; capture exit code without triggering set -e
  if _R --vanilla "$R_SCRIPT"; then
    timer_end
    ok "Level '${LVL}' (${LVL_NAME}) — SUCCESS"
    RESULTS+=("${LVL}:OK")
  else
    timer_end
    warn "Level '${LVL}' (${LVL_NAME}) — FAILED (see output above)"
    warn "  Common causes:"
    warn "    - Too few sequences classified to level '${LVL}' (normal for species)"
    warn "    - Metadata header not valid (must be sample-id or similar)"
    warn "    - level-7.csv missing required taxonomy prefixes"
    RESULTS+=("${LVL}:FAILED")
  fi
  rm -f "$R_SCRIPT"

done

# ─────────────────────────────────────────────────────────────────────────────
step "6/6 — Verify outputs + write mbx_ezclean_info.txt"
# ─────────────────────────────────────────────────────────────────────────────

NOW="$(date '+%Y-%m-%d %H:%M:%S')"
INFO_TXT="${CLEANED_DIR}/mbx_ezclean_info.txt"

# Map level letters to expected output filename prefixes
declare_level_map() {
  case "$1" in
    d) echo "mbX_cleaned_domains_or_kingdom" ;;
    p) echo "mbX_cleaned_phylum" ;;
    c) echo "mbX_cleaned_classes" ;;
    o) echo "mbX_cleaned_orders" ;;
    f) echo "mbX_cleaned_families" ;;
    g) echo "mbX_cleaned_genera" ;;
    s) echo "mbX_cleaned_species" ;;
  esac
}

PASS_COUNT=0
FAIL_COUNT=0
LEVEL_PATHS=""

for i in "${!LEVELS[@]}"; do
  LVL="${LEVELS[$i]}"
  PREFIX="$(declare_level_map "$LVL")"
  MICROBIOME_BASE="$(basename "$LEVEL7_CSV" .csv)"   # level-7
  EXPECTED_DIR="${CLEANED_DIR}/${PREFIX}_${MICROBIOME_BASE}"
  EXPECTED_XLS="${EXPECTED_DIR}/${PREFIX}_${MICROBIOME_BASE}.xlsx"

  if $DRY_RUN; then
    LEVEL_PATHS="${LEVEL_PATHS}LEVEL_${LVL^^}_XLSX=${EXPECTED_XLS}\n"
    continue
  fi

  RESULT="${RESULTS[$i]:-UNKNOWN}"
  if [[ -f "$EXPECTED_XLS" ]]; then
    ok "✔  level '${LVL}': $EXPECTED_XLS"
    LEVEL_PATHS="${LEVEL_PATHS}LEVEL_${LVL^^}_XLSX=${EXPECTED_XLS}\n"
    PASS_COUNT=$(( PASS_COUNT + 1 ))
  elif [[ "$RESULT" == *"FAILED"* ]]; then
    warn "✘  level '${LVL}': output not created (ezclean failed)"
    LEVEL_PATHS="${LEVEL_PATHS}LEVEL_${LVL^^}_XLSX=FAILED\n"
    FAIL_COUNT=$(( FAIL_COUNT + 1 ))
  else
    warn "✘  level '${LVL}': output not found at expected path: $EXPECTED_XLS"
    LEVEL_PATHS="${LEVEL_PATHS}LEVEL_${LVL^^}_XLSX=NOT_FOUND\n"
    FAIL_COUNT=$(( FAIL_COUNT + 1 ))
  fi
done

# Write info file
cat > "$INFO_TXT" << INFO
# ============================================================================
# mbx_ezclean_info.txt
# Generated by mbx_ezclean_all_levels.sh
# Date : $NOW
# ============================================================================

# ── Inputs used ──────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
LEVEL_7_CSV=$LEVEL7_CSV
METADATA_TXT=$METADATA_TXT
R_VERSION=${R_VERSION:-unknown}
MBX_PACKAGE_VERSION=0.2.0 (CRAN)
OUTPUT_DIR=$CLEANED_DIR

# ── Output files (one .xlsx per level) ───────────────────────────────────────
INFO
printf '%b' "$LEVEL_PATHS" >> "$INFO_TXT"
cat >> "$INFO_TXT" << INFO2

# ── Using the cleaned files with ezviz() and ezstat() ────────────────────────
# The cleaned .xlsx files are direct inputs for the next mbX steps.
# Replace <your_group_column> with the column name from your metadata
# that defines your experimental groups (e.g. "Treatment", "SampleType").
#
# Genus level (most common):
#   library(mbX)
#   ezviz(
#     microbiome_data  = "$CLEANED_DIR/mbX_cleaned_genera_level-7/mbX_cleaned_genera_level-7.xlsx",
#     metadata         = "$METADATA_TXT",
#     level            = "g",
#     selected_metadata = "<your_group_column>",
#     top_taxa         = 20
#   )
#   ezstat(
#     microbiome_data  = "$CLEANED_DIR/mbX_cleaned_genera_level-7/mbX_cleaned_genera_level-7.xlsx",
#     metadata         = "$METADATA_TXT",
#     level            = "g",
#     selected_metadata = "<your_group_column>"
#   )
INFO2

ok "Info file written → $INFO_TXT"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete — no R code was executed."
else
  ok "ezclean complete:  ${PASS_COUNT}/7 levels succeeded   ${FAIL_COUNT}/7 failed"
  [[ "$FAIL_COUNT" -gt 0 ]] && \
    warn "Some levels failed. Check output above. Species (s) failures are normal when classification depth is low."
fi
sep
echo ""
echo "  Output structure:"
echo "  $MBX_OUT_DIR/"
echo "  └── 8_cleaned_files/"
for LVL in d p c o f g s; do
  PREFIX="$(declare_level_map "$LVL")"
  MICROBIOME_BASE="$(basename "$LEVEL7_CSV" .csv)"
  printf '      ├── %s_%s/\n' "$PREFIX" "$MICROBIOME_BASE"
  printf '      │   └── %s_%s.xlsx\n' "$PREFIX" "$MICROBIOME_BASE"
done
echo "      └── mbx_ezclean_info.txt"
echo ""
echo "  ── Next steps with mbX R package ────────────────────────────────"
MICROBIOME_BASE="$(basename "$LEVEL7_CSV" .csv)"
printf '  ezviz("%s/mbX_cleaned_genera_%s/mbX_cleaned_genera_%s.xlsx",\n' \
  "$CLEANED_DIR" "$MICROBIOME_BASE" "$MICROBIOME_BASE"
printf '        "%s", "g", "<your_group_column>")\n' "$METADATA_TXT"
echo ""
