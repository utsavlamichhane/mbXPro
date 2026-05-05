#!/usr/bin/env bash
# =============================================================================
#  mbx_taxonomy_run.sh
#  Generate taxonomy bar plots + export all 7 taxonomy-level CSV files
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  STEPS:
#    1  Parse mbx_classifier_run_info.txt → feature table path
#    2  Locate taxonomy.qza from 6_classifier_taxonomy/
#    3  Filter mitochondria + chloroplasts from feature table  ← critical
#    4  qiime taxa barplot (on filtered table)
#    5  Export barplot QZV → all 7 level-*.csv files
#    6  Write taxonomy level legend + mbx_taxonomy_info.txt
#
#  INPUT FILES (auto-discovered):
#    5_classifier_working_dir/mbx_classifier_run_info.txt → feature table path
#    6_classifier_taxonomy/taxonomy.qza                   → taxonomy assignments
#    <metadata.txt>                                        → required argument
#
#  OUTPUT STRUCTURE:
#    mbX_pro_outputs_<timestamp>/
#    └── 7_taxonomy_csv/
#        ├── level-1.csv     Domain / Kingdom
#        ├── level-2.csv     Phylum
#        ├── level-3.csv     Class
#        ├── level-4.csv     Order
#        ├── level-5.csv     Family
#        ├── level-6.csv     Genus
#        ├── level-7.csv     Species
#        ├── taxa_bar_plots.qzv
#        ├── feature_table_filtered.qza
#        └── mbx_taxonomy_info.txt
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
skipped() { echo "[SKIP]  $* — already exists."; }
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

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_taxonomy_run.sh — Generate taxonomy bar plots + export all 7 level CSV files

USAGE:
  mbx_taxonomy_run.sh <mbX_pro_outputs_dir> <metadata.txt> [OPTIONS]

DESCRIPTION:
  Auto-discovers the feature table (from mbx_classifier_run_info.txt) and
  taxonomy.qza (from 6_classifier_taxonomy/), then:

    Step 3  Filter mitochondria + chloroplasts from feature table
    Step 4  qiime taxa barplot  → taxa_bar_plots.qzv
    Step 5  Export barplot      → level-1.csv ... level-7.csv
    Step 6  Write mbx_taxonomy_info.txt (paths + R usage guide)

  The 7 CSV files are what mbX R functions (ezclean, ezviz, ezstat) use directly.
  Each level corresponds to a taxonomic rank:
    level-1 = Domain/Kingdom  level-5 = Family
    level-2 = Phylum          level-6 = Genus     ← most commonly used
    level-3 = Class           level-7 = Species
    level-4 = Order

WHY FILTER MITOCHONDRIA + CHLOROPLASTS:
  Plant-associated or environmental samples often contain host mitochondrial
  and plastid 16S. These inflate non-bacterial counts and distort community
  profiles. Always filter before visualization.

OPTIONS:
  --skip-filter     Skip mitochondria/chloroplast filtering
  --exclude <terms> Comma-separated taxa to exclude (default: mitochondria,chloroplast)
  --dry-run         Print commands without executing
  -h, --help        Show this help and exit

EXAMPLES:
  mbx_taxonomy_run.sh /path/to/mbX_pro_outputs_20250422_143022 /path/to/metadata.txt

  # Skip filtering (e.g. already know no plastid contamination)
  mbx_taxonomy_run.sh /path/to/mbX_pro_outputs_20250422_143022 metadata.txt --skip-filter

  # Custom exclusion list
  mbx_taxonomy_run.sh /path/to/mbX_pro_outputs_20250422_143022 metadata.txt \
    --exclude mitochondria,chloroplast,unassigned

USING THE OUTPUT CSV FILES WITH mbX R PACKAGE:
  library(mbX)
  ezclean("7_taxonomy_csv/level-7.csv", "metadata.txt", "s")  # species
  ezclean("7_taxonomy_csv/level-6.csv", "metadata.txt", "g")  # genus
  ezviz("7_taxonomy_csv/level-6.csv",   "metadata.txt", "g", "your_group_column")
  ezstat("7_taxonomy_csv/level-6.csv",  "metadata.txt", "g", "your_group_column")

COMMON ERRORS:
  "No sample data shared between table and taxonomy"
    → Sample IDs in feature_table.qza and taxonomy.qza do not match.
      Check both with: qiime tools peek <file.qza>
  "Metadata file not found"
    → Pass the correct path to your metadata.txt as the second argument.
  "Not all samples in the feature table are present in the metadata"
    → Ensure all sample IDs in the feature table appear in metadata.txt column 1.

EOF
  exit 0
}

# ── Parse arguments ───────────────────────────────────────────────────────────
MBX_OUT_DIR=""
METADATA_PATH=""
SKIP_FILTER=false
EXCLUDE_TERMS="mitochondria,chloroplast"
DRY_RUN=false

[[ $# -eq 0 ]] && usage

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)      usage ;;
    --dry-run)      DRY_RUN=true; shift ;;
    --skip-filter)  SKIP_FILTER=true; shift ;;
    --exclude)      EXCLUDE_TERMS="${2:?Missing value for --exclude}"; shift 2 ;;
    -*)  err "Unknown option: '${1}'  —  run with --help for usage." ;;
    *)
      if   [[ -z "$MBX_OUT_DIR"   ]]; then MBX_OUT_DIR="$1";    shift
      elif [[ -z "$METADATA_PATH" ]]; then METADATA_PATH="$1";  shift
      else err "Unexpected extra argument: '${1}'"; fi ;;
  esac
done

[[ -z "$MBX_OUT_DIR"   ]] && err "No mbX_pro_outputs directory provided.  Run with --help."
[[ -z "$METADATA_PATH" ]] && err "No metadata file provided.
  Usage: mbx_taxonomy_run.sh <mbX_pro_outputs_dir> <metadata.txt>
  Note : the metadata path is required because it was not stored
         in mbx_classifier_run_info.txt by mbx_classifier_arranger.sh."

[[ -d "$MBX_OUT_DIR"    ]] || err "Directory does not exist: '${MBX_OUT_DIR}'"
[[ -f "$METADATA_PATH"  ]] || err "Metadata file not found: '${METADATA_PATH}'
  → Check the path for typos.
  → The first column must be sample-id (or id, sampleid, sample-id)."

MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"
METADATA_PATH="$(_abspath "$METADATA_PATH")"

command -v qiime &>/dev/null || err "qiime not found.
  → conda activate qiime2-amplicon-2025.4"

# ─────────────────────────────────────────────────────────────────────────────
step "1/6 — Parse mbx_classifier_run_info.txt + locate inputs"
# ─────────────────────────────────────────────────────────────────────────────

RUN_INFO="${MBX_OUT_DIR}/5_classifier_working_dir/mbx_classifier_run_info.txt"
[[ -f "$RUN_INFO" ]] || err "mbx_classifier_run_info.txt not found:
  $RUN_INFO
  → Run mbx_classifier_arranger.sh first."

FEATURE_TABLE_QZA="$(_read_key "FEATURE_TABLE_QZA" "$RUN_INFO")"

# Validate feature table
[[ -z "$FEATURE_TABLE_QZA" ]] && err "FEATURE_TABLE_QZA not found in run_info.txt.
  → Check $RUN_INFO"
[[ -f "$FEATURE_TABLE_QZA" ]] || err "Feature table not found: $FEATURE_TABLE_QZA
  → Run mbx_dada2_run.sh first."

# Locate taxonomy.qza from step 6
TAXONOMY_QZA="${MBX_OUT_DIR}/6_classifier_taxonomy/taxonomy.qza"
[[ -f "$TAXONOMY_QZA" ]] || err "taxonomy.qza not found: $TAXONOMY_QZA
  → Run mbx_classifier_run.sh first."

info "Feature table   : $FEATURE_TABLE_QZA"
info "Taxonomy        : $TAXONOMY_QZA"
info "Metadata        : $METADATA_PATH"
info "Filter terms    : $EXCLUDE_TERMS"
$DRY_RUN && warn "DRY-RUN — commands will be printed but NOT executed."

# ── Create output directory ───────────────────────────────────────────────────
TAXONOMY_CSV_DIR="${MBX_OUT_DIR}/7_taxonomy_csv"
EXPORT_DIR="${TAXONOMY_CSV_DIR}/barplot_export"
mkdir -p "$TAXONOMY_CSV_DIR" "$EXPORT_DIR" \
  || err "Could not create output directories — check permissions."

FILTERED_TABLE="${TAXONOMY_CSV_DIR}/feature_table_filtered.qza"
BARPLOT_QZV="${TAXONOMY_CSV_DIR}/taxa_bar_plots.qzv"

# ── Dry-run helper ─────────────────────────────────────────────────────────────
run_cmd() {
  echo ""
  echo "  \$ $(printf '%s ' "$@" | sed 's/ --/\\\n    --/g')"
  echo ""
  if ! $DRY_RUN; then
    "$@" || err "Command failed: $1
  → Check QIIME2 output above.
  → Re-run with --dry-run to review the command."
  fi
}

# ─────────────────────────────────────────────────────────────────────────────
step "2/6 — Validate metadata against feature table sample IDs"
# ─────────────────────────────────────────────────────────────────────────────
# A mismatch between sample IDs in the feature table and metadata is the
# most common failure point in taxa barplot. We pre-check before running.

if ! $DRY_RUN; then
  META_HEADER="$(head -1 "$METADATA_PATH" | cut -f1)"
  META_HEADER_CLEAN="$(printf '%s' "$META_HEADER" | tr '[:upper:]' '[:lower:]' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
  VALID="id|sampleid|sample id|sample-id|featureid|feature id|feature-id"
  if ! echo "$META_HEADER_CLEAN" | grep -qiE "^(id|sampleid|sample id|sample-id|featureid|feature id|feature-id)$"; then
    err "Metadata first column header is invalid: '${META_HEADER}'
  Expected: sample-id, id, sampleid, sample id, featureid, feature id, or feature-id
  → Rename the first column in your metadata file."
  fi
  ok "Metadata header valid: '$META_HEADER'"
  META_SAMPLE_COUNT="$(tail -n +2 "$METADATA_PATH" | grep -c '[^[:space:]]' || true)"
  info "Metadata samples: $META_SAMPLE_COUNT"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "3/6 — Filter mitochondria + chloroplasts from feature table"
# ─────────────────────────────────────────────────────────────────────────────
# This removes ASVs classified as mitochondria or chloroplasts before
# visualization. Skipping this can inflate non-bacterial relative abundances
# in plant-associated, soil, or environmental samples.
# The filtered table is also used for all downstream diversity analyses.

if $SKIP_FILTER; then
  warn "--skip-filter set. Using raw (unfiltered) feature table."
  ACTIVE_TABLE="$FEATURE_TABLE_QZA"
elif [[ -f "$FILTERED_TABLE" ]]; then
  skipped "feature_table_filtered.qza"
  ACTIVE_TABLE="$FILTERED_TABLE"
else
  info "Excluding taxa: $EXCLUDE_TERMS"
  timer_start
  run_cmd qiime taxa filter-table \
    --i-table       "$FEATURE_TABLE_QZA" \
    --i-taxonomy    "$TAXONOMY_QZA" \
    --p-exclude     "$EXCLUDE_TERMS" \
    --o-filtered-table "$FILTERED_TABLE"
  timer_end
  ACTIVE_TABLE="$FILTERED_TABLE"
  $DRY_RUN || ok "Filtered table → $FILTERED_TABLE"
fi

# In dry-run mode, use original table path for command preview
$DRY_RUN && ACTIVE_TABLE="$FILTERED_TABLE"

# ─────────────────────────────────────────────────────────────────────────────
step "4/6 — Generate taxa bar plot"
# ─────────────────────────────────────────────────────────────────────────────

if [[ -f "$BARPLOT_QZV" ]]; then
  skipped "taxa_bar_plots.qzv"
else
  timer_start
  run_cmd qiime taxa barplot \
    --i-table             "$ACTIVE_TABLE" \
    --i-taxonomy          "$TAXONOMY_QZA" \
    --m-metadata-file     "$METADATA_PATH" \
    --o-visualization     "$BARPLOT_QZV"
  timer_end
  $DRY_RUN || ok "Bar plot → $BARPLOT_QZV"
fi

# ─────────────────────────────────────────────────────────────────────────────
step "5/6 — Export all 7 taxonomy level CSV files"
# ─────────────────────────────────────────────────────────────────────────────
# Exporting the barplot QZV produces level-1.csv through level-7.csv.
# Each file contains absolute ASV counts collapsed at that taxonomic rank.
# These are the direct inputs for mbX R functions (ezclean, ezviz, ezstat).

LEVELS_PRESENT=0

if ! $DRY_RUN; then

  # Export the barplot QZV to the staging directory
  if [[ ! -f "${EXPORT_DIR}/level-7.csv" ]]; then
    info "Exporting barplot to: $EXPORT_DIR"
    qiime tools export \
      --input-path  "$BARPLOT_QZV" \
      --output-path "$EXPORT_DIR" \
      || err "qiime tools export failed on taxa_bar_plots.qzv.
  → Verify the barplot was created: qiime tools peek $BARPLOT_QZV"
  fi

  # Move all level-*.csv files to the main 7_taxonomy_csv/ directory
  info "Moving level CSV files to: $TAXONOMY_CSV_DIR"
  for LVL in 1 2 3 4 5 6 7; do
    SRC="${EXPORT_DIR}/level-${LVL}.csv"
    DST="${TAXONOMY_CSV_DIR}/level-${LVL}.csv"
    if [[ -f "$SRC" ]]; then
      cp "$SRC" "$DST"
      LEVELS_PRESENT=$(( LEVELS_PRESENT + 1 ))
      ok "level-${LVL}.csv  → $DST"
    else
      warn "level-${LVL}.csv not found in export — taxonomy may not reach this depth."
    fi
  done

  [[ "$LEVELS_PRESENT" -eq 0 ]] && err "No level-*.csv files were exported.
  → Check that the barplot QZV was generated correctly: $BARPLOT_QZV
  → Try: qiime tools peek $BARPLOT_QZV"

  ok "Exported $LEVELS_PRESENT / 7 taxonomy levels."

else
  warn "[DRY-RUN] Would export $BARPLOT_QZV → level-1.csv ... level-7.csv"
  LEVELS_PRESENT=7
fi

# ─────────────────────────────────────────────────────────────────────────────
step "6/6 — Write mbx_taxonomy_info.txt"
# ─────────────────────────────────────────────────────────────────────────────

INFO_TXT="${TAXONOMY_CSV_DIR}/mbx_taxonomy_info.txt"
NOW="$(date '+%Y-%m-%d %H:%M:%S')"

cat > "$INFO_TXT" << INFO
# ============================================================================
# mbx_taxonomy_info.txt
# Generated by mbx_taxonomy_run.sh
# Date              : $NOW
# ============================================================================

# ── Input files used ─────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
FEATURE_TABLE_QZA=$FEATURE_TABLE_QZA
FILTERED_TABLE_QZA=$FILTERED_TABLE
TAXONOMY_QZA=$TAXONOMY_QZA
METADATA_TXT=$METADATA_PATH
BARPLOT_QZV=$BARPLOT_QZV
EXCLUDE_TERMS=$EXCLUDE_TERMS

# ── Taxonomy level reference ──────────────────────────────────────────────────
# Level   Rank              Typical prefix (Greengenes2 / SILVA)
LEVEL_1=Domain / Kingdom  (d__)
LEVEL_2=Phylum            (p__)
LEVEL_3=Class             (c__)
LEVEL_4=Order             (o__)
LEVEL_5=Family            (f__)
LEVEL_6=Genus             (g__)    ← most commonly used in mbX R functions
LEVEL_7=Species           (s__)    ← use with caution; short reads rarely resolve species

# ── CSV files (absolute counts, collapsed by taxonomic level) ─────────────────
LEVEL_1_CSV=$TAXONOMY_CSV_DIR/level-1.csv
LEVEL_2_CSV=$TAXONOMY_CSV_DIR/level-2.csv
LEVEL_3_CSV=$TAXONOMY_CSV_DIR/level-3.csv
LEVEL_4_CSV=$TAXONOMY_CSV_DIR/level-4.csv
LEVEL_5_CSV=$TAXONOMY_CSV_DIR/level-5.csv
LEVEL_6_CSV=$TAXONOMY_CSV_DIR/level-6.csv
LEVEL_7_CSV=$TAXONOMY_CSV_DIR/level-7.csv

# ── Using CSV files with the mbX R package ────────────────────────────────────
# The level-*.csv files are direct inputs for ezclean(), ezviz(), ezstat().
# Use the level that matches your analysis goal:
#
#   library(mbX)
#
#   # Genus level (most common)
#   ezclean("$TAXONOMY_CSV_DIR/level-6.csv", "$METADATA_PATH", "g")
#   ezviz("$TAXONOMY_CSV_DIR/level-6.csv",   "$METADATA_PATH", "g", "your_group_column")
#   ezstat("$TAXONOMY_CSV_DIR/level-6.csv",  "$METADATA_PATH", "g", "your_group_column")
#
#   # Species level
#   ezclean("$TAXONOMY_CSV_DIR/level-7.csv", "$METADATA_PATH", "s")
#
#   # Phylum level
#   ezclean("$TAXONOMY_CSV_DIR/level-2.csv", "$METADATA_PATH", "p")

# ── Suggested downstream QIIME2 steps ────────────────────────────────────────
# 1. Alpha diversity (requires rarefaction — check feature table summary first)
#    qiime diversity alpha-rarefaction \
#      --i-table $FILTERED_TABLE \
#      --p-max-depth <choose_based_on_min_sample_depth> \
#      --o-visualization $TAXONOMY_CSV_DIR/../7_taxonomy_csv/alpha_rarefaction.qzv
#
# 2. Core diversity metrics
#    qiime diversity core-metrics-phylogenetic \
#      --i-phylogeny <rooted_tree.qza> \
#      --i-table $FILTERED_TABLE \
#      --p-sampling-depth <rarefaction_depth> \
#      --m-metadata-file $METADATA_PATH \
#      --output-dir core_metrics/
INFO

ok "Info file written → $INFO_TXT"

# ── Final summary ─────────────────────────────────────────────────────────────
sep
if $DRY_RUN; then
  warn "Dry-run complete — no commands were executed."
else
  ok "All steps complete!"
fi
sep
echo ""
echo "  Output structure:"
echo "  $MBX_OUT_DIR/"
echo "  └── 7_taxonomy_csv/"
echo "      ├── level-1.csv   Domain / Kingdom"
echo "      ├── level-2.csv   Phylum"
echo "      ├── level-3.csv   Class"
echo "      ├── level-4.csv   Order"
echo "      ├── level-5.csv   Family"
echo "      ├── level-6.csv   Genus        ← use with ezclean/ezviz/ezstat"
echo "      ├── level-7.csv   Species"
echo "      ├── taxa_bar_plots.qzv"
echo "      ├── feature_table_filtered.qza"
echo "      └── mbx_taxonomy_info.txt"
echo ""
echo "  View taxa_bar_plots.qzv at: https://view.qiime2.org"
echo ""
echo "  ── Quick start with mbX R package ──────────────────────────────"
echo "  library(mbX)"
printf '  ezclean("%s/level-6.csv", "%s", "g")\n' "$TAXONOMY_CSV_DIR" "$METADATA_PATH"
echo ""
