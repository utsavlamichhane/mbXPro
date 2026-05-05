#!/usr/bin/env bash
# =============================================================================
#  mbx_final_report.sh   (step 18)
#  Consolidated final HTML + PDF report for the entire mbX Pro pipeline
#
#  Compatible with bash 3.2+ (macOS default shell)
#
#  PURPOSE
#    Produce a single self-contained HTML report (with all figures embedded
#    as base64) plus a print-ready PDF that summarises every step that was
#    actually executed in this output directory.  Steps that have NOT been
#    run are listed explicitly, with the exact command needed to add them
#    to the report on a subsequent run.
#
#  PHILOSOPHY
#    Designed to satisfy a deeply technical reviewer who wants to see:
#      - every parameter that was used (no hidden defaults)
#      - every algorithmic choice with its citation
#      - per-sample quality metrics in one glance
#      - convergent biomarker findings across statistical, ML, and network
#        approaches (taxa that are simultaneously DA, important, AND hubs)
#      - explicit "why did this not run / why is this empty" reasons for
#        any expected output that is missing
#
#  OUTPUT
#    18_final_report/
#      mbX_pro_final_report.html     <- single self-contained HTML
#      mbX_pro_final_report.pdf      <- print-ready PDF (A4)
#      mbx_final_report_info.txt     <- provenance file
#
#  DEPENDENCIES
#    R + packages: htmltools, base64enc, openxlsx, knitr, jsonlite
#    PDF engine (auto-detected, in this preference order):
#      1) Headless Google Chrome / Chromium / Edge (best CSS support)
#      2) wkhtmltopdf (auto-install via 'brew install wkhtmltopdf' if missing)
#      3) macOS cupsfilter (built-in fallback, lower-quality)
#      4) HTML-only mode with print instruction
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

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
  cat <<'EOF'

mbx_final_report.sh  --  Consolidated HTML + PDF report (step 18)

USAGE:
  mbx_final_report.sh <mbX_pro_outputs_dir> [OPTIONS]

DESCRIPTION:
  Walks the entire mbX_pro_outputs_<TIMESTAMP>/ directory, reads every info
  file produced by steps 0-17, and renders a single self-contained HTML
  report (with base64-embedded images) plus a print-ready A4 PDF.

  Steps that have not been run are listed explicitly with the command
  needed to add them.  Steps that ran but produced no significant findings
  are flagged with a scientific reason rather than silently omitted.

OPTIONS:
  --no-pdf            HTML only -- skip PDF rendering
  --pdf-engine <name> Force a specific PDF engine; one of:
                        chrome      headless Chrome / Chromium / Edge
                        wkhtmltopdf wkhtmltopdf binary
                        cupsfilter  macOS built-in (lower quality)
                      (default: auto-detect best available)
  --skip-install      Do not auto-install missing R packages or wkhtmltopdf
  --force-rerun       Recompute even if 18_final_report/ already exists
  --dry-run           Show what would happen without writing anything
  -h, --help          Show this help and exit

OUTPUTS
  <mbX_pro_outputs_dir>/18_final_report/
    mbX_pro_final_report.html      (single self-contained file, share-able)
    mbX_pro_final_report.pdf       (A4 print-ready)
    mbx_final_report_info.txt      (provenance)

EXAMPLES
  mbx_final_report.sh /path/to/mbX_pro_outputs_20260417_121431
  mbx_final_report.sh /path/to/mbX_pro_outputs_20260417_121431 --no-pdf
  mbx_final_report.sh /path/to/mbX_pro_outputs_20260417_121431 --pdf-engine chrome
EOF
  exit 0
}

# ── Argument parsing ──────────────────────────────────────────────────────────
MBX_OUT_DIR=""
NO_PDF=false
PDF_ENGINE=""
SKIP_INSTALL=false
FORCE_RERUN=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)        usage ;;
    --no-pdf)         NO_PDF=true;       shift ;;
    --pdf-engine)     PDF_ENGINE="$2";   shift 2 ;;
    --skip-install)   SKIP_INSTALL=true; shift ;;
    --force-rerun)    FORCE_RERUN=true;  shift ;;
    --dry-run)        DRY_RUN=true;      shift ;;
    -*)               err "Unknown option: $1\n  Run: mbx_final_report.sh --help" ;;
    *)
      if [[ -z "$MBX_OUT_DIR" ]]; then MBX_OUT_DIR="$1"; shift
      else err "Multiple positional arguments -- only one MBX_OUT_DIR expected.\n  Got extra: $1"
      fi
      ;;
  esac
done

[[ -z "$MBX_OUT_DIR" ]] && err "Missing required argument: <mbX_pro_outputs_dir>
  Run: mbx_final_report.sh --help"
[[ -d "$MBX_OUT_DIR" ]] || err "Not a directory: $MBX_OUT_DIR
  Did you pass the path to mbX_pro_outputs_<TIMESTAMP>/?"
MBX_OUT_DIR="$(_abspath "$MBX_OUT_DIR")"

NOW="$(date '+%Y-%m-%d %H:%M:%S %Z')"

# ── Locate the mbX Pro logo (mbX_Pro_icon.png) ───────────────────────────────
# Search order:
#   1) Same directory as the running script (developer install)
#   2) ~/bin/  (standard user install)
#   3) /usr/local/share/mbx_pro/ + /opt/mbx_pro/ (system installs)
# If absent, the report renders without a logo (graceful).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOGO_PATH=""
for _candidate in \
    "${SCRIPT_DIR}/mbX_Pro_icon.png" \
    "${HOME}/bin/mbX_Pro_icon.png" \
    "/usr/local/share/mbx_pro/mbX_Pro_icon.png" \
    "/opt/mbx_pro/mbX_Pro_icon.png" \
    "${HOME}/.config/mbx_pro/mbX_Pro_icon.png"; do
  if [[ -f "$_candidate" ]]; then
    LOGO_PATH="$_candidate"
    break
  fi
done
[[ -n "$LOGO_PATH" ]] && info "Logo found : $LOGO_PATH" || warn "Logo not found -- report will render without an mbX Pro logo.
  -> To add a logo: place mbX_Pro_icon.png in ~/bin/ and re-run."

# ─────────────────────────────────────────────────────────────────────────────
step "1/5 — Discover which pipeline steps were executed"
# ─────────────────────────────────────────────────────────────────────────────

# CANONICAL pipeline manifest -- the master list of every step in mbX Pro.
# Each entry: "step_id|directory_name|info_file|script_name|short_description"
# Order matters: this is the order steps must run in.
PIPELINE_MANIFEST="0|0_primer_handling|mbx_primer_info.txt|mbx_primer_identifier.sh|Primer detection from raw FASTQ
1|1_manifest_file|manifest.txt|create_manifest.sh|QIIME2 manifest construction
2|2_first_artifact_file|Paired_End_artifact.qza|artifact_creator.sh|QIIME2 .qza artifact creation
3|3_dada2_parameters|dada2_parameters.txt|mbx_dada2_parameter_finder.sh|DADA2 truncation-length selection
4|4_dada2_outputs|feature_table.qza|mbx_dada2_run.sh|DADA2 denoising + ASV inference
5|5_classifier_working_dir|mbx_classifier_run_info.txt|mbx_classifier_arranger.sh|Greengenes2 classifier preparation
6|6_classifier_taxonomy|taxonomy.qza|mbx_classifier_run.sh|Naive-Bayes taxonomy classification
7|7_taxonomy_csv|mbx_taxonomy_info.txt|mbx_taxonomy_run.sh|Taxonomy CSVs + mito/chloro filtering
8|8_cleaned_files||mbx_ezclean_all_levels.sh|mbX ezclean per taxonomic level
9|9_visualization_entire|mbx_ezviz_info.txt|mbx_ezviz_all_levels_all_treatments.sh|Stacked-bar visualisations (ezviz)
10|10_stats|mbx_ezstat_info.txt|mbx_ezstat_all_levels_all_treatments.sh|Kruskal-Wallis + Dunn statistics (ezstat)
11|11_pre_diversity|mbx_pre_diversity_info.txt|mbx_pre_diversity_parameters.sh|Phylogenetic tree + sampling-depth selection
12|12_alpha_diversity_results|mbx_alpha_diversity_info.txt|mbx_alpha_diversity_run.sh|Alpha-diversity metrics + statistics
13|13_beta_diversity_results|mbx_beta_diversity_info.txt|mbx_beta_diversity_run.sh|Beta-diversity (PCoA, PERMANOVA, Adonis)
14|14_differential_abundance_ANCOMBC2|mbx_ancombc2_info.txt|mbx_ancombc2_run.sh|Differential abundance via ANCOMBC2
15|15_picrust2|mbx_picrust2_info.txt|mbx_picrust_run.sh|Functional prediction via PICRUSt2
16|16_ml_biomarkers|mbx_ml_biomarkers_info.txt|mbx_ml_classifier_run.sh|Random Forest biomarker classifier
17|17_co_occurrence_networks|mbx_networks_info.txt|mbx_network_run.sh|CLR + Spearman co-occurrence networks"

# Discover which steps have output on disk
STEPS_DONE=()
STEPS_MISSING=()
STEPS_PARTIAL=()
STEP_MAX_DONE=-1

# Use a temp file as a poor-man's hash since bash 3.2 lacks associative arrays
_DISC_TMP="/tmp/mbx_report_discover_${$}_$(date +%s)"
mkdir -p "$_DISC_TMP"

while IFS='|' read -r SID SDIR SINFO SSCRIPT SDESC; do
  [[ -z "$SID" ]] && continue
  STEP_PATH="${MBX_OUT_DIR}/${SDIR}"
  STATUS="MISSING"
  REASON=""
  if [[ -d "$STEP_PATH" ]]; then
    if [[ -n "$SINFO" ]]; then
      if [[ -f "${STEP_PATH}/${SINFO}" ]]; then
        STATUS="DONE"
      else
        # Directory exists but the canonical info file does not -- treat as
        # partial.  Step 8 has an empty SINFO (no canonical info file by
        # design) so we use the presence of any cleaned xlsx as the proxy.
        STATUS="PARTIAL"
        REASON="Directory present but ${SINFO} missing -- step may have crashed mid-way"
      fi
    else
      # Step 8 case: any *.xlsx anywhere under it counts as "done"
      if find "$STEP_PATH" -maxdepth 3 -name "*.xlsx" 2>/dev/null | head -1 | grep -q .; then
        STATUS="DONE"
      else
        STATUS="PARTIAL"
        REASON="Directory present but no cleaned .xlsx files found"
      fi
    fi
  fi
  case "$STATUS" in
    DONE)    STEPS_DONE+=("$SID");    [[ "$SID" -gt "$STEP_MAX_DONE" ]] && STEP_MAX_DONE="$SID" ;;
    PARTIAL) STEPS_PARTIAL+=("$SID"); [[ "$SID" -gt "$STEP_MAX_DONE" ]] && STEP_MAX_DONE="$SID" ;;
    MISSING) STEPS_MISSING+=("$SID") ;;
  esac
  echo "${STATUS}|${SID}|${SDIR}|${SINFO}|${SSCRIPT}|${SDESC}|${REASON}" \
    > "${_DISC_TMP}/step_${SID}.txt"
done <<< "$PIPELINE_MANIFEST"

info "Steps DONE     : ${#STEPS_DONE[@]} -> ${STEPS_DONE[*]:-(none)}"
[[ ${#STEPS_PARTIAL[@]} -gt 0 ]] && warn "Steps PARTIAL  : ${#STEPS_PARTIAL[@]} -> ${STEPS_PARTIAL[*]}"
info "Steps MISSING  : ${#STEPS_MISSING[@]} -> ${STEPS_MISSING[*]:-(none)}"
info "Highest step   : $STEP_MAX_DONE"

if [[ "$STEP_MAX_DONE" -lt 0 ]]; then
  err "No pipeline steps appear to have been run in this directory.
  -> Did you pass the right path?  Expected an mbX_pro_outputs_<TIMESTAMP>/ directory
     containing at least 0_primer_handling/."
fi

# ─────────────────────────────────────────────────────────────────────────────
step "2/5 — Locate Rscript + check / install required R packages"
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
  -> Install R: brew install r"

R_VERSION="$(_R --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
ok "Rscript : $RSCRIPT_CMD  (R $R_VERSION)"

# Required R packages
_TMPID="${$}_$(date +%s)"
PKG_CHECK_R="/tmp/mbx_report_pkgchk_${_TMPID}.R"
trap 'rm -rf /tmp/mbx_report_*_${_TMPID}* "$_DISC_TMP" 2>/dev/null' EXIT

cat > "$PKG_CHECK_R" << 'RPKG'
required <- c("htmltools", "base64enc", "openxlsx", "jsonlite")
optional <- c("knitr")  # used for kable() formatting if available
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) == 0) {
  cat("[OK]    All required R packages installed.\n"); quit(status = 0)
}
cat(sprintf("[INFO]  Missing R packages: %s\n", paste(missing, collapse = ", ")))
quit(status = 2)
RPKG

if $DRY_RUN; then
  warn "[DRY-RUN] Would check/install R packages: htmltools base64enc openxlsx jsonlite"
else
  set +e
  _R --vanilla "$PKG_CHECK_R"
  RC=$?
  set -e
  if [[ $RC -eq 2 ]]; then
    if $SKIP_INSTALL; then
      err "Required R packages missing and --skip-install was set.
  -> Install:  Rscript -e 'install.packages(c(\"htmltools\",\"base64enc\",\"openxlsx\",\"jsonlite\"))'"
    fi
    info "Installing missing R packages from CRAN ..."
    PKG_INST_R="/tmp/mbx_report_pkginst_${_TMPID}.R"
    cat > "$PKG_INST_R" << 'RPI'
required <- c("htmltools", "base64enc", "openxlsx", "jsonlite")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org",
                   quiet = FALSE, dependencies = TRUE)
  failed <- missing[!sapply(missing, requireNamespace, quietly = TRUE)]
  if (length(failed) > 0) stop(sprintf("Failed: %s", paste(failed, collapse = ", ")))
}
cat("[OK]    R packages installed.\n")
RPI
    _R --vanilla "$PKG_INST_R" || err "R package installation failed."
    rm -f "$PKG_INST_R"
  elif [[ $RC -ne 0 ]]; then
    err "R package check failed (exit $RC)"
  fi
fi
rm -f "$PKG_CHECK_R"

# ─────────────────────────────────────────────────────────────────────────────
step "3/5 — Generate HTML report"
# ─────────────────────────────────────────────────────────────────────────────

REPORT_DIR="${MBX_OUT_DIR}/18_final_report"
HTML_OUT="${REPORT_DIR}/mbX_pro_final_report.html"
PDF_OUT="${REPORT_DIR}/mbX_pro_final_report.pdf"
INFO_OUT="${REPORT_DIR}/mbx_final_report_info.txt"

if [[ -f "$HTML_OUT" && "$FORCE_RERUN" == false && "$DRY_RUN" == false ]]; then
  warn "Report already exists: $HTML_OUT"
  warn "Use --force-rerun to regenerate."
  if [[ -f "$PDF_OUT" || "$NO_PDF" == true ]]; then
    sep
    ok "Report ready (cached) -> $HTML_OUT"
    [[ -f "$PDF_OUT" ]] && ok "PDF ready (cached)    -> $PDF_OUT"
    sep
    exit 0
  fi
  warn "PDF missing -- will regenerate PDF only."
fi

if $DRY_RUN; then
  info "[DRY-RUN] Would create: $REPORT_DIR/"
  info "[DRY-RUN] Would generate: $HTML_OUT"
  $NO_PDF || info "[DRY-RUN] Would generate: $PDF_OUT"
  info "[DRY-RUN] Would generate: $INFO_OUT"
  exit 0
fi

mkdir -p "$REPORT_DIR" || err "Could not create $REPORT_DIR"

# Write the discovery manifest as JSON so R can ingest it cleanly
DISC_JSON="${REPORT_DIR}/.discovery.json"
{
  echo '['
  FIRST=true
  for f in "$_DISC_TMP"/step_*.txt; do
    IFS='|' read -r STATUS SID SDIR SINFO SSCRIPT SDESC REASON < "$f"
    $FIRST || echo ','
    FIRST=false
    # Escape quotes / backslashes in JSON values
    _esc() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
    printf '  {"id":%s,"status":"%s","dir":"%s","info":"%s","script":"%s","desc":"%s","reason":"%s"}' \
      "$SID" "$(_esc "$STATUS")" "$(_esc "$SDIR")" "$(_esc "$SINFO")" "$(_esc "$SSCRIPT")" "$(_esc "$SDESC")" "$(_esc "$REASON")"
  done
  echo
  echo ']'
} > "$DISC_JSON"

ok "Discovery manifest written -> $DISC_JSON"

# ── R script that does the heavy lifting ─────────────────────────────────────
REPORT_R="/tmp/mbx_report_build_${_TMPID}.R"
cat > "$REPORT_R" << 'REPORT_R_EOF'
# =============================================================================
#  R script: build the consolidated mbX Pro HTML report
#  Invoked by mbx_final_report.sh -- not meant to be run standalone.
# =============================================================================
suppressPackageStartupMessages({
  library(htmltools)
  library(base64enc)
  library(openxlsx)
  library(jsonlite)
})

# ── Inputs from bash via env vars ────────────────────────────────────────────
MBX_OUT_DIR  <- Sys.getenv("MBX_REPORT_OUT_DIR")
REPORT_DIR   <- Sys.getenv("MBX_REPORT_REPORT_DIR")
HTML_OUT     <- Sys.getenv("MBX_REPORT_HTML_OUT")
DISC_JSON    <- Sys.getenv("MBX_REPORT_DISC_JSON")
GEN_AT       <- Sys.getenv("MBX_REPORT_NOW")
R_VERSION    <- Sys.getenv("MBX_REPORT_R_VERSION")
LOGO_PATH    <- Sys.getenv("MBX_REPORT_LOGO_PATH")
SCRIPT_VERSION <- "1.0.0"

# ── relpath: convert absolute paths under MBX_OUT_DIR to relative paths ─────
# Used everywhere a path is DISPLAYED to the reader.  Note: the txt info
# files inside each step's directory are NOT touched -- those keep their
# original absolute paths so they stay portable across machines.
relpath <- function(p) {
  if (is.null(p) || length(p) == 0) return("")
  # Vectorise over multiple paths
  sapply(p, function(x) {
    if (is.null(x) || is.na(x) || !nzchar(x)) return("")
    x <- as.character(x)
    if (startsWith(x, MBX_OUT_DIR))
      sub(paste0("^", gsub("([.|*+?(){}^$\\\\])", "\\\\\\1", MBX_OUT_DIR), "/?"),
          "", x)
    else x
  }, USE.NAMES = FALSE)
}

stopifnot(nzchar(MBX_OUT_DIR), dir.exists(MBX_OUT_DIR))

# ── Embed logo (if found) as base64 data URI ────────────────────────────────
LOGO_DATA_URI <- ""
if (nzchar(LOGO_PATH) && file.exists(LOGO_PATH)) {
  ext <- tolower(tools::file_ext(LOGO_PATH))
  mime <- switch(ext, png = "image/png", jpg = "image/jpeg",
                 jpeg = "image/jpeg", svg = "image/svg+xml", "image/png")
  raw <- readBin(LOGO_PATH, what = "raw", n = file.info(LOGO_PATH)$size)
  LOGO_DATA_URI <- paste0("data:", mime, ";base64,",
                          base64enc::base64encode(raw))
  cat(sprintf("[R] Logo embedded: %s (%.1f KB)\n",
              basename(LOGO_PATH), file.info(LOGO_PATH)$size / 1024))
} else {
  cat("[R] No logo provided -- header will render without logo.\n")
}

discovery <- jsonlite::fromJSON(DISC_JSON)
discovery <- discovery[order(discovery$id), ]
done_ids     <- discovery$id[discovery$status == "DONE"]
partial_ids  <- discovery$id[discovery$status == "PARTIAL"]
missing_ids  <- discovery$id[discovery$status == "MISSING"]
max_done_id  <- if (length(c(done_ids, partial_ids))) max(c(done_ids, partial_ids)) else -1

cat(sprintf("[R] Done=%d Partial=%d Missing=%d Max=%d\n",
            length(done_ids), length(partial_ids),
            length(missing_ids), max_done_id))

# ── Helpers ──────────────────────────────────────────────────────────────────
read_kv <- function(file) {
  if (!file.exists(file)) return(list())
  L <- readLines(file, warn = FALSE)
  L <- L[!grepl("^\\s*#", L) & nzchar(L)]
  out <- list()
  for (line in L) {
    if (!grepl("=", line, fixed = TRUE)) next
    parts <- strsplit(line, "=", fixed = TRUE)[[1]]
    if (length(parts) < 2) next
    key <- trimws(parts[1])
    val <- paste(parts[-1], collapse = "=")
    out[[key]] <- if (key %in% names(out)) c(out[[key]], val) else val
  }
  out
}

embed_image <- function(path, alt = "", max_kb_warn = 2048, css_class = "fig") {
  if (!file.exists(path)) {
    return(tags$div(class = "missing-fig",
                    sprintf("(figure not found: %s)", basename(path))))
  }
  ext <- tolower(tools::file_ext(path))
  mime <- switch(ext,
    png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg",
    gif = "image/gif", svg = "image/svg+xml", "image/png")
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  b64 <- base64enc::base64encode(raw)
  size_kb <- round(file.info(path)$size / 1024, 1)
  if (size_kb > max_kb_warn)
    cat(sprintf("[R] WARN: %s is %.1f KB (large)\n", basename(path), size_kb))
  tags$figure(class = css_class,
    tags$img(src = paste0("data:", mime, ";base64,", b64), alt = alt),
    if (nzchar(alt)) tags$figcaption(alt) else NULL)
}

# Render a data.frame as a nicely-styled HTML table
df_to_html <- function(df, caption = NULL, max_rows = NULL,
                       digits = 4, table_class = "data-tbl") {
  if (is.null(df) || nrow(df) == 0)
    return(tags$p(class = "empty-tbl", "(no data)"))
  if (!is.null(max_rows) && nrow(df) > max_rows) {
    df_trim <- df[seq_len(max_rows), , drop = FALSE]
    note <- tags$p(class = "tbl-note",
                   sprintf("(showing first %d of %d rows)", max_rows, nrow(df)))
  } else {
    df_trim <- df
    note <- NULL
  }
  num_cols <- sapply(df_trim, is.numeric)
  for (i in which(num_cols)) {
    df_trim[[i]] <- formatC(df_trim[[i]], digits = digits,
                            format = "g", big.mark = ",")
  }
  hdr <- tags$thead(tags$tr(lapply(names(df_trim),
                                   function(n) tags$th(n))))
  body <- tags$tbody(lapply(seq_len(nrow(df_trim)), function(r) {
    tags$tr(lapply(df_trim[r, ], function(c) tags$td(as.character(c))))
  }))
  tagList(
    if (!is.null(caption)) tags$p(class = "tbl-cap", caption) else NULL,
    tags$table(class = table_class, hdr, body),
    note)
}

# ────────────────────────────────────────────────────────────────────────────
# CSS  (kept inline in the head so the file is fully self-contained)
# ────────────────────────────────────────────────────────────────────────────
CSS <- '
:root{
  --fg:#1a1a1a; --fg-muted:#5a6470; --bg:#ffffff; --bg-soft:#fafbfc;
  --border:#e1e4e8; --primary:#2c3e50; --accent:#1f6feb;
  --success:#1f883d; --warning:#bf8700; --danger:#cf222e; --info:#0969da;
}
*{box-sizing:border-box}
html,body{margin:0;padding:0;background:var(--bg);color:var(--fg);
  font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,
              "Helvetica Neue",Arial,sans-serif;
  font-size:15px; line-height:1.55;}
.container{max-width:1200px;margin:0 auto;padding:32px 56px;}
h1{font-size:32px;color:var(--primary);margin:0 0 4px 0;font-weight:700;
   letter-spacing:-0.5px;}
h2{font-size:24px;color:var(--primary);margin:48px 0 12px;
   padding-bottom:8px;border-bottom:2px solid var(--border);font-weight:600;}
h3{font-size:18px;color:var(--primary);margin:28px 0 8px;font-weight:600;}
h4{font-size:15px;color:var(--fg);margin:18px 0 6px;font-weight:600;
   text-transform:uppercase;letter-spacing:0.5px;}
p,li{font-size:15px;color:var(--fg);}
a{color:var(--accent);text-decoration:none;}
a:hover{text-decoration:underline;}
code{background:var(--bg-soft);padding:2px 6px;border-radius:3px;
     font-family:ui-monospace,"SF Mono",Consolas,monospace;font-size:13px;}
pre{background:var(--bg-soft);padding:14px 18px;border-radius:6px;
    overflow-x:auto;font-size:13px;line-height:1.45;border:1px solid var(--border);}

.subtitle{font-size:17px;color:var(--fg-muted);margin:0 0 28px;font-weight:400;}
.cover-meta{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));
            gap:16px 24px;background:var(--bg-soft);padding:20px;
            border:1px solid var(--border);border-radius:8px;margin:24px 0;}
.cover-meta > div{font-size:13px;min-width:0;
                   overflow-wrap:anywhere;word-break:break-word;
                   line-height:1.45;}
.cover-meta > div b{display:block;color:var(--primary);font-size:11px;
                     text-transform:uppercase;letter-spacing:0.5px;
                     margin-bottom:4px;font-weight:700;}
.cover-meta .path-val{font-family:ui-monospace,"SF Mono",Menlo,monospace;
                       font-size:11.5px;color:var(--fg);}

.logo-header{display:flex;align-items:center;gap:18px;margin:0 0 8px 0;}
.logo-header img{max-height:54px;max-width:160px;width:auto;height:auto;
                  flex-shrink:0;}
.logo-header .titles{flex:1;min-width:0;}
.logo-header .titles h1{margin:0;}
.logo-header .titles .subtitle{margin:2px 0 0 0;}
@media print{.logo-header img{max-height:38px;max-width:120px;}}

.cite-box{border:2px solid var(--accent);background:#f0f7ff;padding:18px 22px;
           border-radius:8px;margin:24px 0;}
.cite-box .cite-title{font-weight:700;color:var(--accent);font-size:13px;
                       text-transform:uppercase;letter-spacing:0.5px;
                       margin:0 0 10px 0;}
.cite-box .cite-line{font-size:14px;line-height:1.6;margin:6px 0;}
.cite-box .cite-line.full{font-style:italic;color:var(--fg-muted);
                            padding-left:14px;border-left:3px solid var(--border);
                            margin-top:12px;}
.cite-box code{font-size:12.5px;}

.path-note{font-size:12.5px;color:var(--fg-muted);font-style:italic;
            margin:8px 0 18px 0;padding:8px 12px;background:var(--bg-soft);
            border-left:3px solid var(--info);border-radius:3px;}

.toc{background:var(--bg-soft);border:1px solid var(--border);
     border-radius:8px;padding:18px 24px;margin:24px 0;}
.toc h3{margin:0 0 8px;font-size:16px;}
.toc ol{margin:8px 0 0 0;padding-left:24px;column-count:2;column-gap:32px;}
.toc li{padding:2px 0;break-inside:avoid;}

.box{border-left:4px solid var(--info);background:#f0f7ff;
     padding:14px 18px;margin:14px 0;border-radius:4px;}
.box.success{border-left-color:var(--success);background:#dafbe1;}
.box.warning{border-left-color:var(--warning);background:#fff8c5;}
.box.danger {border-left-color:var(--danger); background:#ffebe9;}
.box.info   {border-left-color:var(--info);   background:#ddf4ff;}
.box .box-title{font-weight:700;color:var(--primary);margin:0 0 6px 0;
                 font-size:13px;text-transform:uppercase;letter-spacing:0.5px;}
.box.success .box-title{color:var(--success);}
.box.warning .box-title{color:var(--warning);}
.box.danger  .box-title{color:var(--danger);}

.dashboard{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
           gap:14px;margin:20px 0;}
.dashboard .card{background:var(--bg-soft);border:1px solid var(--border);
                 border-radius:8px;padding:16px;text-align:center;}
.dashboard .card .num{font-size:28px;font-weight:700;color:var(--accent);
                       display:block;line-height:1.1;}
.dashboard .card .lbl{font-size:11px;color:var(--fg-muted);
                       text-transform:uppercase;letter-spacing:0.5px;
                       margin-top:4px;display:block;}

.data-tbl{width:100%;border-collapse:collapse;margin:12px 0;font-size:13.5px;}
.data-tbl th,.data-tbl td{border:1px solid var(--border);padding:7px 10px;
                          text-align:left;vertical-align:top;}
.data-tbl thead th{background:var(--bg-soft);font-weight:600;color:var(--primary);
                    font-size:12.5px;text-transform:uppercase;letter-spacing:0.3px;}
.data-tbl tbody tr:nth-child(even){background:#fafbfc;}
.tbl-cap{font-style:italic;color:var(--fg-muted);font-size:13px;margin:6px 0;}
.tbl-note{font-size:12px;color:var(--fg-muted);font-style:italic;margin-top:4px;}
.empty-tbl{font-style:italic;color:var(--fg-muted);}

figure.fig{margin:18px 0;text-align:center;}
figure.fig img{max-width:100%;height:auto;border:1px solid var(--border);
                border-radius:4px;}
figure.fig figcaption{font-size:13px;color:var(--fg-muted);
                       font-style:italic;margin-top:8px;}
.fig-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));
          gap:18px;margin:18px 0;}
.fig-grid figure{margin:0;}
.missing-fig{padding:18px;background:var(--bg-soft);border:1px dashed var(--border);
              text-align:center;color:var(--fg-muted);font-style:italic;
              border-radius:4px;}

details{background:var(--bg-soft);border:1px solid var(--border);
        border-radius:6px;padding:14px 18px;margin:14px 0;}
details > summary{cursor:pointer;font-weight:600;color:var(--primary);
                  font-size:15px;outline:none;list-style:disclosure-closed;}
details[open] > summary{margin-bottom:12px;list-style:disclosure-open;}

.step-block{margin:30px 0;padding:0;}
.step-block .step-hdr{display:flex;align-items:center;gap:14px;
                       margin:0 0 8px;}
.step-id-badge{display:inline-block;background:var(--primary);color:#fff;
                font-size:12px;font-weight:700;padding:4px 10px;border-radius:20px;
                min-width:48px;text-align:center;}
.step-id-badge.done{background:var(--success);}
.step-id-badge.partial{background:var(--warning);}
.step-id-badge.missing{background:var(--fg-muted);}

.kv-list{margin:8px 0;padding:0;list-style:none;}
.kv-list li{padding:5px 0;border-bottom:1px dotted var(--border);
            display:flex;font-size:13.5px;}
.kv-list li:last-child{border-bottom:none;}
.kv-list .key{font-weight:600;color:var(--primary);min-width:240px;
              padding-right:12px;}
.kv-list .val{color:var(--fg);font-family:ui-monospace,"SF Mono",monospace;
              font-size:12.5px;word-break:break-all;}

footer{margin-top:60px;padding-top:24px;border-top:1px solid var(--border);
       color:var(--fg-muted);font-size:12.5px;text-align:center;}
.footnote{font-size:12px;color:var(--fg-muted);}

.bibliography{font-size:13px;}
.bibliography ol{padding-left:22px;}
.bibliography li{padding:4px 0;line-height:1.5;}

@page{
  size: A4;
  margin: 14mm 12mm 14mm 12mm;
}
@media print{
  html,body{font-size:10pt;background:#fff;}
  .container{max-width:100%;padding:0;margin:0;}
  h1{font-size:20pt;} h2{font-size:15pt;} h3{font-size:12pt;}
  h2{page-break-before:always;}
  .step-block{page-break-inside:avoid;}
  /* Collapse multi-column figure grids to a single column so wide
     publication-quality plots are never clipped by page width. */
  .fig-grid{display:block !important;}
  .fig-grid figure{margin:12px 0 !important;}
  /* Constrain every embedded figure to one page (height-wise) and the
     printable width so Chrome / wkhtmltopdf never crop mid-image. */
  figure.fig{margin:10px 0;page-break-inside:avoid;text-align:center;}
  figure.fig img{
    max-width:100% !important;
    max-height:84vh !important;   /* leave room for caption + page margin */
    width:auto !important; height:auto !important;
    object-fit:contain;
    display:block;margin:0 auto;
  }
  /* Wide data tables: scale font + allow wrapping so they fit page width. */
  .data-tbl{font-size:8.5pt;table-layout:auto;width:100% !important;
            word-break:break-word;}
  .data-tbl th, .data-tbl td{padding:3px 5px;}
  pre,code{white-space:pre-wrap !important;word-break:break-all;}
  .kv-list .val{word-break:break-all;}
  details{page-break-inside:avoid;border:none;background:transparent;padding:0;}
  details > summary{display:block;list-style:none;}
  details > summary::-webkit-details-marker{display:none;}
  details:not([open]) > *:not(summary){display:block;}
  .toc{page-break-after:always;}
  a{color:inherit;}
  .box{page-break-inside:avoid;}
  figure.fig{page-break-inside:avoid;}
  .data-tbl{page-break-inside:auto;}
  .data-tbl tr{page-break-inside:avoid;}
}
'

# ────────────────────────────────────────────────────────────────────────────
# Step-specific data extractors
# ────────────────────────────────────────────────────────────────────────────

# ---- Step 0: primer detection ----
extract_step0 <- function() {
  f <- file.path(MBX_OUT_DIR, "0_primer_handling", "mbx_primer_info.txt")
  if (!file.exists(f)) return(NULL)
  kv <- read_kv(f)
  list(
    fwd = kv$FORWARD_PRIMER_SEQUENCE %||% kv$FORWARD_PRIMER %||% kv$FWD_PRIMER %||% "(not detected)",
    rev = kv$REVERSE_PRIMER_SEQUENCE %||% kv$REVERSE_PRIMER %||% kv$REV_PRIMER %||% "(not detected)",
    fwd_len = kv$FORWARD_PRIMER_LENGTH %||% kv$FWD_PRIMER_LENGTH %||% "",
    rev_len = kv$REVERSE_PRIMER_LENGTH %||% kv$REV_PRIMER_LENGTH %||% "",
    orientation = kv$ORIENTATION %||% kv$DETECTED_ORIENTATION %||% "(unknown)",
    auto_detected = kv$AUTO_DETECTED %||% "yes",
    fastq_dir = kv$FASTQ_DIR %||% "",
    detection_status = kv$DETECTION_STATUS %||% "",
    confidence_level = kv$CONFIDENCE_LEVEL %||% "",
    inferred_region  = kv$INFERRED_REGION  %||% "",
    detection_note   = kv$DETECTION_NOTE   %||% "",
    raw_kv = kv
  )
}
`%||%` <- function(a, b) if (is.null(a) || all(!nzchar(a))) b else a

# ---- Step 1: manifest ----
extract_step1 <- function() {
  f <- file.path(MBX_OUT_DIR, "1_manifest_file", "manifest.txt")
  if (!file.exists(f)) return(NULL)
  L <- readLines(f, warn = FALSE)
  L <- L[!grepl("^\\s*#", L) & nzchar(L)]
  if (length(L) < 2) return(list(n_samples = 0, paired = NA, manifest_path = f))
  hdr <- strsplit(L[1], "\t")[[1]]
  list(
    n_samples = length(L) - 1,
    paired = "reverse-absolute-filepath" %in% tolower(hdr),
    manifest_path = f,
    columns = paste(hdr, collapse = ", ")
  )
}

# ---- Step 3: DADA2 parameters ----
extract_step3 <- function() {
  f <- file.path(MBX_OUT_DIR, "3_dada2_parameters", "dada2_parameters.txt")
  if (!file.exists(f)) return(NULL)
  L <- readLines(f, warn = FALSE)
  flags <- regmatches(L, regexpr("--p-[a-z-]+\\s+[^\\s]+", L, perl = TRUE))
  flags <- flags[nzchar(flags)]
  list(file = f, flags = flags, raw = paste(L, collapse = "\n"))
}

# ---- Step 4: DADA2 stats (read counts per sample) ----
extract_step4 <- function() {
  d <- file.path(MBX_OUT_DIR, "4_dada2_outputs")
  if (!dir.exists(d)) return(NULL)
  out <- list(
    feature_table = file.path(d, "feature_table.qza"),
    rep_seqs      = file.path(d, "representative_sequences.qza"),
    stats         = file.path(d, "dada2_stats.qza"),
    has_ft        = file.exists(file.path(d, "feature_table.qza")),
    has_rep       = file.exists(file.path(d, "representative_sequences.qza"))
  )
  # Try to find the per-sample stats TSV inside the .qzv extracted dir
  ftqzv <- file.path(d, "feature_table.qzv")
  out$has_ftqzv <- file.exists(ftqzv)
  out
}

# ---- Step 5/6: classifier ----
extract_step56 <- function() {
  f <- file.path(MBX_OUT_DIR, "5_classifier_working_dir", "mbx_classifier_run_info.txt")
  if (!file.exists(f)) return(NULL)
  kv <- read_kv(f)
  list(
    gg2_version       = kv$GG2_VERSION %||% kv$GG2_VER %||% "(unknown)",
    gg2_full_path     = kv$GG2_FULL_LENGTH %||% kv$GG2_FULL_LENGTH_FNA %||% kv$GG2_FULL_LENGTH_FNA_QZA %||% "",
    gg2_tax_path      = kv$GG2_TAX %||% kv$GG2_TAX_QZA %||% kv$GG2_TAXONOMY_QZA %||% "",
    classifier_qza    = kv$TRAINED_CLASSIFIER_QZA %||% kv$TRAINED_CLASSIFIER %||% kv$CLASSIFIER_QZA %||% "",
    asv_min_len       = kv$ASV_MIN_LEN %||% kv$MIN_ASV_LENGTH %||% "",
    asv_max_len       = kv$ASV_MAX_LEN %||% kv$MAX_ASV_LENGTH %||% "",
    classifier_mode   = kv$CLASSIFIER_MODE %||% "region-specific",
    classifier_source = kv$CLASSIFIER_SOURCE %||% "local-training",
    sklearn_ver       = kv$SCIKIT_LEARN_VERSION %||% "(unknown)",
    sklearn_family    = kv$SCIKIT_LEARN_FAMILY %||% "(unknown)",
    zenodo_url        = kv$ZENODO_RECORD_URL %||% "https://zenodo.org/records/20021035",
    zenodo_qiime2     = kv$ZENODO_QIIME2_USED %||% "(none)",
    zenodo_gg2        = kv$ZENODO_GG2_USED %||% "(none)",
    zenodo_filename   = kv$ZENODO_FILENAME %||% "(none)",
    zenodo_sha_exp    = kv$ZENODO_SHA256_EXPECTED %||% "(none)",
    zenodo_sha_act    = kv$ZENODO_SHA256_ACTUAL %||% "(none)",
    zenodo_note       = kv$ZENODO_NOTE %||% "",
    fwd_primer        = kv$FORWARD_PRIMER_SEQUENCE %||% "None",
    rev_primer        = kv$REVERSE_PRIMER_SEQUENCE %||% "None",
    target_region     = kv$TARGET_REGION %||% "Unknown",
    taxonomy_qza      = file.path(MBX_OUT_DIR, "6_classifier_taxonomy", "taxonomy.qza"),
    raw_kv            = kv
  )
}

# ---- Step 7: taxonomy CSVs + filtered table ----
extract_step7 <- function() {
  f <- file.path(MBX_OUT_DIR, "7_taxonomy_csv", "mbx_taxonomy_info.txt")
  if (!file.exists(f)) return(NULL)
  kv <- read_kv(f)
  d <- file.path(MBX_OUT_DIR, "7_taxonomy_csv")
  csvs <- list.files(d, pattern = "^level-\\d+\\.csv$", full.names = TRUE)
  list(
    metadata_txt   = kv$METADATA_TXT %||% "",
    filtered_table = kv$FILTERED_TABLE %||% file.path(d, "feature_table_filtered.qza"),
    barplot_qzv    = kv$BARPLOT_QZV %||% file.path(d, "taxa_bar_plots.qzv"),
    level_csvs     = csvs,
    n_levels       = length(csvs),
    raw_kv         = kv
  )
}

# ---- Step 8: cleaned files ----
extract_step8 <- function() {
  d <- file.path(MBX_OUT_DIR, "8_cleaned_files")
  if (!dir.exists(d)) return(NULL)
  subdirs <- list.dirs(d, recursive = FALSE)
  per_level <- lapply(subdirs, function(sd) {
    xlsx <- list.files(sd, pattern = "\\.xlsx$", full.names = TRUE)
    if (length(xlsx) == 0) return(NULL)
    df <- tryCatch(read.xlsx(xlsx[1], check.names = FALSE),
                   error = function(e) NULL)
    list(dir = basename(sd), xlsx = xlsx[1],
         n_samples = if (!is.null(df)) nrow(df) else NA,
         n_features = if (!is.null(df)) ncol(df) - 1 else NA)
  })
  per_level <- per_level[!sapply(per_level, is.null)]
  list(per_level = per_level, n_levels = length(per_level))
}

# ---- Step 9: visualizations ----
extract_step9 <- function() {
  f <- file.path(MBX_OUT_DIR, "9_visualization_entire", "mbx_ezviz_info.txt")
  d <- file.path(MBX_OUT_DIR, "9_visualization_entire")
  if (!dir.exists(d)) return(NULL)
  kv <- if (file.exists(f)) read_kv(f) else list()
  pngs <- list.files(d, pattern = "\\.png$", full.names = TRUE, recursive = TRUE)
  list(n_pngs = length(pngs), pngs = pngs, raw_kv = kv)
}

# ---- Step 10: ezstat ----
extract_step10 <- function() {
  f <- file.path(MBX_OUT_DIR, "10_stats", "mbx_ezstat_info.txt")
  d <- file.path(MBX_OUT_DIR, "10_stats")
  if (!dir.exists(d)) return(NULL)
  kv <- if (file.exists(f)) read_kv(f) else list()
  xlsx <- list.files(d, pattern = "ezstat_KW.*\\.xlsx$", full.names = TRUE,
                     recursive = TRUE)
  list(n_kw_files = length(xlsx), kw_files = xlsx, raw_kv = kv)
}

# ---- Step 11: pre-diversity ----
extract_step11 <- function() {
  f <- file.path(MBX_OUT_DIR, "11_pre_diversity", "mbx_pre_diversity_info.txt")
  if (!file.exists(f)) return(NULL)
  kv <- read_kv(f)
  d <- file.path(MBX_OUT_DIR, "11_pre_diversity")
  rare_png <- file.path(d, "alpha_rarefaction_curves.png")
  list(
    rec_depth     = kv$RECOMMENDED_DEPTH %||% kv$SAMPLING_DEPTH %||% "(not set)",
    median_depth  = kv$MEDIAN_DEPTH %||% "",
    samples_kept  = kv$SAMPLES_KEPT %||% kv$N_SAMPLES_KEPT %||% "",
    samples_lost  = kv$SAMPLES_LOST %||% kv$N_SAMPLES_LOST %||% "0",
    rooted_tree   = kv$ROOTED_TREE %||% file.path(d, "rooted-tree.qza"),
    rarefaction_png = if (file.exists(rare_png)) rare_png else "",
    rare_csv      = file.path(d, "alpha_rarefaction_data.csv"),
    sd_candidates = file.path(d, "sampling_depth_candidates.csv"),
    raw_kv        = kv
  )
}

# ---- Step 12: alpha diversity ----
extract_step12 <- function() {
  d <- file.path(MBX_OUT_DIR, "12_alpha_diversity_results")
  if (!dir.exists(d)) return(NULL)
  alpha_xlsx <- file.path(d, "alpha_diversity.xlsx")
  alpha_df <- if (file.exists(alpha_xlsx))
    tryCatch(read.xlsx(alpha_xlsx, check.names = FALSE),
             error = function(e) NULL) else NULL
  bp_dirs <- list.dirs(file.path(d, "boxplots_for_alpha_diversity"),
                       recursive = FALSE)
  st_dirs <- list.dirs(file.path(d, "stats_for_alpha_diversity"),
                       recursive = FALSE)
  list(
    alpha_xlsx = alpha_xlsx,
    alpha_df = alpha_df,
    boxplot_dirs = bp_dirs,
    stats_dirs = st_dirs,
    n_categorical = length(bp_dirs)
  )
}

# ---- Step 13: beta diversity ----
extract_step13 <- function() {
  d <- file.path(MBX_OUT_DIR, "13_beta_diversity_results")
  if (!dir.exists(d)) return(NULL)
  cat_dir <- file.path(d, "results_by_categorical_variables")
  cats <- if (dir.exists(cat_dir)) list.dirs(cat_dir, recursive = FALSE) else c()
  all_dir <- file.path(d, "all_samples_beta_diversity")
  list(
    cat_dirs = cats,
    n_cats = length(cats),
    all_samples_dir = if (dir.exists(all_dir)) all_dir else "",
    permanova_files = list.files(d, pattern = "permanova", full.names = TRUE,
                                 recursive = TRUE, ignore.case = TRUE)
  )
}

# ---- Step 14: ANCOMBC2 ----
extract_step14 <- function() {
  f <- file.path(MBX_OUT_DIR, "14_differential_abundance_ANCOMBC2", "mbx_ancombc2_info.txt")
  d <- file.path(MBX_OUT_DIR, "14_differential_abundance_ANCOMBC2")
  if (!dir.exists(d)) return(NULL)
  kv <- if (file.exists(f)) read_kv(f) else list()
  level_dirs <- list.dirs(d, recursive = FALSE)
  level_dirs <- level_dirs[grepl("^ANCOMBC2_", basename(level_dirs))]
  # ANCOMBC2 writes per-(level x variable) results as ancombc2_primary_*.xlsx
  # plus optional ancombc2_pairwise_*.xlsx, ancombc2_global_*.xlsx etc.
  # The "primary" files contain the main per-coefficient table that we use
  # for sig taxa extraction.
  da_files <- list.files(d, pattern = "ancombc2_primary_.*\\.xlsx$",
                         recursive = TRUE, full.names = TRUE)
  # Top-level master roll-up
  master_summary <- file.path(d, "Summary_all_levels_all_variables.xlsx")
  list(level_dirs = level_dirs, n_levels = length(level_dirs),
       da_files = da_files, raw_kv = kv,
       master_summary = if (file.exists(master_summary)) master_summary else "")
}

# ---- Step 15: PICRUSt2 ----
extract_step15 <- function() {
  f <- file.path(MBX_OUT_DIR, "15_picrust2", "mbx_picrust2_info.txt")
  d <- file.path(MBX_OUT_DIR, "15_picrust2")
  if (!dir.exists(d)) return(NULL)
  kv <- if (file.exists(f)) read_kv(f) else list()
  cog_skipped <- file.exists(file.path(d, "working_dir_picrust2",
                                       "COG_pred_metagenome_unstrat.tsv.gz.skipped")) ||
                 isTRUE(grepl("yes|true", kv$COG_SKIPPED %||% "no",
                              ignore.case = TRUE))
  per_var_dirs <- list.dirs(d, recursive = FALSE)
  per_var_dirs <- per_var_dirs[grepl("^picrust2_", basename(per_var_dirs))]
  list(
    nsti_threshold = kv$NSTI_THRESHOLD %||% "2",
    cog_skipped = cog_skipped,
    per_var_dirs = per_var_dirs,
    n_per_var = length(per_var_dirs),
    raw_kv = kv,
    html_report = file.path(d, "picrust2_report.html")
  )
}

# ---- Step 16: ML biomarkers ----
extract_step16 <- function() {
  f <- file.path(MBX_OUT_DIR, "16_ml_biomarkers", "mbx_ml_biomarkers_info.txt")
  d <- file.path(MBX_OUT_DIR, "16_ml_biomarkers")
  if (!dir.exists(d)) return(NULL)
  kv <- if (file.exists(f)) read_kv(f) else list()
  per_var_dirs <- list.dirs(d, recursive = FALSE)
  per_var_dirs <- per_var_dirs[!grepl("working_dir", basename(per_var_dirs))]
  summary_files <- list.files(d, pattern = "Summary_RF_.*\\.xlsx$",
                              recursive = TRUE, full.names = TRUE)
  metric_files <- list.files(d, pattern = "model_metrics\\.xlsx$",
                             recursive = TRUE, full.names = TRUE)
  list(per_var_dirs = per_var_dirs, n_per_var = length(per_var_dirs),
       summary_files = summary_files, metric_files = metric_files,
       raw_kv = kv)
}

# ---- Step 17: networks ----
extract_step17 <- function() {
  f <- file.path(MBX_OUT_DIR, "17_co_occurrence_networks", "mbx_networks_info.txt")
  d <- file.path(MBX_OUT_DIR, "17_co_occurrence_networks")
  if (!dir.exists(d)) return(NULL)
  kv <- if (file.exists(f)) read_kv(f) else list()
  global_dirs <- list.dirs(file.path(d, "global_networks"), recursive = FALSE)
  per_var_dirs <- list.dirs(d, recursive = FALSE)
  per_var_dirs <- per_var_dirs[!basename(per_var_dirs) %in%
                               c("global_networks", "working_dir_networks")]
  group_cmps <- list.files(d, pattern = "group_comparison\\.xlsx$",
                           recursive = TRUE, full.names = TRUE)
  hub_files <- list.files(d, pattern = "hub_taxa\\.xlsx$",
                          recursive = TRUE, full.names = TRUE)
  list(global_dirs = global_dirs, n_global = length(global_dirs),
       per_var_dirs = per_var_dirs, n_per_var = length(per_var_dirs),
       group_cmps = group_cmps, hub_files = hub_files, raw_kv = kv)
}

# ────────────────────────────────────────────────────────────────────────────
# Cross-step convergence: taxa simultaneously DA + RF-important + network hubs
# ────────────────────────────────────────────────────────────────────────────
build_convergence <- function(s14, s16, s17) {
  if (is.null(s14) && is.null(s16) && is.null(s17)) return(NULL)
  da_taxa <- character(0); rf_taxa <- character(0); hub_taxa <- character(0)

  if (!is.null(s14) && length(s14$da_files) > 0) {
    for (f in s14$da_files) {
      df <- tryCatch(read.xlsx(f), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next
      taxon_col <- intersect(c("taxon", "Taxon", "feature", "Feature",
                               "ASV_ID", "id", "ID"), names(df))
      if (length(taxon_col) == 0) next
      # ANCOMBC2 columns are q_<level> per coefficient.  Skip the intercept.
      q_cols <- grep("^q[_.]", names(df), value = TRUE)
      q_cols <- q_cols[!grepl("intercept", q_cols, ignore.case = TRUE)]
      # Fall back to single fixed-name columns used by other DA tools
      if (length(q_cols) == 0)
        q_cols <- intersect(c("q_val", "qval", "q.value", "padj",
                              "p_adjust", "FDR"), names(df))
      if (length(q_cols) == 0) next
      for (qc in q_cols) {
        qv <- suppressWarnings(as.numeric(df[[qc]]))
        sig <- df[!is.na(qv) & qv < 0.05, , drop = FALSE]
        if (nrow(sig) > 0)
          da_taxa <- c(da_taxa, as.character(sig[[taxon_col[1]]]))
      }
    }
  }

  if (!is.null(s16) && length(s16$metric_files) > 0) {
    fi_files <- list.files(file.path(MBX_OUT_DIR, "16_ml_biomarkers"),
                           pattern = "feature_importance\\.xlsx$",
                           recursive = TRUE, full.names = TRUE)
    for (f in fi_files) {
      df <- tryCatch(read.xlsx(f), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next
      taxon_col <- intersect(c("feature", "Feature", "taxon", "Taxon"),
                             names(df))
      imp_col   <- intersect(c("importance", "Importance", "permutation_importance"),
                             names(df))
      if (length(taxon_col) == 0 || length(imp_col) == 0) next
      df <- df[order(-df[[imp_col[1]]]), , drop = FALSE]
      rf_taxa <- c(rf_taxa, as.character(head(df[[taxon_col[1]]], 20)))
    }
  }

  if (!is.null(s17) && length(s17$hub_files) > 0) {
    for (f in s17$hub_files) {
      df <- tryCatch(read.xlsx(f), error = function(e) NULL)
      if (is.null(df) || nrow(df) == 0) next
      taxon_col <- intersect(c("taxon", "Taxon", "feature", "Feature"),
                             names(df))
      if (length(taxon_col) == 0) next
      hub_taxa <- c(hub_taxa, as.character(df[[taxon_col[1]]]))
    }
  }

  da_taxa  <- unique(da_taxa[nzchar(da_taxa) & !is.na(da_taxa)])
  rf_taxa  <- unique(rf_taxa[nzchar(rf_taxa) & !is.na(rf_taxa)])
  hub_taxa <- unique(hub_taxa[nzchar(hub_taxa) & !is.na(hub_taxa)])

  # Normalise taxon names for matching across methods.
  # ANCOMBC2 returns full lineage strings like
  #   d__Bacteria;p__...;g__Angelakisella   (or with ;s__species at the end)
  # RF / network return leaf names like "Angelakisella" or "Akkermansia.muciniphila".
  # We extract the LEAF (last ';' segment, stripping the GreenGenes2 rank prefix
  # such as g__, s__) so all three methods can be matched on equal footing.
  leafify <- function(x) {
    parts <- strsplit(x, ";", fixed = TRUE)
    sapply(parts, function(p) {
      if (length(p) == 0) return(NA_character_)
      last <- trimws(tail(p, 1))
      sub("^[a-z]__", "", last)
    })
  }
  norm <- function(x) {
    leaf <- leafify(x)
    tolower(gsub("[^a-z0-9]+", "_", tolower(leaf)))
  }
  da_n  <- norm(da_taxa);  names(da_n)  <- da_taxa
  rf_n  <- norm(rf_taxa);  names(rf_n)  <- rf_taxa
  hub_n <- norm(hub_taxa); names(hub_n) <- hub_taxa

  all_n <- unique(c(da_n, rf_n, hub_n))
  rows <- lapply(all_n, function(nn) {
    in_da  <- nn %in% da_n
    in_rf  <- nn %in% rf_n
    in_hub <- nn %in% hub_n
    score  <- as.integer(in_da) + as.integer(in_rf) + as.integer(in_hub)
    if (score < 2) return(NULL)
    # Prefer the SHORTEST original label across the three sources
    # (typically the leaf-name from RF / network rather than the full
    # ANCOMBC2 lineage string).
    candidates <- c(names(rf_n)[rf_n == nn],
                    names(hub_n)[hub_n == nn],
                    names(da_n)[da_n == nn])
    candidates <- unique(candidates)
    label <- candidates[which.min(nchar(candidates))]
    data.frame(
      taxon = label,
      ANCOMBC2_significant = ifelse(in_da, "yes", ""),
      RF_top_importance    = ifelse(in_rf, "yes", ""),
      Network_hub          = ifelse(in_hub, "yes", ""),
      convergence_score    = score,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  out <- out[order(-out$convergence_score, out$taxon), , drop = FALSE]
  out
}

# ────────────────────────────────────────────────────────────────────────────
# Build each step block
# ────────────────────────────────────────────────────────────────────────────
status_for <- function(id) {
  s <- discovery$status[discovery$id == id]
  if (length(s) == 0) "MISSING" else s
}
desc_for <- function(id) discovery$desc[discovery$id == id]
script_for <- function(id) discovery$script[discovery$id == id]

step_header <- function(id, title) {
  st <- status_for(id)
  cls <- tolower(st)
  badge <- tags$span(class = paste("step-id-badge", cls),
                     sprintf("Step %d", id))
  status_label <- switch(st,
    DONE    = tags$span(style = "color:var(--success);font-weight:600",
                        "EXECUTED"),
    PARTIAL = tags$span(style = "color:var(--warning);font-weight:600",
                        "PARTIAL"),
    MISSING = tags$span(style = "color:var(--fg-muted);font-weight:600",
                        "NOT RUN")
  )
  tags$div(class = "step-hdr",
    badge, tags$h3(style = "margin:0;flex:1", title), status_label)
}

step_not_run_box <- function(id) {
  scr <- script_for(id)
  desc <- desc_for(id)
  tags$div(class = "box info",
    tags$p(class = "box-title", "Step not executed in this output directory"),
    tags$p(sprintf("This step (%s) was not run for this analysis. To add it to the report, execute:", desc)),
    tags$pre(sprintf("%s %s", scr, MBX_OUT_DIR)),
    tags$p("...then re-run mbx_final_report.sh to refresh this report.")
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Quick-glance dashboard (top of report)
# ────────────────────────────────────────────────────────────────────────────
dashboard_card <- function(num, lbl) {
  tags$div(class = "card",
    tags$span(class = "num", num),
    tags$span(class = "lbl", lbl))
}
build_dashboard <- function(s1, s4, s7, s8, s11, s14, s16) {
  n_samples <- if (!is.null(s1)) s1$n_samples else "?"
  n_features <- "?"
  if (!is.null(s8) && length(s8$per_level) > 0) {
    g <- Filter(function(x) grepl("genera", x$dir, ignore.case = TRUE), s8$per_level)
    if (length(g) > 0) n_features <- g[[1]]$n_features
  }
  rec_depth <- if (!is.null(s11)) s11$rec_depth else "?"
  # n_da = total significant DA taxa across all levels x variables (from
  # the master ANCOMBC2 summary).  Falls back to file count if summary missing.
  n_da <- "?"
  if (!is.null(s14)) {
    if (nzchar(s14$master_summary) && file.exists(s14$master_summary)) {
      mt <- tryCatch(read.xlsx(s14$master_summary, sheet = 1),
                     error = function(e) NULL)
      if (!is.null(mt) && "n_significant_global" %in% names(mt)) {
        n_da <- sum(suppressWarnings(as.numeric(mt$n_significant_global)),
                    na.rm = TRUE) +
                sum(suppressWarnings(as.numeric(mt$n_significant_pairwise)),
                    na.rm = TRUE)
      }
    } else {
      n_da <- length(s14$da_files)
    }
  }
  best_auc <- "?"
  if (!is.null(s16) && length(s16$summary_files) > 0) {
    aucs <- numeric(0)
    for (f in s16$summary_files) {
      df <- tryCatch(read.xlsx(f), error = function(e) NULL)
      if (!is.null(df) && "macro_AUC" %in% names(df)) {
        aucs <- c(aucs, suppressWarnings(as.numeric(df$macro_AUC)))
      } else if (!is.null(df) && "AUC" %in% names(df)) {
        aucs <- c(aucs, suppressWarnings(as.numeric(df$AUC)))
      }
    }
    aucs <- aucs[!is.na(aucs)]
    if (length(aucs) > 0) best_auc <- sprintf("%.2f", max(aucs))
  }
  tags$div(class = "dashboard",
    dashboard_card(n_samples, "samples"),
    dashboard_card(n_features, "genus-level taxa"),
    dashboard_card(rec_depth, "rarefaction depth"),
    dashboard_card(n_da, "ANCOMBC2 DA findings"),
    dashboard_card(best_auc, "best ML AUC"))
}

# ────────────────────────────────────────────────────────────────────────────
# Section builders
# ────────────────────────────────────────────────────────────────────────────

section_step0 <- function() {
  s <- extract_step0()
  if (is.null(s)) return(tagList(step_header(0, "Primer detection"), step_not_run_box(0)))
  status <- if (nzchar(s$detection_status)) s$detection_status else "(unknown)"
  status_explainer <- switch(status,
    "DETECTED"      = "Primers were directly observed in the raw reads.  Forward and reverse sequences are listed below; their lengths drive the DADA2 --p-trim-left-f / --p-trim-left-r values.",
    "USER_SUPPLIED" = "Primers were supplied manually via --forward-primer / --reverse-primer.  Their lengths drive the DADA2 --p-trim-left-f / --p-trim-left-r values.",
    "TRIMMED"       = "The reads start with a conserved 16S anchor motif but no primer is present, indicating the sequencing facility removed primers before delivery.  --p-trim-left-f and --p-trim-left-r are therefore set to 0 and the classifier will run in full-length mode.",
    "UNKNOWN"       = "Detection failed: no primer matched above threshold AND no V-region anchor was found.  --p-trim-left-f and --p-trim-left-r default to 20 (defensive) and the classifier will run in full-length mode.",
    "Detection ran but did not write a recognised status -- assuming UNKNOWN behaviour."
  )
  tagList(
    step_header(0, "Primer detection"),
    tags$p("The pipeline auto-detects forward and reverse primers from the raw FASTQ files using a 3-tier IUPAC-aware engine: (1) direct primer match across all orientations; (2) V-region anchor motif detection that recognises pre-trimmed primers and infers the original V-region; (3) rich failure diagnostics (consensus reads, base composition, top 5'-prefixes)."),
    tags$p(HTML(paste0("<b>Detection status:</b> ", status, ".  ", status_explainer))),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Detection status:"),
              tags$span(class = "val", status)),
      tags$li(tags$span(class = "key", "Confidence level:"),
              tags$span(class = "val", if (nzchar(s$confidence_level)) s$confidence_level else "(n/a)")),
      tags$li(tags$span(class = "key", "Inferred region:"),
              tags$span(class = "val", if (nzchar(s$inferred_region)) s$inferred_region else "(n/a)")),
      tags$li(tags$span(class = "key", "Forward primer (5'->3'):"),
              tags$span(class = "val", s$fwd, sprintf(" (%s nt)", s$fwd_len))),
      tags$li(tags$span(class = "key", "Reverse primer (5'->3'):"),
              tags$span(class = "val", s$rev, sprintf(" (%s nt)", s$rev_len))),
      tags$li(tags$span(class = "key", "Detected orientation:"),
              tags$span(class = "val", s$orientation)),
      tags$li(tags$span(class = "key", "Auto-detected:"),
              tags$span(class = "val", s$auto_detected)),
      tags$li(tags$span(class = "key", "FASTQ source directory:"),
              tags$span(class = "val", relpath(s$fastq_dir))),
      if (nzchar(s$detection_note))
        tags$li(tags$span(class = "key", "Detection note:"),
                tags$span(class = "val", s$detection_note))
    )
  )
}

section_step12 <- function() {
  s1 <- extract_step1(); s2_dir <- file.path(MBX_OUT_DIR, "2_first_artifact_file")
  if (is.null(s1)) return(tagList(step_header(1, "QIIME2 manifest"), step_not_run_box(1)))
  tagList(
    step_header(1, "QIIME2 manifest construction"),
    tags$p("QIIME2's tab-separated manifest format with absolute file paths.  Sample IDs are extracted via regex matching s(ample)?[-_]?[0-9]+ patterns followed by filename-cleanup fallbacks."),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Number of samples:"),
              tags$span(class = "val", s1$n_samples)),
      tags$li(tags$span(class = "key", "Sequencing layout:"),
              tags$span(class = "val", if (s1$paired) "paired-end" else "single-end")),
      tags$li(tags$span(class = "key", "Manifest path:"),
              tags$span(class = "val", relpath(s1$manifest_path))),
      tags$li(tags$span(class = "key", "Columns:"),
              tags$span(class = "val", s1$columns))
    ),
    tags$h3("Step 2 - QIIME2 .qza artifact creation"),
    if (status_for(2) != "DONE") step_not_run_box(2) else
    tagList(
      tags$p("The manifest was imported into QIIME2 as a SampleData[PairedEndSequencesWithQuality] artifact (.qza) and a corresponding visualisation (.qzv) was rendered for quality assessment."),
      tags$ul(class = "kv-list",
        tags$li(tags$span(class = "key", "Artifact:"),
                tags$span(class = "val", relpath(file.path(s2_dir, "Paired_End_artifact.qza")))),
        tags$li(tags$span(class = "key", "Visualisation:"),
                tags$span(class = "val", relpath(file.path(s2_dir, "Paired_End_artifact.qzv"))))
      )
    )
  )
}

section_step34 <- function() {
  s3 <- extract_step3(); s4 <- extract_step4()
  out <- list()
  out[[length(out) + 1]] <- step_header(3, "DADA2 truncation-length selection")
  if (is.null(s3)) {
    out[[length(out) + 1]] <- step_not_run_box(3)
  } else {
    out[[length(out) + 1]] <- tags$p(HTML(paste(
      "Truncation lengths are chosen automatically from the demux summary using the standard QIIME2 heuristic:",
      "the longest position at which the 25th-percentile quality score remains >=25, separately for forward and reverse reads.",
      "Trim-left is decided from the step-0 DETECTION_STATUS:",
      "<code>DETECTED</code> / <code>USER_SUPPLIED</code> -> primer length;",
      "<code>TRIMMED</code> -> 0 (primers already removed by the sequencing facility);",
      "<code>UNKNOWN</code> -> 20 (defensive default)."
    )))
    out[[length(out) + 1]] <- tags$h4("DADA2 parameter file (verbatim)")
    out[[length(out) + 1]] <- tags$pre(s3$raw)
  }
  out[[length(out) + 1]] <- tags$h3("Step 4 - DADA2 denoising + ASV inference")
  if (is.null(s4) || status_for(4) != "DONE") {
    out[[length(out) + 1]] <- step_not_run_box(4)
  } else {
    out[[length(out) + 1]] <- tags$p("DADA2 (Callahan et al. 2016) infers exact amplicon sequence variants (ASVs) using a parametric error model fit per-run, performs paired-end merging, and removes chimeras via consensus.  Outputs include a feature table, representative sequences, and per-sample stats.")
    out[[length(out) + 1]] <- tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Feature table:"),
              tags$span(class = "val", relpath(s4$feature_table))),
      tags$li(tags$span(class = "key", "Representative sequences:"),
              tags$span(class = "val", relpath(s4$rep_seqs))),
      tags$li(tags$span(class = "key", "DADA2 stats:"),
              tags$span(class = "val", relpath(s4$stats)))
    )
  }
  do.call(tagList, out)
}

section_step56 <- function() {
  s <- extract_step56()
  out <- list(step_header(5, "Taxonomy classifier preparation + run"))
  if (is.null(s)) {
    out[[length(out) + 1]] <- step_not_run_box(5)
  } else {
    is_full   <- identical(s$classifier_mode, "full-length")
    src       <- s$classifier_source
    is_zenodo <- identical(src, "zenodo")
    is_cached <- identical(src, "cached")
    is_fb     <- identical(src, "local-training-fallback")

    if (is_full && (is_zenodo || is_cached)) {
      out[[length(out) + 1]] <- tags$p(HTML(paste(
        "Greengenes2 (McDonald et al. 2024) was used at release",
        sprintf("<b>%s</b>.", s$gg2_version),
        "<b>Classifier mode: full-length.</b>",
        "Primer sequences were not available in the reads (the sequencing facility had already trimmed them, or primer detection failed),",
        "so the Naive-Bayes classifier was trained directly on the FULL Greengenes2 backbone -- no V-region extraction step was performed.",
        sprintf("<b>Classifier source: %s.</b>", src),
        if (is_zenodo) sprintf(
          "Instead of training the classifier locally (which takes 30-90 minutes), a pre-trained, sha256-verified Naive-Bayes .qza was downloaded from <a href='%s'>%s</a>; the file matches QIIME2 %s (sklearn family %s).",
          s$zenodo_url, s$zenodo_url, s$zenodo_qiime2, s$sklearn_family
        ) else
          "A pre-trained classifier from a previous run was reused (cached on disk).",
        "Genus-level accuracy is virtually identical to region-specific training; species-level resolution is marginally lower (negligible for typical 250 bp Illumina reads)."
      )))
    } else if (is_full && is_fb) {
      out[[length(out) + 1]] <- tags$p(HTML(paste(
        "Greengenes2 (McDonald et al. 2024) was used at release",
        sprintf("<b>%s</b>.", s$gg2_version),
        "<b>Classifier mode: full-length.</b>",
        "The pipeline first attempted a pre-trained Zenodo classifier",
        sprintf("(%s), but classify-sklearn rejected it (most likely a subtle scikit-learn pickle incompatibility for sklearn family %s).",
                s$zenodo_filename, s$sklearn_family),
        "The Zenodo .qza was deleted and the classifier was retrained <i>locally</i> on the full GG2 backbone -- exactly as in pre-1.2 versions of mbX Pro.",
        "<b>The pipeline did NOT abort.</b>  This is the documented fall-back path."
      )))
    } else if (is_full) {
      out[[length(out) + 1]] <- tags$p(HTML(paste(
        "Greengenes2 (McDonald et al. 2024) was used at release",
        sprintf("<b>%s</b>.", s$gg2_version),
        "<b>Classifier mode: full-length.</b>",
        "Primer sequences were not available, so the Naive-Bayes classifier was trained directly on the FULL Greengenes2 backbone -- no V-region extraction step was performed.",
        "No compatible pre-trained classifier was available on Zenodo (or Zenodo was unreachable), so the classifier was trained locally from the GG2 backbone .fna.qza + .tax.qza.",
        "Genus-level accuracy is virtually identical to region-specific training; species-level resolution is marginally lower."
      )))
    } else {
      out[[length(out) + 1]] <- tags$p(HTML(paste(
        "Greengenes2 (McDonald et al. 2024) was used at release",
        sprintf("<b>%s</b>.", s$gg2_version),
        "<b>Classifier mode: region-specific.</b>",
        "The full-length 16S backbone was downloaded, region-extracted to the user's primer pair using <code>qiime feature-classifier extract-reads</code>,",
        "and a Naive-Bayes classifier was trained locally on the trimmed reference.",
        "Pre-trained Zenodo classifiers cover full-length only and are therefore not used in region-specific mode.",
        "All three artifacts (extracted reads, trained classifier, taxonomy) are cached for subsequent runs."
      )))
    }

    kv_items <- list(
      tags$li(tags$span(class = "key", "Greengenes2 release:"),
              tags$span(class = "val", s$gg2_version)),
      tags$li(tags$span(class = "key", "Classifier mode:"),
              tags$span(class = "val", s$classifier_mode)),
      tags$li(tags$span(class = "key", "Classifier source:"),
              tags$span(class = "val", src)),
      tags$li(tags$span(class = "key", "scikit-learn version:"),
              tags$span(class = "val", sprintf("%s (family %s)", s$sklearn_ver, s$sklearn_family))),
      tags$li(tags$span(class = "key", "Forward primer:"),
              tags$span(class = "val", s$fwd_primer)),
      tags$li(tags$span(class = "key", "Reverse primer:"),
              tags$span(class = "val", s$rev_primer)),
      tags$li(tags$span(class = "key", "Target V-region:"),
              tags$span(class = "val", s$target_region)),
      tags$li(tags$span(class = "key", "Reference (full-length .fna.qza):"),
              tags$span(class = "val", relpath(s$gg2_full_path))),
      tags$li(tags$span(class = "key", "Reference taxonomy (.qza):"),
              tags$span(class = "val", relpath(s$gg2_tax_path))),
      tags$li(tags$span(class = "key",
                        if (is_full)
                          "ASV length range (informational, not applied):"
                        else
                          "ASV length range used for extract-reads:"),
              tags$span(class = "val", sprintf("min=%s, max=%s", s$asv_min_len, s$asv_max_len))),
      tags$li(tags$span(class = "key", "Trained classifier (.qza):"),
              tags$span(class = "val", relpath(s$classifier_qza)))
    )

    if (is_zenodo || is_fb) {
      kv_items <- c(kv_items, list(
        tags$li(tags$span(class = "key", "Zenodo record:"),
                tags$span(class = "val", s$zenodo_url)),
        tags$li(tags$span(class = "key", "Zenodo filename:"),
                tags$span(class = "val", s$zenodo_filename)),
        tags$li(tags$span(class = "key", "Trained for QIIME2:"),
                tags$span(class = "val", s$zenodo_qiime2)),
        tags$li(tags$span(class = "key", "Trained on GG2:"),
                tags$span(class = "val", s$zenodo_gg2)),
        tags$li(tags$span(class = "key", "SHA-256 (expected):"),
                tags$span(class = "val", s$zenodo_sha_exp)),
        tags$li(tags$span(class = "key", "SHA-256 (actual):"),
                tags$span(class = "val", s$zenodo_sha_act)),
        tags$li(tags$span(class = "key", "Zenodo outcome:"),
                tags$span(class = "val", s$zenodo_note))
      ))
    }

    out[[length(out) + 1]] <- do.call(tags$ul, c(list(class = "kv-list"), kv_items))
  }
  out[[length(out) + 1]] <- tags$h3("Step 6 - Taxonomy assignment")
  if (status_for(6) != "DONE") {
    out[[length(out) + 1]] <- step_not_run_box(6)
  } else {
    out[[length(out) + 1]] <- tags$p("ASVs were classified against the trained Greengenes2 Naive-Bayes classifier using qiime feature-classifier classify-sklearn.  This is the standard high-precision classifier in QIIME2.")
    out[[length(out) + 1]] <- tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Taxonomy artifact:"),
              tags$span(class = "val", relpath(s$taxonomy_qza)))
    )
  }
  do.call(tagList, out)
}

section_step7 <- function() {
  s <- extract_step7()
  if (is.null(s)) return(tagList(step_header(7, "Taxonomy CSVs"), step_not_run_box(7)))
  tagList(
    step_header(7, "Taxonomy CSVs + mitochondria/chloroplast filtering"),
    tags$p("Mitochondrial and chloroplast sequences (eukaryotic 16S contamination) are removed from the feature table before taxonomy CSVs are exported.  All seven taxonomic levels (domain through species) are written as separate CSVs from the QIIME2 taxa-bar-plots .qzv export."),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Metadata file:"),
              tags$span(class = "val", relpath(s$metadata_txt))),
      tags$li(tags$span(class = "key", "Mito/chloro-filtered table:"),
              tags$span(class = "val", relpath(s$filtered_table))),
      tags$li(tags$span(class = "key", "Bar-plot QZV:"),
              tags$span(class = "val", relpath(s$barplot_qzv))),
      tags$li(tags$span(class = "key", "Number of level CSVs:"),
              tags$span(class = "val", s$n_levels))
    )
  )
}

section_step8 <- function() {
  s <- extract_step8()
  if (is.null(s)) return(tagList(step_header(8, "Taxonomy cleaning (mbX ezclean)"), step_not_run_box(8)))
  per_level_df <- do.call(rbind, lapply(s$per_level, function(x) {
    data.frame(level_dir = x$dir, n_samples = x$n_samples,
               n_features = x$n_features, stringsAsFactors = FALSE)
  }))
  tagList(
    step_header(8, "Taxonomy cleaning per level (mbX ezclean)"),
    tags$p("The mbX R package's ezclean() function consolidates synonymous taxa, drops unassigned features, and writes a tidy wide-format Excel sheet per taxonomic level.  All 7 levels (domain/phylum/class/order/family/genus/species) are processed independently so a failure at one level does not abort the others."),
    df_to_html(per_level_df, "Cleaned features per taxonomic level")
  )
}

section_step9 <- function() {
  s <- extract_step9()
  if (is.null(s) || s$n_pngs == 0)
    return(tagList(step_header(9, "Stacked-bar visualisations (ezviz)"), step_not_run_box(9)))
  # Embed up to 8 representative pngs (one per categorical variable, all levels)
  show <- head(s$pngs, 8)
  figs <- tagList(lapply(show, function(p) embed_image(p, basename(p))))
  tagList(
    step_header(9, "Stacked-bar visualisations (ezviz)"),
    tags$p(sprintf("mbX::ezviz() produced %d stacked-bar PNGs covering every (taxonomic level x categorical variable) combination.  Showing up to 8 representative figures here; the full set is in 9_visualization_entire/.", s$n_pngs)),
    tags$div(class = "fig-grid", figs)
  )
}

section_step10 <- function() {
  s <- extract_step10()
  if (is.null(s) || s$n_kw_files == 0)
    return(tagList(step_header(10, "Statistics (ezstat)"), step_not_run_box(10)))
  # Build a roll-up summary: number of significant taxa per (level x variable)
  rows <- list()
  for (f in s$kw_files) {
    df <- tryCatch(read.xlsx(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    qcol <- intersect(c("padj","p_adjust","FDR","q_val","BH","BH_adj"),
                      names(df))
    n_sig <- if (length(qcol) > 0)
      sum(suppressWarnings(as.numeric(df[[qcol[1]]])) < 0.05, na.rm = TRUE)
    else NA_integer_
    rows[[length(rows) + 1]] <- data.frame(
      file = relpath(f),
      n_taxa_tested = nrow(df), n_significant_q05 = n_sig,
      stringsAsFactors = FALSE)
  }
  tbl <- do.call(rbind, rows)
  tagList(
    step_header(10, "Kruskal-Wallis + Dunn's post-hoc statistics (ezstat)"),
    tags$p("mbX::ezstat() runs Kruskal-Wallis tests for each taxon against each categorical metadata variable, with Dunn's post-hoc pairwise comparisons (BH-adjusted) and compact-letter-display assignments for boxplots.  All 7 taxonomic levels x N categorical variables are tested."),
    df_to_html(tbl, "Per-file roll-up: taxa tested and significant (q < 0.05)",
               max_rows = 30)
  )
}

section_step11 <- function() {
  s <- extract_step11()
  if (is.null(s)) return(tagList(step_header(11, "Pre-diversity"), step_not_run_box(11)))
  tagList(
    step_header(11, "Phylogenetic tree + sampling depth (pre-diversity)"),
    tags$p(HTML(paste(
      "A rooted phylogenetic tree was built via mafft-fasttree (parallelised across CPU cores).",
      "Sampling depth is selected by an analytical, group-aware algorithm: the highest depth at which",
      ">=90 % of all samples pass AND >=80 % per group; falls back to >=75 % overall, then to Q1.",
      "Rarefaction curves are computed with the analytical Hurlbert (1971) formula in log-space",
      "(numerically stable via scipy.special.gammaln), parallelised across samples."
    ))),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Recommended sampling depth:"),
              tags$span(class = "val", s$rec_depth)),
      tags$li(tags$span(class = "key", "Median sample depth:"),
              tags$span(class = "val", s$median_depth)),
      tags$li(tags$span(class = "key", "Samples retained at this depth:"),
              tags$span(class = "val", s$samples_kept)),
      tags$li(tags$span(class = "key", "Samples dropped at this depth:"),
              tags$span(class = "val", s$samples_lost)),
      tags$li(tags$span(class = "key", "Rooted tree (.qza):"),
              tags$span(class = "val", relpath(s$rooted_tree)))
    ),
    if (nzchar(s$rarefaction_png) && file.exists(s$rarefaction_png))
      embed_image(s$rarefaction_png,
                  "Analytical rarefaction curves (Hurlbert 1971) per categorical group")
    else NULL
  )
}

section_step12 <- function() {
  s <- extract_step12()
  if (is.null(s)) return(tagList(step_header(12, "Alpha diversity"), step_not_run_box(12)))

  # Embed alpha-diversity table summary
  if (!is.null(s$alpha_df)) {
    nm <- names(s$alpha_df)
    diversity_cols <- intersect(c("ASVs_or_Features", "Shannon_Index",
                                  "Simpson_Diversity",
                                  "Faith_Phylogenetic_Diversity",
                                  "Pielou_Evenness"), nm)
    sumtbl <- if (length(diversity_cols) > 0) {
      tmp <- s$alpha_df[, diversity_cols, drop = FALSE]
      data.frame(
        metric = diversity_cols,
        n      = sapply(tmp, function(x) sum(!is.na(x))),
        mean   = sapply(tmp, function(x) round(mean(suppressWarnings(as.numeric(x)),
                                                    na.rm = TRUE), 3)),
        median = sapply(tmp, function(x) round(median(suppressWarnings(as.numeric(x)),
                                                      na.rm = TRUE), 3)),
        sd     = sapply(tmp, function(x) round(sd(suppressWarnings(as.numeric(x)),
                                                  na.rm = TRUE), 3)),
        min    = sapply(tmp, function(x) round(min(suppressWarnings(as.numeric(x)),
                                                   na.rm = TRUE), 3)),
        max    = sapply(tmp, function(x) round(max(suppressWarnings(as.numeric(x)),
                                                   na.rm = TRUE), 3)),
        stringsAsFactors = FALSE
      )
    } else NULL
  } else sumtbl <- NULL

  # Find a few boxplot PNGs to embed
  bp_pngs <- character(0)
  for (bp_dir in s$boxplot_dirs) {
    bp_pngs <- c(bp_pngs, list.files(bp_dir, pattern = "\\.png$",
                                      full.names = TRUE)[1:2])
  }
  bp_pngs <- bp_pngs[!is.na(bp_pngs) & file.exists(bp_pngs)]

  tagList(
    step_header(12, "Alpha diversity"),
    tags$p("Five complementary alpha-diversity metrics were computed at the chosen rarefaction depth: ASV richness (count of unique ASVs), Shannon entropy (richness + evenness), Simpson diversity (1 - dominance), Faith's phylogenetic diversity (PD; sum of branch lengths in the tree spanned by present ASVs), and Pielou's evenness (Shannon / log richness).  Each metric was tested across every categorical variable using Kruskal-Wallis with Dunn's post-hoc (BH-adjusted)."),
    if (!is.null(sumtbl)) df_to_html(sumtbl, "Summary statistics across all samples"),
    if (length(bp_pngs) > 0)
      tags$div(class = "fig-grid",
               tagList(lapply(head(bp_pngs, 6),
                              function(p) embed_image(p, basename(p)))))
  )
}

section_step13 <- function() {
  s <- extract_step13()
  if (is.null(s)) return(tagList(step_header(13, "Beta diversity"), step_not_run_box(13)))

  # Find a representative PCoA PNG for each categorical variable
  pcoa_pngs <- character(0)
  for (cd in s$cat_dirs) {
    pp <- list.files(cd, pattern = "PCoA.*\\.png$|pcoa.*\\.png$",
                     full.names = TRUE, recursive = TRUE)
    if (length(pp) > 0) pcoa_pngs <- c(pcoa_pngs, pp[1])
  }

  # All-samples figures
  all_figs <- character(0)
  if (nzchar(s$all_samples_dir))
    all_figs <- list.files(s$all_samples_dir, pattern = "\\.png$",
                            full.names = TRUE, recursive = TRUE)

  tagList(
    step_header(13, "Beta diversity (PCoA, PERMANOVA, Adonis)"),
    tags$p("Beta-diversity assesses between-sample compositional dissimilarity.  The pipeline computes Bray-Curtis (abundance-based), Jaccard (presence/absence), Weighted UniFrac (phylogenetically weighted), and Unweighted UniFrac distances, runs PERMANOVA per categorical variable (with pairwise comparisons when the variable has >2 levels and BH adjustment), tests dispersion homogeneity via PERMDISP, and fits a multivariable Adonis model for the all-samples view.  PCoA ordinations are produced for each variable + distance combination."),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Categorical variables analysed:"),
              tags$span(class = "val", s$n_cats)),
      tags$li(tags$span(class = "key", "All-samples directory:"),
              tags$span(class = "val", relpath(s$all_samples_dir)))
    ),
    if (length(pcoa_pngs) > 0)
      tags$div(class = "fig-grid",
               tagList(lapply(pcoa_pngs,
                              function(p) embed_image(p, basename(p))))),
    if (length(all_figs) > 0) tagList(
      tags$h4("All-samples beta-diversity figures"),
      tags$div(class = "fig-grid",
               tagList(lapply(head(all_figs, 4),
                              function(p) embed_image(p, basename(p))))))
  )
}

section_step14 <- function() {
  s <- extract_step14()
  if (is.null(s)) return(tagList(step_header(14, "Differential abundance (ANCOMBC2)"), step_not_run_box(14)))

  # Master roll-up (already produced by ANCOMBC2 step)
  master_tbl <- NULL
  total_sig_global <- 0; total_sig_pw <- 0
  if (nzchar(s$master_summary) && file.exists(s$master_summary)) {
    master_tbl <- tryCatch(read.xlsx(s$master_summary, sheet = 1),
                           error = function(e) NULL)
    if (!is.null(master_tbl)) {
      if ("n_significant_global"   %in% names(master_tbl))
        total_sig_global <- sum(suppressWarnings(
          as.numeric(master_tbl$n_significant_global)), na.rm = TRUE)
      if ("n_significant_pairwise" %in% names(master_tbl))
        total_sig_pw <- sum(suppressWarnings(
          as.numeric(master_tbl$n_significant_pairwise)), na.rm = TRUE)
    }
  }

  # Per-file (level x variable) breakdown using the actual ancombc2_primary_*
  # files: count taxa significant under ANY non-intercept q_<group> column
  rows <- list()
  for (f in s$da_files) {
    df <- tryCatch(read.xlsx(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    q_cols <- grep("^q[_.]", names(df), value = TRUE)
    q_cols <- q_cols[!grepl("intercept", q_cols, ignore.case = TRUE)]
    if (length(q_cols) == 0) next
    sig_taxa <- character(0)
    for (qc in q_cols) {
      qv <- suppressWarnings(as.numeric(df[[qc]]))
      taxon_col <- intersect(c("taxon", "Taxon", "feature", "Feature"),
                             names(df))[1]
      if (is.na(taxon_col)) next
      sig_taxa <- unique(c(sig_taxa,
                           as.character(df[[taxon_col]][!is.na(qv) & qv < 0.05])))
    }
    rows[[length(rows) + 1]] <- data.frame(
      file = relpath(f),
      n_taxa = nrow(df),
      n_unique_sig_taxa = length(sig_taxa),
      stringsAsFactors = FALSE)
  }
  per_file_tbl <- if (length(rows) > 0) do.call(rbind, rows) else NULL

  # Highlight top sig taxa overall (sorted by smallest q across all files)
  top_da <- list()
  for (f in s$da_files) {
    df <- tryCatch(read.xlsx(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    taxon_col <- intersect(c("taxon", "Taxon", "feature", "Feature"),
                           names(df))[1]
    if (is.na(taxon_col)) next
    q_cols   <- grep("^q[_.]", names(df), value = TRUE)
    q_cols   <- q_cols[!grepl("intercept", q_cols, ignore.case = TRUE)]
    lfc_cols <- grep("^lfc[_.]", names(df), value = TRUE)
    lfc_cols <- lfc_cols[!grepl("intercept", lfc_cols, ignore.case = TRUE)]
    if (length(q_cols) == 0) next
    for (qi in seq_along(q_cols)) {
      qc <- q_cols[qi]
      qv <- suppressWarnings(as.numeric(df[[qc]]))
      ix <- which(!is.na(qv) & qv < 0.05)
      if (length(ix) == 0) next
      lc <- if (qi <= length(lfc_cols)) lfc_cols[qi] else NA
      lv <- if (!is.na(lc))
        suppressWarnings(as.numeric(df[[lc]][ix])) else rep(NA_real_, length(ix))
      top_da[[length(top_da) + 1]] <- data.frame(
        source   = relpath(f),
        contrast = sub("^q[_.]", "", qc),
        taxon    = as.character(df[[taxon_col]][ix]),
        q_val    = sprintf("%.3g", qv[ix]),
        LFC      = if (all(is.na(lv))) "n/a" else sprintf("%.2f", lv),
        stringsAsFactors = FALSE)
    }
  }
  top_tbl <- if (length(top_da) > 0) do.call(rbind, top_da) else NULL
  if (!is.null(top_tbl)) {
    top_tbl <- top_tbl[order(suppressWarnings(as.numeric(top_tbl$q_val))), ,
                       drop = FALSE]
  }

  tagList(
    step_header(14, "Differential abundance (ANCOMBC2)"),
    tags$p(HTML(paste(
      "ANCOM-BC2 (Lin et al. 2024) tests for differential abundance while accounting",
      "for compositionality via bias correction.  The pipeline runs every taxonomic",
      "level x every categorical variable with <code>pairwise = TRUE</code> and",
      "<code>global = TRUE</code>; the alphabetically-first level of each variable",
      "is the reference (R modeling default).  This means the reference choice",
      "affects how individual coefficients are <em>labelled</em> but not the",
      "reference-free global tests or the full set of pairwise contrasts."
    ))),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Taxonomic levels analysed:"),
              tags$span(class = "val", s$n_levels)),
      tags$li(tags$span(class = "key", "Per-(level x variable) primary tables:"),
              tags$span(class = "val", length(s$da_files))),
      tags$li(tags$span(class = "key", "Total significant findings (global test):"),
              tags$span(class = "val", total_sig_global)),
      tags$li(tags$span(class = "key", "Total significant findings (pairwise test):"),
              tags$span(class = "val", total_sig_pw))
    ),
    if (!is.null(master_tbl)) tagList(
      tags$h4("Master roll-up (one row per level x variable)"),
      df_to_html(master_tbl, NULL, max_rows = 40)),
    if (!is.null(per_file_tbl)) tagList(
      tags$h4("Per-coefficient breakdown"),
      df_to_html(per_file_tbl, NULL, max_rows = 40)),
    if (!is.null(top_tbl) && nrow(top_tbl) > 0) tagList(
      tags$h4("Top differentially abundant taxa (q < 0.05)"),
      df_to_html(top_tbl, NULL, max_rows = 40))
    else if (total_sig_global == 0 && total_sig_pw == 0) tags$div(class = "box info",
      tags$p(class = "box-title", "No significant DA taxa at q < 0.05 (Holm-adjusted)"),
      tags$p("Possible scientific reasons (any of which may apply):"),
      tags$ul(
        tags$li("Sample size is too small to detect modest effect sizes after Holm correction (Holm is conservative; consider re-running with --p-adjust BH for a less stringent control)"),
        tags$li("True compositional differences are below the detection threshold of bias-corrected ANCOMBC2"),
        tags$li("Prevalence filter (default 10%) removed taxa that may have been DA but were rare"),
        tags$li("Groups are biologically very similar at this taxonomic level (try genus/species levels for more resolution)")
      ))
  )
}

section_step15 <- function() {
  s <- extract_step15()
  if (is.null(s)) return(tagList(step_header(15, "Functional prediction (PICRUSt2)"), step_not_run_box(15)))
  tagList(
    step_header(15, "Functional prediction (PICRUSt2)"),
    tags$p(HTML(paste(
      "PICRUSt2 (Douglas et al. 2020) predicts gene-family abundances from 16S rRNA placements onto a reference tree.",
      "ASVs with NSTI > <code>", s$nsti_threshold, "</code> were dropped (high NSTI = unreliable phylogenetic placement).",
      "Per-sample mean NSTI is reported in the per-variable HTML report; samples with mean NSTI > 1.0 are flagged as having unreliable functional predictions."
    ))),
    if (s$cog_skipped) tags$div(class = "box warning",
      tags$p(class = "box-title", "Known limitation: COG database skipped"),
      tags$p("Three of four functional databases (KO, EC, MetaCyc) were generated successfully.  COG was NOT generated due to a known incompatibility in PICRUSt2 v2.6.x where the COG mapping file uses an HDF5 format incompatible with the current pathway-pipeline implementation.  This affects only the COG-category breakdown; KEGG (KO), EC, and MetaCyc-pathway predictions are unaffected.")),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "NSTI threshold:"),
              tags$span(class = "val", s$nsti_threshold)),
      tags$li(tags$span(class = "key", "Per-variable directories:"),
              tags$span(class = "val", s$n_per_var)),
      tags$li(tags$span(class = "key", "Per-variable HTML report (with NSTI per-sample warnings):"),
              tags$span(class = "val", if (file.exists(s$html_report)) relpath(s$html_report) else "(not generated)"))
    )
  )
}

section_step16 <- function() {
  s <- extract_step16()
  if (is.null(s)) return(tagList(step_header(16, "Random Forest biomarker classifier"), step_not_run_box(16)))
  # Build per-variable summary
  rows <- list()
  for (f in s$summary_files) {
    df <- tryCatch(read.xlsx(f), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    df$source_file <- basename(f)
    rows[[length(rows) + 1]] <- df
  }
  tbl <- if (length(rows) > 0) do.call(rbind, rows) else NULL

  # Embed best ROC curve per variable (highest accuracy)
  roc_pngs <- character(0)
  for (vd in s$per_var_dirs) {
    rp <- list.files(vd, pattern = "roc_curves\\.png$",
                     full.names = TRUE, recursive = TRUE)
    if (length(rp) > 0) roc_pngs <- c(roc_pngs, rp[1])
  }

  tagList(
    step_header(16, "Random Forest biomarker classifier"),
    tags$p(HTML(paste(
      "A Random Forest classifier (Breiman 2001) was trained at every taxonomic level x every categorical variable",
      "using the <code>ranger</code> implementation (Wright & Ziegler 2017).",
      "Cross-validation strategy is auto-selected: 5-fold stratified CV when N >= 20 samples, leave-one-out CV when N < 20.",
      "Class imbalance is handled via <code>case.weights</code>.",
      "Outputs include accuracy, macro-averaged AUC (one-vs-rest), per-class sensitivity/specificity/F1,",
      "permutation feature importance, and SHAP-style per-sample local-importance heatmaps."
    ))),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Variables modelled:"),
              tags$span(class = "val", s$n_per_var)),
      tags$li(tags$span(class = "key", "Per-variable summary files:"),
              tags$span(class = "val", length(s$summary_files)))
    ),
    if (!is.null(tbl)) tagList(
      tags$h4("Per-variable summary (one row per taxonomic level)"),
      df_to_html(tbl, NULL, max_rows = 50)),
    if (length(roc_pngs) > 0) tagList(
      tags$h4("ROC curves (one panel per variable)"),
      tags$div(class = "fig-grid",
               tagList(lapply(roc_pngs,
                              function(p) embed_image(p, basename(dirname(p)))))))
  )
}

section_step17 <- function() {
  s <- extract_step17()
  if (is.null(s)) return(tagList(step_header(17, "Co-occurrence networks"), step_not_run_box(17)))

  # Global network plots
  global_pngs <- character(0)
  for (gd in s$global_dirs) {
    p <- file.path(gd, "network_plot.png")
    if (file.exists(p)) global_pngs <- c(global_pngs, p)
  }
  # Per-variable multi-panel plots
  multi_pngs <- list.files(file.path(MBX_OUT_DIR, "17_co_occurrence_networks"),
                           pattern = "multi_panel_plot\\.png$",
                           recursive = TRUE, full.names = TRUE)

  # Roll up group_comparison.xlsx files
  cmp_rows <- list()
  for (f in s$group_cmps) {
    df <- tryCatch(read.xlsx(f, sheet = 1), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    df$source_file <- relpath(f)
    cmp_rows[[length(cmp_rows) + 1]] <- df
  }
  cmp_tbl <- if (length(cmp_rows) > 0) {
    # Reduce columns from possibly heterogenous cmp dfs by union
    cn <- unique(unlist(lapply(cmp_rows, names)))
    do.call(rbind, lapply(cmp_rows, function(d) {
      missing <- setdiff(cn, names(d))
      for (m in missing) d[[m]] <- NA
      d[, cn, drop = FALSE]
    }))
  } else NULL

  tagList(
    step_header(17, "Co-occurrence / correlation networks"),
    tags$p(HTML(paste(
      "Compositional-aware correlation networks were built using a CLR-transform + Spearman + BH-FDR pipeline (Gloor et al. 2017).",
      "Taxa present in <30 % of samples were dropped before CLR transformation.",
      "Edges with |Spearman rho| >= 0.6 AND BH-q <= 0.05 are retained.",
      "Modules are detected via Louvain (Blondel et al. 2008) using |rho| as edge weight.",
      "Hub taxa are the top 10 % by combined degree + betweenness + hub_score rank."
    ))),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Global networks:"),
              tags$span(class = "val", s$n_global)),
      tags$li(tags$span(class = "key", "Per-variable comparisons:"),
              tags$span(class = "val", s$n_per_var))
    ),
    if (length(global_pngs) > 0) tagList(
      tags$h4("Global networks (the ecological backbone)"),
      tags$div(class = "fig-grid",
               tagList(lapply(global_pngs,
                              function(p) embed_image(p, basename(dirname(p))))))),
    if (length(multi_pngs) > 0) tagList(
      tags$h4("Per-variable multi-panel networks"),
      tags$div(class = "fig-grid",
               tagList(lapply(multi_pngs,
                              function(p) embed_image(p,
                                sprintf("%s / %s",
                                        basename(dirname(dirname(p))),
                                        basename(dirname(p))))))))
    ,
    if (!is.null(cmp_tbl)) tagList(
      tags$h4("Per-group network metrics roll-up"),
      df_to_html(cmp_tbl, NULL, max_rows = 60))
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Cross-step section: convergence
# ────────────────────────────────────────────────────────────────────────────
section_convergence <- function(s14, s16, s17) {
  conv <- build_convergence(s14, s16, s17)
  tagList(
    tags$h2(id = "convergence", "Cross-step convergence: strongest biomarker candidates"),
    tags$p("This table identifies taxa that are simultaneously flagged by two or more independent methods: differential abundance (ANCOMBC2, q < 0.05), Random Forest top-importance (top 20 features per model), and co-occurrence network hub status (top 10 % combined-rank).  Convergent findings are far more likely to reflect true biology than any single method's hits in isolation."),
    if (is.null(conv)) tags$div(class = "box info",
      tags$p(class = "box-title", "No multi-method convergent taxa were identified"),
      tags$p("Possible reasons (any of which may apply):"),
      tags$ul(
        tags$li("Sample size too low for robust statistical signal"),
        tags$li("Effect sizes below detection threshold for at least one of the three methods"),
        tags$li("Methods are detecting genuinely orthogonal aspects of variation (variance- vs hub- vs DA-driven)"),
        tags$li("Steps 14, 16, or 17 were not executed in this run")
      ))
    else tagList(
      tags$div(class = "box success",
        tags$p(class = "box-title", sprintf("%d convergent taxa identified", nrow(conv))),
        tags$p("Taxa with convergence_score = 3 are the most reliable biomarker candidates -- significant in ANCOMBC2 AND high-importance in Random Forest AND a hub in the co-occurrence network.")),
      df_to_html(conv, NULL, max_rows = 60))
  )
}

# ────────────────────────────────────────────────────────────────────────────
# What's still available section (graceful handling per user request)
# ────────────────────────────────────────────────────────────────────────────
section_available <- function() {
  not_run <- discovery[discovery$status == "MISSING", , drop = FALSE]
  if (nrow(not_run) == 0) return(tagList(
    tags$h2(id = "available", "Pipeline coverage"),
    tags$div(class = "box success",
      tags$p(class = "box-title", "Complete pipeline executed"),
      tags$p("All 18 steps of the mbX Pro pipeline were run for this analysis.  No additional analyses are available."))
  ))
  # Build a list of still-available steps with the exact command to run
  rows <- lapply(seq_len(nrow(not_run)), function(i) {
    tags$li(
      tags$b(sprintf("Step %d -- %s", not_run$id[i], not_run$desc[i])), tags$br(),
      tags$code(sprintf("%s %s", not_run$script[i], MBX_OUT_DIR))
    )
  })
  tagList(
    tags$h2(id = "available", "Pipeline coverage"),
    tags$p(HTML(sprintf(
      "<b>The user selected to run the pipeline up to step %d.</b>  The mbX Pro pipeline covers 18 steps in total (steps 0-17).  The following analyses are available but were NOT executed for this output directory.  To add any of them to a future report, run the listed command and then re-run <code>mbx_final_report.sh</code>.",
      max_done_id))),
    tags$ul(rows)
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Methods + Citations
# ────────────────────────────────────────────────────────────────────────────
section_methods <- function() {
  tagList(
    tags$h2(id = "methods", "Methods (paper-ready)"),
    tags$p("Raw paired-end FASTQ files were processed through the mbX Pro 16S rRNA pipeline.  Primers were detected automatically using an IUPAC-aware sliding-window algorithm, then trimmed via the QIIME2 import + DADA2 trim-left parameter.  Amplicon sequence variants (ASVs) were inferred with DADA2 (Callahan et al. 2016) using truncation lengths determined from the demux quality summary (longest position with median Q >= 25).  Taxonomy was assigned via a Naive-Bayes classifier trained on the Greengenes2 backbone (McDonald et al. 2024), region-extracted to the user's primer pair.  Mitochondrial and chloroplast 16S sequences were filtered before downstream analysis.  Tabular outputs at all seven taxonomic levels were tidied with mbX::ezclean()."),
    tags$p("Within-sample diversity was assessed at the analytical-rarefaction depth (selected via a group-aware 90 %/80 % retention rule) using ASV richness, Shannon entropy, Simpson, Faith's PD (Faith 1992), and Pielou's evenness.  Between-sample diversity used Bray-Curtis, Jaccard, weighted UniFrac, and unweighted UniFrac distances visualised via PCoA, with significance assessed by PERMANOVA (with pairwise BH-adjusted comparisons when groups > 2) and dispersion homogeneity by PERMDISP.  A multivariable Adonis model captured the joint contribution of all metadata factors."),
    tags$p("Differential abundance was tested with ANCOM-BC2 (Lin et al. 2024) using <code>pairwise = TRUE</code> and <code>global = TRUE</code>; the alphabetically-first level of each categorical variable was used as the (reference-free) baseline.  Functional gene-family profiles were predicted with PICRUSt2 (Douglas et al. 2020) for the KEGG-Ortholog (KO), Enzyme Commission (EC), and MetaCyc-pathway databases, with NSTI filtering at threshold 2.  A Random Forest classifier (Breiman 2001; ranger implementation: Wright & Ziegler 2017) was trained at every taxonomic level for every categorical variable, with auto-selected cross-validation (5-fold stratified for N >= 20, leave-one-out otherwise) and SHAP-style per-sample local-importance via <code>ranger</code>'s permutation framework.  Co-occurrence networks were built via centred-log-ratio transform + Spearman correlation + Benjamini-Hochberg FDR (Aitchison 1982; Gloor et al. 2017), with edge inclusion threshold |rho| >= 0.6 and q <= 0.05 and module detection by Louvain (Blondel et al. 2008).  All analyses were orchestrated by the mbX Pro shell pipeline running on QIIME2 amplicon-2025.4 and R 4.4."),
    tags$h3(id = "bibliography", "Citations"),
    tags$div(class = "bibliography",
      tags$ol(
        tags$li(HTML("Lamichhane U., Lourenco J. (2025). <b>mbX: An R Package for Streamlined Microbiome Analysis.</b> Stats 8(2):44.  DOI: <a href='https://doi.org/10.3390/stats8020044' target='_blank'>10.3390/stats8020044</a>.  <i>(Cite this paper for the mbX backbone of mbX Pro until the dedicated mbX Pro paper is released.)</i>")),
        tags$li("Aitchison J. (1982). The statistical analysis of compositional data. Journal of the Royal Statistical Society Series B 44:139-177."),
        tags$li("Blondel V.D., Guillaume J.-L., Lambiotte R., Lefebvre E. (2008). Fast unfolding of communities in large networks. Journal of Statistical Mechanics 2008:P10008."),
        tags$li("Bolyen E. et al. (2019). Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37:852-857."),
        tags$li("Breiman L. (2001). Random forests. Machine Learning 45:5-32."),
        tags$li("Callahan B.J., McMurdie P.J., Rosen M.J., Han A.W., Johnson A.J.A., Holmes S.P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature Methods 13:581-583."),
        tags$li("Douglas G.M. et al. (2020). PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38:685-688."),
        tags$li("Faith D.P. (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation 61:1-10."),
        tags$li("Gloor G.B., Macklaim J.M., Pawlowsky-Glahn V., Egozcue J.J. (2017). Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology 8:2224."),
        tags$li("Hurlbert S.H. (1971). The nonconcept of species diversity: a critique and alternative parameters. Ecology 52:577-586."),
        tags$li("Lin H., Peddada S.D. (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods 21:83-91. (ANCOM-BC2)"),
        tags$li("McDonald D., Jiang Y., Balaban M., Cantrell K. et al. (2024). Greengenes2 unifies microbial data in a single reference tree. Nature Biotechnology 42:715-718."),
        tags$li("Spearman C. (1904). The proof and measurement of association between two things. American Journal of Psychology 15:72-101."),
        tags$li("Wright M.N., Ziegler A. (2017). ranger: a fast implementation of random forests for high dimensional data in C++ and R. Journal of Statistical Software 77:1-17.")
      ))
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Caveats / limitations
# ────────────────────────────────────────────────────────────────────────────
section_caveats <- function(s11, s15) {
  items <- list()
  # Sample size note
  items[[length(items)+1]] <- tags$li(tags$b("Sample size."), " The statistical power of every step (DA, RF, networks, PERMANOVA) scales with N.  All p-values and q-values reported here are conditional on the sample sizes shown in the inventory above.  Pairwise PERMANOVA, ANCOMBC2, and Random Forest are particularly sensitive to small per-group N.")
  # PICRUSt2 caveat (only if step 15 ran)
  if (!is.null(s15)) {
    items[[length(items)+1]] <- tags$li(tags$b("Functional predictions are inferences."), " PICRUSt2 outputs are gene-family abundances PREDICTED from 16S rRNA placements onto a reference tree; they are NOT direct measurements of expressed function or metagenomic content.  Samples with mean NSTI > 1.0 should be treated with extra caution (see step 15 per-variable HTML report).")
    if (isTRUE(s15$cog_skipped))
      items[[length(items)+1]] <- tags$li(tags$b("COG database not generated."), " Three of four PICRUSt2 functional databases (KO, EC, MetaCyc) were generated successfully.  COG was skipped due to a known PICRUSt2 v2.6.x format incompatibility -- this is a tool-side issue, not a data issue.")
  }
  # Compositional data caveat
  items[[length(items)+1]] <- tags$li(tags$b("Compositionality."), " 16S rRNA-derived abundances are compositional (constrained to a simplex).  Conventional Pearson correlation and absolute-abundance assumptions are inappropriate.  This pipeline uses ANCOM-BC2 (which models the bias term) for DA and CLR + Spearman for networks (Gloor et al. 2017).")
  # Reference DB note
  items[[length(items)+1]] <- tags$li(tags$b("Reference-database dependency."), " Taxonomy is only as good as the reference (Greengenes2 in this pipeline).  Taxa absent from the reference cannot be classified; fast-evolving environmental clades may be under-represented.")
  # Rarefaction note (only if step 11 ran)
  if (!is.null(s11) && nzchar(s11$samples_lost) && s11$samples_lost != "0")
    items[[length(items)+1]] <- tags$li(tags$b(sprintf("%s sample(s) dropped at rarefaction depth.", s11$samples_lost)), " The rarefaction depth was selected to balance retention and richness; this trade-off explicitly excluded samples below the depth threshold.  See step 11 for the per-sample retention table and step 12 for which samples are absent from alpha-diversity outputs.")

  tagList(
    tags$h2(id = "caveats", "Caveats and limitations"),
    tags$ul(items)
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Reproducibility manifest
# ────────────────────────────────────────────────────────────────────────────
section_reproducibility <- function() {
  tagList(
    tags$h2(id = "reproducibility", "Reproducibility manifest"),
    tags$p("Every component listed below contributes to the figures and numbers in this report.  To reproduce the analysis from raw FASTQ to this PDF, run the commands in order on a machine with QIIME2 amplicon-2025.4 and R installed system-wide."),
    tags$h3("Software versions"),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Operating system:"),
              tags$span(class = "val", sprintf("%s %s",
                Sys.info()[["sysname"]], Sys.info()[["release"]]))),
      tags$li(tags$span(class = "key", "Machine:"),
              tags$span(class = "val", Sys.info()[["nodename"]])),
      tags$li(tags$span(class = "key", "R:"),
              tags$span(class = "val", R_VERSION)),
      tags$li(tags$span(class = "key", "mbx_final_report.sh version:"),
              tags$span(class = "val", SCRIPT_VERSION))
    ),
    tags$h3("Pipeline reproduction sequence"),
    tags$p(HTML(sprintf("Set <code>OUTDIR</code> to the absolute path of your mbX Pro output directory (this run: <code>%s</code>), then execute the commands below in order.  The pipeline is idempotent -- already-completed steps are skipped automatically.", MBX_OUT_DIR))),
    tags$pre(paste(c(
      sprintf("export OUTDIR=\"%s\"", MBX_OUT_DIR),
      "",
      "# 1. Detect primers",
      "mbx_primer_identifier.sh \"$OUTDIR\"",
      "# 2. Build manifest",
      "create_manifest.sh \"$OUTDIR\"",
      "# 3. Make QIIME2 artifact",
      "artifact_creator.sh \"$OUTDIR/1_manifest_file/manifest.txt\"",
      "# 4. Find DADA2 truncation",
      "mbx_dada2_parameter_finder.sh \"$OUTDIR/2_first_artifact_file/Paired_End_artifact.qza\"",
      "# 5. Run DADA2",
      "mbx_dada2_run.sh \"$OUTDIR\"",
      "# 6-7. Train classifier + assign taxonomy",
      "mbx_classifier_arranger.sh \"$OUTDIR\"",
      "mbx_classifier_run.sh \"$OUTDIR\"",
      "# 8. Export taxonomy CSVs",
      "mbx_taxonomy_run.sh \"$OUTDIR\" <metadata.txt>",
      "# 9-11. Clean / visualise / pre-diversity",
      "mbx_ezclean_all_levels.sh \"$OUTDIR\"",
      "mbx_ezviz_all_levels_all_treatments.sh \"$OUTDIR\"",
      "mbx_ezstat_all_levels_all_treatments.sh \"$OUTDIR\"",
      "mbx_pre_diversity_parameters.sh \"$OUTDIR\"",
      "# 12-17. Diversity / DA / functional / ML / networks",
      "mbx_alpha_diversity_run.sh \"$OUTDIR\"",
      "mbx_beta_diversity_run.sh \"$OUTDIR\"",
      "mbx_ancombc2_run.sh \"$OUTDIR\"",
      "mbx_picrust_run.sh \"$OUTDIR\"",
      "mbx_ml_classifier_run.sh \"$OUTDIR\"",
      "mbx_network_run.sh \"$OUTDIR\"",
      "# 18. This report",
      "mbx_final_report.sh \"$OUTDIR\""
    ), collapse = "\n"))
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Sample inventory + per-sample quality
# ────────────────────────────────────────────────────────────────────────────
section_inventory <- function(s8, s11, s12) {
  if (is.null(s8) || length(s8$per_level) == 0)
    return(tagList(
      tags$h2(id = "inventory", "Sample inventory"),
      tags$p("(insufficient information to build inventory; step 8 must be complete)")
    ))
  # Use the genus-level cleaned xlsx as the inventory anchor
  g <- Filter(function(x) grepl("genera", x$dir, ignore.case = TRUE), s8$per_level)
  src <- if (length(g) > 0) g[[1]]$xlsx else s8$per_level[[1]]$xlsx
  df <- tryCatch(read.xlsx(src, check.names = FALSE), error = function(e) NULL)
  if (is.null(df)) return(tagList(
    tags$h2(id = "inventory", "Sample inventory"),
    tags$p(sprintf("(could not read inventory anchor: %s)", src))
  ))
  sid_col <- names(df)[1]
  inv <- data.frame(sample_id = as.character(df[[sid_col]]),
                    stringsAsFactors = FALSE)
  # Optional: alpha diversity columns
  if (!is.null(s12) && !is.null(s12$alpha_df)) {
    a <- s12$alpha_df
    a_sid <- names(a)[1]
    keep_cols <- intersect(c("ASVs_or_Features", "Shannon_Index"), names(a))
    if (length(keep_cols) > 0)
      inv <- merge(inv,
                   a[, c(a_sid, keep_cols), drop = FALSE],
                   by.x = "sample_id", by.y = a_sid, all.x = TRUE)
  }
  tagList(
    tags$h2(id = "inventory", "Sample inventory"),
    tags$p("Per-sample manifest derived from the genus-level cleaned table.  Alpha diversity (ASVs and Shannon) added when step 12 has been run."),
    df_to_html(inv, NULL, max_rows = 100)
  )
}

# ────────────────────────────────────────────────────────────────────────────
# Build report
# ────────────────────────────────────────────────────────────────────────────
cat("[R] Extracting per-step data ...\n")
s0  <- if ("0"  %in% done_ids || "0"  %in% partial_ids) extract_step0()  else NULL
s1  <- if ("1"  %in% done_ids || "1"  %in% partial_ids) extract_step1()  else NULL
s4  <- if ("4"  %in% done_ids || "4"  %in% partial_ids) extract_step4()  else NULL
s7  <- if ("7"  %in% done_ids || "7"  %in% partial_ids) extract_step7()  else NULL
s8  <- if ("8"  %in% done_ids || "8"  %in% partial_ids) extract_step8()  else NULL
s11 <- if ("11" %in% done_ids || "11" %in% partial_ids) extract_step11() else NULL
s12 <- if ("12" %in% done_ids || "12" %in% partial_ids) extract_step12() else NULL
s14 <- if ("14" %in% done_ids || "14" %in% partial_ids) extract_step14() else NULL
s15 <- if ("15" %in% done_ids || "15" %in% partial_ids) extract_step15() else NULL
s16 <- if ("16" %in% done_ids || "16" %in% partial_ids) extract_step16() else NULL
s17 <- if ("17" %in% done_ids || "17" %in% partial_ids) extract_step17() else NULL

cat("[R] Building HTML sections ...\n")

# ── Cover header (logo + title + subtitle), reused in both doc builds ───────
cover_header <- if (nzchar(LOGO_DATA_URI)) {
  tags$div(class = "logo-header",
    tags$img(src = LOGO_DATA_URI, alt = "mbX Pro logo"),
    tags$div(class = "titles",
      tags$h1("mbX Pro 16S rRNA Microbiome Analysis Report"),
      tags$p(class = "subtitle",
             "Single-step automated pipeline -- raw FASTQ to functional insights")
    )
  )
} else {
  tagList(
    tags$h1("mbX Pro 16S rRNA Microbiome Analysis Report"),
    tags$p(class = "subtitle",
           "Single-step automated pipeline -- raw FASTQ to functional insights")
  )
}

# ── How-to-cite box (prominent, near the top of the report) ────────────────
cite_box <- tags$div(class = "cite-box",
  tags$p(class = "cite-title", "How to cite this report"),
  tags$p(class = "cite-line",
    HTML("The peer-reviewed paper describing the full <b>mbX Pro</b> pipeline is currently in preparation.  For now, please cite the original <b>mbX</b> R-package paper (which forms the analytic backbone of this pipeline):")),
  tags$p(class = "cite-line full",
    "Lamichhane U., Lourenco J. (2025).  ",
    HTML("<b>mbX: An R Package for Streamlined Microbiome Analysis</b>.  "),
    HTML("<i>Stats</i> 8(2): 44.  "),
    "DOI: ",
    tags$a(href = "https://doi.org/10.3390/stats8020044", target = "_blank",
           "10.3390/stats8020044"),
    "  |  URL: ",
    tags$a(href = "https://www.mdpi.com/2571-905X/8/2/44", target = "_blank",
           "https://www.mdpi.com/2571-905X/8/2/44"))
)

toc_items <- c(
  "Executive summary", "Sample inventory", "Reproducibility manifest",
  "Step 0: Primer detection", "Steps 1-2: Manifest + artifact",
  "Steps 3-4: DADA2", "Steps 5-6: Classifier + taxonomy",
  "Step 7: Taxonomy CSVs", "Step 8: Cleaned files",
  "Step 9: Visualisations (ezviz)", "Step 10: Statistics (ezstat)",
  "Step 11: Pre-diversity", "Step 12: Alpha diversity",
  "Step 13: Beta diversity", "Step 14: Differential abundance (ANCOMBC2)",
  "Step 15: Functional prediction (PICRUSt2)",
  "Step 16: ML biomarkers", "Step 17: Co-occurrence networks",
  "Cross-step convergence", "Caveats and limitations",
  "Pipeline coverage", "Methods + citations"
)
toc <- tags$div(class = "toc",
  tags$h3("Contents"),
  tags$ol(lapply(toc_items, function(x) tags$li(x))))

doc <- tags$html(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$title("mbX Pro Final Report"),
    tags$style(HTML(CSS))
  ),
  tags$body(
    tags$div(class = "container",
      cover_header,
      tags$div(class = "cover-meta",
        tags$div(tags$b("Generated"), GEN_AT),
        tags$div(tags$b("Output directory"),
                 tags$span(class = "path-val", MBX_OUT_DIR)),
        tags$div(tags$b("Steps executed"), sprintf("%d of 18 (max = step %d)",
                 length(done_ids), max_done_id)),
        tags$div(tags$b("Pipeline version"),
                 sprintf("mbx_final_report.sh %s | R %s", SCRIPT_VERSION, R_VERSION))
      ),
      cite_box,
      build_dashboard(s1, s4, s7, s8, s11, s14, s16),
      toc,

      tags$h2(id = "summary", "Executive summary"),
      tags$p(HTML(sprintf(
        "This report consolidates every step of the mbX Pro 16S rRNA microbiome pipeline that has been executed for output directory <code>%s</code>.  Out of 18 total pipeline steps, %d were executed (highest step = %d).  Each section below documents the parameters used, key outputs, representative figures, and -- where applicable -- a scientific reason if a particular sub-result is missing or empty.  A dedicated cross-step convergence analysis identifies taxa flagged as biologically meaningful by two or more independent methods (differential abundance, machine learning, network hub status); these are the strongest biomarker candidates worth follow-up validation.",
        MBX_OUT_DIR, length(done_ids), max_done_id))),
      tags$p(class = "path-note",
             HTML(sprintf("All file paths shown in this report are <b>relative to the output directory</b> above (<code>%s</code>).  The provenance <code>.txt</code> files inside each step directory keep their original absolute paths.",
                          MBX_OUT_DIR))),

      section_inventory(s8, s11, s12),
      section_reproducibility(),

      tags$div(class = "step-block", section_step0()),
      tags$div(class = "step-block", section_step12()),
      tags$div(class = "step-block", section_step34()),
      tags$div(class = "step-block", section_step56()),
      tags$div(class = "step-block", section_step7()),
      tags$div(class = "step-block", section_step8()),
      tags$div(class = "step-block", section_step9()),
      tags$div(class = "step-block", section_step10()),
      tags$div(class = "step-block", section_step11()),
      tags$div(class = "step-block", section_step12()),
      tags$div(class = "step-block", section_step13()),
      tags$div(class = "step-block", section_step14()),
      tags$div(class = "step-block", section_step15()),
      tags$div(class = "step-block", section_step16()),
      tags$div(class = "step-block", section_step17()),

      section_convergence(s14, s16, s17),
      section_caveats(s11, s15),
      section_available(),
      section_methods(),

      tags$footer(
        sprintf("mbX Pro Final Report v%s | generated %s | %d of 18 pipeline steps executed",
                SCRIPT_VERSION, GEN_AT, length(done_ids)))
    )
  )
)

# Note: `section_step12()` defined above is used twice in the document --
# once as the early-pipeline "Steps 1-2: Manifest + artifact" wrapper
# (yes, the function name collides because step 1+2 are reported together),
# and once as alpha diversity.  These are two DIFFERENT functions defined
# by the same name in this script -- the SECOND definition wins.  Fix by
# renaming the alpha section call:
# (Below: monkey-patch -- we already passed section_step12 as alpha above.)
# Actually we have a name collision; we will rebuild the doc with explicit calls.

cat("[R] Saving HTML ...\n")
# Re-build with corrected section name resolution: we need the early-pipeline
# steps 1-2 wrapper to be a SEPARATE function from alpha-diversity.  Above
# I defined two different functions with the same name `section_step12`, which
# means only the second (alpha) survives.  Let's rename and redo:

section_steps12 <- function() {
  s1 <- extract_step1(); s2_dir <- file.path(MBX_OUT_DIR, "2_first_artifact_file")
  if (is.null(s1)) return(tagList(step_header(1, "QIIME2 manifest"), step_not_run_box(1)))
  tagList(
    step_header(1, "QIIME2 manifest construction"),
    tags$p("QIIME2's tab-separated manifest format with absolute file paths.  Sample IDs are extracted via regex matching s(ample)?[-_]?[0-9]+ patterns followed by filename-cleanup fallbacks."),
    tags$ul(class = "kv-list",
      tags$li(tags$span(class = "key", "Number of samples:"),
              tags$span(class = "val", s1$n_samples)),
      tags$li(tags$span(class = "key", "Sequencing layout:"),
              tags$span(class = "val", if (s1$paired) "paired-end" else "single-end")),
      tags$li(tags$span(class = "key", "Manifest path:"),
              tags$span(class = "val", relpath(s1$manifest_path))),
      tags$li(tags$span(class = "key", "Columns:"),
              tags$span(class = "val", s1$columns))
    ),
    tags$h3("Step 2 - QIIME2 .qza artifact creation"),
    if (status_for(2) != "DONE") step_not_run_box(2) else
    tagList(
      tags$p("The manifest was imported into QIIME2 as a SampleData[PairedEndSequencesWithQuality] artifact (.qza) and a corresponding visualisation (.qzv) was rendered for quality assessment."),
      tags$ul(class = "kv-list",
        tags$li(tags$span(class = "key", "Artifact:"),
                tags$span(class = "val", relpath(file.path(s2_dir, "Paired_End_artifact.qza")))),
        tags$li(tags$span(class = "key", "Visualisation:"),
                tags$span(class = "val", relpath(file.path(s2_dir, "Paired_End_artifact.qzv"))))
      )
    )
  )
}

# Rebuild doc with the correct sections
doc <- tags$html(
  tags$head(
    tags$meta(charset = "utf-8"),
    tags$title("mbX Pro Final Report"),
    tags$style(HTML(CSS))
  ),
  tags$body(
    tags$div(class = "container",
      cover_header,
      tags$div(class = "cover-meta",
        tags$div(tags$b("Generated"), GEN_AT),
        tags$div(tags$b("Output directory"),
                 tags$span(class = "path-val", MBX_OUT_DIR)),
        tags$div(tags$b("Steps executed"), sprintf("%d of 18 (max = step %d)",
                 length(done_ids), max_done_id)),
        tags$div(tags$b("Pipeline version"),
                 sprintf("mbx_final_report.sh %s | R %s", SCRIPT_VERSION, R_VERSION))
      ),
      cite_box,
      build_dashboard(s1, s4, s7, s8, s11, s14, s16),
      toc,

      tags$h2(id = "summary", "Executive summary"),
      tags$p(HTML(sprintf(
        "This report consolidates every step of the mbX Pro 16S rRNA microbiome pipeline that has been executed for output directory <code>%s</code>.  Out of 18 total pipeline steps, %d were executed (highest step = %d).  Each section below documents the parameters used, key outputs, representative figures, and -- where applicable -- a scientific reason if a particular sub-result is missing or empty.  A dedicated cross-step convergence analysis identifies taxa flagged as biologically meaningful by two or more independent methods (differential abundance, machine learning, network hub status); these are the strongest biomarker candidates worth follow-up validation.",
        MBX_OUT_DIR, length(done_ids), max_done_id))),
      tags$p(class = "path-note",
             HTML(sprintf("All file paths shown in this report are <b>relative to the output directory</b> above (<code>%s</code>).  The provenance <code>.txt</code> files inside each step directory keep their original absolute paths.",
                          MBX_OUT_DIR))),

      section_inventory(s8, s11, s12),
      section_reproducibility(),

      tags$div(class = "step-block", section_step0()),
      tags$div(class = "step-block", section_steps12()),
      tags$div(class = "step-block", section_step34()),
      tags$div(class = "step-block", section_step56()),
      tags$div(class = "step-block", section_step7()),
      tags$div(class = "step-block", section_step8()),
      tags$div(class = "step-block", section_step9()),
      tags$div(class = "step-block", section_step10()),
      tags$div(class = "step-block", section_step11()),
      tags$div(class = "step-block", section_step12()),
      tags$div(class = "step-block", section_step13()),
      tags$div(class = "step-block", section_step14()),
      tags$div(class = "step-block", section_step15()),
      tags$div(class = "step-block", section_step16()),
      tags$div(class = "step-block", section_step17()),

      section_convergence(s14, s16, s17),
      section_caveats(s11, s15),
      section_available(),
      section_methods(),

      tags$footer(
        sprintf("mbX Pro Final Report v%s | generated %s | %d of 18 pipeline steps executed",
                SCRIPT_VERSION, GEN_AT, length(done_ids)))
    )
  )
)

# save_html() correctly handles tags$head() hoisting + writes <!DOCTYPE>;
# as.character() silently drops <head> (htmltools quirk).
htmltools::save_html(doc, file = HTML_OUT, libdir = "_assets",
                     background = "white")
# Remove the auto-created _assets/ dir (we have nothing to put in it; everything
# is already inlined as base64) so the report stays a single self-contained file.
.assets <- file.path(dirname(HTML_OUT), "_assets")
if (dir.exists(.assets) && length(list.files(.assets)) == 0) {
  unlink(.assets, recursive = TRUE)
}
cat(sprintf("[OK]    HTML written: %s (%.1f MB)\n",
            HTML_OUT, file.info(HTML_OUT)$size / 1024 / 1024))
REPORT_R_EOF

# Invoke R
ok "R script written ($(wc -l < "$REPORT_R" | tr -d ' ') lines)"
info "Generating HTML report (this may take 1-3 minutes for image embedding)..."

timer_start
MBX_REPORT_OUT_DIR="$MBX_OUT_DIR" \
MBX_REPORT_REPORT_DIR="$REPORT_DIR" \
MBX_REPORT_HTML_OUT="$HTML_OUT" \
MBX_REPORT_DISC_JSON="$DISC_JSON" \
MBX_REPORT_NOW="$NOW" \
MBX_REPORT_R_VERSION="$R_VERSION" \
MBX_REPORT_LOGO_PATH="$LOGO_PATH" \
  _R --vanilla "$REPORT_R" \
    || err "HTML generation failed.  See R output above."
timer_end

if [[ ! -s "$HTML_OUT" ]]; then
  err "HTML file was not produced or is empty."
fi
HTML_SIZE_MB="$(du -m "$HTML_OUT" | cut -f1)"
ok "HTML report ready: $HTML_OUT  (${HTML_SIZE_MB} MB)"

# ─────────────────────────────────────────────────────────────────────────────
step "4/5 — Convert HTML to PDF"
# ─────────────────────────────────────────────────────────────────────────────

if $NO_PDF; then
  info "PDF rendering skipped (--no-pdf was set)."
  PDF_OUT=""
else
  # Auto-detect PDF engine if not forced
  CHOSEN_ENGINE=""
  if [[ -n "$PDF_ENGINE" ]]; then
    CHOSEN_ENGINE="$PDF_ENGINE"
  else
    # 1) Headless Chrome / Chromium / Edge
    for _b in \
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
        "/Applications/Chromium.app/Contents/MacOS/Chromium" \
        "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge" \
        "$(command -v google-chrome 2>/dev/null || true)" \
        "$(command -v chromium-browser 2>/dev/null || true)" \
        "$(command -v chromium 2>/dev/null || true)"; do
      if [[ -n "$_b" && -x "$_b" ]]; then
        CHOSEN_ENGINE="chrome"; CHROME_BIN="$_b"; break
      fi
    done
    # 2) wkhtmltopdf
    if [[ -z "$CHOSEN_ENGINE" ]] && command -v wkhtmltopdf &>/dev/null; then
      CHOSEN_ENGINE="wkhtmltopdf"
    fi
    # 3) macOS cupsfilter
    if [[ -z "$CHOSEN_ENGINE" ]] && command -v cupsfilter &>/dev/null; then
      CHOSEN_ENGINE="cupsfilter"
    fi
  fi

  if [[ -z "$CHOSEN_ENGINE" ]]; then
    if $SKIP_INSTALL; then
      warn "No PDF engine found and --skip-install set.  Skipping PDF generation."
      warn "  -> Install Google Chrome OR run: brew install wkhtmltopdf"
      warn "  -> Or open the HTML in any browser and File > Print > Save as PDF"
      PDF_OUT=""
    else
      info "No PDF engine found.  Attempting to install wkhtmltopdf via Homebrew..."
      if command -v brew &>/dev/null; then
        if brew install --cask wkhtmltopdf 2>&1 | tail -5; then
          CHOSEN_ENGINE="wkhtmltopdf"
        else
          warn "wkhtmltopdf install failed.  Skipping PDF generation."
          warn "  -> Open the HTML in any browser and File > Print > Save as PDF"
          PDF_OUT=""
        fi
      else
        warn "Homebrew not found.  Skipping PDF generation."
        warn "  -> Install Chrome or wkhtmltopdf manually, then re-run this step."
        PDF_OUT=""
      fi
    fi
  fi

  if [[ -n "$CHOSEN_ENGINE" ]]; then
    info "PDF engine: $CHOSEN_ENGINE"
    case "$CHOSEN_ENGINE" in
      chrome)
        # Resolve binary if user passed --pdf-engine chrome (re-run lookup)
        if [[ -z "${CHROME_BIN:-}" ]]; then
          for _b in \
              "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
              "/Applications/Chromium.app/Contents/MacOS/Chromium" \
              "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge" \
              "$(command -v google-chrome 2>/dev/null || true)" \
              "$(command -v chromium 2>/dev/null || true)"; do
            [[ -n "$_b" && -x "$_b" ]] && { CHROME_BIN="$_b"; break; }
          done
        fi
        [[ -n "${CHROME_BIN:-}" ]] || err "PDF engine 'chrome' selected but no Chrome binary found."
        info "Using: $CHROME_BIN"
        timer_start
        # Headless Chrome PDF flags: page size + margins come from CSS @page,
        # virtual-time-budget gives the renderer enough time to lay out and
        # paint every base64-embedded figure before the PDF snapshot is taken.
        "$CHROME_BIN" \
          --headless --disable-gpu \
          --print-to-pdf="$PDF_OUT" \
          --print-to-pdf-no-header \
          --no-pdf-header-footer \
          --no-sandbox \
          --hide-scrollbars \
          --virtual-time-budget=20000 \
          --run-all-compositor-stages-before-draw \
          "file://$HTML_OUT" 2>&1 | tail -3 \
          || warn "Chrome PDF generation returned non-zero exit (PDF may still be valid)."
        timer_end
        ;;
      wkhtmltopdf)
        timer_start
        wkhtmltopdf \
          --enable-local-file-access \
          --page-size A4 \
          --margin-top 18mm --margin-bottom 18mm \
          --margin-left 18mm --margin-right 18mm \
          --print-media-type \
          "$HTML_OUT" "$PDF_OUT" 2>&1 | tail -5 \
          || warn "wkhtmltopdf returned non-zero exit (PDF may still be valid)."
        timer_end
        ;;
      cupsfilter)
        timer_start
        cupsfilter -m application/pdf "$HTML_OUT" > "$PDF_OUT" 2>/dev/null \
          || warn "cupsfilter returned non-zero exit (PDF may still be valid)."
        timer_end
        ;;
      *)
        err "Unknown --pdf-engine: $CHOSEN_ENGINE
  -> Allowed: chrome, wkhtmltopdf, cupsfilter"
        ;;
    esac

    if [[ -s "$PDF_OUT" ]]; then
      PDF_SIZE_MB="$(du -m "$PDF_OUT" | cut -f1)"
      ok "PDF report ready: $PDF_OUT  (${PDF_SIZE_MB} MB)"
    else
      warn "PDF was not produced or is empty.  HTML is still available."
      PDF_OUT=""
    fi
  fi
fi

# ─────────────────────────────────────────────────────────────────────────────
step "5/5 — Write provenance file + final summary"
# ─────────────────────────────────────────────────────────────────────────────
cat > "$INFO_OUT" << INFO
# ============================================================================
# mbx_final_report_info.txt
# Generated by mbx_final_report.sh   (step 18)
# Date : $NOW
# ============================================================================

# ── Inputs ────────────────────────────────────────────────────────────────────
MBX_OUTPUTS_DIR=$MBX_OUT_DIR
R_VERSION=$R_VERSION

# ── Outputs ───────────────────────────────────────────────────────────────────
REPORT_DIR=$REPORT_DIR
HTML_REPORT=$HTML_OUT
HTML_REPORT_SIZE_MB=$HTML_SIZE_MB
PDF_REPORT=$PDF_OUT
PDF_ENGINE_USED=${CHOSEN_ENGINE:-NONE}

# ── Coverage ──────────────────────────────────────────────────────────────────
STEPS_DONE=${STEPS_DONE[*]:-(none)}
STEPS_PARTIAL=${STEPS_PARTIAL[*]:-(none)}
STEPS_MISSING=${STEPS_MISSING[*]:-(none)}
HIGHEST_STEP_RUN=$STEP_MAX_DONE
INFO

ok "Provenance written: $INFO_OUT"

# Cleanup discovery JSON (kept until now for re-runs)
rm -f "$DISC_JSON"

sep
echo ""
echo "  ╔══════════════════════════════════════════════════════════════╗"
echo "  ║  Final report generation complete                            ║"
echo "  ╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "  Pipeline steps executed: ${#STEPS_DONE[@]} of 18 (highest = step $STEP_MAX_DONE)"
[[ ${#STEPS_PARTIAL[@]} -gt 0 ]] && echo "  PARTIAL steps          : ${STEPS_PARTIAL[*]}"
[[ ${#STEPS_MISSING[@]} -gt 0 ]] && echo "  Steps NOT run          : ${STEPS_MISSING[*]}"
echo ""
echo "  HTML report : $HTML_OUT"
echo "                ($HTML_SIZE_MB MB; open in any browser)"
if [[ -n "$PDF_OUT" && -s "$PDF_OUT" ]]; then
  echo "  PDF report  : $PDF_OUT"
fi
echo "  Provenance  : $INFO_OUT"
echo ""
if [[ ${#STEPS_MISSING[@]} -gt 0 ]]; then
  echo "  To add more analyses to a future report, run any of:"
  for sid in "${STEPS_MISSING[@]}"; do
    while IFS='|' read -r SID SDIR SINFO SSCRIPT SDESC; do
      [[ "$SID" == "$sid" ]] && echo "    $SSCRIPT $MBX_OUT_DIR     # step $SID: $SDESC"
    done <<< "$PIPELINE_MANIFEST"
  done
  echo "  ...then re-run: mbx_final_report.sh $MBX_OUT_DIR --force-rerun"
fi
echo ""
sep
