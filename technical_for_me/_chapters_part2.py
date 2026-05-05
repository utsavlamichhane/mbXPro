"""
Chapters 4-12 (DADA2 parameter finder + DADA2 run + classifiers + taxonomy + ezclean/viz/stat + pre-diversity).
"""

from pathlib import Path
import sys
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from _build_technical import (
    DIAGRAMS_DIR, HERE,
    make_doc_with_cover,
    section_purpose, section_inputs_outputs, section_parameters,
    section_algorithm, section_flow, section_edge_cases,
    section_impl_notes, section_testing,
    add_h1, add_h2, add_h3, add_para, add_bullet, add_code, add_kv, add_table,
)
from _flowchart import render_flowchart


def build_04_dada2_param_finder():
    NAME = "4_mbx_dada2_parameter_finder.docx"
    doc = make_doc_with_cover(
        "QC the demux artifact and compute optimal DADA2 parameters",
        "Pipeline step 3 -- the most parameter-intensive auto-tuning stage.",
        "mbx_dada2_parameter_finder.sh + create_dada2_parameters_txt.sh",
        step_id="3",
    )

    section_purpose(doc, (
        "This step is split across two scripts that run as a unit. "
        "mbx_dada2_parameter_finder.sh is the orchestrator: it runs seven "
        "QIIME2 QC commands then delegates to create_dada2_parameters_txt.sh, "
        "an embedded Python algorithm that derives the four DADA2 parameters "
        "(--p-trunc-len-f, --p-trunc-len-r, --p-trim-left-f, --p-trim-left-r) "
        "from quality scores in the demux artifact. "
        "Both scripts are documented together because they are inseparable in "
        "practice."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<artifact.qza>", "From step 2", "QIIME2 SampleData[PairedEndSequencesWithQuality] artifact", "yes"],
            ["0_primer_handling/mbx_primer_info.txt", "From step 0 (auto-discovered)", "Drives DETECTION_STATUS-based trim-left rule (NEW in v1.2.0)", "no (auto)"],
            ["--forward-primer SEQ", "User flag", "Manual override of detected primer", "no"],
            ["--reverse-primer SEQ", "User flag", "Manual override of detected primer", "no"],
            ["--assume-primer-length N", "User flag (NEW in v1.2.0)", "Force --p-trim-left-f / --p-trim-left-r to N bp", "no"],
            ["--primer-info FILE", "User flag (NEW in v1.2.0)", "Explicit path to mbx_primer_info.txt", "no"],
            ["--amplicon-length N", "User flag", "Integer (insert length hint)", "no"],
            ["--min-overlap N", "User flag", "Integer (default 12)", "no"],
            ["--max-ee-f F", "User flag", "Float (default 2.0)", "no"],
            ["--max-ee-r F", "User flag", "Float (default 2.0)", "no"],
        ],
        outputs_table=[
            ["3_dada2_parameters/exported_demux_summary/", "qiime tools export of the demux QZV (forward-seven-number-summaries.csv etc.)", "step 4 / report"],
            ["3_dada2_parameters/extracted_qza/", "qiime tools extract of the QZA (raw FASTQ access)", "audit/debug"],
            ["3_dada2_parameters/extracted_qzv/", "qiime tools extract of the QZV", "audit/debug"],
            ["3_dada2_parameters/demux_read_counts.qza", "Per-sample read count artifact", "step 4 / report"],
            ["3_dada2_parameters/demux_read_counts.qzv", "Visualization of the above", "report"],
            ["3_dada2_parameters/dada2_parameters.txt", "Recommended DADA2 parameters + ranked alternatives + Input artifact path + provenance", "step 4 (mbx_dada2_run.sh)"],
        ],
        upstream_consumers=[
            "Step 4 (mbx_dada2_run.sh) reads dada2_parameters.txt to get trunc-len-f/r, trim-left-f/r, and the path to the input artifact.",
            "Step 18 (mbx_final_report.sh) embeds dada2_parameters.txt verbatim plus the per-sample read-count summary.",
        ],
    )

    section_parameters(doc,
        intro=(
            "Two layers of parameters: the bash wrapper's (mostly pass-through) "
            "flags, and the Python algorithm's environment variables. The Python "
            "internals are ALL exposed via env vars so power users can tune "
            "without editing source."
        ),
        params_table=[
            ["--forward-primer SEQ", "(empty)",
             "Explicit override.  If provided, --p-trim-left-f = len(primer)."],
            ["--reverse-primer SEQ", "(empty)",
             "Explicit override.  If provided, --p-trim-left-r = len(primer)."],
            ["--assume-primer-length N", "(empty)",
             "NEW in v1.2.0.  Force --p-trim-left-f and --p-trim-left-r to exactly N.  Use --assume-primer-length 0 to disable trimming when primers were already removed by the sequencing facility."],
            ["DETECTION_STATUS auto-rule",
             "(applied when no explicit override is given)",
             "NEW in v1.2.0.  The script reads 0_primer_handling/mbx_primer_info.txt and applies: DETECTED/USER_SUPPLIED -> primer-length trim; TRIMMED -> 0 (was incorrectly 20 in <=1.1.x); UNKNOWN/missing-file -> 20 (defensive)."],
            ["--amplicon-length N", "(unset)",
             "If supplied, the Python engine SKIPS k-mer alignment-based insert estimation and uses N directly as median/p99/requirement. Useful when the user knows the exact 16S region length (e.g. V4 = 253 nt)."],
            ["--min-overlap N (MIN_OVERLAP)", "12",
             "Minimum bp the post-truncation forward and reverse reads must overlap to merge. The candidate ranker enforces (trunc_f + trunc_r - insert_requirement) >= MIN_OVERLAP. 12 is DADA2's recommended default."],
            ["--max-ee-f F (MAX_EE_F)", "2.0",
             "Per-read maximum expected errors after truncation. EE = sum_i 10^(-Q_i/10). The cumulative EE walk identifies the longest length where EE <= max-ee. 2.0 is DADA2 default; raise to 5.0 for very low-quality data."],
            ["--max-ee-r F (MAX_EE_R)", "2.0", "Same logic for reverse reads."],
            ["READS_PER_SAMPLE (env)", "1000",
             "Number of read pairs sampled PER SAMPLE for the parameter search. 1000 x 20 samples = 20k pairs is empirically enough to converge on the same answer as scanning the full FASTQ."],
            ["INSERT_SAMPLE_PAIRS (env)", "2500",
             "Subset of the sampled pairs that go through k-mer-based pair alignment for insert-length estimation. Capped because the alignment is the slowest single step."],
            ["MIN_PAIR_PASS_RATE (env)", "0.80",
             "After candidate enumeration, candidates whose suffix-matrix pair-pass rate is below 0.80 are dropped first. Relaxed automatically if no candidate passes."],
            ["MIN_FORWARD_PASS_RATE (env)", "0.90",
             "Forward-only pass-rate floor."],
            ["MIN_REVERSE_PASS_RATE (env)", "0.84",
             "Reverse-only pass-rate floor (lower than forward because Illumina R2 is intrinsically lower-quality on the 3' end)."],
            ["LENGTH_SLACK (env)", "1",
             "Within the eligible candidates, the script considers all candidates whose effective_f + effective_r is within LENGTH_SLACK of the maximum. Among those it tie-breaks on pair-pass rate and balance. SLACK=1 means literally 'longest combined length OR one bp shorter'."],
            ["TOP_CANDIDATES (env)", "20",
             "Number of ranked alternatives to print at the bottom of dada2_parameters.txt for the auditor."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 The bash wrapper -- 7 QIIME2 calls", [
            ("p", "Step 1 (peek) verifies the artifact type. Steps 2-3 export and extract the QZA + QZV (two different operations: export gives human-readable CSV/HTML inside the QZV, extract gives raw FASTQ inside the QZA). Step 4 calls `qiime demux tabulate-read-counts` to produce a per-sample read-count artifact, step 5 visualizes it, step 6 exports + extracts that artifact, step 7 hands off to the Python parameter recommender."),
            ("p", "All 7 outputs land in subdirectories of 3_dada2_parameters/. The user can audit any of them post-hoc. The QZVs are dragged into https://view.qiime2.org for interactive inspection."),
        ]),
        ("4.2 The Python algorithm -- conceptual overview", [
            ("p", "The recommender solves an optimization problem:"),
            ("code",
             "Given paired reads with quality scores, find\n"
             "    (trunc_f, trunc_r, trim_f, trim_r)\n"
             "such that:\n"
             "    1. trunc_f + trunc_r - estimated_insert_length >= min_overlap\n"
             "    2. fraction(reads passing EE<=max_ee at trunc_f) >= min_forward_pass\n"
             "    3. fraction(reads passing EE<=max_ee at trunc_r) >= min_reverse_pass\n"
             "    4. fraction(pairs where BOTH pass) >= min_pair_pass\n"
             "and we PREFER the longest effective length (= trunc - trim_left)."),
        ]),
        ("4.3 Stage A -- Manifest scan + read sampling", [
            ("p", "QzaFastqSource() opens the .qza zipfile, locates data/MANIFEST inside it, and parses the per-sample direction (forward/reverse) and filename. sample_pairs() then opens R1 and R2 in lockstep via gzip.open and yields ReadPair(forward_seq, reverse_seq, forward_qual, reverse_qual). Quality strings are converted to integer arrays (q = ord(ch) - 33) up-front so the inner loop is pure integer math."),
            ("p", "Sampling stops at READS_PER_SAMPLE per sample (default 1000). With 20 samples that is 20k pairs -- enough to estimate quantiles with negligible Monte Carlo error."),
        ]),
        ("4.4 Stage B -- Insert length estimation (k-mer seed-and-extend)", [
            ("p", "best_pair_alignment(forward, reverse, kmer=13, min_overlap=40) implements a fast k-mer chain alignment:"),
            ("bullets", [
                "Build positions: dict mapping every 13-mer of the forward read to its starting indices.",
                "For each 13-mer in the reverse-complemented R2, look up matching positions and increment a Counter[shift] where shift = i - j.",
                "Take the 12 most common shift values; for each, compute matches/comparable identity over the implied overlap region.",
                "Score = identity * comparable. Best score wins.",
                "Reject overlaps shorter than min_overlap (40 bp).",
                "Return (insert_length, identity, overlap).",
            ]),
            ("p", "estimate_insert_length() runs alignment on up to INSERT_SAMPLE_PAIRS pairs, keeps those with identity >= 0.90, applies an IQR outlier filter (then a fallback ±50bp filter if too many were dropped), and returns median, p95, p99. The 99th percentile is used as the 'requirement' downstream -- a conservative choice that protects against unusually long inserts."),
            ("p", "If --amplicon-length is provided, the entire alignment stage is bypassed: median/p95/p99/requirement are all set to that value."),
        ]),
        ("4.5 Stage C -- Trim-left decision", [
            ("p", "decide_trim() chooses --p-trim-left-f and --p-trim-left-r:"),
            ("kv", [
                ("Primer provided", "trim_left = len(primer). primer_match_counts() additionally walks the sample reads and counts how many have <= 15% IUPAC-mismatches in the primer-length prefix -- this match count is reported to the user as a sanity check."),
                ("Primer NOT provided", "trim_left = ASSUME_PRIMER_LENGTH (default 20). The 'unknown-primer' source string is recorded in the output."),
            ]),
        ]),
        ("4.6 Stage D -- Per-read max-trunc walk", [
            ("p", "max_trunc_from_read(seq, qual, trim_left, max_ee) walks position by position from trim_left to len(qual). At each position, it accumulates expected error EE += 10^(-q/10). The largest position where EE <= max_ee is returned. The walk also short-circuits at the first 'N' base."),
            ("p", "This single 1D function is the entire quality model. It is conservative (Q40 contributes 0.0001, Q20 contributes 0.01, Q10 contributes 0.1) and matches DADA2's own filterAndTrim() expectation-based filter exactly."),
        ]),
        ("4.7 Stage E -- The suffix matrix", [
            ("p", "build_suffix_matrix() is the heart of the algorithm. It pre-computes, for every (max_f_len, max_r_len) pair, how many sampled pairs have BOTH max_trunc_f >= forward_index AND max_trunc_r >= reverse_index."),
            ("p", "Implementation:"),
            ("bullets", [
                "Create matrix[i][j] = count of pairs whose (f_max, r_max) == (i, j).",
                "Convert to suffix sums: matrix[i][j] -> sum over (i' >= i, j' >= j).",
                "This is done in two passes: first column-suffix (bottom-up) then row-suffix (right-to-left). Both passes are O(M*N).",
                "Also compute f_pass_counts[i] = number of pairs with f_max >= i, and similarly r_pass_counts[j].",
            ]),
            ("p", "Now ranking N x M candidate (trunc_f, trunc_r) combinations is O(N*M) instead of O(N*M*pairs). For typical Illumina data (max_f=300, max_r=250, 20k pairs) this is ~10^4 candidates rather than ~10^9."),
        ]),
        ("4.8 Stage F -- Candidate enumeration + ranking", [
            ("p", "rank_candidates() iterates trunc_f from trim_f+50 to max_f_len, trunc_r from trim_r+50 to max_r_len. For each (trunc_f, trunc_r):"),
            ("bullets", [
                "Compute overlap = trunc_f + trunc_r - insert_requirement.",
                "Skip if overlap < MIN_OVERLAP.",
                "Compute pair_pass = matrix[trunc_f][trunc_r] / N.",
                "Skip if pair_pass == 0.",
                "Record forward_pass = f_pass_counts[trunc_f] / N, reverse_pass = r_pass_counts[trunc_r] / N.",
                "Append a Candidate dataclass entry.",
            ]),
            ("p", "Ranking proceeds in three filtered tiers:"),
            ("kv", [
                ("Tier 1", "Eligible = candidates passing all three pass-rate floors (0.80 / 0.90 / 0.84). If empty, move to Tier 2."),
                ("Tier 2", "Eligible = candidates passing only MIN_PAIR_PASS_RATE. If empty, Tier 3."),
                ("Tier 3", "Eligible = ALL candidates (graceful degradation -- something is better than nothing)."),
            ]),
            ("p", "Within the eligible set, find max(effective_f + effective_r). Take the 'near-longest' subset = candidates within LENGTH_SLACK of that max. Sort the subset by (pair_pass DESC, |effective_f - effective_r| ASC, overlap DESC, length DESC) and pick #1."),
            ("p", "The 'balance' tiebreaker (smaller |F-R| difference preferred) is the subtle but important step: it prevents pathological recommendations like (trunc_f=300, trunc_r=50) when (200, 200) would also pass."),
        ]),
        ("4.9 Output emission", [
            ("p", "write_output() produces a self-documenting .txt with these blocks: header (input QZA, sample count, total pair count), the four parameters one-per-line in QIIME2-compatible form, primer-handling provenance, quality+overlap estimates, the ready-to-run `qiime dada2 denoise-paired` command, and a top-N candidate table for auditing."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "04_dada2_param.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.6, h=0.6, label="START\nbash wrapper", kind="term"),
        dict(id="qiime", x=5, y=10.0, w=4.6, h=0.7, label="7 QIIME2 calls\n(peek/export/extract/tabulate)"),
        dict(id="open", x=5, y=8.9, w=4.6, h=0.7, label="Python: open .qza as zipfile\nparse data/MANIFEST"),
        dict(id="sample", x=5, y=7.9, w=5.4, h=0.7, label="Sample 1000 pairs/sample\n(quality bytes converted up-front)"),
        dict(id="trim", x=2.0, y=6.5, w=3.4, h=0.9, label="decide_trim()\nprimer? len(primer)\nelse: 20", kind="decision"),
        dict(id="insert", x=7.0, y=6.5, w=3.6, h=0.9, label="estimate_insert_length()\nk-mer alignment p99\nor --amplicon-length"),
        dict(id="walk", x=5, y=5.0, w=5.6, h=0.8, label="max_trunc_from_read()\nfor every R1 + every R2\nEE = sum 10^(-q/10)"),
        dict(id="matrix", x=5, y=3.7, w=6.0, h=0.9, label="build_suffix_matrix()\n2-pass O(M*N) cumulative sum"),
        dict(id="rank", x=5, y=2.4, w=6.0, h=1.0, label="rank_candidates()\nTier 1 (strict) -> Tier 2 -> Tier 3\nlongest within LENGTH_SLACK + balance tiebreak", kind="decision"),
        dict(id="write", x=5, y=0.7, w=5.0, h=0.6, label="dada2_parameters.txt", kind="io"),
    ]
    edges = [
        dict(**{"from": "start", "to": "qiime"}),
        dict(**{"from": "qiime", "to": "open"}),
        dict(**{"from": "open", "to": "sample"}),
        dict(**{"from": "sample", "to": "trim"}),
        dict(**{"from": "sample", "to": "insert"}),
        dict(**{"from": "trim", "to": "walk"}),
        dict(**{"from": "insert", "to": "walk"}),
        dict(**{"from": "walk", "to": "matrix"}),
        dict(**{"from": "matrix", "to": "rank"}),
        dict(**{"from": "rank", "to": "write"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="DADA2 parameter recommender -- algorithm")
    section_flow(doc, flow_png, "Two-stage flow: bash QIIME2 calls (top) feed the Python recommender (middle), which writes the parameter file.")

    section_edge_cases(doc, [
        ("Insert length cannot be inferred",
         "If <8 pairs align with identity >= 0.90, estimate_insert_length() falls back to median ±50 bp filter. If still too few, it dies with a request for --amplicon-length. Common cause: reads do not actually overlap (e.g. V1-V9 long-read data)."),
        ("All candidates fail strict tier",
         "Tiers 2 and 3 progressively relax constraints. The third tier accepts ANY candidate. The user is still notified via the pass-rate columns in the candidate table."),
        ("Amplicon shorter than 2 * trim_left + 50",
         "The candidate-enumeration `for trunc_f in range(trim_f+50, max_f_len+1)` will yield no candidates. Will trigger 'No candidate truncation lengths' error -- the user must drop `--min-overlap` and/or use a shorter primer."),
        ("Mixed-length sample pools",
         "If individual samples have very different read lengths (e.g. one re-run sample), max_f_len and max_r_len reflect the global maximum. Candidates beyond an individual sample's max are still considered for that sample but yield max_trunc capped at the sample's actual length, naturally bringing the score down."),
        ("--p-trim-left equal across all samples",
         "DADA2 enforces that trim-left is identical for every sample in a run. This is a QIIME2 limitation, not ours. If samples have different primer-trim conventions, pre-process them separately."),
    ])

    section_impl_notes(doc, [
        ("dataclasses for clarity",
         "ManifestEntry, ReadPair, InsertEstimate, TrimInfo, Candidate are all @dataclass-typed. This is purely for readability: the algorithmic intent is much clearer when each function takes/returns named fields rather than tuples."),
        ("Why pre-convert quality strings to bytes",
         "Python's str -> int conversion on Illumina quality strings is expensive when done in the inner loop. We pay the cost once per pair when reading the FASTQ and use byte-array indexing in max_trunc_from_read."),
        ("The 'effective length' metric",
         "effective_f = trunc_f - trim_f. This is the actual amplicon length contributing to the merged read after primer trimming. We optimize EFFECTIVE length, not raw trunc-len, because two recommendations with the same trunc-len but different trim-left have different merged-read lengths."),
        ("Why two scripts and not one",
         "create_dada2_parameters_txt.sh is intentionally standalone -- power users can run it on any QIIME2 paired-end .qza without going through the bash wrapper's QC steps. The wrapper exists to bundle the QC steps with the parameter computation for the lay user."),
    ])

    section_testing(doc, [
        "Run with a known-good demux artifact and verify the recommended trunc-len matches DADA2's own forward/reverse FilterAndTrim recommendation (typically within 1-2 bp).",
        "Run with --amplicon-length 253 -> verify the candidate table reflects insert_requirement = 253.",
        "Run on artificially low-quality data (max-ee 2 fails for all candidates) -> verify Tier 3 fallback and that warnings appear in the report.",
        "Run with primers explicitly provided (--forward-primer/--reverse-primer) -> verify trim-left equals primer length.",
        "Run with --dry-run from the wrapper -> verify the wrapper prints all 7 commands without invoking them.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_05_dada2_run():
    NAME = "5_mbx_dada2_run.docx"
    doc = make_doc_with_cover(
        "Run DADA2 denoising and produce all summary visualizations",
        "Pipeline step 4 -- the slowest and most CPU-bound stage.",
        "mbx_dada2_run.sh",
        step_id="4",
    )

    section_purpose(doc, (
        "Read the four DADA2 parameters chosen by step 3, validate the user's "
        "metadata file (mirroring mbX R-package validation), and run "
        "`qiime dada2 denoise-paired` plus three QZV-summary commands. "
        "The result is the canonical feature table + representative-sequences + "
        "denoising-stats triple that drives every downstream step."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<dada2_parameters.txt>", "From step 3", "Plain-text k=v + ranked-candidate table", "yes"],
            ["<metadata.txt>", "User-provided", "QIIME2 metadata file (.txt/.tsv/.csv)", "yes"],
            ["--threads N", "User flag", "Integer (default = sysctl/nproc)", "no"],
            ["--dry-run", "User flag", "Boolean", "no"],
        ],
        outputs_table=[
            ["4_dada2_outputs/feature_table.qza", "FeatureTable[Frequency] (samples x ASVs)", "every later step"],
            ["4_dada2_outputs/feature_table.qzv", "Per-sample/per-feature visualization", "report"],
            ["4_dada2_outputs/representative_sequences.qza", "FeatureData[Sequence] (one DNA seq per ASV)", "step 5/6 (classifier), step 11 (tree)"],
            ["4_dada2_outputs/representative_sequences.qzv", "Tabulated seqs (BLAST-able)", "report"],
            ["4_dada2_outputs/dada2_stats.qza", "Per-sample reads-in / reads-out / merged / non-chim counts", "report"],
            ["4_dada2_outputs/dada2_stats.qzv", "Visualization of the above", "report"],
        ],
        upstream_consumers=[
            "Step 5/6 read representative_sequences.qza for classifier extract-reads + classify-sklearn.",
            "Step 7 reads feature_table.qza to build feature_table_filtered.qza (mito/chloro removed).",
            "Step 11 reads representative_sequences.qza to build the phylogenetic tree.",
            "Step 18 embeds dada2_stats.qzv-derived numbers into the report's QC table.",
        ],
    )

    section_parameters(doc,
        intro="The script is parameter-light: trunc-len/trim-left are inherited from step 3.",
        params_table=[
            ["--threads N", "auto",
             "Passed verbatim to --p-n-threads. auto resolves via `nproc` (Linux) or `sysctl -n hw.logicalcpu` (macOS). DADA2 scales sub-linearly past 8 threads -- typical sweet spot is 4-8."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Metadata validation (mirrors mbX R-package logic)", [
            ("p", "The script applies six validation rules, designed to catch every failure mode we have seen in real datasets BEFORE running the 30-90 minute DADA2 step:"),
            ("bullets", [
                "Extension must be .txt/.tsv/.csv. .xls/.xlsx is explicitly refused with an Excel-export tip.",
                "First column header (after BOM stripping + lowercasing + trim) must match one of: id, sampleid, sample id, sample-id, featureid, feature id, feature-id.",
                "If second line starts with #q2:types it is recognised as a QIIME2 directive and DATA_START moves to row 3.",
                "Sample IDs must be non-empty, unique, and have no leading/trailing whitespace.",
                "Sample IDs containing characters outside [A-Za-z0-9._-] trigger a warning (QIIME2 may reject them later).",
            ]),
            ("p", "Empty/duplicate detection uses bash sort + sort -u and comm -23 to enumerate the duplicates for the error message."),
        ]),
        ("4.2 Parameter parsing", [
            ("p", "The four DADA2 parameters and the input artifact path are extracted from dada2_parameters.txt via grep + cut. Each parameter is then sanity-checked: integer type and trunc > trim_left (DADA2 enforces this internally but a clearer error message here saves debugging time)."),
        ]),
        ("4.3 DADA2 denoising", [
            ("p", "Six QIIME2 calls are run sequentially:"),
            ("kv", [
                ("Step 3", "qiime dada2 denoise-paired with the four parameters + --p-n-threads. This is the long-running step (CPU-bound, scales with sample count and read length)."),
                ("Step 4", "qiime metadata tabulate on dada2_stats.qza -> dada2_stats.qzv."),
                ("Step 5", "qiime feature-table summarize on feature_table.qza + metadata -> feature_table.qzv (the per-sample frequency table the user looks at to choose rarefaction depth)."),
                ("Step 6", "qiime feature-table tabulate-seqs on representative_sequences.qza -> representative_sequences.qzv (BLAST-clickable)."),
            ]),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "05_dada2_run.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.4, h=0.6, label="START", kind="term"),
        dict(id="parse", x=5, y=10.0, w=4.4, h=0.7, label="Parse args (params.txt + metadata)"),
        dict(id="metaval", x=5, y=8.9, w=5.4, h=0.9, label="Validate metadata\nheader / empty / dup / whitespace / chars", kind="decision"),
        dict(id="paramparse", x=5, y=7.6, w=5.4, h=0.7, label="Extract trunc/trim-left + input artifact path"),
        dict(id="sanity", x=5, y=6.6, w=4.4, h=0.7, label="trunc > trim_left ?", kind="decision"),
        dict(id="dada2", x=5, y=5.4, w=5.6, h=0.9, label="qiime dada2 denoise-paired\n(longest step, --p-n-threads)", kind="io"),
        dict(id="stats", x=5, y=4.0, w=4.6, h=0.6, label="metadata tabulate -> stats.qzv"),
        dict(id="ftsumm", x=5, y=3.0, w=5.4, h=0.7, label="feature-table summarize -> ft.qzv\n(uses --m-sample-metadata-file)"),
        dict(id="seqs", x=5, y=1.9, w=5.4, h=0.7, label="feature-table tabulate-seqs -> rep_seqs.qzv"),
        dict(id="end", x=5, y=0.8, w=2.0, h=0.6, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "parse"}),
        dict(**{"from": "parse", "to": "metaval"}),
        dict(**{"from": "metaval", "to": "paramparse"}),
        dict(**{"from": "paramparse", "to": "sanity"}),
        dict(**{"from": "sanity", "to": "dada2"}),
        dict(**{"from": "dada2", "to": "stats"}),
        dict(**{"from": "stats", "to": "ftsumm"}),
        dict(**{"from": "ftsumm", "to": "seqs"}),
        dict(**{"from": "seqs", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_dada2_run.sh -- control flow")
    section_flow(doc, flow_png, "Validation gates run BEFORE the slow DADA2 step.")

    section_edge_cases(doc, [
        ("UTF-8 BOM in metadata",
         "Excel-exported files start with U+FEFF. The header validator strips it via `sed 's/^\\xef\\xbb\\xbf//'`. Without this stripping, the BOM is part of the first column header, which fails to match any allowed alias."),
        ("Sample IDs in metadata that are not in the artifact",
         "DADA2 itself does not check this -- step 5 (feature-table summarize) will silently drop those samples. We pass --m-sample-metadata-file so the QZV exposes the mismatch visually."),
        ("Threads = 0",
         "Treated identically to 'auto'. Resolves via nproc/sysctl."),
    ])

    section_impl_notes(doc, [
        ("Why pre-validate metadata?",
         "DADA2 is the slowest step (30+ min). A metadata typo discovered AFTER DADA2 finishes is hours of wasted time. Pre-validation catches every metadata-format error in seconds."),
        ("BOM stripping",
         "We use the literal byte sequence `\\xef\\xbb\\xbf` (UTF-8 BOM) rather than relying on iconv. This works on every macOS/Linux without extra dependencies."),
        ("Reading parameters with grep + cut",
         "Faster and more portable than parsing via awk/python. Each parameter line in dada2_parameters.txt has a fixed `--p-key value` format -- grep -m1 'key ' + awk '{print $2}' is sufficient."),
    ])

    section_testing(doc, [
        "Run with valid metadata + valid params -> verify all 6 outputs exist.",
        "Run with metadata containing duplicate sample IDs -> verify error names the duplicates.",
        "Run with metadata first-column header 'SampleName' (invalid) -> verify clear actionable error.",
        "Run with --threads 4 -> verify --p-n-threads 4 appears in the dada2 command.",
        "Run with a manually edited dada2_parameters.txt where trunc < trim -> verify pre-flight catches it.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_06_classifier_arranger():
    NAME = "6_mbx_classifier_arranger.docx"
    doc = make_doc_with_cover(
        "Auto-detect QIIME2/GG2 versions, download references, prepare classifier inputs",
        "Pipeline step 5 -- the orchestration prep for taxonomic classification.",
        "mbx_classifier_arranger.sh",
        step_id="5",
    )

    section_purpose(doc, (
        "Solve the practical questions of taxonomic classification BEFORE the "
        "expensive train+classify step: which Greengenes2 version to use, "
        "where to download it from, what the per-ASV length distribution is, "
        "and how to construct the three QIIME2 commands that step 6 will run."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI arg", "Directory containing 0_primer_handling/, 4_dada2_outputs/, etc.", "yes"],
            ["--gg2-version V", "User flag", "Override GG2 version (e.g. 2024.09)", "no"],
            ["--skip-download", "User flag", "Boolean", "no"],
            ["--skip-zenodo", "User flag (NEW in v1.2.0)", "Boolean -- forces local training even when a compatible pre-trained classifier exists on Zenodo", "no"],
            ["--dry-run", "User flag", "Boolean", "no"],
            ["0_primer_handling/mbx_primer_info.txt", "Auto-discover", "Read FORWARD_PRIMER_*, REVERSE_PRIMER_*, DETECTION_STATUS", "no (graceful degradation)"],
            ["4_dada2_outputs/representative_sequences.qza", "Auto-discover", "QIIME2 artifact", "yes"],
            ["https://zenodo.org/records/20021035 (Zenodo)", "Network (NEW in v1.2.0)", "Pre-trained sha256-verified Naive-Bayes classifiers", "no (graceful degradation)"],
        ],
        outputs_table=[
            ["5_classifier_working_dir/exported_rep_seqs/dna-sequences.fasta", "FASTA export of rep_seqs", "step 6 (length stats)"],
            ["5_classifier_working_dir/<ver>.backbone.full-length.fna.qza", "Downloaded GG2 reference sequences (1.5-2.5 GB) -- skipped if Zenodo classifier was downloaded", "step 6"],
            ["5_classifier_working_dir/<ver>.backbone.tax.qza", "Downloaded GG2 reference taxonomy (0.5-1 GB) -- skipped if Zenodo classifier was downloaded", "step 6"],
            ["5_classifier_working_dir/gg2_full_length_trained_classifier.qza", "NEW in v1.2.0 -- pre-trained classifier from Zenodo (when CLASSIFIER_SOURCE=zenodo) OR locally trained .qza", "step 6"],
            ["5_classifier_working_dir/MANIFEST.zenodo.tsv", "NEW in v1.2.0 -- copy of the Zenodo manifest used for the version-matching decision (only when Zenodo download was attempted)", "audit"],
            ["5_classifier_working_dir/length_summary.txt", "Min / max / total ASV lengths", "step 6"],
            ["5_classifier_working_dir/mbx_classifier_run_info.txt", "Machine-readable: every path + ready-to-run commands + CLASSIFIER_SOURCE + Zenodo provenance", "step 6"],
        ],
        upstream_consumers=[
            "Step 6 (mbx_classifier_run.sh) reads mbx_classifier_run_info.txt for everything it needs.",
            "Step 7 (mbx_taxonomy_run.sh) reads FEATURE_TABLE_QZA path from the same info file.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--gg2-version V", "(auto)",
             "Override the auto-detected GG2 version. Auto rule: QIIME2 >= 2024.5 -> GG2 2024.09 (sklearn 1.4.2 compatible). QIIME2 < 2024.5 -> GG2 2022.10 (sklearn 0.24.1)."],
            ["--skip-download", "off",
             "Assume the .qza files are already on disk in the expected paths (e.g. shared install)."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 QIIME2 version detection", [
            ("p", "Run `qiime info | grep 'QIIME 2 release' | grep -oE '[0-9]{4}\\.[0-9]+'`. Fallback to grep 'version' if the format changes. Default to '2025.4' if all parsing fails (a defensive guess that picks the modern GG2)."),
        ]),
        ("4.2 GG2 version selection", [
            ("p", "Compare the major.minor of the detected QIIME2 version. The threshold (2024.5) is significant because that release shipped the sklearn 1.4.2 stack -- mismatched sklearn versions between training and classification trigger 'incompatible classifier' errors."),
        ]),
        ("4.3 Disk-space guard", [
            ("p", "Run `df -k <target>` and parse `awk 'NR==2 {print $4}'` (kilobytes). Convert to GB. Warn (but do not fail) if < 4 GB available -- the user might know they have an external drive elsewhere."),
        ]),
        ("4.4 Download with resume support", [
            ("p", "Prefer wget -c (continues partial downloads). Fall back to `curl -L --retry 3 --continue-at -`. Both support HTTP range requests, which matters for ~2 GB downloads on flaky networks."),
        ]),
        ("4.5 Length statistics", [
            ("p", "Single-pipeline awk extracts per-ASV sequence lengths from dna-sequences.fasta:"),
            ("code",
             "awk '/^>/ {if (seq) print length(seq); seq=\"\"; next} {seq=seq$0} END {print length(seq)}' \\\n"
             "    \"$FASTA_PATH\" | sort -n | \\\n"
             "awk 'NR==1 {min=$1} {max=$1; count++} END {\n"
             "  printf \"Min Length: %d\\nMax Length: %d\\nTotal ASVs: %d\\n\", min, max, count\n"
             "}' > length_summary.txt"),
            ("p", "min/max are passed to `qiime feature-classifier extract-reads --p-min-length / --p-max-length`. This narrows the GG2 reference to roughly the same length distribution as the user's ASVs, which improves classification accuracy at genus/species level."),
        ]),
        ("4.6 CLASSIFIER_MODE auto-fallback", [
            ("p", "The arranger decides one of two modes BEFORE writing the run-info file, and embeds that decision as CLASSIFIER_MODE=region-specific or CLASSIFIER_MODE=full-length. As of v1.2.0 it ALSO honours DETECTION_STATUS=TRIMMED as a strict force-to-full-length signal (a TRIMMED status means the primers are not in the reads even if they happen to be listed in the info file):"),
            ("code",
             "PRIMERS_AVAILABLE=false\n"
             "[[ \"$FWD_PRIMER\" != \"None\" && \"$REV_PRIMER\" != \"None\" \\\n"
             "   && -n \"$FWD_PRIMER\" && -n \"$REV_PRIMER\" ]] && PRIMERS_AVAILABLE=true\n"
             "\n"
             "# v1.2.0 addition: a TRIMMED/UNKNOWN status forces full-length\n"
             "case \"$DETECTION_STATUS\" in\n"
             "  TRIMMED|UNKNOWN|\"\") PRIMERS_AVAILABLE=false ;;\n"
             "esac\n"
             "\n"
             "if $PRIMERS_AVAILABLE; then\n"
             "  CLASSIFIER_MODE=\"region-specific\"\n"
             "  TRAINED_CLASSIFIER_QZA=\"${CLASSIFIER_DIR}/gg2_trained_classifier.qza\"\n"
             "else\n"
             "  CLASSIFIER_MODE=\"full-length\"\n"
             "  TRAINED_CLASSIFIER_QZA=\"${CLASSIFIER_DIR}/gg2_full_length_trained_classifier.qza\"\n"
             "fi"),
            ("p", "Different filenames so both classifiers can co-exist on disk and be reused independently across projects (region-specific is keyed to a primer set, full-length is reusable for ANY primer set / V-region)."),
        ]),
        ("4.7 NEW in v1.2.0 -- Zenodo classifier registry + sha256-verified download", [
            ("p", "When CLASSIFIER_MODE=full-length, the arranger tries to download a pre-trained classifier from https://zenodo.org/records/20021035 BEFORE downloading the 2 GB GG2 backbone. The Zenodo bundle contains 8 pre-trained .qza files, one per QIIME2 release: 2023.2, 2023.5, 2023.7, 2023.9, 2024.2 (sklearn 0.24.x family) and 2024.5, 2024.10, 2025.4 (sklearn 1.4.x family)."),
            ("p", "Selection algorithm: (1) try EXACT QIIME2 version match for guaranteed pickle compatibility; (2) if no exact match, pick the closest release in the SAME sklearn major.minor family; (3) if no compatible release exists, do not attempt download."),
            ("p", "The script downloads MANIFEST.tsv from Zenodo at run-time so newly-uploaded versions are automatically discovered (with the embedded ZENODO_REGISTRY in mbx_classifier_arranger.sh as a hard-coded fallback when the network is unavailable)."),
            ("p", "Integrity: every downloaded .qza is sha256-verified against the manifest entry; mismatches trigger immediate fall-through to local training. The artifact type is also confirmed via `qiime tools peek` -- it must be TaxonomicClassifier."),
            ("p", "If ANY step fails (network down, no compatible match, sha256 mismatch, wrong artifact type), CLASSIFIER_SOURCE stays at 'local-training' and the script proceeds to download GG2 + write a 'train locally' run-info file. The pipeline never aborts because of Zenodo issues."),
        ]),
        ("4.8 Run-info file emission", [
            ("p", "The output run-info file has TWO sections: a machine-readable key=value block (read by step 6/7), and a human-readable 'ready-to-run commands' block at the bottom (Step A: extract-reads, Step B: fit-classifier, Step C: classify-sklearn, Step D: tabulate, Step E: barplot)."),
            ("p", "Three new keys in v1.2.0: CLASSIFIER_SOURCE (zenodo|cached|local-training|local-training-fallback), SCIKIT_LEARN_VERSION/SCIKIT_LEARN_FAMILY (recorded for audit), ZENODO_QIIME2_USED + ZENODO_FILENAME + ZENODO_SHA256_EXPECTED + ZENODO_SHA256_ACTUAL + ZENODO_NOTE (full provenance of the Zenodo download attempt)."),
            ("p", "When CLASSIFIER_SOURCE=zenodo|cached, both Step A and Step B in the ready-to-run command block are replaced with explanatory 'AUTO-SKIPPED' notes -- the trained classifier already exists. Step 6 reads CLASSIFIER_SOURCE and skips both extract-reads and fit-classifier accordingly."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "06_classifier_arranger.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.6, h=0.6, label="START", kind="term"),
        dict(id="qver", x=5, y=10.1, w=4.6, h=0.7, label="qiime info -> QIIME2 version"),
        dict(id="ggver", x=5, y=9.0, w=5.6, h=0.9, label="GG2 version select\n>=2024.5 -> 2024.09\n<2024.5 -> 2022.10", kind="decision"),
        dict(id="primer", x=2.0, y=7.6, w=3.4, h=0.9, label="Read primer_info.txt\n(graceful if missing)"),
        dict(id="fa", x=7.0, y=7.6, w=3.6, h=0.9, label="Verify rep_seqs.qza\n(REQUIRED)"),
        dict(id="disk", x=5, y=6.4, w=4.6, h=0.7, label="Disk space check"),
        dict(id="dl", x=5, y=5.3, w=5.6, h=0.8, label="wget -c (or curl -C -)\nGG2 fna + tax", kind="io"),
        dict(id="export", x=5, y=4.1, w=5.0, h=0.7, label="qiime tools export rep_seqs"),
        dict(id="awk", x=5, y=3.0, w=5.4, h=0.7, label="awk -> min/max/count -> length_summary.txt"),
        dict(id="info", x=5, y=1.8, w=5.4, h=0.8, label="Write mbx_classifier_run_info.txt\n(KV + ready commands)", kind="io"),
        dict(id="end", x=5, y=0.7, w=2.0, h=0.5, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "qver"}),
        dict(**{"from": "qver", "to": "ggver"}),
        dict(**{"from": "ggver", "to": "primer"}),
        dict(**{"from": "ggver", "to": "fa"}),
        dict(**{"from": "primer", "to": "disk"}),
        dict(**{"from": "fa", "to": "disk"}),
        dict(**{"from": "disk", "to": "dl"}),
        dict(**{"from": "dl", "to": "export"}),
        dict(**{"from": "export", "to": "awk"}),
        dict(**{"from": "awk", "to": "info"}),
        dict(**{"from": "info", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_classifier_arranger.sh -- control flow")
    section_flow(doc, flow_png, "Decisions are made up-front (versions, paths) so step 6 has a deterministic plan.")

    section_edge_cases(doc, [
        ("QIIME2 version cannot be parsed",
         "Defaults to 2025.4 -> picks GG2 2024.09. The user can always override with --gg2-version."),
        ("Primer info file missing or DETECTION_STATUS=TRIMMED/UNKNOWN",
         "FWD_PRIMER/REV_PRIMER are set to 'None' and CLASSIFIER_MODE is set to 'full-length'. The Step A extract-reads block in run_info.txt becomes an 'AUTO-SKIPPED' note. The script then attempts a Zenodo download (v1.2.0+); if that succeeds Step B is also auto-skipped. If it fails, Step 6 trains Naive-Bayes directly on the full GG2 backbone -- no user intervention. The pipeline runs end-to-end without error."),
        ("Zenodo unreachable / DNS down / firewall blocks zenodo.org",
         "_zenodo_fetch_manifest returns non-zero, the script falls back to the embedded ZENODO_REGISTRY hard-coded in mbx_classifier_arranger.sh, then attempts the actual file download. If the download itself fails, CLASSIFIER_SOURCE stays at 'local-training' and the GG2 backbone is downloaded for local training -- exactly as in mbX Pro <= 1.1.x."),
        ("Zenodo file downloaded but sha256 does not match manifest",
         "The temp file is deleted, ZENODO_NOTE is set to 'sha256 mismatch', and the script falls through to local training. Never silently uses a corrupted artifact."),
        ("Existing GG2 download in target dir",
         "The download function checks `[[ -f \"$dest\" ]]` first. Idempotent re-run: re-uses the previous download."),
        ("Existing trained classifier .qza in target dir",
         "CLASSIFIER_SOURCE is set to 'cached' and the Zenodo step is skipped. Delete the file to force a fresh Zenodo download or local training."),
        ("`df` output format differs across BSD/GNU",
         "We use -k for portability (kilobytes) and parse `awk 'NR==2 {print $4}'`. Tested on macOS BSD df and Linux GNU df."),
    ])

    section_impl_notes(doc, [
        ("Why this script and not just bake commands into step 6?",
         "Separating 'plan' from 'execute' makes the pipeline auditable: a researcher can read mbx_classifier_run_info.txt and verify every parameter without running the slow training step. It also makes the trained classifier reusable across projects with the same primers + GG2 version."),
        ("Resumable downloads",
         "wget -c re-uses partial downloads byte-for-byte via HTTP Range requests. Curl --continue-at - has the same effect. Both honour the GG2 server's Last-Modified+If-Range semantics."),
    ])

    section_testing(doc, [
        "Run with QIIME2 2025.4 active -> verify GG2 2024.09 is selected.",
        "Pre-create empty GG2 .qza files -> verify --skip-download mode picks them up.",
        "Hand-delete primer_info.txt -> verify run_info.txt has CLASSIFIER_MODE=full-length and an 'AUTO-SKIPPED' Step A note pointing at the full backbone for Step B.",
        "Run with --dry-run -> verify no downloads, no FASTA export.",
        "Verify length_summary.txt min/max match the actual rep_seqs.fasta lengths.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_07_classifier_run():
    NAME = "7_mbx_classifier_run.docx"
    doc = make_doc_with_cover(
        "Train a Naive Bayes classifier on a primer-extracted GG2 region, then classify ASVs",
        "Pipeline step 6 -- one of the longest (20-90 min) but fully cacheable.",
        "mbx_classifier_run.sh",
        step_id="6",
    )

    section_purpose(doc, (
        "Execute the taxonomy classification plan that step 5 wrote into "
        "mbx_classifier_run_info.txt. The script honours TWO independent "
        "selectors, both auto-decided by step 5:\n\n"
        "CLASSIFIER_MODE -- the science decision:\n"
        "* region-specific : primers are known -> extract amplicon-region GG2 "
        "reads -> train a Naive Bayes classifier on those reads -> classify the "
        "user's ASVs. Slightly more precise at species level.\n"
        "* full-length : primers were absent -> extract-reads is SKIPPED -> "
        "Naive-Bayes is trained directly on the full GG2 backbone (or "
        "downloaded from Zenodo, see below) -> classify the user's ASVs.\n\n"
        "CLASSIFIER_SOURCE -- the supply decision (NEW in v1.2.0):\n"
        "* zenodo : a pre-trained, sha256-verified classifier was downloaded "
        "from https://zenodo.org/records/20021035 in step 5. BOTH extract-reads "
        "AND fit-classifier are SKIPPED -- the script jumps straight to "
        "classify-sklearn, saving 30-90 minutes per run.\n"
        "* cached : a previously produced classifier is being reused.\n"
        "* local-training : the classifier will be trained on the fly here.\n"
        "* local-training-fallback : Zenodo classifier was tried but rejected "
        "by classify-sklearn (rare sklearn pickle incompatibility); the "
        "downloaded file was deleted, GG2 was re-fetched, and the classifier "
        "was retrained locally -- all without aborting the pipeline.\n\n"
        "The script is fully idempotent in EVERY combination of the above. "
        "The dual-selector design is what makes the pipeline robust to "
        "(a) pre-trimmed primers and (b) Zenodo failures."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Directory containing 5_classifier_working_dir/", "yes"],
            ["--classifier <path>", "User flag", "Pre-built .qza classifier (skip extract+train)", "no"],
            ["--confidence F", "User flag", "Naive Bayes posterior probability cutoff (default 0.7)", "no"],
            ["--force-extract / --force-train / --force-classify", "User flags", "Boolean toggles to overwrite existing outputs", "no"],
            ["--skip-extract", "User flag", "Boolean", "no"],
        ],
        outputs_table=[
            ["5_classifier_working_dir/gg2_extracted_reads.qza", "Region-trimmed GG2 reference (region-specific mode only; REUSABLE)", "step 6 fit-classifier"],
            ["5_classifier_working_dir/gg2_trained_classifier.qza", "V-region NB classifier (region-specific mode; REUSABLE)", "step 6 classify-sklearn"],
            ["5_classifier_working_dir/gg2_full_length_trained_classifier.qza", "Full-length NB classifier (full-length mode; REUSABLE across ANY primer set)", "step 6 classify-sklearn"],
            ["6_classifier_taxonomy/taxonomy.qza", "FeatureData[Taxonomy] -- per-ASV taxonomy + confidence", "step 7 (tax csv export)"],
            ["6_classifier_taxonomy/taxonomy.qzv", "Tabulated taxonomy", "report"],
        ],
        upstream_consumers=[
            "Step 7 reads taxonomy.qza for `qiime taxa filter-table` (mito/chloro removal) and `qiime taxa barplot`.",
            "Step 14 (ANCOMBC2), 15 (PICRUSt2), 17 (networks) read the level-*.csv files derived from this taxonomy.",
        ],
    )

    section_parameters(doc,
        intro="Idempotency flags are the heart of this script -- each major step has a per-step on-disk sentinel.",
        params_table=[
            ["--classifier PATH", "(unset)",
             "Use a pre-built classifier .qza directly. Skips both extract-reads and fit-classifier-naive-bayes. Useful when reusing a classifier from a previous project (same primers + GG2 version)."],
            ["--confidence F", "0.7",
             "qiime feature-classifier classify-sklearn --p-confidence. Posterior probability threshold above which a taxonomic label is assigned. 0.7 = QIIME2 default. Lower (0.5-0.6) yields more assignments at the cost of more false positives. Higher (0.8-0.9) yields fewer but more reliable assignments."],
            ["--force-extract", "off", "Re-run extract-reads even if gg2_extracted_reads.qza exists (e.g. primer changed)."],
            ["--force-train",   "off", "Re-run fit-classifier-naive-bayes even if gg2_trained_classifier.qza exists."],
            ["--force-classify","off", "Re-run classify-sklearn even if taxonomy.qza exists."],
            ["--skip-extract", "off",
             "Skip extract-reads (use existing gg2_extracted_reads.qza or pre-built --classifier)."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Read run-info file", [
            ("p", "Every parameter (GG2 paths, primer sequences, ASV length range, output paths) is read from `5_classifier_working_dir/mbx_classifier_run_info.txt` via the _read_key() helper. This is the contract that lets step 5 plan and step 6 execute."),
        ]),
        ("4.2 Idempotency engine", [
            ("p", "For each of the three major steps, the rule is the same:"),
            ("code",
             "if [[ -f \"$OUTPUT_QZA\" ]] && ! $FORCE_FLAG; then\n"
             "  skipped \"step name\"\n"
             "else\n"
             "  qiime ... --o-... \"$OUTPUT_QZA\"\n"
             "fi"),
            ("p", "Each skip prints a [SKIP] line with the existing file path. Power users can `rm` a single artifact to re-run only that step."),
        ]),
        ("4.3 Extract-reads (region-specific mode only)", [
            ("p", "When CLASSIFIER_MODE=region-specific, qiime feature-classifier extract-reads trims the full-length GG2 backbone to amplicons matching --p-f-primer + --p-r-primer, with min/max length filters."),
            ("p", "When CLASSIFIER_MODE=full-length the script branches early: SKIP_EXTRACT is forced to true, this entire step is bypassed, and Step 4.4 below points at the full backbone instead."),
        ]),
        ("4.4 Fit Naive Bayes (mode-aware)", [
            ("p", "qiime feature-classifier fit-classifier-naive-bayes trains a multinomial NB on k-mer features."),
            ("p", "Region-specific mode: --i-reference-reads is gg2_extracted_reads.qza (the V-region trim from step 4.3). Runtime 20-60 min, classifier .qza ~500 MB - 1 GB."),
            ("p", "Full-length mode: --i-reference-reads is the full GG2 backbone .fna.qza directly. Runtime 30-90 min, classifier .qza ~1.5-3 GB. Both modes use 8-16 GB RAM. Both classifiers are reusable across projects -- the full-length one is reusable across ANY primer set; the region-specific one is reusable only for matching primers + GG2 version."),
            ("code",
             "if [[ \"$CLASSIFIER_MODE\" == \"full-length\" ]]; then\n"
             "  REF_QZA=\"$GG2_FNA_QZA\"   # full backbone, no extract step\n"
             "else\n"
             "  REF_QZA=\"$EXTRACTED_QZA\"  # V-region trimmed\n"
             "fi\n"
             "qiime feature-classifier fit-classifier-naive-bayes \\\n"
             "  --i-reference-reads    \"$REF_QZA\" \\\n"
             "  --i-reference-taxonomy \"$GG2_TAX_QZA\" \\\n"
             "  --o-classifier         \"$TRAINED_QZA\""),
        ]),
        ("4.5 Classify-sklearn", [
            ("p", "qiime feature-classifier classify-sklearn applies the trained model to the user's representative_sequences.qza. --p-n-jobs is set to the auto-detected core count. --p-confidence defaults to 0.7."),
        ]),
        ("4.6 Tabulate taxonomy", [
            ("p", "qiime metadata tabulate produces a clickable QZV. This is the only one of the four artifacts that ALWAYS gets re-tabulated when --force-classify fires (cheap to regenerate)."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "07_classifier_run.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.0, h=0.6, label="START", kind="term"),
        dict(id="info", x=5, y=10.0, w=4.6, h=0.7, label="Read run_info.txt"),
        dict(id="ext_ext", x=5, y=8.8, w=5.0, h=0.9, label="external classifier provided?\nor --skip-extract?", kind="decision"),
        dict(id="ext1", x=5, y=7.5, w=4.4, h=0.7, label="extracted_reads.qza exists?", kind="decision"),
        dict(id="extract", x=5, y=6.4, w=5.0, h=0.7, label="qiime feature-classifier\nextract-reads", kind="io"),
        dict(id="t1", x=5, y=5.2, w=4.4, h=0.7, label="trained_classifier.qza exists?", kind="decision"),
        dict(id="train", x=5, y=4.0, w=5.0, h=0.8, label="qiime feature-classifier\nfit-classifier-naive-bayes\n(8-16 GB RAM)", kind="io"),
        dict(id="c1", x=5, y=2.7, w=4.4, h=0.7, label="taxonomy.qza exists?", kind="decision"),
        dict(id="classify", x=5, y=1.5, w=5.0, h=0.7, label="qiime feature-classifier\nclassify-sklearn", kind="io"),
        dict(id="qzv", x=5, y=0.4, w=4.6, h=0.6, label="metadata tabulate -> taxonomy.qzv"),
    ]
    edges = [
        dict(**{"from": "start", "to": "info"}),
        dict(**{"from": "info", "to": "ext_ext"}),
        dict(**{"from": "ext_ext", "to": "ext1"}),
        dict(**{"from": "ext1", "to": "extract"}),
        dict(**{"from": "extract", "to": "t1"}),
        dict(**{"from": "t1", "to": "train"}),
        dict(**{"from": "train", "to": "c1"}),
        dict(**{"from": "c1", "to": "classify"}),
        dict(**{"from": "classify", "to": "qzv"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_classifier_run.sh -- idempotent 3-step flow")
    section_flow(doc, flow_png, "Each artifact is checked on disk before being recomputed. Force flags override.")

    section_edge_cases(doc, [
        ("Primers are 'None' but no --classifier provided (mbX Pro 1.1.0+)",
         "The script auto-detects CLASSIFIER_MODE=full-length from the run-info file and runs in full-length mode (instead of aborting as in <=1.0.x). v1.2.0+ additionally tries Zenodo first."),
        ("CLASSIFIER_SOURCE=zenodo (NEW in v1.2.0)",
         "Steps 2 and 3 (extract-reads, fit-classifier) are both skipped. Step 4 runs classify-sklearn against the downloaded pre-trained classifier. If classify-sklearn errors out (e.g. unanticipated sklearn pickle incompatibility) the script DELETES the Zenodo .qza, calls `mbx_classifier_arranger.sh --skip-zenodo` to re-fetch the GG2 backbone if it was previously skipped, retrains the classifier locally, and retries classify-sklearn. CLASSIFIER_SOURCE is rewritten to 'local-training-fallback' in run_info.txt for the report. The pipeline never aborts on a Zenodo problem."),
        ("Memory error during fit-classifier",
         "Naive Bayes feature counting needs >> the size of the reference alone. We document the 8-16 GB requirement up-front. If the user's machine cannot accommodate, the workaround is to pass --classifier <prebuilt.qza> or rely on the v1.2.0 Zenodo path (which avoids local training entirely)."),
        ("Stale extracted_reads.qza vs new primers",
         "If the user re-runs step 5 with different primers but doesn't pass --force-extract here, the old extracted reads are reused. Mitigation: --force-extract is documented prominently."),
    ])

    section_impl_notes(doc, [
        ("Reusability of intermediates",
         "extracted_reads.qza and trained_classifier.qza live in 5_classifier_working_dir/, NOT in the per-run output. This intentionally invites cross-project reuse: copy that directory to a new project's mbX_pro_outputs and step 6 will skip extract+train automatically."),
        ("Confidence threshold rationale",
         "0.7 is QIIME2's default and matches Bokulich et al. 2018 Mol Ecol Resour benchmarks. We expose --confidence so users running comparison studies can match their lab's previous setting."),
    ])

    section_testing(doc, [
        "Run with all artifacts absent -> verify all 3 commands run.",
        "Re-run -> verify all 3 are skipped.",
        "Re-run with --force-classify -> verify only classify-sklearn runs.",
        "Run with --classifier <prebuilt.qza> -> verify extract+train both skipped.",
        "Run with primers='None' and no --classifier -> verify the 3-option error fires.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_08_taxonomy_run():
    NAME = "8_mbx_taxonomy_run.docx"
    doc = make_doc_with_cover(
        "Filter mitochondria/chloroplasts, build bar plot, export 7 level CSVs",
        "Pipeline step 7 -- the taxonomy table-formatting bridge.",
        "mbx_taxonomy_run.sh",
        step_id="7",
    )

    section_purpose(doc, (
        "Convert the QIIME2 .qza taxonomy + feature table into the seven "
        "level-*.csv files that the mbX R package consumes. Filter "
        "mitochondria + chloroplasts before bar-plotting (critical for "
        "plant/soil samples) and emit the canonical filtered feature table "
        "that every downstream step uses."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Directory containing 6_classifier_taxonomy/", "yes"],
            ["<metadata.txt>", "User CLI", "QIIME2 metadata", "yes"],
            ["--exclude TERMS", "User flag", "Comma-separated taxa names (default: mitochondria,chloroplast)", "no"],
            ["--skip-filter", "User flag", "Boolean (do NOT filter)", "no"],
        ],
        outputs_table=[
            ["7_taxonomy_csv/feature_table_filtered.qza", "Mito/chloro-removed FeatureTable[Frequency]", "every later step (alpha/beta/ANCOMBC/PICRUSt/...)"],
            ["7_taxonomy_csv/taxa_bar_plots.qzv", "Multi-level bar plot QZV", "report + 8/9/10 R steps"],
            ["7_taxonomy_csv/level-1.csv ... level-7.csv", "Wide-format absolute counts (sample-id rows, taxon columns)", "step 8 (ezclean), 9 (ezviz), 10 (ezstat), 16 (ML), 17 (networks)"],
            ["7_taxonomy_csv/mbx_taxonomy_info.txt", "K=V file with all paths + R usage examples", "step 8 onwards"],
        ],
        upstream_consumers=[
            "Steps 8/9/10 use the level CSVs as direct input to mbX::ezclean/ezviz/ezstat.",
            "Step 11 reads feature_table_filtered.qza for tree+depth selection.",
            "Step 14 collapses the filtered table at multiple levels for ANCOMBC2.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--exclude TERMS", "mitochondria,chloroplast",
             "Pass-through to qiime taxa filter-table --p-exclude. Comma-separated, partial-match (any taxon string containing this substring at any level is dropped)."],
            ["--skip-filter", "off",
             "Bypass filtering -- only useful when the user has confirmed there is no host or plastid contamination."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Auto-discover inputs", [
            ("p", "FEATURE_TABLE_QZA is read from 5_classifier_working_dir/mbx_classifier_run_info.txt. taxonomy.qza is taken from 6_classifier_taxonomy/. The user only supplies the metadata.txt path explicitly because it was never stored in the run_info file (intentional -- metadata path is project-specific)."),
        ]),
        ("4.2 Metadata header validation", [
            ("p", "Same logic as step 4: first column header must match the QIIME2-valid set. We pre-check before invoking qiime to surface the error in seconds."),
        ]),
        ("4.3 Filter mitochondria + chloroplasts", [
            ("p", "qiime taxa filter-table is called with --p-exclude. This drops any feature whose taxonomy string contains 'mitochondria' or 'chloroplast' at ANY level. The result is feature_table_filtered.qza. Idempotent: if the file exists we skip."),
            ("p", "WHY this matters: in plant/soil samples, mitochondrial and plastid 16S can constitute 30-70% of reads. Without filtering, relative abundance plots are dominated by these eukaryotic sequences, and downstream stats become uninterpretable."),
        ]),
        ("4.4 Bar plot generation", [
            ("p", "qiime taxa barplot consumes the filtered table + taxonomy + metadata. Output is a single QZV that contains relative-abundance bar charts at every taxonomic level. The same QZV exports cleanly to 7 level-*.csv files."),
        ]),
        ("4.5 Level-CSV export", [
            ("p", "qiime tools export on the bar-plot QZV produces level-1.csv ... level-7.csv (named by the QIIME2 taxa-barplot internals). Each CSV is a wide table: rows = samples, columns = taxa at that level + a column per metadata variable. The `cp` step moves them out of the staging subdir into the canonical `7_taxonomy_csv/` directory."),
            ("p", "If a level cannot be exported (e.g. taxonomy did not reach species level for any feature), a warning fires but the pipeline proceeds. LEVELS_PRESENT is reported."),
        ]),
        ("4.6 Info-file emission", [
            ("p", "mbx_taxonomy_info.txt records the rank of each level (1=Domain, 2=Phylum, ..., 7=Species), every CSV path, and inline R usage examples. Steps 8+ rely on this file to find inputs."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "08_taxonomy_run.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.0, h=0.6, label="START", kind="term"),
        dict(id="parse", x=5, y=10.0, w=5.4, h=0.7, label="Read run_info.txt + verify metadata"),
        dict(id="metaval", x=5, y=8.9, w=5.4, h=0.7, label="Validate metadata header"),
        dict(id="filt", x=5, y=7.7, w=5.4, h=0.9, label="qiime taxa filter-table\n--p-exclude mitochondria,chloroplast", kind="io"),
        dict(id="bar", x=5, y=6.3, w=5.0, h=0.8, label="qiime taxa barplot\n-> taxa_bar_plots.qzv", kind="io"),
        dict(id="export", x=5, y=4.9, w=5.0, h=0.8, label="qiime tools export QZV\n-> staging/level-N.csv"),
        dict(id="cp", x=5, y=3.6, w=5.0, h=0.7, label="cp staging/level-*.csv -> 7_taxonomy_csv/"),
        dict(id="report", x=5, y=2.5, w=5.4, h=0.7, label="Report LEVELS_PRESENT (1..7)"),
        dict(id="info", x=5, y=1.4, w=5.0, h=0.7, label="Write mbx_taxonomy_info.txt", kind="io"),
        dict(id="end", x=5, y=0.4, w=2.0, h=0.5, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "parse"}),
        dict(**{"from": "parse", "to": "metaval"}),
        dict(**{"from": "metaval", "to": "filt"}),
        dict(**{"from": "filt", "to": "bar"}),
        dict(**{"from": "bar", "to": "export"}),
        dict(**{"from": "export", "to": "cp"}),
        dict(**{"from": "cp", "to": "report"}),
        dict(**{"from": "report", "to": "info"}),
        dict(**{"from": "info", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_taxonomy_run.sh -- control flow")
    section_flow(doc, flow_png, "The barplot QZV is the universal source from which all 7 level CSVs derive.")

    section_edge_cases(doc, [
        ("Sample IDs in feature table absent from metadata",
         "qiime taxa barplot will fail with 'No sample data shared'. We pre-validate metadata header but do NOT do a full sample-ID intersection -- that would be redundant since QIIME2's error is clear."),
        ("Custom taxa exclusion list",
         "--exclude lets users add 'unassigned' or 'eukaryota' in addition to the default. Each term is partial-matched against the full taxonomy string."),
        ("Less than 7 levels exported",
         "Some classifiers (e.g. SILVA) commonly stop at level 6. The script tolerates missing higher levels with a per-level warning. Subsequent scripts (8/9/10) loop over [1..7] but skip missing files."),
    ])

    section_impl_notes(doc, [
        ("Why the bar plot is the source of truth for level CSVs",
         "QIIME2 stores per-level CSVs as files INSIDE the .qzv. There is no `qiime taxa collapse` -> `feature-table summarize` -> CSV path that produces the exact same format. Going via the barplot QZV is the supported, idiomatic route."),
        ("filtered table as the canonical input downstream",
         "Every later step preferentially reads `7_taxonomy_csv/feature_table_filtered.qza` if it exists. This is captured explicitly in step 11's auto-discovery logic."),
    ])

    section_testing(doc, [
        "Run with default exclusion -> verify mitochondria/chloroplast features dropped (compare counts in feature_table.qza vs feature_table_filtered.qza).",
        "Run with --skip-filter -> verify the unfiltered table is used for the bar plot.",
        "Run with metadata containing fewer samples than the feature table -> verify QIIME2's actual error message is surfaced.",
        "Verify level-6.csv (genus) opens cleanly in mbX::ezclean(\"...\", \"metadata.txt\", \"g\").",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_09_ezclean():
    NAME = "9_mbx_ezclean_all_levels.docx"
    doc = make_doc_with_cover(
        "Install mbX R package, run ezclean() for all 7 taxonomic levels",
        "Pipeline step 8 -- the first R-based step. Establishes the system-R environment.",
        "mbx_ezclean_all_levels.sh",
        step_id="8",
    )

    section_purpose(doc, (
        "Run mbX::ezclean() once per taxonomic level (1-7). ezclean() takes "
        "the QIIME2-exported level-7.csv (which contains the full taxonomy "
        "string) and produces one cleaned, taxonomically-deduplicated XLSX "
        "per level. This script also bootstraps the entire R environment for "
        "subsequent R-based steps (9, 10, 12, 13, 14, 15, 16, 17, 18)."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Directory containing 7_taxonomy_csv/", "yes"],
            ["7_taxonomy_csv/level-7.csv", "Auto-discover", "Wide-format CSV with full taxonomy strings", "yes"],
            ["metadata.txt", "Auto-discover (path stored in mbx_taxonomy_info.txt)", "QIIME2 metadata", "yes"],
            ["--levels LIST", "User flag", "Subset to clean (default: d,p,c,o,f,g,s -- all)", "no"],
            ["--skip-install", "User flag", "Boolean (assume mbX already installed)", "no"],
        ],
        outputs_table=[
            ["8_cleaned_files/mbX_cleaned_<plural>_level-7/<plural>_level-7.xlsx", "ezclean output for each level (one per level)", "step 9, 10, 16, 17"],
            ["8_cleaned_files/mbx_ezclean_info.txt", "K=V mapping level letter -> XLSX path", "step 9, 10, 16, 17, 18"],
        ],
        upstream_consumers=[
            "Steps 9 (ezviz) and 10 (ezstat) call mbX::ezviz() / mbX::ezstat() which internally re-run ezclean() on the same input.",
            "Steps 16 (ML) and 17 (networks) read these XLSX files directly via openxlsx.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--levels LIST", "d,p,c,o,f,g,s",
             "Comma-separated single-letter codes. Each letter maps to one level (d=Domain, p=Phylum, c=Class, o=Order, f=Family, g=Genus, s=Species)."],
            ["--skip-install", "off",
             "Skip the R-package install step. Useful for HPC where root cannot write to /usr/local/R."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 The system-R discovery + env-stripping pattern", [
            ("p", "mbX is installed via Homebrew R (system-wide), NOT inside the QIIME2 conda env. Installing R inside conda corrupts QIIME2's Python ABI. But running this script INSIDE the conda env (typical user pattern) means conda's R_LIBS_USER and R_HOME poison the system R."),
            ("p", "The mitigation, used by every R-based step in the pipeline, is the `_R` wrapper:"),
            ("code",
             "_R() {\n"
             "  env -u R_LIBS_USER -u R_LIBS -u R_HOME -u R_INCLUDE_DIR \\\n"
             "      -u R_DOC_DIR -u R_SHARE_DIR -u R_PROFILE -u R_PROFILE_USER \\\n"
             "      \"$RSCRIPT_CMD\" \"$@\"\n"
             "}"),
            ("p", "This strips every R_* variable. The system R then searches its own .libPaths() correctly. Errors like \"shared object 'methods.dylib' not found\" or \"Rcpp.so not found\" disappear."),
        ]),
        ("4.2 Rscript discovery preference order", [
            ("p", "1) /opt/homebrew/bin/Rscript (Apple Silicon Homebrew). 2) /usr/local/bin/Rscript (Intel macOS Homebrew). 3) /usr/bin/Rscript (system R on Linux). 4) `command -v Rscript` as a last-resort fallback. Absolute paths come first specifically to avoid `command -v` returning conda's R inside an active env."),
        ]),
        ("4.3 R package installation", [
            ("p", "Required packages: mbX (CRAN), openxlsx, ggplot2, dplyr, tidyr, stringr, dunn.test. The script writes a temp R script (PID-based filename to avoid macOS mktemp .R-suffix bug), runs `_R \"$TMP_R\"`. Inside the R script, install.packages() is called for any missing package via setdiff(required, installed.packages()[,'Package'])."),
        ]),
        ("4.4 ezclean() invocation", [
            ("p", "For each level letter L:"),
            ("bullets", [
                "Compute the matching plural name (d->domains_or_kingdom, p->phylum, c->classes, o->orders, f->families, g->genera, s->species).",
                "Compute the output sub-dir: `8_cleaned_files/mbX_cleaned_<plural>_level-7/`.",
                "Write a temp R script that setwd()s to `8_cleaned_files/` (CRITICAL -- ezclean writes intermediate files to cwd) and calls `mbX::ezclean(level_7_csv, metadata_txt, L)`.",
                "Invoke `_R` with the temp script.",
                "Verify the expected XLSX appeared.",
            ]),
            ("p", "Each level runs in an INDEPENDENT R session. This means a failure at level 's' (species) does not abort levels 'd' through 'g'."),
        ]),
        ("4.5 Info file emission", [
            ("p", "After all 7 levels (or the user-selected subset), the script writes `mbx_ezclean_info.txt` with one CLEANED_<LETTER>=<absolute_path> line per level. Steps 9, 10, 16, 17 read these to locate XLSX inputs without globbing the directory."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "09_ezclean.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.0, h=0.6, label="START", kind="term"),
        dict(id="info", x=5, y=10.0, w=5.4, h=0.7, label="Read mbx_taxonomy_info.txt"),
        dict(id="rdisc", x=5, y=8.9, w=5.6, h=0.9, label="Discover system Rscript\n(absolute paths first)"),
        dict(id="env", x=5, y=7.5, w=5.6, h=0.9, label="Define _R() wrapper\nstrip R_LIBS_USER / R_HOME / etc."),
        dict(id="install", x=5, y=6.1, w=5.4, h=0.9, label="Install missing R packages\n(mbX, openxlsx, ...)"),
        dict(id="loop", x=5, y=4.7, w=5.4, h=0.9, label="for L in d,p,c,o,f,g,s:\n  setwd(8_cleaned_files); ezclean(level7, meta, L)", kind="decision"),
        dict(id="check", x=5, y=3.3, w=5.4, h=0.7, label="verify expected XLSX appeared"),
        dict(id="info2", x=5, y=2.2, w=5.0, h=0.7, label="Write mbx_ezclean_info.txt", kind="io"),
        dict(id="end", x=5, y=1.0, w=2.0, h=0.5, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "info"}),
        dict(**{"from": "info", "to": "rdisc"}),
        dict(**{"from": "rdisc", "to": "env"}),
        dict(**{"from": "env", "to": "install"}),
        dict(**{"from": "install", "to": "loop"}),
        dict(**{"from": "loop", "to": "check"}),
        dict(**{"from": "check", "to": "info2"}),
        dict(**{"from": "info2", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_ezclean_all_levels.sh -- control flow")
    section_flow(doc, flow_png, "Each level runs in its own R session so failures isolate.")

    section_edge_cases(doc, [
        ("level-7.csv has fewer than 7 levels populated",
         "ezclean()'s level argument controls how the taxonomy string is parsed; even if level-1.csv is the only file present, we ALWAYS pass level-7.csv because it carries the full string. ezclean parses based on the requested rank, not the file."),
        ("R installation needed but no admin access",
         "The script tries `Rscript -e 'install.packages(...)'` first. install.packages defaults to a user-writable lib path on macOS Homebrew. If that fails (HPC), the user must `--skip-install` and pre-install packages via the system administrator."),
        ("conda's Rscript shadows system Rscript",
         "Mitigated by the absolute-path preference order. Even if conda's Rscript appears first in PATH, our discovery picks the system one."),
    ])

    section_impl_notes(doc, [
        ("Why level-7.csv is the input for every level",
         "level-7.csv is the only CSV that preserves the FULL taxonomy string (Domain;Phylum;...;Species) per ASV. ezclean parses this string to whatever depth the user requested. level-2.csv (Phylum-only) does NOT have enough information to produce a Genus-level cleaned output."),
        ("setwd() before ezclean()",
         "ezclean writes ~12 intermediate XLSX files relative to the current working directory, then deletes them. Running it from the wrong cwd litters the user's home directory. We always setwd() to 8_cleaned_files/ first."),
    ])

    section_testing(doc, [
        "Run with all 7 levels -> verify 7 XLSX files in 7 sub-directories.",
        "Run with --levels g -> verify only the genus XLSX is produced.",
        "Hand-corrupt one level (rename the file) -> verify only that level's R session fails, others succeed.",
        "Inspect mbx_ezclean_info.txt and verify all reported paths exist.",
        "Run inside an active conda env -> verify _R wrapper successfully avoids 'methods.dylib not found'.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_10_ezviz():
    NAME = "10_mbx_ezviz_all_levels_all_treatments.docx"
    doc = make_doc_with_cover(
        "Auto-detect categorical columns, run mbX::ezviz() for every level x variable",
        "Pipeline step 9 -- the first multi-variable cross-product step.",
        "mbx_ezviz_all_levels_all_treatments.sh",
        step_id="9",
    )

    section_purpose(doc, (
        "For every categorical metadata column (auto-detected), generate one "
        "stacked-bar visualization per taxonomic level (1-7) using "
        "mbX::ezviz(). The script handles QIIME2 #q2:types directives, "
        "skips ID-like columns, and sanitizes filesystem-unfriendly variable "
        "names."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Directory containing 7_taxonomy_csv/", "yes"],
            ["--variables LIST", "User flag", "Override auto-detection with explicit list", "no"],
            ["--threshold F", "User flag", "ezviz aggregation threshold (default 0.5%)", "no"],
            ["--levels LIST", "User flag", "Subset (default all 7)", "no"],
        ],
        outputs_table=[
            ["9_visualization_entire/<Variable>/ezviz_<level>_<Variable>.png", "Stacked-bar PNG (one per variable x level)", "report"],
            ["9_visualization_entire/mbx_ezviz_info.txt", "K=V manifest of variables x levels generated", "report"],
        ],
        upstream_consumers=[
            "Step 18 (final report) embeds these PNGs as base64.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--variables LIST", "(auto)",
             "Override the categorical-detection step. Useful when the user wants to plot only ONE variable."],
            ["--threshold F", "0.5",
             "Pass-through to ezviz threshold argument. Taxa with relative abundance below this percentage are aggregated into 'Others'. 0.5% is mbX's default."],
            ["--levels LIST", "all 7", "Same as ezclean."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Categorical column detection (R-based)", [
            ("p", "An R script reads the metadata file (handling tab/comma delimiters), skips the optional #q2:types directive row if present, then evaluates each column with the following rules:"),
            ("bullets", [
                "Skip the first column (sample-id).",
                "Skip the column if every value is numeric (suppressWarnings(as.numeric()) returns no NAs).",
                "Skip if every value is unique (free-text / IDs / dates).",
                "Skip if there is only one unique value (no variation).",
                "Require at least one group with >= 2 samples (single-sample groups give degenerate plots).",
            ]),
            ("p", "The R script prints surviving column names between MBX_CATS_BEGIN / MBX_CATS_END sentinels on stdout. Bash captures them and iterates."),
        ]),
        ("4.2 Variable name sanitisation", [
            ("p", "Some legitimate categorical column names contain spaces or special characters that filesystems handle poorly (`Sample Type`, `Site/Region`). A sanitiser replaces spaces with underscores, strips parens/brackets/slashes, and collapses repeated underscores."),
        ]),
        ("4.3 ezviz invocation", [
            ("p", "For every (variable, level) pair, write a temp R script that:"),
            ("code",
             "library(mbX); library(ggplot2)\n"
             "setwd(\"9_visualization_entire/<sanitized_variable>\")\n"
             "p <- ezviz(\"<level-7.csv>\", \"<metadata.txt>\", \"<level_letter>\",\n"
             "           \"<original_variable>\", threshold = 0.5)\n"
             "ggsave(\"ezviz_<level_name>_<sanitized_variable>.png\",\n"
             "       plot = p, width = 10, height = 7, dpi = 300)"),
            ("p", "ezviz returns a ggplot object rather than auto-saving. We must ggsave() with explicit dimensions to ensure publication-ready output. Idempotency: if the PNG already exists, the temp R script is not even generated."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "10_ezviz.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.0, h=0.6, label="START", kind="term"),
        dict(id="info", x=5, y=10.0, w=5.4, h=0.7, label="Read mbx_taxonomy_info.txt"),
        dict(id="rd", x=5, y=8.9, w=4.6, h=0.7, label="Locate Rscript + _R wrapper"),
        dict(id="cats", x=5, y=7.7, w=5.6, h=1.0, label="Categorical detection (R)\nskip num / unique / single / 1-sample groups", kind="decision"),
        dict(id="loop", x=5, y=6.0, w=6.0, h=1.4, label="for var in CATS:\n  for L in d,p,c,o,f,g,s:\n    if PNG missing: ggsave(ezviz(...))", kind="decision"),
        dict(id="info2", x=5, y=4.0, w=5.4, h=0.7, label="Write mbx_ezviz_info.txt", kind="io"),
        dict(id="end", x=5, y=2.8, w=2.0, h=0.5, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "info"}),
        dict(**{"from": "info", "to": "rd"}),
        dict(**{"from": "rd", "to": "cats"}),
        dict(**{"from": "cats", "to": "loop"}),
        dict(**{"from": "loop", "to": "info2"}),
        dict(**{"from": "info2", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_ezviz_all_levels_all_treatments.sh -- control flow")
    section_flow(doc, flow_png, "Categorical detection runs once; the level x variable cross-product runs ezviz once per cell.")

    section_edge_cases(doc, [
        ("All-numeric metadata column",
         "Skipped automatically (numeric variables would yield gradient bar plots that are misleading). User can force inclusion via --variables but ezviz will still treat the values as factor levels."),
        ("Variable names with slashes",
         "'Site/Region' would create a sub-sub-directory. The sanitizer replaces / with _ to produce 'Site_Region'."),
    ])

    section_impl_notes(doc, [
        ("Why R returns ggplot, not a saved file",
         "ezviz is designed for interactive use first; auto-saving in R packages is considered a side-effect. Our script makes the side-effect explicit via ggsave with controlled dimensions."),
        ("Sentinel-based stdout parsing",
         "We do not parse R output by line number. Instead R prints `MBX_CATS_BEGIN` then column names then `MBX_CATS_END`. Bash uses awk between these markers. This is robust to extra messages from ggplot or library load."),
    ])

    section_testing(doc, [
        "Run with metadata containing 1 numeric + 2 categorical columns -> verify only the 2 categorical columns get plot dirs.",
        "Run with --variables 'sample-id' (intentionally broken) -> verify the script refuses or warns.",
        "Run twice -> verify second run skips all (PNGs already exist).",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_11_ezstat():
    NAME = "11_mbx_ezstat_all_levels_all_treatments.docx"
    doc = make_doc_with_cover(
        "Per-taxon statistical tests + box-plots, every level x every variable",
        "Pipeline step 10 -- the descriptive-statistics step.",
        "mbx_ezstat_all_levels_all_treatments.sh",
        step_id="10",
    )

    section_purpose(doc, (
        "Run mbX::ezstat() once per (categorical variable x taxonomic level) "
        "combination. ezstat computes Kruskal-Wallis (omnibus), Dunn pairwise "
        "(BH-adjusted), and CLD (Compact Letter Display) summaries, and plots "
        "a box plot per significant taxon (p <= 0.05). Output is XLSX + PNG "
        "directories suitable for direct insertion into supplementary "
        "material."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Directory containing 7_taxonomy_csv/", "yes"],
            ["--variables LIST / --levels LIST", "User flags", "Same as ezviz", "no"],
        ],
        outputs_table=[
            ["10_stats/<Var>/ezstat_KW_<level>_by_<Var>.xlsx", "Kruskal-Wallis omnibus per taxon", "report"],
            ["10_stats/<Var>/ezstat_pairwise_<level>_by_<Var>.xlsx", "Dunn pairwise comparisons (BH-adjusted)", "report"],
            ["10_stats/<Var>/ezstat_CLD_Summary_<level>_by_<Var>.xlsx", "Compact letter display per taxon", "report"],
            ["10_stats/<Var>/Boxplots_<Var>/", "PNG per significant taxon (one per p<=0.05 hit)", "report"],
            ["10_stats/mbx_ezstat_info.txt", "K=V manifest", "report"],
        ],
        upstream_consumers=[
            "Step 18 (final report) embeds the KW xlsx as a stats table and a sampling of significant boxplots.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--variables / --levels", "(auto / all)", "Identical to ezviz."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 ezstat() vs ezviz()", [
            ("p", "ezstat is similar to ezviz in input shape but writes ALL files internally (no ggsave needed). Output files are deterministic-named, so the bash wrapper only needs to setwd() into the output directory and call ezstat()."),
        ]),
        ("4.2 The skip-check", [
            ("p", "The wrapper checks for `ezstat_KW_<level>_by_<var>.xlsx` (the FIRST file ezstat produces). If present, the entire (var, level) cell is skipped. This is conservative -- a partially-completed run that produced KW but not pairwise/CLD will be incorrectly skipped. Mitigation: delete the entire variable directory to force regeneration."),
        ]),
        ("4.3 Statistics performed (inside ezstat)", [
            ("bullets", [
                "Kruskal-Wallis: non-parametric omnibus test for >= 3 groups; collapses to Mann-Whitney U for 2 groups.",
                "Dunn pairwise: post-hoc test for KW. p-values BH-adjusted.",
                "CLD (Compact Letter Display): groups receive shared letters when no pairwise contrast separates them at p<=0.05.",
                "Box plot per significant taxon: ggplot with per-group quartiles + outliers + sample-count overlay.",
            ]),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "11_ezstat.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.0, h=0.6, label="START", kind="term"),
        dict(id="info", x=5, y=10.0, w=5.4, h=0.7, label="Read mbx_taxonomy_info.txt"),
        dict(id="rd", x=5, y=8.9, w=4.6, h=0.7, label="_R wrapper + Rscript"),
        dict(id="cats", x=5, y=7.7, w=5.6, h=0.9, label="Categorical detection (same as step 9)", kind="decision"),
        dict(id="loop", x=5, y=6.0, w=6.0, h=1.4, label="for var in CATS:\n  for L in 1..7:\n    if KW xlsx missing: ezstat(level7, meta, L, var)", kind="decision"),
        dict(id="emit", x=5, y=4.0, w=5.4, h=0.9, label="ezstat writes:\n  KW.xlsx + pairwise.xlsx + CLD.xlsx + Boxplots/", kind="io"),
        dict(id="info2", x=5, y=2.5, w=5.0, h=0.7, label="Write mbx_ezstat_info.txt", kind="io"),
        dict(id="end", x=5, y=1.4, w=2.0, h=0.5, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "info"}),
        dict(**{"from": "info", "to": "rd"}),
        dict(**{"from": "rd", "to": "cats"}),
        dict(**{"from": "cats", "to": "loop"}),
        dict(**{"from": "loop", "to": "emit"}),
        dict(**{"from": "emit", "to": "info2"}),
        dict(**{"from": "info2", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_ezstat_all_levels_all_treatments.sh -- control flow")
    section_flow(doc, flow_png, "Per-cell control flow mirrors ezviz; ezstat handles all I/O internally.")

    section_edge_cases(doc, [
        ("Single-group variable",
         "Already excluded by categorical detection. KW requires >= 2 groups."),
        ("Group with 1 sample",
         "Categorical detection requires >= 2 samples per group. ezstat would otherwise produce NaN p-values."),
        ("Partial-completion skip false-positive",
         "If KW exists but pairwise was never written, the cell is skipped. Workaround: rm the entire variable directory and re-run."),
    ])

    section_impl_notes(doc, [
        ("Per-cell R session isolation",
         "Same pattern as ezviz: each cell runs in its own _R invocation. A failure in one cell does not break the others."),
        ("ezstat handles file output internally",
         "Unlike ezviz, ezstat is opinionated about its outputs. We do not need to capture and save anything from R's stdout -- just check existence after."),
    ])

    section_testing(doc, [
        "Run with 3 categorical variables x 7 levels -> verify 21 cell-directories exist.",
        "Confirm KW + pairwise + CLD + Boxplots/ all present in one selected cell.",
        "Confirm Boxplots/ directory contains one PNG per significant taxon (p<=0.05).",
        "Re-run -> verify all cells skipped.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_12_pre_diversity():
    NAME = "12_mbx_pre_diversity_parameters.docx"
    doc = make_doc_with_cover(
        "Build the rooted tree, choose rarefaction depth, run analytical rarefaction",
        "Pipeline step 11 -- the most algorithmically dense script. The diversity gatekeeper.",
        "mbx_pre_diversity_parameters.sh",
        step_id="11",
    )

    section_purpose(doc, (
        "Pre-diversity is the gateway to all alpha/beta/diversity steps. It "
        "(1) builds the rooted phylogenetic tree, (2) selects a "
        "scientifically-defensible rarefaction depth using a multi-criterion "
        "decision rule, (3) runs an analytical (Hurlbert + Chao-Shen) "
        "rarefaction in seconds (instead of QIIME2's hours-long Monte Carlo "
        "alternative), and (4) writes a STATUS field that gates every "
        "subsequent step."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing earlier step outputs", "yes"],
            ["--group-col COL", "User flag", "Force per-group retention check on a specific column", "no"],
            ["--min-overall F", "User flag", "Minimum overall sample retention (default 0.90)", "no"],
            ["--min-group F", "User flag", "Minimum per-group retention (default 0.80)", "no"],
            ["--rare-steps N / --rare-iter N", "User flags", "QIIME alpha-rarefaction sampling resolution", "no"],
            ["--skip-tree / --skip-qc", "User flags", "Boolean toggles", "no"],
        ],
        outputs_table=[
            ["11_pre_diversity/rooted-tree.qza", "Phylogeny[Rooted]", "step 12-15 (diversity, ANCOMBC)"],
            ["11_pre_diversity/sample-frequencies.tsv", "Per-sample read counts (TSV-exported)", "depth selection + report"],
            ["11_pre_diversity/sampling_depth_candidates.csv", "Per-depth metrics that drive selection", "report"],
            ["11_pre_diversity/alpha_rarefaction_data.csv + curves.png", "Analytical rarefaction (Hurlbert+ChaoShen)", "report"],
            ["11_pre_diversity/alpha-rarefaction-qiime.qzv", "Official QIIME alpha-rarefaction visualization", "report"],
            ["11_pre_diversity/mbx_pre_diversity_info.txt", "K=V machine-readable + STATUS + READY_FOR_DIVERSITY", "step 12, 13, 14, 15 (gate)"],
            ["11_pre_diversity/mbx_pre_diversity_summary.txt", "Plain-language explanation for the lay user", "user, report"],
        ],
        upstream_consumers=[
            "Steps 12-15 read STATUS / READY_FOR_DIVERSITY and refuse to run on FAIL unless --force.",
            "Steps 12 (alpha) and 13 (beta) read RECOMMENDED_DEPTH for core-metrics rarefaction.",
            "Step 18 embeds the rarefaction curves and depth-decision summary verbatim.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--min-overall F", "0.90",
             "Minimum fraction of samples that must be retained at the chosen depth. Below this, the candidate is rejected. 0.90 is the conservative microbiome-publication norm."],
            ["--min-group F", "0.80",
             "Per-group minimum retention when --group-col is used (or auto-selected). Defaults to 0.80 to allow some uneven sampling."],
            ["--group-col COL", "(auto)",
             "Force the per-group check to use this metadata column. Default: pick the categorical column with the MOST groups (the most stringent constraint)."],
            ["--rare-steps N", "20", "QIIME alpha-rarefaction --p-steps."],
            ["--rare-iter N", "10", "QIIME alpha-rarefaction --p-iterations."],
            ["--skip-tree", "off", "Reuse an existing rooted-tree.qza (e.g. from an aborted run)."],
            ["--skip-qc", "off", "Skip metadata + feature-table QZV regeneration."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Auto-discovery (graceful chain of fallbacks)", [
            ("p", "The script walks the previous step info files in order:"),
            ("bullets", [
                "metadata path: from 7_taxonomy_csv/mbx_taxonomy_info.txt (canonical).",
                "feature_table_filtered.qza: from 7_taxonomy_csv/ -- REQUIRED. Refuses to fall back to the unfiltered DADA2 table because mito/chloro contamination materially distorts diversity metrics.",
                "rep_seqs.qza: from 4_dada2_outputs/.",
                "taxonomy.qza: from 6_classifier_taxonomy/ (used only for sanity check + report).",
            ]),
        ]),
        ("4.2 QC visualizations", [
            ("p", "qiime metadata tabulate -> metadata_summary.qzv. qiime feature-table summarize -> feature_table_summary.qzv. Both are user-facing 'have you got what you think you have' sanity checks."),
        ]),
        ("4.3 Phylogenetic tree construction", [
            ("p", "qiime phylogeny align-to-tree-mafft-fasttree wraps four operations in a single call: mafft alignment, alignment masking, FastTree2 unrooted phylogeny, midpoint-rooting. --p-n-threads is auto-detected. Output: aligned-rep-seqs.qza, masked-aligned-rep-seqs.qza, unrooted-tree.qza, rooted-tree.qza."),
        ]),
        ("4.4 Sample-frequencies tabulation + the comma-thousands bug", [
            ("p", "qiime feature-table tabulate-sample-frequencies -> sample-frequencies.qza. Exporting that QZA gives a TSV with TWO non-obvious format quirks:"),
            ("bullets", [
                "First two rows often start with `#q2:types` and `#q2:tabulate-` directives -- they MUST be skipped before parsing.",
                "The frequency column contains commas as thousands separators (e.g. '28,005.0') -- these MUST be stripped before as.numeric(), or every value parses as NA.",
            ]),
            ("p", "Both quirks are handled in R: skip lines starting with #, then `gsub(',', '', val)` before `as.numeric(val)`."),
        ]),
        ("4.5 The depth-selection algorithm (the science)", [
            ("p", "Three concurrent criteria evaluated at every candidate depth d:"),
            ("kv", [
                ("(a) Overall retention",     "fraction of samples with N_i >= d must be >= MIN_OVERALL (default 0.90)."),
                ("(b) Good's coverage",        "mean Good's coverage of retained samples at depth d must be >= GOOD_COV_MIN (default 0.98). Good's coverage = 1 - f1/N where f1 is singletons. This catches under-sampling that retention alone misses."),
                ("(c) Plateau slope",          "slope of the mean observed_features curve at d must be < PLATEAU_SLOPE_MAX (0.5 features per 1000 reads). Communities below the asymptote contribute non-comparable diversity estimates."),
            ]),
            ("p", "If --group-col is provided OR auto-selected, an ADDITIONAL constraint is enforced: every group must retain >= MIN_GROUP fraction of samples at the candidate depth."),
            ("p", "The script picks the HIGHEST depth d that satisfies ALL active criteria. This maximises the depth (-> more diversity captured) without sacrificing the retention/coverage/plateau guarantees."),
        ]),
        ("4.6 Group-column auto-selection rule", [
            ("p", "Of all categorical variables, pick the one with the MOST groups (most stringent constraint). Rationale: if the user is going to compare across treatment x site x time, the most-divided variable will have the smallest per-group cells; satisfying retention there guarantees retention for coarser groupings."),
        ]),
        ("4.7 Analytical rarefaction (key innovation)", [
            ("p", "QIIME2's `qiime diversity alpha-rarefaction` runs Monte Carlo: at each depth, sample reads N times, compute observed_features each time, average. With 10 iterations x 20 steps x N samples, on the test dataset this took ~30 minutes. Hours on real datasets."),
            ("p", "Our analytical alternative uses the closed-form Hurlbert (1971) formula:"),
            ("code",
             "E[S(n)] = sum_i [ 1 - C(N - n_i, n) / C(N, n) ]\n"
             "       = sum_i [ 1 - exp( log_binom(N - n_i, n) - log_binom(N, n) ) ]"),
            ("p", "Where N = sample's total read count, n = subsampling depth, n_i = count of feature i. We use scipy.special.gammaln for numerically-stable log-space computation. This is EXACT (not estimated) and runs in seconds."),
            ("p", "Shannon entropy uses Chao-Shen (2003) coverage-adjusted analytical estimator. Faith PD is intentionally NOT included in the curve because the tree topology is fixed -- including it adds 70% runtime for no insight (it gets computed at the chosen depth in step 12 anyway)."),
            ("p", "Parallelization: ProcessPoolExecutor distributes per-sample computation across CPU cores. Result: same science, smoother curves (exact vs Monte Carlo), seconds instead of hours."),
        ]),
        ("4.8 ID validation", [
            ("p", "Three ID intersections are verified:"),
            ("bullets", [
                "metadata sample-IDs intersect with feature-table sample-IDs (>= 95% overlap or FAIL).",
                "feature-table feature-IDs intersect with rep-seqs IDs.",
                "rep-seqs IDs intersect with rooted-tree leaf names.",
            ]),
            ("p", "Failures here set OVERALL_STATUS=FAIL. Subsequent steps refuse to run unless --force."),
        ]),
        ("4.9 Decision-supporting visualizations", [
            ("p", "depth_vs_retention.png shows fraction-retained vs depth with the selected depth annotated. sequencing_depth_distribution.png shows per-sample histogram with median + selected-depth lines. Both are referenced by the lay-user summary."),
        ]),
        ("4.10 STATUS field semantics", [
            ("kv", [
                ("PASS",                "All three criteria met cleanly with the auto-selected depth. Diversity steps run without warnings."),
                ("PASS_WITH_WARNINGS",  "Criteria met but borderline (e.g. retention at 90.5% rather than the 95% comfort zone). Steps run but show a yellow banner."),
                ("REVIEW_REQUIRED",     "A fallback rule was used (e.g. 75% retention) or one criterion failed but the others compensated. User should read the summary file before continuing."),
                ("FAIL",                "Pipeline-stopping problem (no metadata-table overlap, group fully wiped, etc.). Subsequent steps MUST refuse without --force."),
            ]),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "12_pre_diversity.png"
    nodes = [
        dict(id="start", x=5, y=11.6, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="disc", x=5, y=10.8, w=5.4, h=0.6, label="Auto-discover inputs (chain of info files)"),
        dict(id="qc", x=5, y=10.0, w=5.0, h=0.5, label="QC visualizations (metadata + ft summary)"),
        dict(id="tree", x=5, y=9.0, w=5.4, h=0.7, label="Phylogenetic tree\nmafft -> fasttree -> midpoint-root", kind="io"),
        dict(id="freq", x=5, y=7.9, w=5.4, h=0.7, label="Tabulate sample-frequencies\nexport TSV (skip #, strip commas)"),
        dict(id="depth", x=5, y=6.5, w=6.4, h=1.2, label="Depth selection (R)\n(a) overall retention >= MIN_OVERALL\n(b) Good's coverage >= GOOD_COV_MIN\n(c) plateau slope < SLOPE_MAX", kind="decision"),
        dict(id="rare", x=5, y=4.8, w=6.4, h=1.0, label="Analytical rarefaction (Python)\nHurlbert 1971 + Chao-Shen 2003\nProcessPoolExecutor across samples"),
        dict(id="qrare", x=5, y=3.5, w=5.4, h=0.6, label="Official QIIME alpha-rarefaction (audit)"),
        dict(id="ids", x=5, y=2.6, w=5.4, h=0.7, label="ID intersections (metadata/table/seqs/tree)", kind="decision"),
        dict(id="status", x=5, y=1.4, w=5.4, h=0.7, label="Compute STATUS\nPASS / WARN / REVIEW / FAIL", kind="decision"),
        dict(id="info", x=5, y=0.4, w=5.4, h=0.5, label="Write info + summary", kind="io"),
    ]
    edges = [
        dict(**{"from": "start", "to": "disc"}),
        dict(**{"from": "disc", "to": "qc"}),
        dict(**{"from": "qc", "to": "tree"}),
        dict(**{"from": "tree", "to": "freq"}),
        dict(**{"from": "freq", "to": "depth"}),
        dict(**{"from": "depth", "to": "rare"}),
        dict(**{"from": "rare", "to": "qrare"}),
        dict(**{"from": "qrare", "to": "ids"}),
        dict(**{"from": "ids", "to": "status"}),
        dict(**{"from": "status", "to": "info"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_pre_diversity_parameters.sh -- 12-stage flow", figsize=(8.5, 11.0))
    section_flow(doc, flow_png, "12 stages, but every artifact is reusable on rerun.")

    section_edge_cases(doc, [
        ("Comma thousands separator in TSV",
         "QIIME2 tabulate-sample-frequencies emits '28,005.0'. R's as.numeric() returns NA without prior gsub(',',''). This was a real production bug; fix is in stage 4."),
        ("Group fully wiped at chosen depth",
         "If --group-col enforcement fails, STATUS drops to REVIEW_REQUIRED with the binding criterion named. The user can lower --min-group, change --group-col, or drop the constraint."),
        ("FAIL means downstream steps refuse",
         "STATUS=FAIL writes READY_FOR_DIVERSITY=no. Steps 12-15 read this field and abort with an actionable message."),
        ("Hurlbert formula cancellation for very small N-n_i",
         "When n_i ~ N (singleton features in low-N samples), C(N-n_i, n) can underflow. The log-space computation via gammaln avoids this entirely; raw factorials would have failed."),
    ])

    section_impl_notes(doc, [
        ("Why analytical rarefaction beats Monte Carlo",
         "Hurlbert 1971 gives the EXPECTED number of features at depth n -- the same quantity Monte Carlo estimates. Computing it directly is exact, deterministic, and parallelisable across samples (each sample's curve is independent)."),
        ("Why Faith PD is excluded from the curve",
         "Faith PD = sum of branch lengths spanned by present features. The branch-length sum changes very little once we have ~80% of the features (most of the tree is sampled by then). Including it in the rarefaction curve adds 70% runtime for ~3% additional information. We compute it ONCE at the chosen depth in step 12."),
        ("STATUS as the diversity gate",
         "Every diversity-related step (12, 13, 14, 15) reads STATUS at start and aborts on FAIL. --force bypass is intentionally noisy."),
        ("Resilience to partial reruns",
         "Tree, sample-frequencies, official QIIME alpha-rarefaction are all idempotent file checks. A user re-running after a network blip skips the slow stages."),
    ])

    section_testing(doc, [
        "Run on test data with PASS expected -> verify STATUS=PASS, READY_FOR_DIVERSITY=yes.",
        "Hand-corrupt the metadata to break sample-ID intersection -> verify STATUS=FAIL.",
        "Set --min-overall 0.99 (artificially strict) -> verify the script falls to REVIEW_REQUIRED with retention named as the binding criterion.",
        "Run with --skip-tree after deleting rooted-tree.qza -> verify the script aborts with a clear message.",
        "Verify analytical rarefaction matches QIIME's Monte Carlo within Monte Carlo error (compare alpha_rarefaction_data.csv to the QIIME QZV).",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def main():
    builders = [
        build_04_dada2_param_finder,
        build_05_dada2_run,
        build_06_classifier_arranger,
        build_07_classifier_run,
        build_08_taxonomy_run,
        build_09_ezclean,
        build_10_ezviz,
        build_11_ezstat,
        build_12_pre_diversity,
    ]
    for b in builders:
        b()


if __name__ == "__main__":
    main()
