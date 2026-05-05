"""
Chapters 13-20: alpha, beta, ANCOMBC, PICRUSt, ML, networks, final_report, orchestrator.
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


def build_13_alpha():
    NAME = "13_mbx_alpha_diversity_run.docx"
    doc = make_doc_with_cover(
        "Compute 5 alpha-diversity metrics + per-variable Kruskal-Wallis stats + boxplots",
        "Pipeline step 12 -- the first proper diversity step.",
        "mbx_alpha_diversity_run.sh",
        step_id="12",
    )

    section_purpose(doc, (
        "Rarefy the filtered feature table at the depth chosen by step 11, "
        "compute five canonical alpha-diversity metrics (ASVs/Features, "
        "Shannon, Simpson, Faith PD, Pielou Evenness), join them with the "
        "metadata into a single tidy alpha_diversity.xlsx, and run "
        "Kruskal-Wallis + Dunn pairwise + CLD per categorical variable "
        "(plus boxplots) for every metric."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing 11_pre_diversity/", "yes"],
            ["--metrics LIST", "User flag", "Subset of metrics (default all 5)", "no"],
            ["--depth N", "User flag", "Override RECOMMENDED_DEPTH from step 11", "no"],
            ["--force", "User flag", "Bypass STATUS=FAIL gate", "no"],
        ],
        outputs_table=[
            ["12_alpha_diversity_results/all_alpha_outputs/<metric>_vector.qza", "Raw QIIME2 alpha-vector artifacts (5 of them)", "audit/report"],
            ["12_alpha_diversity_results/all_alpha_outputs/exported/<metric>/alpha-diversity.tsv", "TSV exports for joining", "report"],
            ["12_alpha_diversity_results/alpha_diversity.xlsx", "Tidy table: sample-id + 5 metrics + every metadata column", "report"],
            ["12_alpha_diversity_results/stats_for_alpha_diversity/<Var>/{KW,Pairwise,CLD}_<metric>_by_<Var>.xlsx", "Per-variable per-metric stats", "report"],
            ["12_alpha_diversity_results/boxplots_for_alpha_diversity/<Var>/boxplot_<metric>_by_<Var>.png + panel", "Per-variable boxplots, one per metric + a 5-panel summary", "report"],
            ["12_alpha_diversity_results/mbx_alpha_diversity_info.txt", "K=V manifest", "step 13 (beta), 18 (report)"],
        ],
        upstream_consumers=[
            "Step 13 (beta) chains from this info file for metadata path + STATUS check.",
            "Step 18 embeds alpha_diversity.xlsx + a sample of the boxplots.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--metrics LIST", "ASVs,Shannon,Simpson,Faith,Pielou", "Comma-separated metric names. Subsetting saves time on Faith PD which requires the tree."],
            ["--depth N", "(read from step 11)", "Override RECOMMENDED_DEPTH. Required for sensitivity analyses."],
            ["--force", "off", "Bypass STATUS=FAIL gate from step 11. Use only when you understand WHY the gate fired."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Step-11 gate", [
            ("p", "Read OVERALL_STATUS and READY_FOR_DIVERSITY from `11_pre_diversity/mbx_pre_diversity_info.txt`. If STATUS=FAIL or READY_FOR_DIVERSITY=no, abort with the relevant message UNLESS --force is provided. PASS_WITH_WARNINGS is allowed (the user has been warned)."),
        ]),
        ("4.2 Rarefaction", [
            ("p", "qiime feature-table rarefy on the filtered table at RECOMMENDED_DEPTH. Output: `rarefied_table.qza`. Idempotent (skip if exists)."),
        ]),
        ("4.3 Per-metric computation", [
            ("p", "Five QIIME2 calls run in sequence:"),
            ("kv", [
                ("ASVs / Features",        "qiime diversity alpha --p-metric observed_features"),
                ("Shannon",                "qiime diversity alpha --p-metric shannon"),
                ("Simpson (1 - Simpson)",  "qiime diversity alpha --p-metric simpson"),
                ("Pielou's Evenness",      "qiime diversity alpha --p-metric pielou_e"),
                ("Faith Phylogenetic D",   "qiime diversity alpha-phylogenetic --p-metric faith_pd (requires tree)"),
            ]),
            ("p", "Each output is a SampleData[AlphaDiversity] artifact. We then `qiime tools export` each to TSV (sample-id\\tvalue). Idempotent: each metric's QZA AND its TSV existence is checked separately."),
        ]),
        ("4.4 Tidy joining (R)", [
            ("p", "An R helper reads the metadata file (handles delimiter detection + #q2:types skip), reads each of the 5 TSVs, joins them by sample-id, and emits alpha_diversity.xlsx with a uniform column ordering: sample-id, ASVs_or_Features, Shannon_Index, Simpson_Diversity, Faith_Phylogenetic_Diversity, Pielou_Evenness, then every metadata column."),
        ]),
        ("4.5 Per-variable statistics", [
            ("p", "For every categorical variable (auto-detected with the same R logic as step 9/10), we run, for each metric:"),
            ("bullets", [
                "Kruskal-Wallis (omnibus across groups). Holm correction across metrics so the user does not have to worry about the multiple-metrics inflation.",
                "Dunn pairwise (BH-adjusted within each variable).",
                "CLD computation: groups receive shared letters when no contrast separates them at q<=0.05.",
            ]),
            ("p", "Output: KW_<metric>_by_<Var>.xlsx, Pairwise_<metric>_by_<Var>.xlsx, CLD_Summary_<metric>_by_<Var>.xlsx, plus a Summary_KW_all_metrics_by_<Var>.xlsx that stacks all 5 KW results into one sheet."),
        ]),
        ("4.6 Boxplots", [
            ("p", "ggplot per (variable, metric) -> PNG. Plus a 5-panel patchwork combining all metrics for a single variable, useful for slides."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "13_alpha.png"
    nodes = [
        dict(id="start", x=5, y=11.0, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="gate", x=5, y=10.0, w=5.4, h=0.8, label="STATUS gate (step 11)\nFAIL -> abort unless --force", kind="decision"),
        dict(id="rarefy", x=5, y=8.7, w=5.0, h=0.7, label="qiime feature-table rarefy", kind="io"),
        dict(id="metrics", x=5, y=7.5, w=5.6, h=1.0, label="5 metrics (qiime diversity alpha[-phylogenetic])\nObserved / Shannon / Simpson / Pielou / Faith PD", kind="io"),
        dict(id="export", x=5, y=6.1, w=5.0, h=0.7, label="qiime tools export per metric -> TSV"),
        dict(id="join", x=5, y=5.0, w=5.6, h=0.8, label="R: join metadata + 5 TSVs -> alpha_diversity.xlsx"),
        dict(id="cats", x=5, y=3.7, w=5.6, h=0.8, label="Categorical detection\n(skip first col / numeric / unique / single)", kind="decision"),
        dict(id="stats", x=5, y=2.5, w=6.0, h=0.9, label="Per-(var, metric):\n  KW + Pairwise + CLD + boxplot"),
        dict(id="info", x=5, y=1.2, w=5.0, h=0.6, label="Write mbx_alpha_diversity_info.txt", kind="io"),
        dict(id="end", x=5, y=0.4, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "gate"}),
        dict(**{"from": "gate", "to": "rarefy"}),
        dict(**{"from": "rarefy", "to": "metrics"}),
        dict(**{"from": "metrics", "to": "export"}),
        dict(**{"from": "export", "to": "join"}),
        dict(**{"from": "join", "to": "cats"}),
        dict(**{"from": "cats", "to": "stats"}),
        dict(**{"from": "stats", "to": "info"}),
        dict(**{"from": "info", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_alpha_diversity_run.sh -- control flow", figsize=(8.0, 11.0))
    section_flow(doc, flow_png, "5 metrics x N variables x stats + plots, all driven by step 11's depth.")

    section_edge_cases(doc, [
        ("Faith PD requires the tree",
         "If rooted-tree.qza is missing or has fewer leaves than the rarefied table, Faith PD fails. The script reports the metric as failed but continues with the other 4."),
        ("Missing samples after rarefaction",
         "qiime feature-table rarefy drops samples below depth. The TSV exports therefore have FEWER samples than metadata. The R join uses left-join on metadata so dropped samples appear with NA metric values -- visually obvious in the XLSX."),
        ("Single-group variable",
         "Already excluded by categorical detection."),
    ])

    section_impl_notes(doc, [
        ("Why Holm correction across metrics",
         "5 metrics x N variables = many KW tests. Holm is conservative but uniformly more powerful than Bonferroni and safe to use without distributional assumptions. Holm is applied PER variable across the 5 metrics, not across all variables."),
        ("Why a 5-panel patchwork",
         "Reviewers commonly ask 'show me everything for this variable on one figure'. The patchwork combines all 5 boxplots with a shared variable label."),
    ])

    section_testing(doc, [
        "Run with PASS status -> verify all 5 metrics + alpha_diversity.xlsx produced.",
        "Run with --metrics ASVs,Shannon -> verify only those 2 metrics computed.",
        "Run with --depth 2000 -> verify rarefaction at 2000 reads (ignoring step 11's recommendation).",
        "Hand-edit STATUS=FAIL in step 11 info file -> verify abort.",
        "Re-run -> verify rarefied_table.qza skip and per-metric skip both fire.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_14_beta():
    NAME = "14_mbx_beta_diversity_run.docx"
    doc = make_doc_with_cover(
        "Beta diversity, PCoA, PERMANOVA, PERMDISP, Adonis, distance heatmaps, UPGMA",
        "Pipeline step 13 -- the most reviewer-targeted output stage.",
        "mbx_beta_diversity_run.sh",
        step_id="13",
    )

    section_purpose(doc, (
        "Compute four canonical beta-diversity distance matrices "
        "(Bray-Curtis, Jaccard, Weighted UniFrac, Unweighted UniFrac), then "
        "for every categorical variable produce: PCoA plots (per metric + "
        "4-panel summary), PERMANOVA + PERMDISP + pairwise PERMANOVA tables, "
        "distance-to-centroid boxplots, plus a global Adonis table, distance "
        "heatmaps, and UPGMA dendrograms across all samples."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing 12_alpha_diversity_results/", "yes"],
            ["--metrics LIST", "User flag", "Subset of 4 distance metrics", "no"],
            ["--permutations N", "User flag", "PERMANOVA/PERMDISP permutation count (default 999)", "no"],
        ],
        outputs_table=[
            ["13_beta_diversity_results/results_by_categorical_variables/<Var>/PCoA_*.png", "Per-variable PCoA per metric + 4-panel", "report"],
            ["13_beta_diversity_results/results_by_categorical_variables/<Var>/PERMANOVA_results_<Var>.xlsx", "1 sheet, 4 rows (one per metric)", "report"],
            ["13_beta_diversity_results/results_by_categorical_variables/<Var>/Pairwise_PERMANOVA_<Var>.xlsx", "Only when k>2 groups", "report"],
            ["13_beta_diversity_results/results_by_categorical_variables/<Var>/PERMDISP_results_<Var>.xlsx", "Per-metric dispersion test", "report"],
            ["13_beta_diversity_results/results_by_categorical_variables/<Var>/Boxplot_DistanceToCentroid_*.png", "Visualizes PERMDISP", "report"],
            ["13_beta_diversity_results/all_samples_beta_diversity/Adonis_multivariable_PERMANOVA.xlsx", "Type III + univariate Adonis with all metadata columns", "report"],
            ["13_beta_diversity_results/all_samples_beta_diversity/Distance_heatmap_<Metric>.png", "Clustered heatmap, all samples", "report"],
            ["13_beta_diversity_results/all_samples_beta_diversity/UPGMA_dendrogram_<Metric>.png", "All-sample dendrogram", "report"],
            ["13_beta_diversity_results/working_dir_beta_diversity/", "Distance matrices + PCoA QZAs + exports", "audit"],
            ["13_beta_diversity_results/mbx_beta_diversity_info.txt", "K=V manifest", "step 14 (ANCOMBC), 18 (report)"],
        ],
        upstream_consumers=[
            "Step 14 chains via this file. Step 18 embeds many of the figures and tables.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--metrics LIST", "BC,Jaccard,WUF,UUF",
             "Comma-separated metric tokens. WUF = Weighted UniFrac, UUF = Unweighted UniFrac. Subsetting halves runtime when the user only wants Bray-Curtis."],
            ["--permutations N", "999",
             "PERMANOVA / PERMDISP permutation count. 999 is QIIME2 default; 9999 for publication-grade p-values."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 The dual-R conda-env interaction (CRITICAL)", [
            ("p", "This script needs TWO different R installations to coexist:"),
            ("kv", [
                ("System R", "/opt/homebrew/bin/Rscript or /usr/local/bin/Rscript -- has openxlsx, ggplot2, vegan, pheatmap installed by step 8. Compiled extensions are .so on Apple Silicon."),
                ("Conda R",  "$CONDA_PREFIX/lib/R -- used by QIIME2's q2-composition via rpy2. Compiled extensions are .dylib."),
            ]),
            ("p", "Two wrappers solve the cross-env pollution:"),
            ("code",
             "_R()     { env -u R_LIBS_USER -u R_HOME ... \"$RSCRIPT\" \"$@\"; }\n"
             "_QIIME() { R_HOME=\"$CONDA_R_HOME\" qiime \"$@\"; }"),
            ("p", "_R is used for our R scripts (vegan, plotting). _QIIME ensures qiime's embedded rpy2 finds the matching .dylib R extensions. Without this dual-wrapper, q2-composition fails with 'methods.dylib not found'."),
        ]),
        ("4.2 Distance matrices", [
            ("p", "qiime diversity beta on the rarefied feature table -> Bray-Curtis and Jaccard. qiime diversity beta-phylogenetic with the rooted tree -> Weighted/Unweighted UniFrac. Each output is a DistanceMatrix QZA."),
        ]),
        ("4.3 PCoA + Emperor", [
            ("p", "qiime diversity pcoa per matrix; qiime emperor plot per matrix. The Emperor QZV is kept as supplementary; we ALSO export the PCoA coordinates for our own ggplot rendering (so the report has a uniform ggplot-styled figure rather than a raw Emperor screenshot)."),
        ]),
        ("4.4 Per-variable hypothesis tests", [
            ("p", "For each categorical variable:"),
            ("bullets", [
                "qiime diversity beta-group-significance (PERMANOVA, default 999 permutations).",
                "If groups > 2: qiime diversity beta-group-significance --p-pairwise -> Pairwise_PERMANOVA xlsx (BH adjusted).",
                "qiime diversity beta-group-significance --p-method permdisp.",
                "ggplot distance-to-centroid boxplot from PERMDISP-exported distances.",
                "ggplot PCoA: 4 plots (one per metric) + a 4-panel patchwork.",
            ]),
            ("p", "All four distance metrics are tested, and the script records which metric gave the strongest signal (lowest q-value) per variable -- helpful for the report's 'biological summary'."),
        ]),
        ("4.5 Global Adonis (multivariable PERMANOVA)", [
            ("p", "vegan::adonis2 in R, called with formula `distance ~ var1 + var2 + ...` for ALL detected categorical variables. Type III sequential SS gives variable-controlling-for-others. Univariate Adonis is also run per variable. Results combined into one XLSX with 4 sheets (one per metric)."),
            ("p", "This single output answers the most common reviewer question: 'when you control for covariate X, does Y still drive the microbiome?'"),
        ]),
        ("4.6 Distance heatmaps + UPGMA", [
            ("p", "pheatmap on the symmetric distance matrix with hierarchical clustering on both axes. Samples are auto-coloured by the most-divisive categorical variable (annotation_col). UPGMA dendrogram via hclust(method='average') -- complementary view of the same data."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "14_beta.png"
    nodes = [
        dict(id="start", x=5, y=11.4, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="gate", x=5, y=10.6, w=5.0, h=0.6, label="STATUS gate (step 11)", kind="decision"),
        dict(id="dual", x=5, y=9.7, w=5.6, h=0.7, label="Setup _R + _QIIME wrappers"),
        dict(id="dist", x=5, y=8.6, w=5.4, h=0.9, label="4 distance matrices\n(BC, Jaccard, WUF, UUF)", kind="io"),
        dict(id="pcoa", x=5, y=7.3, w=5.0, h=0.7, label="PCoA + Emperor per metric", kind="io"),
        dict(id="cats", x=5, y=6.3, w=5.4, h=0.7, label="Categorical detection", kind="decision"),
        dict(id="perm", x=5, y=5.0, w=6.0, h=1.0, label="Per variable:\n PERMANOVA + Pairwise + PERMDISP\n + DistanceToCentroid + PCoA panels"),
        dict(id="adonis", x=5, y=3.4, w=5.6, h=0.8, label="Adonis multivariable\nType III + univariate"),
        dict(id="heat", x=5, y=2.2, w=5.6, h=0.7, label="Distance heatmaps + UPGMA"),
        dict(id="info", x=5, y=1.2, w=5.0, h=0.6, label="Write info file", kind="io"),
        dict(id="end", x=5, y=0.4, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "gate"}),
        dict(**{"from": "gate", "to": "dual"}),
        dict(**{"from": "dual", "to": "dist"}),
        dict(**{"from": "dist", "to": "pcoa"}),
        dict(**{"from": "pcoa", "to": "cats"}),
        dict(**{"from": "cats", "to": "perm"}),
        dict(**{"from": "perm", "to": "adonis"}),
        dict(**{"from": "adonis", "to": "heat"}),
        dict(**{"from": "heat", "to": "info"}),
        dict(**{"from": "info", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_beta_diversity_run.sh -- control flow", figsize=(8.0, 12.0))
    section_flow(doc, flow_png, "Distance matrices feed every test; per-variable + global tests run in parallel.")

    section_edge_cases(doc, [
        ("Two-group variable",
         "Pairwise PERMANOVA is skipped (no pairs to compare beyond the omnibus). The PERMANOVA result IS the only contrast."),
        ("PERMDISP fails when one group has < 3 samples",
         "vegan::betadisper requires >= 3 samples per group. Captured and reported as 'PERMDISP not computed for this variable -- one or more groups too small'."),
        ("Adonis with collinear variables",
         "Adonis handles this internally (degenerate terms get NA). The XLSX flags those rows."),
    ])

    section_impl_notes(doc, [
        ("Per-metric output mirroring",
         "Every per-variable directory has parallel files for all 4 metrics. This pattern simplifies the report builder -- it can iterate metrics blindly."),
        ("Why not just rely on Emperor",
         "Emperor's plots are interactive but cannot be embedded in a static report. We re-render via ggplot to get colored, captioned, dimension-labelled (% variance) PCoA plots."),
    ])

    section_testing(doc, [
        "Verify 4 distance matrices produced.",
        "Verify per-variable PERMANOVA xlsx has 4 rows when all 4 metrics requested.",
        "Verify Pairwise_PERMANOVA absent when variable has only 2 groups, present when >= 3.",
        "Verify Adonis xlsx has 4 sheets named after the metrics.",
        "Verify all_samples_beta_diversity/ contains exactly 4 heatmap PNGs and 4 dendrogram PNGs.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_15_ancombc2():
    NAME = "15_mbx_ancombc2_run.docx"
    doc = make_doc_with_cover(
        "ANCOMBC2 differential abundance at multiple taxonomic levels",
        "Pipeline step 14 -- the rigorous compositional-aware DA step.",
        "mbx_ancombc2_run.sh",
        step_id="14",
    )

    section_purpose(doc, (
        "Differential abundance via ANCOMBC2 (Lin & Peddada 2024). The script "
        "tests every categorical variable at every requested taxonomic level "
        "(phylum -> species), running with pairwise=TRUE AND global=TRUE so "
        "users get reference-free omnibus results PLUS all pairwise contrasts "
        "with mdFDR correction. Crucially, it does NOT modify the user's "
        "QIIME2 conda env."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing 13_beta_diversity_results/", "yes"],
            ["--levels LIST", "User flag", "Default: phylum,class,order,family,genus,species (skip domain)", "no"],
            ["--variables LIST", "User flag", "Override auto-detection", "no"],
            ["--alpha F", "User flag", "Significance level (default 0.05)", "no"],
            ["--prv-cut F", "User flag", "Prevalence cutoff (default 0.10)", "no"],
            ["--lib-cut N", "User flag", "Library-size cutoff (default 1000)", "no"],
        ],
        outputs_table=[
            ["14_differential_abundance_ANCOMBC2/ANCOMBC2_<level>/ancombc2_<Var>/ancombc2_primary_<Var>.xlsx", "Per-taxon LFC + se + W + p + q for each category contrast", "report"],
            ["14_differential_abundance_ANCOMBC2/ANCOMBC2_<level>/ancombc2_<Var>/ancombc2_pairwise_<Var>.xlsx", "Pairwise contrasts with mdFDR", "report"],
            ["14_differential_abundance_ANCOMBC2/ANCOMBC2_<level>/ancombc2_<Var>/ancombc2_global_<Var>.xlsx", "Reference-free global test", "report"],
            ["14_differential_abundance_ANCOMBC2/ANCOMBC2_<level>/ancombc2_<Var>/ancombc2_structural_zeros_<Var>.xlsx", "Taxa structurally zero in some groups", "report"],
            ["14_differential_abundance_ANCOMBC2/ANCOMBC2_<level>/ancombc2_<Var>/{volcano,heatmap}*.png", "Volcano + significant-taxa heatmap", "report"],
            ["14_differential_abundance_ANCOMBC2/Summary_all_levels_all_variables.xlsx", "Cross-level summary (count of significant taxa)", "report"],
            ["14_differential_abundance_ANCOMBC2/working_dir_differential_abundance/", "Collapsed tables, RDS objects, run logs", "audit"],
        ],
        upstream_consumers=[
            "Step 15 (PICRUSt) optionally reads ANCOMBC2 genus-level outputs to compute genus x pathway correlations.",
            "Step 16 (ML) reads to highlight DA taxa with high feature importance.",
            "Step 18 embeds the per-variable summaries.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--alpha F", "0.05", "Significance threshold for the q-value column. Lower = stricter."],
            ["--prv-cut F", "0.10",
             "Prevalence cutoff. Taxa present in <10% of samples are filtered before testing. ANCOMBC2 default."],
            ["--lib-cut N", "1000",
             "Library-size cutoff. Samples with <1000 reads are excluded from the test (standard for ANCOMBC2)."],
            ["--levels LIST", "p,c,o,f,g,s",
             "Domain (d) is excluded by default because it has too few groups for meaningful DA."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 The smart-Rscript locator (cost-weighted)", [
            ("p", "ANCOMBC2 needs ANCOMBC + microbiome + phyloseq -- all heavy Bioconductor packages. Compiling them on Apple Silicon takes 10-20 minutes apiece. The script picks the BEST Rscript based on installation cost rather than mere presence:"),
            ("kv", [
                ("Cost rule",          "Missing CRAN package = 1 point. Missing Bioconductor package = 100 points."),
                ("Tie-breaker",        "Earlier in the candidate list wins. Conda's R is listed FIRST because the qiime2-amplicon-2025.4 env ships ANCOMBC pre-installed."),
                ("Short-circuit",      "A candidate with 0 missing wins immediately."),
            ]),
            ("p", "This avoids the worst case of compiling 5 Bioc packages from source on system R when conda's R already has them. Saves 30-90 minutes on first run."),
        ]),
        ("4.2 Reference-level convention", [
            ("p", "For every categorical variable, the reference level is set to the alphabetically-first level (R's default for factor()). Running with pairwise=TRUE and global=TRUE means the reference choice only affects how coefficients are LABELLED -- the global test is reference-free, and every pairwise contrast is computed regardless of which level was the reference."),
            ("p", "Documented prominently because reviewers often ask 'why is treatment X the reference?' and the answer is 'it does not matter -- you have all contrasts plus the global test'."),
        ]),
        ("4.3 Table collapse", [
            ("p", "qiime taxa collapse --p-level <N> on the filtered feature table for each requested level. Output: feature_table_collapsed_L<N>.qza. ANCOMBC2 needs RAW counts per taxonomic level."),
        ]),
        ("4.4 Phyloseq construction", [
            ("p", "An R helper builds a phyloseq object from each collapsed table + taxonomy + metadata. ANCOMBC2 takes phyloseq directly."),
        ]),
        ("4.5 ANCOMBC2 invocation", [
            ("p", "ANCOMBC2 args used:"),
            ("kv", [
                ("formula",     "~ <variable> (one variable at a time)."),
                ("p_adj_method","holm (per-variable Holm correction)."),
                ("prv_cut",     "0.10 by default."),
                ("lib_cut",     "1000 by default."),
                ("group",       "<variable>."),
                ("struc_zero",  "TRUE -- detect structural zeros."),
                ("neg_lb",      "TRUE."),
                ("global",      "TRUE -- omnibus test."),
                ("pairwise",    "TRUE -- all pairwise contrasts."),
                ("alpha",       "0.05."),
            ]),
        ]),
        ("4.6 Output rendering", [
            ("p", "Five XLSX files per (variable, level): primary (model output table), pairwise, global, structural_zeros, summary. Two PNGs: volcano (LFC vs -log10(q)), heatmap (DA taxa across groups). The full ANCOMBC2 result object is also pickled to RDS in the working dir for audit."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "15_ancombc2.png"
    nodes = [
        dict(id="start", x=5, y=11.5, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="gate", x=5, y=10.7, w=5.4, h=0.6, label="STATUS gate (step 11)", kind="decision"),
        dict(id="rdisc", x=5, y=9.6, w=5.6, h=0.9, label="Smart Rscript locator\n(CRAN=1pt, Bioc=100pt)", kind="decision"),
        dict(id="install", x=5, y=8.3, w=5.4, h=0.9, label="Install missing R packages\nBiocManager::install for Bioc"),
        dict(id="cats", x=5, y=7.0, w=5.4, h=0.7, label="Categorical detection", kind="decision"),
        dict(id="loop", x=5, y=5.7, w=6.0, h=1.0, label="for L in p,c,o,f,g,s:\n  qiime taxa collapse + export\n  for var in CATS: ANCOMBC2(formula = ~ var)"),
        dict(id="emit", x=5, y=4.0, w=5.6, h=1.0, label="primary + pairwise + global + struc_zeros\n+ volcano.png + heatmap.png + RDS", kind="io"),
        dict(id="summ", x=5, y=2.5, w=5.4, h=0.7, label="Cross-level summary xlsx"),
        dict(id="info", x=5, y=1.4, w=5.0, h=0.6, label="Write info file", kind="io"),
        dict(id="end", x=5, y=0.5, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "gate"}),
        dict(**{"from": "gate", "to": "rdisc"}),
        dict(**{"from": "rdisc", "to": "install"}),
        dict(**{"from": "install", "to": "cats"}),
        dict(**{"from": "cats", "to": "loop"}),
        dict(**{"from": "loop", "to": "emit"}),
        dict(**{"from": "emit", "to": "summ"}),
        dict(**{"from": "summ", "to": "info"}),
        dict(**{"from": "info", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_ancombc2_run.sh -- control flow", figsize=(8.0, 12.0))
    section_flow(doc, flow_png, "Smart Rscript discovery + collapse + per-cell ANCOMBC2 + cross-cell summary.")

    section_edge_cases(doc, [
        ("ANCOMBC2 fails on small group",
         "When a group has < 3 samples after lib_cut filtering, ANCOMBC2 returns NA for that contrast. We report it in the run log without aborting."),
        ("Structural zeros vs sampling zeros",
         "ANCOMBC2 distinguishes them. Structural zeros (taxon completely absent in a group) are reported in the dedicated XLSX so users know which 'absences' are biologically meaningful vs sequencing-depth artifacts."),
        ("Cross-platform R compilation",
         "On Apple Silicon, igraph and ANCOMBC need C++ compilation. The smart locator picks conda's R when system R lacks them, sidestepping the compile cost."),
    ])

    section_impl_notes(doc, [
        ("Why not run all variables in a single ANCOMBC2 call",
         "ANCOMBC2 does support multi-covariate formulas, but the per-variable runs give cleaner pairwise contrasts and easier-to-interpret tables. Multi-covariate analysis is captured by step 13's Adonis instead."),
        ("Why exclude domain",
         "Domain has typically 1-3 levels (Bacteria/Archaea/Eukaryota). DA testing on this level is statistically uninformative for 16S microbiome studies."),
    ])

    section_testing(doc, [
        "Verify all 6 levels (p..s) produce ANCOMBC2 outputs.",
        "Verify pairwise xlsx has C(k,2) row groups for k-group variables.",
        "Verify global xlsx exists for every (level, variable).",
        "Verify volcano + heatmap PNGs for at least one cell.",
        "Re-run -> verify the cross-level summary xlsx is regenerated but per-cell RDS is reused.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_16_picrust():
    NAME = "16_mbx_picrust_run.docx"
    doc = make_doc_with_cover(
        "PICRUSt2 functional prediction with NSTI filtering and DA + correlation",
        "Pipeline step 15 -- the most complex environment-management script.",
        "mbx_picrust_run.sh",
        step_id="15",
    )

    section_purpose(doc, (
        "Predict metagenomic content from 16S marker genes via PICRUSt2 "
        "(Douglas et al. 2020). Generates KO, EC, COG, MetaCyc pathway "
        "abundance tables, applies NSTI filtering with audit reporting, runs "
        "Kruskal-Wallis differential abundance per categorical variable, "
        "produces stacked-bar / heatmap / functional PCoA plots, computes "
        "DA-genera x DA-pathway correlation matrices, and emits a "
        "self-contained HTML report. PICRUSt2 lives in its OWN conda env to "
        "protect QIIME2."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing 14_differential_abundance_ANCOMBC2/ (optional)", "yes"],
            ["--nsti F", "User flag", "NSTI filter threshold (default 2.0)", "no"],
            ["--picrust-env-path PATH", "User flag", "Override auto-discovery of PICRUSt2 conda env", "no"],
            ["--threads N", "User flag", "Default = auto-detected", "no"],
            ["--variables LIST", "User flag", "Override auto-detection", "no"],
        ],
        outputs_table=[
            ["15_picrust2/all_picrust2_outputs/{KO,EC,COG}_metagenome/", "Functional metagenome predictions per database", "report"],
            ["15_picrust2/all_picrust2_outputs/pathways_metacyc/", "MetaCyc pathway abundances", "report"],
            ["15_picrust2/all_picrust2_outputs/nsti/", "Per-ASV NSTI tables + filtering summary", "report"],
            ["15_picrust2/all_picrust2_outputs/place_seqs/", "Phylogenetic placement results + log", "audit"],
            ["15_picrust2/picrust2_<Var>/DA_*_KW_<Var>.xlsx + pairwise.xlsx", "KW + Dunn DA per database per variable", "report"],
            ["15_picrust2/picrust2_<Var>/stacked_bar_top20_metacyc_<Var>.{png,pdf}", "Top-20 pathway stacked bar", "report"],
            ["15_picrust2/picrust2_<Var>/heatmap_DA_pathways_<Var>.{png,pdf}", "LFC heatmap with q-stars", "report"],
            ["15_picrust2/picrust2_<Var>/PCoA_BrayCurtis_FUNCTIONAL_<Var>.{png,pdf}", "Functional beta diversity", "report"],
            ["15_picrust2/picrust2_<Var>/correlation_DAgenera_x_DApathways_<Var>.{tsv,png,pdf}", "Spearman correlation matrix", "report"],
            ["15_picrust2/picrust2_report.html", "Self-contained HTML report", "user"],
            ["15_picrust2/Summary_picrust2_NSTI.xlsx", "Per-sample NSTI + sample reliability flags", "report"],
        ],
        upstream_consumers=[
            "Step 18 embeds the picrust2_report.html links + the NSTI table.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--nsti F", "2.0",
             "Drop ASVs with NSTI > F. NSTI = Nearest Sequenced Taxon Index = phylogenetic distance to nearest reference. Default 2.0 is appropriate for gut/soil. Rumen/termite-gut may need higher (3.0) because of fewer reference genomes."],
            ["--picrust-env-path PATH", "(auto)", "Bypass env discovery. Useful on shared HPC."],
            ["--threads N", "auto", "Pass-through to picrust2_pipeline.py --processes."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Conda env management (independent of QIIME2)", [
            ("p", "PICRUSt2 has dependency conflicts with QIIME2 (different scikit-learn). To avoid corrupting the user's QIIME2 env, this script does the following:"),
            ("bullets", [
                "Search every conda/mamba install for an env containing place_seqs.py in bin/.",
                "Order: env named exactly $PICRUST_ENV (default 'mbx_picrust2'), then 'picrust2', then any env name ending in 'picrust2' or 'picrust'.",
                "If none found, create env via mamba create -n mbx_picrust2 -c bioconda picrust2. Falls back to conda create on machines without mamba.",
                "Invoke picrust2 commands by PREPENDING the env's bin/ to PATH for that single command (no `conda activate`, no `conda run`).",
            ]),
            ("p", "PATH-prefixing works because picrust2's shebangs are absolute env-python paths. This sidesteps both `conda run`'s --no-capture-output awkwardness and mamba run's flag-rejection issues."),
        ]),
        ("4.2 PICRUSt2 pipeline", [
            ("p", "picrust2_pipeline.py runs in 4 internal stages: place_seqs (ASV placement onto reference tree), hsp.py (hidden state prediction for 16S copy number + KO/EC/COG/MetaCyc), metagenome_pipeline.py (per-sample functional abundance), pathway_pipeline.py (MetaCyc pathway from EC). Each stage emits its own log."),
        ]),
        ("4.3 NSTI filtering with audit", [
            ("p", "After placement, every ASV has a NSTI score. The script:"),
            ("bullets", [
                "Logs per-sample mean NSTI. Flags samples with mean NSTI > 1.0 (predictions unreliable).",
                "Drops ASVs with NSTI > --nsti.",
                "Records: how many ASVs dropped, what % of total reads they represented, which samples lost most reads.",
                "Flags samples that lost > 50% of reads after filtering -- these become ineligible for downstream DA.",
            ]),
            ("p", "All audit numbers are written to nsti_filtering_summary.txt and Summary_picrust2_NSTI.xlsx -- the latter has a per-sample reliability flag column."),
        ]),
        ("4.4 DA per variable per database", [
            ("p", "Kruskal-Wallis + Dunn pairwise per (variable, database) on functional abundances. BH-corrected within each database. 4 databases x N variables = 4N tests per variable group."),
        ]),
        ("4.5 Visualizations", [
            ("p", "Stacked bars of top-20 MetaCyc pathways grouped by variable (mirrors the taxonomy bar plots). LFC heatmaps with q-value star overlay (visual consistency with ANCOMBC2 heatmaps). Functional PCoA on Bray-Curtis distance between MetaCyc pathway samples -- LABELLED 'functional beta diversity' so users do not confuse it with the taxonomic PCoA."),
        ]),
        ("4.6 Genus x pathway correlation", [
            ("p", "If step 14 (ANCOMBC2) produced genus-level DA results, the script computes Spearman correlation between every DA genus and every DA pathway. Output: TSV + heatmap. This is the most reviewer-requested cross-step output: it ties microbial taxa to predicted function."),
        ]),
        ("4.7 Self-contained HTML report", [
            ("p", "A standalone picrust2_report.html with embedded plots (base64), provenance, parameter tables, and the documented LIMITATION of 16S-based functional prediction (printed prominently as a caveat)."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "16_picrust.png"
    nodes = [
        dict(id="start", x=5, y=11.5, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="gate", x=5, y=10.6, w=5.0, h=0.6, label="STATUS gate (step 11)", kind="decision"),
        dict(id="env", x=5, y=9.5, w=6.0, h=1.0, label="Detect / install PICRUSt2 conda env\n(mamba create -n mbx_picrust2 ...)\nPATH-prefix at run time", kind="decision"),
        dict(id="run", x=5, y=8.0, w=5.6, h=0.9, label="picrust2_pipeline.py\nplace + hsp + metagenome + pathway", kind="io"),
        dict(id="nsti", x=5, y=6.6, w=5.6, h=1.0, label="NSTI filtering + audit\nflag samples with mean NSTI > 1.0\ndrop ASVs > --nsti", kind="decision"),
        dict(id="da", x=5, y=5.1, w=5.6, h=0.9, label="DA per (variable, database)\nKW + Dunn BH"),
        dict(id="plots", x=5, y=3.7, w=5.6, h=0.9, label="stacked bars + heatmaps + functional PCoA"),
        dict(id="corr", x=5, y=2.4, w=5.6, h=0.7, label="ANCOMBC2 outputs available?\n-> genus x pathway correlation", kind="decision"),
        dict(id="html", x=5, y=1.3, w=5.0, h=0.6, label="Build picrust2_report.html"),
        dict(id="end", x=5, y=0.4, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "gate"}),
        dict(**{"from": "gate", "to": "env"}),
        dict(**{"from": "env", "to": "run"}),
        dict(**{"from": "run", "to": "nsti"}),
        dict(**{"from": "nsti", "to": "da"}),
        dict(**{"from": "da", "to": "plots"}),
        dict(**{"from": "plots", "to": "corr"}),
        dict(**{"from": "corr", "to": "html"}),
        dict(**{"from": "html", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_picrust_run.sh -- control flow", figsize=(8.0, 12.5))
    section_flow(doc, flow_png, "Isolated conda env + 4-stage pipeline + per-variable DA + correlation.")

    section_edge_cases(doc, [
        ("PICRUSt2 placement fails for an ASV",
         "Some ASVs fail phylogenetic placement (very divergent sequences). They are reported in the place_seqs log and excluded from downstream calculations."),
        ("Sample with mean NSTI > 1.0",
         "Flagged in Summary_picrust2_NSTI.xlsx with a 'reliability=low' indicator. The user is warned in the HTML report and the sample is highlighted in red on functional PCoA."),
        ("ANCOMBC2 outputs absent",
         "Genus x pathway correlation step is skipped with a friendly note in the report. Other PICRUSt2 outputs are unaffected."),
        ("Apple Silicon mamba post-link issues",
         "Some PICRUSt2 dependencies trigger post-link script failures on M-series. The script tries mamba first (which handles these better) before falling back to conda."),
    ])

    section_impl_notes(doc, [
        ("Functional prediction caveat",
         "All visualizations are subtitled 'Predicted functional profiles -- inferred from 16S rRNA data'. The HTML report has a prominent caveat box explaining that PICRUSt2 inference is NOT a substitute for shotgun metagenomics."),
        ("Why correlation only on DA-genera x DA-pathways",
         "Computing the full N_taxa x N_pathways correlation is O(taxa * pathways) and produces a giant uninterpretable matrix. Restricting to differentially-abundant items (from steps 14 and 15) yields a focused, biologically interpretable matrix."),
    ])

    section_testing(doc, [
        "Verify mbx_picrust2 conda env exists after first run.",
        "Verify NSTI filter actually drops ASVs and the count is logged.",
        "Verify all 4 databases have output files.",
        "Verify HTML report opens in browser and shows embedded plots.",
        "Verify correlation file is empty when ANCOMBC2 outputs are missing (graceful degradation).",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_17_ml():
    NAME = "17_mbx_ml_classifier_run.docx"
    doc = make_doc_with_cover(
        "Random Forest biomarker classifier for every level x categorical variable",
        "Pipeline step 16 -- the predictive (vs descriptive) biomarker step.",
        "mbx_ml_classifier_run.sh",
        step_id="16",
    )

    section_purpose(doc, (
        "For every (taxonomic level x categorical variable) combination, "
        "train a Random Forest with auto cross-validation and report: "
        "accuracy, AUC, F1, sensitivity/specificity, confusion matrix, ROC "
        "curves (one-vs-rest for multi-class), permutation feature "
        "importance (top-20 + full XLSX), SHAP-style per-sample "
        "contributions (Saabas-equivalent), and a per-variable summary "
        "table. Complements ANCOMBC2: predictive vs descriptive."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing 7_taxonomy_csv/ + 8_cleaned_files/", "yes"],
            ["--levels LIST", "User flag", "Default: all 7", "no"],
            ["--variables LIST", "User flag", "Override auto-detection", "no"],
            ["--ntree N", "User flag", "RF tree count (default: 1000)", "no"],
        ],
        outputs_table=[
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/model_metrics.xlsx", "accuracy / AUC / F1 / sens / spec", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/confusion_matrix.{png,pdf}", "Heatmap with cell counts", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/roc_curves.{png,pdf}", "OvR ROC for multi-class", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/feature_importance.xlsx", "All taxa ranked by permutation importance", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/top20_importance.{png,pdf}", "Top-20 feature bar chart", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/shap_per_sample.{png,pdf}", "Per-sample feature contribution heatmap", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/predicted_vs_actual.xlsx", "Per-sample prediction + truth + class probabilities", "report"],
            ["16_ml_biomarkers/<Var>/RF_<level>_by_<Var>/model.rds", "Pickled ranger model for re-use", "audit"],
            ["16_ml_biomarkers/<Var>/Summary_RF_<Var>.xlsx", "Cross-level summary (one row per level)", "report"],
        ],
        upstream_consumers=[
            "Step 18 embeds top-20 importance bars + the cross-variable summary table.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--ntree N", "1000", "Number of trees in the random forest. Higher = more stable importance scores."],
            ["--ncv K", "auto",
             "K-fold CV count. Auto-rule: K = min(5, n_min_class). With small classes K is reduced so each fold has at least one example of every class."],
            ["--levels / --variables", "(default all)", "Same as ezviz."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Why Random Forest after ANCOMBC2", [
            ("p", "ANCOMBC2 answers 'which taxa differ between groups?' (descriptive). Random Forest answers 'can we PREDICT the group from the taxa?' (predictive). Reviewers commonly ask BOTH questions:"),
            ("bullets", [
                "A taxon that is statistically different (ANCOMBC2 q < 0.05) but does not help the classifier may be a false-positive driven by a single sample.",
                "A taxon that helps the classifier but is not flagged by ANCOMBC2 is likely a non-linear / interaction effect (which RF naturally captures).",
            ]),
            ("p", "The two analyses are complementary; reporting both elevates the rigor of the biomarker findings."),
        ]),
        ("4.2 Inputs", [
            ("p", "Reads CLEANED_<L>= paths from 8_cleaned_files/mbx_ezclean_info.txt. Falls back to globbing the 8_cleaned_files/mbX_cleaned_<plural>_level-7/*.xlsx pattern if the info file is absent. Each XLSX has rows = samples, columns = taxa + metadata."),
        ]),
        ("4.3 ranger Random Forest", [
            ("p", "ranger is used for speed (10-100x faster than randomForest on the same data). Args:"),
            ("kv", [
                ("formula",      "<variable> ~ . (using the metadata column as response)."),
                ("num.trees",    "1000."),
                ("importance",   "permutation (gives unbiased importance).") ,
                ("classification","TRUE."),
                ("seed",         "42 (reproducibility)."),
            ]),
        ]),
        ("4.4 Cross-validation strategy", [
            ("p", "K-fold CV with K = min(5, n_min_class). For each fold, ranger is fit on K-1 folds and predictions are accumulated for the held-out fold. After all folds, we have one out-of-fold prediction per sample. Metrics computed on the OOF predictions (accuracy, F1, AUC via pROC)."),
        ]),
        ("4.5 ROC curves (one-vs-rest)", [
            ("p", "For multi-class, we build OvR ROC: for each class, treat it as positive vs all-others. pROC::multiclass.roc gives per-class AUCs. Plot them on shared axes for the ROC PNG."),
        ]),
        ("4.6 Permutation feature importance", [
            ("p", "ranger's importance='permutation' computes, for each feature: the drop in OOB accuracy when that feature's values are shuffled. Higher = more important. We sort and plot the top 20."),
        ]),
        ("4.7 SHAP-style per-sample contributions", [
            ("p", "True SHAP for RF is computationally expensive. We use the Saabas-equivalent: for each sample, walk every tree and accumulate per-feature contributions to the predicted class probability. Saabas approximation matches SHAP within ~5% for tree models and is 10-50x faster."),
            ("p", "Output: a heatmap with samples on Y, top-N features on X, color = signed contribution. Lets the reader see WHICH features mattered for WHICH samples."),
        ]),
        ("4.8 Cross-level summary", [
            ("p", "After all 7 levels for one variable are done, a Summary_RF_<Var>.xlsx is built with one row per level: level, accuracy, AUC, F1, top-1 feature, top-5 features. Lets the reader pick the most-informative taxonomic resolution at a glance."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "17_ml.png"
    nodes = [
        dict(id="start", x=5, y=11.4, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="info", x=5, y=10.5, w=5.4, h=0.7, label="Read mbx_ezclean_info.txt"),
        dict(id="rd", x=5, y=9.4, w=4.6, h=0.7, label="Locate Rscript + _R wrapper"),
        dict(id="pkg", x=5, y=8.3, w=5.6, h=0.9, label="Install: ranger, pROC, openxlsx, ggplot2, pheatmap"),
        dict(id="cats", x=5, y=6.9, w=5.4, h=0.7, label="Categorical detection", kind="decision"),
        dict(id="loop", x=5, y=5.5, w=6.0, h=1.4, label="for var in CATS:\n  for L in 1..7:\n    K-fold CV ranger\n    importance + Saabas SHAP\n    metrics + plots"),
        dict(id="summ", x=5, y=3.5, w=5.0, h=0.7, label="Per-variable cross-level summary xlsx"),
        dict(id="info2", x=5, y=2.4, w=5.0, h=0.6, label="Write mbx_ml_biomarkers_info.txt", kind="io"),
        dict(id="end", x=5, y=1.4, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "info"}),
        dict(**{"from": "info", "to": "rd"}),
        dict(**{"from": "rd", "to": "pkg"}),
        dict(**{"from": "pkg", "to": "cats"}),
        dict(**{"from": "cats", "to": "loop"}),
        dict(**{"from": "loop", "to": "summ"}),
        dict(**{"from": "summ", "to": "info2"}),
        dict(**{"from": "info2", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_ml_classifier_run.sh -- control flow", figsize=(8.0, 11.0))
    section_flow(doc, flow_png, "All 7 levels x N variables x ranger CV + plots in a single per-cell pipeline.")

    section_edge_cases(doc, [
        ("Tiny class size",
         "K-fold breaks down if a class has fewer than K samples. Auto K = min(5, n_min_class) handles this. Below 3 samples per class, the cell is skipped with a warning."),
        ("Multi-class AUC ill-defined",
         "Macro-AUC via pROC::multiclass.roc handles >=3 classes. For 2-class, the standard binary AUC is reported."),
        ("Imbalanced classes",
         "ranger does not auto-balance. We report per-class precision/recall in metrics XLSX; if the user wants class weighting, they can re-run with the model.rds and a custom weight vector."),
    ])

    section_impl_notes(doc, [
        ("Why ranger over randomForest",
         "ranger is multi-threaded, written in C++, and 10-100x faster on datasets >500 features. For our typical 1000-feature genus-level run on 20 samples, ranger finishes in <1 s vs randomForest's 10-30 s."),
        ("Saabas vs full SHAP",
         "Saabas gives per-feature contributions by walking each tree's decision path. It is an approximation of SHAP for tree models and matches within ~5%. Computing exact SHAP for RF is O(2^F) per sample -- impractical."),
    ])

    section_testing(doc, [
        "Verify ROC PNG axes are labelled with per-class AUC values.",
        "Verify top20_importance.png matches the top 20 rows of feature_importance.xlsx.",
        "Verify SHAP heatmap dimensions = (samples, top-N features).",
        "Run with --ntree 100 -> verify a noticeably faster run with similar accuracy.",
        "Re-run with model.rds present -> verify cells are skipped (idempotency).",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_18_networks():
    NAME = "18_mbx_network_run.docx"
    doc = make_doc_with_cover(
        "Co-occurrence networks (CLR + Spearman + BH-FDR + Louvain)",
        "Pipeline step 17 -- the ecological-systems view.",
        "mbx_network_run.sh",
        step_id="17",
    )

    section_purpose(doc, (
        "For every (taxonomic level x categorical variable) combination "
        "build a global network (all samples = ecological backbone) and "
        "per-group networks. Use CLR-transform + Spearman + BH-FDR + Louvain "
        "modularity (Gloor et al. 2017's lightweight compositional-aware "
        "alternative to SparCC and SpiecEasi). Output: edge / node tables, "
        "GraphML for Cytoscape, network plots, hub-taxa lists, and a "
        "side-by-side group_comparison.xlsx."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing 8_cleaned_files/", "yes"],
            ["--levels LIST", "User flag", "Default: g,f,s (most ecologically meaningful)", "no"],
            ["--prevalence-threshold F", "User flag", "Filter (default 0.20)", "no"],
            ["--rho-threshold F", "User flag", "Edge inclusion (default 0.30)", "no"],
            ["--q-threshold F", "User flag", "Edge BH-q ceiling (default 0.05)", "no"],
            ["--min-group-n N", "User flag", "Minimum N per per-group network (default 8)", "no"],
        ],
        outputs_table=[
            ["17_co_occurrence_networks/global_networks/network_<level>/edges.tsv + nodes.tsv", "Edge / node tables", "report, Cytoscape"],
            ["17_co_occurrence_networks/global_networks/network_<level>/network.graphml", "Cytoscape/Gephi import", "user"],
            ["17_co_occurrence_networks/global_networks/network_<level>/network_plot.{png,pdf}", "Force-directed visualization with module coloring", "report"],
            ["17_co_occurrence_networks/global_networks/network_<level>/network_summary.xlsx", "Edge count, density, modularity, etc.", "report"],
            ["17_co_occurrence_networks/global_networks/network_<level>/hub_taxa.xlsx", "Top 10% by combined degree+betweenness+hub_score", "report"],
            ["17_co_occurrence_networks/global_networks/network_<level>/modules.tsv", "Louvain module assignments", "report"],
            ["17_co_occurrence_networks/<Var>/network_<level>_by_<Var>/per_group/<group>/...", "Per-group equivalents", "report"],
            ["17_co_occurrence_networks/<Var>/network_<level>_by_<Var>/group_comparison.xlsx", "Side-by-side metrics across groups", "report"],
            ["17_co_occurrence_networks/<Var>/network_<level>_by_<Var>/multi_panel_plot.{png,pdf}", "Group networks side-by-side", "report"],
        ],
        upstream_consumers=[
            "Step 18 embeds the network plots + hub-taxa tables.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--prevalence-threshold F", "0.20",
             "Drop taxa present in fewer than this fraction of samples. Lower -> more taxa, more spurious edges. 0.20 is a community-standard floor."],
            ["--rho-threshold F", "0.30",
             "Include edge only when |Spearman rho| >= F. 0.30 = 'moderate' correlation."],
            ["--q-threshold F", "0.05",
             "Include edge only when BH-adjusted p < F. Strict significance gate."],
            ["--min-group-n N", "8",
             "Skip per-group networks where the group has < N samples. Below 8 the inferred network has too few edges to be reliable."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Why CLR + Spearman, not SparCC / SpiecEasi", [
            ("p", "Microbiome data is COMPOSITIONAL: counts are constrained to a simplex, so naive Pearson correlation gives spurious 'negative correlation' artifacts. Three solutions exist:"),
            ("kv", [
                ("SparCC",     "iterative bootstrap on log-ratios. Python-only, slow, no R port."),
                ("SpiecEasi",  "graphical-lasso on CLR-transformed data. Bioconductor, heavy compile, often fails on Apple Silicon."),
                ("CLR + Spearman", "Gloor et al. 2017 'Microbiome datasets are compositional'. Light, CRAN-only, fast, scientifically defensible."),
            ]),
            ("p", "We chose CLR + Spearman because (a) it runs entirely on CRAN, (b) Spearman is rank-based -> robust to outliers, (c) it composes cleanly with BH-FDR adjustment."),
        ]),
        ("4.2 Pipeline per cell", [
            ("p", "For each (level, scope) where scope = 'global' or 'per_group':"),
            ("bullets", [
                "Read cleaned XLSX. Identify metadata vs taxa columns by intersecting with metadata header.",
                "Filter to scope: global = all samples; per_group = samples in one group of one variable.",
                "Prevalence filter: drop taxa present in < threshold of samples.",
                "Pseudocount: add 0.5 of the smallest non-zero value (or absolute 0.5 if smaller).",
                "CLR transform: x_clr = log(x_i) - mean_log(x_row).",
                "psych::corr.test: returns r and p matrices over CLR values.",
                "BH-adjust ALL pairwise raw p-values (not just per-row).",
                "Filter edges: |rho| >= rho_threshold AND q <= q_threshold.",
                "Build igraph object; drop isolated nodes.",
                "Compute degree, betweenness, hub_score per node.",
                "Louvain max-modularity for module detection.",
                "Hub taxa: top 10% by combined rank of (degree + betweenness + hub_score).",
                "Export edges.tsv, nodes.tsv, modules.tsv, GraphML, plot PNG/PDF, hub_taxa.xlsx, network_summary.xlsx.",
            ]),
        ]),
        ("4.3 Group-comparison output", [
            ("p", "After all per-group networks for one (level, variable) are built, group_comparison.xlsx aggregates: edge count, density, modularity, # modules, # hub taxa for each group, side-by-side. multi_panel_plot.png renders all per-group networks on a single canvas."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "18_networks.png"
    nodes = [
        dict(id="start", x=5, y=11.4, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="info", x=5, y=10.5, w=5.4, h=0.7, label="Read mbx_ezclean_info.txt"),
        dict(id="pkg", x=5, y=9.4, w=5.4, h=0.7, label="igraph + psych + openxlsx"),
        dict(id="cats", x=5, y=8.3, w=5.4, h=0.7, label="Categorical detection", kind="decision"),
        dict(id="loop", x=5, y=6.6, w=6.4, h=1.5, label="for L in (g,f,s):\n  global network (all samples)\n  for var in CATS:\n    for grp in groups(var):\n      if N(grp) >= min_group_n: per-group network"),
        dict(id="pipe", x=5, y=4.7, w=6.0, h=1.6, label="Per cell:\n prevalence -> CLR -> Spearman + BH\n |rho|>=, q<= filter\n graph + Louvain + hubs"),
        dict(id="export", x=5, y=2.8, w=5.6, h=1.0, label="edges.tsv + nodes.tsv + modules.tsv\nnetwork.graphml + plot + summary + hubs", kind="io"),
        dict(id="cmp", x=5, y=1.3, w=5.0, h=0.7, label="group_comparison.xlsx + multi_panel_plot"),
        dict(id="end", x=5, y=0.4, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "info"}),
        dict(**{"from": "info", "to": "pkg"}),
        dict(**{"from": "pkg", "to": "cats"}),
        dict(**{"from": "cats", "to": "loop"}),
        dict(**{"from": "loop", "to": "pipe"}),
        dict(**{"from": "pipe", "to": "export"}),
        dict(**{"from": "export", "to": "cmp"}),
        dict(**{"from": "cmp", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_network_run.sh -- algorithm flow", figsize=(8.0, 12.5))
    section_flow(doc, flow_png, "CLR + Spearman + BH + Louvain pipeline applied per cell.")

    section_edge_cases(doc, [
        ("All edges fail the q-threshold",
         "An empty edge table is still written for transparency, but no plot or hub table is produced. Logged to network_summary.xlsx as 'no significant edges'."),
        ("Group with N < min_group_n",
         "Skipped. Logged in group_comparison.xlsx with status 'too few samples'."),
        ("Pseudocount choice",
         "0.5 * min_nonzero is empirically robust. If the table has zero non-zero counts (pathological), absolute 0.5 is used."),
    ])

    section_impl_notes(doc, [
        ("Why default to g,f,s (genus, family, species)",
         "These three levels are the most ecologically informative for co-occurrence: domain/phylum aggregate too aggressively, class/order too coarsely. Family/Genus/Species capture functional cohorts."),
        ("BH adjustment is global, not per-row",
         "Each cell tests N*(N-1)/2 pairs. Adjusting BH across ALL pairs gives proper FDR control; per-row adjustment would over-reject."),
    ])

    section_testing(doc, [
        "Verify GraphML opens in Cytoscape.",
        "Verify hub_taxa.xlsx contains exactly 10% of nodes.",
        "Verify modules.tsv assigns every node to a module ID.",
        "Verify per-group network exists only for groups with N >= min_group_n.",
        "Verify multi_panel_plot has one subplot per group.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_19_final_report():
    NAME = "19_mbx_final_report.docx"
    doc = make_doc_with_cover(
        "Self-contained HTML + PDF report consolidating every executed step",
        "Pipeline step 18 -- the user-facing summary.",
        "mbx_final_report.sh",
        step_id="18",
    )

    section_purpose(doc, (
        "Produce ONE self-contained HTML report (with all figures embedded "
        "as base64) plus a print-ready A4 PDF that summarises every step "
        "actually run in this output directory. Steps that were not run are "
        "listed explicitly with the exact command needed to add them on a "
        "subsequent run. Designed to satisfy a deeply technical reviewer "
        "(every parameter, every citation, every algorithmic choice)."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<mbX_pro_outputs_dir>", "User CLI", "Containing every previous step's directory", "yes"],
            ["--no-pdf", "User flag", "Skip PDF rendering (HTML only)", "no"],
            ["--logo PATH", "User flag", "Override logo path", "no"],
        ],
        outputs_table=[
            ["18_final_report/mbX_pro_final_report.html", "Self-contained HTML, all images base64-embedded", "user"],
            ["18_final_report/mbX_pro_final_report.pdf", "Print-ready A4 PDF (single column, no cropping)", "user"],
            ["18_final_report/mbx_final_report_info.txt", "Provenance file (timestamp, R version, PDF engine used)", "audit"],
        ],
        upstream_consumers=[
            "End user (the only consumer).",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--no-pdf", "off", "Skip PDF generation. Useful when no headless browser is installed."],
            ["--logo PATH", "(auto)", "Override logo discovery."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Step discovery (canonical manifest)", [
            ("p", "The script defines a hard-coded manifest of every pipeline step:"),
            ("code",
             "step_id|directory_name|info_file|script_name|short_description"),
            ("p", "For each entry, it checks if the directory exists AND the info file exists. The result is a JSON-encoded ran/skipped table that R consumes. This is what makes the report 'partial-run aware': only ran steps render full sections; skipped ones render explanatory placeholders."),
        ]),
        ("4.2 Logo + PDF engine discovery", [
            ("p", "Logo: search current dir, ~/bin/, /usr/local/share/mbx_pro/, /opt/mbx_pro/. If none found, render without logo (graceful)."),
            ("p", "PDF engine preference: (1) Headless Chrome / Chromium / Edge -- best CSS support. (2) wkhtmltopdf (auto-install via brew if missing). (3) macOS cupsfilter (built-in fallback, lower quality). (4) HTML-only with print-instruction page."),
        ]),
        ("4.3 R-rendered HTML", [
            ("p", "The R script uses htmltools to build the HTML. Critical implementation note: we use `htmltools::save_html()`, NOT `writeLines(as.character(doc))`, because the latter silently drops the entire <head> section -- causing CSS to vanish from PDFs."),
            ("p", "All images are read via base64enc::base64encode() and inlined as `data:image/png;base64,...` URIs. The result is a single self-contained HTML file the user can attach to email."),
            ("p", "Path display: every absolute path under MBX_OUT_DIR is converted to a relative path by stripping the MBX_OUT_DIR prefix. Info files inside each step's directory keep absolute paths -- they stay portable across machines."),
        ]),
        ("4.4 CSS strategy for both screen and print", [
            ("p", "Inline CSS in <head> with media queries:"),
            ("kv", [
                ("@page", "size: A4; margin: 14mm 12mm. Critical for PDF rendering."),
                ("@media print", ".fig-grid -> display: block !important. Forces single-column layout in PDF (figures grid on screen, stacked in print)."),
                ("figure.fig img", "max-width: 100%; max-height: 84vh; object-fit: contain. Prevents oversize figures from being cropped."),
                (".cover-meta > div", "min-width: 0; overflow-wrap: anywhere; word-break: break-word. Prevents long path strings from overflowing the cover-meta grid."),
            ]),
        ]),
        ("4.5 Headless Chrome PDF rendering", [
            ("p", "The Chrome flags used to ensure proper PDF generation:"),
            ("bullets", [
                "--headless=new",
                "--disable-gpu",
                "--no-sandbox",
                "--virtual-time-budget=15000  (let JS/CSS settle)",
                "--run-all-compositor-stages-before-draw",
                "--hide-scrollbars",
                "--print-to-pdf-no-header",
                "--print-to-pdf=<output.pdf>",
                "<input.html>",
            ]),
            ("p", "These flags collectively prevent the most common PDF rendering bugs: missing CSS (settle time), partial figures (compositor sync), scrollbar artifacts, browser-default header/footer."),
        ]),
        ("4.6 What goes in each section", [
            ("p", "For each step that ran, the report shows:"),
            ("bullets", [
                "Section heading + one-line scientific rationale.",
                "Citation block (paper + DOI for the algorithm).",
                "Parameter table (every flag's value, including defaults that were inherited).",
                "Verbatim info file content (collapsible).",
                "Embedded figures + tables.",
                "Caveats specific to that step.",
            ]),
            ("p", "For each step that did NOT run, the report shows:"),
            ("bullets", [
                "Section heading + 'NOT RUN' badge.",
                "Plain-language explanation of what would have been computed.",
                "Exact copy-pasteable command to add it.",
            ]),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "19_final_report.png"
    nodes = [
        dict(id="start", x=5, y=11.5, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="disc", x=5, y=10.5, w=5.4, h=0.7, label="Discover ran/skipped steps\n(canonical manifest x disk check)"),
        dict(id="pdf", x=5, y=9.4, w=5.4, h=0.7, label="Find PDF engine\n(Chrome > wkhtmltopdf > cups > none)", kind="decision"),
        dict(id="logo", x=5, y=8.5, w=4.6, h=0.6, label="Find logo (4 search paths)"),
        dict(id="rd", x=5, y=7.6, w=4.6, h=0.6, label="_R wrapper + Rscript"),
        dict(id="json", x=5, y=6.7, w=5.6, h=0.7, label="Write run-manifest JSON"),
        dict(id="html", x=5, y=5.5, w=5.6, h=0.9, label="R: htmltools::save_html()\nbase64-embed every image", kind="io"),
        dict(id="css", x=5, y=4.0, w=5.6, h=0.9, label="Apply @page + print CSS\n(single-col print, max-height clip)"),
        dict(id="pdfgen", x=5, y=2.6, w=5.6, h=0.9, label="Headless Chrome\n--virtual-time-budget --print-to-pdf", kind="io"),
        dict(id="prov", x=5, y=1.3, w=5.0, h=0.6, label="Write info file (provenance)", kind="io"),
        dict(id="end", x=5, y=0.4, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "disc"}),
        dict(**{"from": "disc", "to": "pdf"}),
        dict(**{"from": "pdf", "to": "logo"}),
        dict(**{"from": "logo", "to": "rd"}),
        dict(**{"from": "rd", "to": "json"}),
        dict(**{"from": "json", "to": "html"}),
        dict(**{"from": "html", "to": "css"}),
        dict(**{"from": "css", "to": "pdfgen"}),
        dict(**{"from": "pdfgen", "to": "prov"}),
        dict(**{"from": "prov", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbx_final_report.sh -- control flow", figsize=(8.0, 12.5))
    section_flow(doc, flow_png, "Discovery + R-rendered HTML + headless-Chrome PDF.")

    section_edge_cases(doc, [
        ("htmltools::as.character() drops <head>",
         "Fixed by switching to htmltools::save_html(doc, file=path, libdir='_assets'). Without this, the printed PDF has none of the inline CSS, leading to oversized logo + cropped figures."),
        ("Cover-meta box overflow",
         "Long absolute paths could break out of the grid layout. Solved with min-width:0 + overflow-wrap:anywhere on .cover-meta > div, plus a `.path-val` monospace span for the output path."),
        ("Figures cropped in PDF but full in HTML",
         "The print media query was missing the single-column override. Adding `@media print { .fig-grid { display: block !important; } figure.fig img { max-height: 84vh; object-fit: contain; }}` resolved both issues."),
        ("No headless browser available",
         "Falls through wkhtmltopdf -> cupsfilter -> HTML-only. The HTML-only mode prints a banner asking the user to open the file in a browser and use 'Save as PDF'."),
    ])

    section_impl_notes(doc, [
        ("Self-contained HTML",
         "Embedding every image as base64 makes the file 5-15 MB but lets the user email it without losing figures. Trade-off accepted."),
        ("Relative paths in the report, absolute in info files",
         "Reasoning: the report is shared across machines (relative paths survive); the info files are run-local provenance (absolute paths are unambiguous)."),
        ("Why Chrome over wkhtmltopdf",
         "wkhtmltopdf does not render modern CSS Grid / Flexbox correctly. Chrome's print engine matches what the user sees on screen byte-for-byte."),
    ])

    section_testing(doc, [
        "Run after a complete pipeline -> verify all 18 sections present.",
        "Run after a partial pipeline (skip step 16) -> verify ML section shows 'NOT RUN' placeholder with command.",
        "Open PDF in Adobe Reader and Preview -> verify figures NOT cropped.",
        "Print HTML from Chrome via cmd-P -> verify same single-column layout.",
        "Run --no-pdf -> verify HTML produced and PDF skipped silently.",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def build_20_orchestrator():
    NAME = "20_mbXPro_orchestrator.docx"
    doc = make_doc_with_cover(
        "End-to-end single-command orchestration of every pipeline step",
        "The mbXPro command -- the lay-user-facing front door.",
        "mbXPro",
        step_id="orchestrator",
    )

    section_purpose(doc, (
        "Provide a single command (mbXPro <fastq_dir> <metadata_file>) that "
        "runs every pipeline step (0 -> 18) in the correct order, into a "
        "SINGLE shared output directory, with full logging, idempotent "
        "resume, optional step-skipping, and a final summary. The user does "
        "not need to know any internal command, parameter, or file path."
    ))

    section_inputs_outputs(doc,
        inputs_table=[
            ["<fastq_dir>", "User CLI", "Directory containing raw FASTQ", "yes"],
            ["<metadata_file>", "User CLI", "QIIME2 metadata", "yes"],
            ["--resume", "User flag", "Reuse the most-recent existing output dir", "no"],
            ["--out PATH", "User flag", "Force a specific output dir", "no"],
            ["--from N / --to N", "User flags", "Step range", "no"],
            ["--skip LIST", "User flag", "Comma-separated step IDs to skip", "no"],
            ["--forward-primer / --reverse-primer", "User flags", "Manual primer override (skip step 0)", "no"],
            ["--no-pdf", "User flag", "HTML-only step 18", "no"],
            ["--threads N", "User flag", "Override auto-detection", "no"],
            ["--keep-going", "User flag", "Continue past failed steps", "no"],
        ],
        outputs_table=[
            ["mbX_pro_outputs_<TS>/", "Single shared output dir, all 18 step subdirectories", "user"],
            ["mbX_pro_outputs_<TS>/_pipeline_log/pipeline_master.log", "Full run log (tee'd from every step)", "audit"],
            ["mbX_pro_outputs_<TS>/_pipeline_log/<step>.log", "Per-step log files", "audit"],
            ["mbX_pro_outputs_<TS>/mbXPro_run_info.txt", "Provenance: invocation, timestamps, versions", "audit"],
        ],
        upstream_consumers=[
            "Step 0/1 reuse MBX_OUT_DIR (set by this orchestrator).",
            "Step 18 reads the entire output tree to build the final report.",
        ],
    )

    section_parameters(doc,
        params_table=[
            ["--resume", "off",
             "Find the most-recent mbX_pro_outputs_* dir alongside <fastq_dir>. Use it as MBX_OUT_DIR. Idempotent steps skip themselves automatically."],
            ["--out PATH", "(auto)",
             "Force a specific output directory. Overrides both auto-create and --resume."],
            ["--from N / --to N", "0/18",
             "Run only steps in the inclusive range [from, to]. Earlier steps must already exist if from > 0."],
            ["--skip LIST", "(empty)",
             "e.g. --skip 15,16,17 to skip PICRUSt + ML + networks (the three slowest optional steps)."],
            ["--forward-primer / --reverse-primer", "(unset)",
             "When set, step 0 is skipped and these primers are passed directly to step 3 / step 5."],
            ["--keep-going", "off",
             "Continue running subsequent steps even if one fails. Default is abort-on-first-failure."],
        ],
    )

    section_algorithm(doc, [
        ("4.1 Single shared output directory (the invariant)", [
            ("p", "The orchestrator creates EXACTLY ONE output directory per run. It exports MBX_OUT_DIR as an environment variable. Steps 0 and 1 (which would otherwise generate their own timestamps) check this variable and reuse it. Result: every step's output lands in the same directory."),
        ]),
        ("4.2 Pre-flight checks", [
            ("bullets", [
                "Verify <fastq_dir> exists and contains *.fastq.gz files.",
                "Verify <metadata_file> exists and is readable.",
                "Verify QIIME2 conda env is active (qiime in PATH).",
                "Verify Rscript is on PATH.",
                "Verify all 19 step scripts are installed in PATH.",
                "Auto-detect CPU thread count.",
            ]),
        ]),
        ("4.3 Output dir resolution", [
            ("p", "Three modes:"),
            ("kv", [
                ("Default",  "create <parent_of_FASTQ>/mbX_pro_outputs_<TS>/."),
                ("--resume", "find the latest mbX_pro_outputs_* matching the parent. Reuse it."),
                ("--out X",  "use X regardless. Create if missing."),
            ]),
        ]),
        ("4.4 Step-level idempotency", [
            ("p", "step_already_done(id) checks for a per-step sentinel file (typically the info file written at the END of that step). When sentinel exists, the step is skipped with a [SKIP] line. The user can `rm` a single sentinel to force re-run of just that step."),
        ]),
        ("4.5 run_step()", [
            ("p", "For each step:"),
            ("bullets", [
                "Print a banner with step id + title.",
                "Open _pipeline_log/<id>_<name>.log via `tee -a`.",
                "Time the run (date +%s).",
                "Invoke the step script with the full command.",
                "On non-zero exit, abort (or continue with --keep-going).",
                "Append per-step status to mbXPro_run_info.txt.",
            ]),
        ]),
        ("4.6 STEP_*_CMD definitions", [
            ("p", "Each step's command is built once with the appropriate arguments. For example:"),
            ("code",
             "STEP3_CMD=\"mbx_dada2_parameter_finder.sh \\\n"
             "    \\\"$ARTIFACT_QZA\\\" \\\n"
             "    --forward-primer \\\"$FWD_PRIMER\\\" \\\n"
             "    --reverse-primer \\\"$REV_PRIMER\\\"\""),
            ("p", "Every command is recorded in mbXPro_run_info.txt for auditing."),
        ]),
        ("4.7 Final summary", [
            ("p", "After all steps, print a wall-time summary, per-step status, and the path to the final report. If step 18 ran successfully, print a clickable hyperlink (file:// URL) to mbX_pro_final_report.html."),
        ]),
    ])

    flow_png = DIAGRAMS_DIR / "20_orchestrator.png"
    nodes = [
        dict(id="start", x=5, y=11.6, w=2.0, h=0.5, label="START", kind="term"),
        dict(id="parse", x=5, y=10.7, w=5.0, h=0.6, label="Parse args (fastq + meta + flags)"),
        dict(id="pre", x=5, y=9.7, w=5.6, h=0.8, label="Pre-flight: env, scripts, qiime, R", kind="decision"),
        dict(id="outdir", x=5, y=8.5, w=5.6, h=0.9, label="Resolve / create MBX_OUT_DIR\n(default | --resume | --out)"),
        dict(id="export", x=5, y=7.2, w=5.0, h=0.6, label="export MBX_OUT_DIR"),
        dict(id="log", x=5, y=6.4, w=5.4, h=0.6, label="Init pipeline_master.log (tee)"),
        dict(id="loop", x=5, y=4.7, w=6.0, h=2.0, label="for step in 0..18:\n  step_already_done? -> [SKIP]\n  if --skip / --from / --to: skip\n  run_step + per-step log\n  fail -> abort (or --keep-going)", kind="decision"),
        dict(id="final", x=5, y=2.4, w=5.4, h=0.7, label="Final summary + report URL"),
        dict(id="info", x=5, y=1.4, w=5.0, h=0.6, label="Write mbXPro_run_info.txt", kind="io"),
        dict(id="end", x=5, y=0.5, w=2.0, h=0.4, label="END", kind="term"),
    ]
    edges = [
        dict(**{"from": "start", "to": "parse"}),
        dict(**{"from": "parse", "to": "pre"}),
        dict(**{"from": "pre", "to": "outdir"}),
        dict(**{"from": "outdir", "to": "export"}),
        dict(**{"from": "export", "to": "log"}),
        dict(**{"from": "log", "to": "loop"}),
        dict(**{"from": "loop", "to": "final"}),
        dict(**{"from": "final", "to": "info"}),
        dict(**{"from": "info", "to": "end"}),
    ]
    render_flowchart(nodes, edges, flow_png, title="mbXPro orchestrator -- control flow", figsize=(8.0, 12.5))
    section_flow(doc, flow_png, "Single shared output dir + per-step idempotency + clean failure semantics.")

    section_edge_cases(doc, [
        ("--resume picks up wrong directory",
         "If multiple mbX_pro_outputs_* exist alongside <fastq_dir>, --resume picks the most-recently-modified. The user can force a specific one via --out instead."),
        ("Step succeeds but sentinel missing",
         "step_already_done relies on the per-step info file existing. If a step finishes everything but crashes during info-file write, the next --resume re-runs the entire step. Acceptable trade-off: safer than running on incomplete state."),
        ("--from 5 with steps 0-4 missing",
         "Each step has its own auto-discovery of inputs (e.g. step 5 reads step 0/4 outputs). Running --from 5 with prior outputs missing causes that step's pre-flight to abort, not the orchestrator's. Pre-flight messages are descriptive ('feature_table.qza not found')."),
        ("--keep-going + late failure",
         "If step 4 fails but --keep-going is set, the orchestrator runs step 5 anyway. Step 5 will likely fail at its own pre-flight. The final summary clearly reports which steps ran vs failed."),
    ])

    section_impl_notes(doc, [
        ("Logging via exec + tee",
         "exec > >(tee -a logfile) 2>&1 redirects ALL stdout AND stderr to both terminal and log file. This works on bash 3.2 (no `&>>` ambiguity) and captures every echo / printf without per-call decoration."),
        ("env -u for child env-var scrubbing",
         "Step scripts that need `_R` (system R, not conda R) do their own env stripping. The orchestrator does NOT strip globally because qiime needs the conda env intact. Each step is responsible for its own env hygiene."),
        ("Why not use Snakemake / Nextflow",
         "Both are excellent but require their own runtime + DSL. mbX Pro's audience is lay users who run one command. A bash orchestrator with idempotent step-level checks gives the same correctness with zero extra dependency."),
    ])

    section_testing(doc, [
        "Run mbXPro <fastq> <meta> -> verify all 18 steps produce outputs in a SINGLE timestamped directory.",
        "Re-run with --resume -> verify every step is skipped (idempotent).",
        "Hand-delete one step's info file -> verify --resume re-runs ONLY that step.",
        "Run with --from 12 --to 15 -> verify only steps 12-15 ran.",
        "Run with --skip 15,16,17 -> verify those three steps not in the log.",
        "Run with --keep-going + a deliberately broken step -> verify subsequent steps still attempt.",
        "Verify mbXPro_run_info.txt has per-step exit codes + timestamps.",
        "Verify --help output is not truncated (the awk-script fix from past iteration).",
    ])

    doc.save(HERE / NAME)
    print(f"[OK] {NAME}  ({(HERE / NAME).stat().st_size/1024:.1f} KB)")


def main():
    builders = [
        build_13_alpha,
        build_14_beta,
        build_15_ancombc2,
        build_16_picrust,
        build_17_ml,
        build_18_networks,
        build_19_final_report,
        build_20_orchestrator,
    ]
    for b in builders:
        b()


if __name__ == "__main__":
    main()
