<p align="center">
  <img src="assets/mbX_Pro_icon.png" alt="mbX Pro logo" width="280">
</p>

<h1 align="center">mbX Pro</h1>

<p align="center">
  <b>Single-command, end-to-end 16S rRNA microbiome pipeline.</b><br>
  Raw FASTQ &rarr; ASVs &rarr; Taxonomy &rarr; Diversity &rarr; Differential abundance &rarr; Functional inference &rarr; ML biomarkers &rarr; Networks &rarr; Publication-ready HTML + PDF report.
</p>

<p align="center">
  <a href="#one-command-quickstart"><b>Quickstart</b></a> &middot;
  <a href="documentation/how_to_run_mbX_Pro.docx"><b>End-user manual</b></a> &middot;
  <a href="documentation/mbXPro_documentation.docx"><b>Technical reference</b></a> &middot;
  <a href="#citation"><b>Citation</b></a>
</p>

---

## One-command quickstart

```bash
# 1. Clone this repo
git clone https://github.com/utsavlamichhane/mbXPro.git
cd mbXPro

# 2. Install (copies scripts to ~/bin and adds it to PATH)
bash install/install_mbXPro.sh

# 3. Open a NEW terminal and activate QIIME2
conda activate qiime2-amplicon-2025.4

# 4. Run the entire pipeline with ONE command
mbXPro /path/to/FASTQ /path/to/metadata.txt
```

That's it. mbXPro will create exactly **one** output directory next to your FASTQ folder (`mbX_pro_outputs_<timestamp>/`) and run all 18 analytical steps in order, ending in an HTML and PDF report.

---

## What you get

| Step | Output sub-dir | What it produces |
|---|---|---|
| 0 | `0_primer_handling/` | Auto-detected forward + reverse primers (or manual override) |
| 1 | `1_manifest_file/` | QIIME2-format sample manifest |
| 2 | `2_first_artifact_file/` | QIIME2 `.qza` artifact + summary `.qzv` |
| 3 | `3_dada2_parameters/` | Truncation lengths picked from quality profile |
| 4 | `4_dada2_outputs/` | Feature table, representative sequences, DADA2 stats |
| 5 | `5_classifier_working_dir/` | Greengenes2 reference + trained Naive-Bayes classifier |
| 6 | `6_classifier_taxonomy/` | Per-ASV taxonomic assignment |
| 7 | `7_taxonomy_csv/` | level-1 to level-7 CSVs + mito/chloro-filtered table |
| 8 | `8_cleaned_files/` | mbX::ezclean outputs per taxonomic level |
| 9 | `9_visualization_entire/` | Stacked-bar plots (ezviz) per categorical variable |
| 10 | `10_stats/` | Kruskal-Wallis + Dunn statistics (ezstat) per variable |
| 11 | `11_pre_diversity/` | Phylogenetic tree + sampling-depth selection |
| 12 | `12_alpha_diversity_results/` | ASV richness, Shannon, Simpson, Faith PD, Pielou |
| 13 | `13_beta_diversity_results/` | PCoA + PERMANOVA + Adonis + dispersion tests |
| 14 | `14_differential_abundance_ANCOMBC2/` | Compositional-aware DA (global + pairwise) |
| 15 | `15_picrust2/` | KO + EC + MetaCyc functional inference + heatmaps |
| 16 | `16_ml_biomarkers/` | Random Forest classifier + SHAP-style importance |
| 17 | `17_co_occurrence_networks/` | CLR + Spearman + Louvain modules + hub taxa |
| 18 | `18_final_report/` | Self-contained HTML report + A4 PDF |

---

## Repository layout

```
mbXPro/
├── README.md                          (this file)
├── LICENSE                            (MIT)
├── VERSION
├── CITATION.cff
├── CHANGELOG.md
├── assets/
│   └── mbX_Pro_icon.png               (logo, embedded in HTML report)
├── scripts/
│   ├── mbXPro                         (the one-command orchestrator)
│   └── mbx_*.sh, create_*.sh, ...     (20 step scripts)
├── install/
│   ├── install_mbXPro.sh
│   └── uninstall_mbXPro.sh
├── documentation/
│   ├── mbXPro_documentation.docx      (full technical reference)
│   ├── how_to_run_mbX_Pro.docx        (end-user manual)
│   └── how_to_upload_mbX_Pro.docx     (publishing guide)
└── examples/
    └── example_metadata.txt
```

---

## Requirements

| Tool | Version | Notes |
|---|---|---|
| OS | macOS 12+ or Linux (any modern distro) | Bash 3.2+ is enough |
| Conda | Miniconda or Anaconda | for the QIIME2 env |
| QIIME2 | `qiime2-amplicon-2025.4` | activate before running |
| R | 4.2+ | install **system-wide**, NOT inside conda |
| Internet | yes (first run) | downloads Greengenes2 + installs PICRUSt2 conda env |
| Disk | 5-15 GB | per typical project |
| RAM | 8 GB minimum, 16+ recommended | larger projects need more |

---

## Citation

The dedicated **mbX Pro** paper is in preparation. While you wait, please cite the underlying **mbX** R-package paper that powers steps 8 (ezclean), 9 (ezviz), and 10 (ezstat):

> Lamichhane U., Lourenco J. (2025). **mbX: An R Package for Streamlined Microbiome Analysis.** *Stats* 8(2): 44. doi:[10.3390/stats8020044](https://doi.org/10.3390/stats8020044)

Every step that uses an external method also cites its underlying paper in the final HTML report.

---

## License

MIT -- see [LICENSE](LICENSE).

---

## Authors & support

- **Utsav Lamichhane** -- Department of Animal and Dairy Science, University of Georgia
- **Jeferson Lourenco** -- Department of Animal and Dairy Science, University of Georgia

For issues, open one on this repo's GitHub issue tracker.
