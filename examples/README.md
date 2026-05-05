# Example metadata file

`example_metadata.txt` shows the **minimum acceptable format** that mbX Pro
expects for the metadata file you pass as the second argument to `mbXPro`.

## Format rules (enforced)

1. Tab-separated values (`.txt` or `.tsv`) -- comma-separated `.csv` is also
   accepted.
2. The **first column header** must be one of:
   `sample-id`, `id`, `sampleid`, `sample id`, `featureid`, `feature id`, `feature-id`
   (case-insensitive).
3. Sample IDs must match the sample IDs that mbX Pro extracts from your
   FASTQ filenames (typically `S1`, `S2`, ..., `Sample-01`, etc.).
4. No duplicate sample IDs.
5. No empty sample IDs.
6. No leading or trailing whitespace inside any cell.
7. (Optional) An immediate row 2 starting with `#q2:types` containing per-column
   QIIME2 type hints (`categorical` / `numeric`) is detected and skipped
   automatically by every step.

## What the categorical columns are used for

mbX Pro **auto-detects** all categorical metadata columns (anything that's not
numeric and not all-unique and not a constant) and runs every diversity test,
visualisation, statistical test, ANCOMBC2 contrast, PICRUSt2 functional comparison,
ML classifier, and co-occurrence network on **every categorical column you provide**.

So the more thoughtfully you label your samples, the more useful the report.

## Where to put it

Anywhere. You pass its absolute path to `mbXPro`:

```bash
mbXPro /path/to/FASTQ /path/to/metadata.txt
```
