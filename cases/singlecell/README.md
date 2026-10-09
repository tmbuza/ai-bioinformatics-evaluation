# ABE-006: original synthetic single-cell teaching case

Four distinct donors contribute one sample each: D1 and D2 are controls; D3 and D4 are treated. This is an observational comparison, with no randomization or clinical outcomes supplied. Cells are nested within donors. Counts and cell labels are constructed, not derived from sequencing.

## Data dictionary and scope
- donors.tsv: unique donor_id and condition.
- cells.tsv: 2,000 unique cell_id values, donor_id, and an assigned cluster (C1 or C2). Rows represent retained cells; raw counts and pre-filter cells are absent.
- marker-evidence.tsv: invented percentages of C1 cells detecting each fictional marker. These are supplied summaries, not independently recoverable from cells.tsv.
- teaching-marker-reference.tsv: fictional candidate types and expected high marker combinations. For this exercise, high detection means at least 50%. This threshold is a teaching convention, not a biological annotation standard. Neither markers nor types refer to real genes or cell types.

No expression matrix, QC metrics, doublet calls, clustering settings, embedding coordinates, batch record, or annotation validation is supplied. Do not treat missing information as proof that QC failed or succeeded.

Reference tables were calculated independently with Python. The base-R verification script has not been executed in the authoring environment. It summarizes donor and pooled proportions and compares fictional marker patterns. It does not perform differential expression, differential abundance, clustering, or cell-type validation.
