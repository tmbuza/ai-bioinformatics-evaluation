# ABE-005: synthetic RNA-seq evaluation

This is an original teaching case with two deliberately separate evidence components.

## Study-design component
sample-metadata.tsv describes six distinct biological samples: three controls processed in B1 and three treated samples processed in B2. There are no repeated measurements. The intended model is ~ batch + condition. No sample in either condition was processed in the other batch.

## Arithmetic component
illustrative-results.tsv contains six invented gene IDs, log2 fold changes (treated/control), and p-values. These numbers were hand-specified for checking interpretation and multiple-testing arithmetic. They were NOT estimated from the six samples, produced by DESeq2, or derived from any count matrix. No count matrix is supplied. They cannot establish biological differential expression or validate the confounded model.

For this arithmetic exercise only, all six rows form the complete testing family; there are no missing p-values or filtered hypotheses. Use BH adjustment across all six rows and a strict adjusted p-value < 0.05 rule. Real workflows require their own documented testing and filtering rules.

## Files and reproducibility
sample_id is unique; condition and batch are categorical variables. gene_id is unique; log2FoldChange is on a base-2 scale; pvalue lies between zero and one.

reference/ contains expected values calculated independently in Python. The R script uses base R, not DESeq2, and has not been executed in the authoring environment. The script checks design rank and arithmetic; it does not fit a gene-expression model. Source documentation: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html (section “Model matrix not full rank”).
