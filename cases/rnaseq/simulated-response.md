# Simulated AI response

Written for teaching; includes sound statements, errors, and unsupported claims.

1. There are three biological samples per condition.
2. Including batch in ~ batch + condition separates the treatment effect from batch, despite their complete overlap.
3. Dropping batch from the formula fixes the design and recovers an unconfounded treatment effect.
4. Four genes have raw p-values below 0.05, so all four pass the specified BH-adjusted threshold.
5. G1 has log2FoldChange = 1, corresponding to a twofold treated/control ratio.
6. G2 has log2FoldChange = -1, meaning expression is negative.
7. These results demonstrate that treatment changes G1 and G2 independently of batch.
