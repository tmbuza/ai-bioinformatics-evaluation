# ABE-004 evidence package

Original deterministic CDI teaching data. No real participants, organisms, sequences, or biological findings are represented. Taxonomic names are fictional genus labels; no species assignments are supplied.

Scenario: four participants each have one before and one after sample following a dietary change. There is no concurrent untreated comparison group, no randomization, no absolute microbial-load measurement, and no clinical outcome. Each participant has two observations. No exclusions or filtering have been applied to the supplied three-feature table. No extraction blanks, read-level QC, or laboratory validation are supplied.

## Data dictionary
- feature-counts.tsv: feature_id plus one column per sample; nonnegative integer read counts, not cell counts.
- sample-metadata.tsv: sample_id uniquely identifies an observation; participant_id identifies the repeated-measure unit; visit is before or after. Rows are intentionally reversed to require identifier-based alignment.
- taxonomy.tsv: one fictional genus label per feature; not a validated classification.

Read the simulated response before the worked review. Run scripts from the repository root. Reference tables were independently calculated with Python from these inputs; they are expected results, not evidence of an R run. The R script has not been executed in the authoring environment. No inferential test is supplied or claimed.
