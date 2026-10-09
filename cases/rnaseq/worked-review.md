# Worked review

| Claim | Judgment and priority | Evidence and action |
|---|---|---|
| 1 | Supported; none | Metadata and scenario identify three distinct samples per condition. Replication alone does not resolve confounding. |
| 2 | Incorrect; blocking | The design matrix has three columns but rank two. Batch B2 and treated condition generate identical indicator columns. Their separate effects cannot be estimated from this design. |
| 3 | Incorrect; blocking | Removing batch produces a fit-able condition-only design, but its comparison still combines condition and batch. A different formula cannot recover the absent comparison. |
| 4 | Incorrect; major | Four raw p-values are < 0.05. BH values are 0.006, 0.03, 0.06, 0.06, 0.24, and 0.8. Only G1 and G2 pass the arithmetic rule. |
| 5 | Supported as arithmetic; none | 2^1 = 2. This interprets an invented effect estimate, not an established biological effect. |
| 6 | Incorrect; major | 2^-1 = 0.5: the illustrative treated/control ratio is one half. A negative log ratio does not mean negative expression. |
| 7 | Not established; blocking | The result rows are invented, not model outputs, and the metadata cannot separate batch from condition. No biological differential-expression conclusion is supported. |

Corrected summary: The six-sample design completely confounds batch and condition, preventing estimation of their separate effects. In the separate six-row arithmetic exercise, only G1 and G2 meet the BH-adjusted < 0.05 rule, with illustrative fold changes of 2 and 0.5. These calculations do not establish treatment effects in the described study.

Request a design with condition represented across batches and appropriate biological replication, plus actual counts, provenance, and a documented analysis. Do not manufacture sample labels or assume batch correction restores missing information. For existing data, clearly state which combined contrasts are estimable and restrict interpretation accordingly.
