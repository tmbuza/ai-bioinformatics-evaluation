# ABE-004: base-R descriptive verification of synthetic teaching data.
# Run from repository root. Outputs are overwritten on rerun.
in_dir <- "cases/microbiome/data"
out_dir <- "results/ch04"
read_tsv <- function(name) read.delim(file.path(in_dir, name), check.names = FALSE,
                                     stringsAsFactors = FALSE)
x <- read_tsv("feature-counts.tsv")
m <- read_tsv("sample-metadata.tsv")
tax <- read_tsv("taxonomy.tsv")
valid_ids <- function(z) !anyNA(z) && all(nzchar(trimws(z))) && !anyDuplicated(z)
stopifnot(valid_ids(x$feature_id), valid_ids(m$sample_id), valid_ids(tax$feature_id))
ids <- names(x)[-1]
stopifnot(valid_ids(ids), setequal(ids, m$sample_id), setequal(x$feature_id, tax$feature_id))
cts <- as.matrix(x[-1]); rownames(cts) <- x$feature_id
stopifnot(is.numeric(cts), all(is.finite(cts)), all(cts >= 0), all(cts == floor(cts)))
m <- m[match(ids, m$sample_id), , drop = FALSE]
stopifnot(identical(m$sample_id, ids), !anyNA(m$participant_id),
          all(nzchar(m$participant_id)), !anyNA(m$visit),
          all(m$visit %in% c("before", "after")))
pair_table <- table(m$participant_id, factor(m$visit, levels = c("before", "after")))
stopifnot(all(pair_table == 1))
totals <- colSums(cts)
stopifnot(all(totals > 0), "F1" %in% rownames(cts))
summary <- data.frame(m, total_reads = totals, F1_reads = as.numeric(cts["F1", ]),
                      F1_percent = 100 * as.numeric(cts["F1", ]) / totals,
                      row.names = NULL)
before <- summary[summary$visit == "before", ]
after <- summary[summary$visit == "after", ]
before <- before[order(before$participant_id), ]
after <- after[match(before$participant_id, after$participant_id), ]
paired <- data.frame(participant_id = before$participant_id,
                     before_percent = before$F1_percent, after_percent = after$F1_percent,
                     change_percentage_points = after$F1_percent - before$F1_percent)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(z, name) write.table(z, file.path(out_dir, name), sep = "\t",
                                          row.names = FALSE, quote = FALSE)
write_tsv(summary, "sample-summary.tsv")
write_tsv(paired, "paired-changes.tsv")
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
write_tsv(data.frame(file = list.files(in_dir, full.names = TRUE),
                    md5 = unname(tools::md5sum(list.files(in_dir, full.names = TRUE)))),
          "input-checksums.tsv")
print(aggregate(cbind(F1_reads, F1_percent) ~ visit, summary, mean))
print(paired)
message("Descriptive checks completed; no significance test or causal analysis performed.")
