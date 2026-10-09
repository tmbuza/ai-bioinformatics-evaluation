# ABE-005: design and arithmetic verification only; no DESeq2 model is fitted.
# Run from repository root; outputs are overwritten on rerun.
in_dir <- "cases/rnaseq/data"
out_dir <- "results/ch05"
read_tsv <- function(name) read.delim(file.path(in_dir, name), stringsAsFactors = FALSE)
m <- read_tsv("sample-metadata.tsv")
x <- read_tsv("illustrative-results.tsv")
valid_ids <- function(z) !anyNA(z) && all(nzchar(trimws(z))) && !anyDuplicated(z)
stopifnot(valid_ids(m$sample_id), valid_ids(x$gene_id),
          !anyNA(m$condition), !anyNA(m$batch),
          setequal(m$condition, c("control", "treated")),
          setequal(m$batch, c("B1", "B2")),
          is.numeric(x$pvalue), is.numeric(x$log2FoldChange),
          all(is.finite(x$pvalue)), all(x$pvalue >= 0 & x$pvalue <= 1),
          all(is.finite(x$log2FoldChange)))
m$condition <- factor(m$condition, levels = c("control", "treated"))
m$batch <- factor(m$batch, levels = c("B1", "B2"))
design <- model.matrix(~ batch + condition, m)
design_check <- data.frame(n_samples = nrow(m), n_columns = ncol(design),
                           rank = qr(design)$rank,
                           full_rank = qr(design)$rank == ncol(design))
x$BH_adjusted <- p.adjust(x$pvalue, method = "BH")
x$passes_BH_0.05 <- x$BH_adjusted < 0.05
x$fold_change <- 2^x$log2FoldChange
stopifnot(all(is.finite(x$fold_change)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(z, name) write.table(z, file.path(out_dir, name), sep = "\t",
                                          row.names = FALSE, quote = FALSE)
write_tsv(design_check, "design-check.tsv")
write_tsv(data.frame(sample_id = m$sample_id, design, check.names = FALSE), "design-matrix.tsv")
write_tsv(as.data.frame(table(m$batch, m$condition)), "batch-condition-table.tsv")
write_tsv(x, "result-checks.tsv")
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
files <- list.files(in_dir, full.names = TRUE)
write_tsv(data.frame(file = files, md5 = unname(tools::md5sum(files))), "input-checksums.tsv")
print(design_check)
print(x)
if (!design_check$full_rank) message("BLOCKING: separate batch and condition effects are not estimable.")
message("Arithmetic checked on invented results only; no biological inference performed.")
