# Run from the repository root. Base R only; original synthetic teaching inputs.
options(stringsAsFactors = FALSE)
read_tsv <- function(p) read.delim(p, check.names = FALSE, na.strings = "NA")
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
out <- "results/ch07"
files <- file.path("cases/proteomics/data", c("sample-metadata.tsv", "log2-intensities.tsv"))
stopifnot(all(file.exists(files)))
m <- read_tsv(files[1]); d <- read_tsv(files[2])
stopifnot(identical(names(m), c("sample_id", "condition")), names(d)[1] == "protein_id")
stopifnot(!anyNA(m), !anyDuplicated(m$sample_id), all(nzchar(m$sample_id)))
stopifnot(!anyNA(d$protein_id), !anyDuplicated(d$protein_id), all(nzchar(d$protein_id)))
stopifnot(setequal(names(d)[-1], m$sample_id), nrow(m) == 6)
stopifnot(setequal(m$condition, c("control", "treated")), all(table(m$condition) == 3))
x <- as.matrix(d[, m$sample_id, drop = FALSE])
stopifnot(is.numeric(x), all(is.finite(x) | is.na(x)))
rownames(x) <- d$protein_id
missingness <- do.call(rbind, lapply(c("control", "treated"), function(g) {
 z <- x[, m$condition == g, drop = FALSE]
 data.frame(protein_id = rownames(x), condition = g, observed = rowSums(!is.na(z)), missing = rowSums(is.na(z)))
}))
sensitivity <- do.call(rbind, lapply(c("observed_only", "fill_18", "fill_22"), function(method) {
 z <- x
 if (method != "observed_only") z[is.na(z)] <- if (method == "fill_18") 18 else 22
 a <- rowMeans(z[, m$condition == "control", drop = FALSE], na.rm = TRUE)
 b <- rowMeans(z[, m$condition == "treated", drop = FALSE], na.rm = TRUE)
 a[is.nan(a)] <- NA_real_; b[is.nan(b)] <- NA_real_
 data.frame(protein_id = rownames(x), method = method, control_mean_log2 = a, treated_mean_log2 = b, difference = b-a, geometric_mean_ratio = 2^(b-a), row.names = NULL)
}))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv(missingness, "missingness.tsv")
write_tsv(sensitivity, "sensitivity.tsv")
write_tsv(data.frame(file = files, md5 = unname(tools::md5sum(files))), "input-checksums.tsv")
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
message("Verification outputs written to ", out)
