# Run from the repository root. Base R only; original synthetic teaching inputs.
options(stringsAsFactors = FALSE)
read_tsv <- function(p) read.delim(p, check.names = FALSE, na.strings = "NA")
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
out <- "results/ch08"
files <- "cases/gwas/data/carrier-counts.tsv"
stopifnot(file.exists(files))
d <- read_tsv(files)
stopifnot(identical(names(d), c("stratum", "status", "carriers", "noncarriers")), !anyNA(d))
stopifnot(setequal(d$stratum, c("A", "B")), setequal(d$status, c("case", "control")))
stopifnot(nrow(d) == 4, !anyDuplicated(paste(d$stratum, d$status)))
stopifnot(is.numeric(d$carriers), is.numeric(d$noncarriers))
n <- as.matrix(d[, c("carriers", "noncarriers")])
stopifnot(all(is.finite(n)), all(n > 0), all(n == floor(n)))
checks <- do.call(rbind, lapply(c("A", "B", "pooled"), function(g) {
 z <- if (g == "pooled") d else d[d$stratum == g, ]
 a <- colSums(z[z$status == "case", c("carriers", "noncarriers"), drop = FALSE])
 b <- colSums(z[z$status == "control", c("carriers", "noncarriers"), drop = FALSE])
 data.frame(comparison = g, case_n = sum(a), control_n = sum(b), case_carrier_proportion = unname(a[1]/sum(a)), control_carrier_proportion = unname(b[1]/sum(b)), odds_ratio = unname(a[1]*b[2]/(a[2]*b[1])))
}))
d$total <- d$carriers + d$noncarriers
d$within_status_proportion <- d$total / ave(d$total, d$status, FUN = sum)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tsv(checks, "association-checks.tsv")
write_tsv(d[, c("stratum", "status", "total", "within_status_proportion")], "composition.tsv")
write_tsv(data.frame(file = files, md5 = unname(tools::md5sum(files))), "input-checksums.tsv")
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
message("Verification outputs written to ", out)
