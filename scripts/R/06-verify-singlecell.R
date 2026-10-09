# ABE-006: descriptive checks only; base R. Run from repository root.
in_dir <- "cases/singlecell/data"
out_dir <- "results/ch06"
read_tsv <- function(n) read.delim(file.path(in_dir,n), stringsAsFactors=FALSE)
d <- read_tsv("donors.tsv"); x <- read_tsv("cells.tsv")
m <- read_tsv("marker-evidence.tsv"); ref <- read_tsv("teaching-marker-reference.tsv")
valid_ids <- function(z) !anyNA(z) && all(nzchar(trimws(z))) && !anyDuplicated(z)
stopifnot(valid_ids(d$donor_id), valid_ids(x$cell_id), !anyNA(x$donor_id),
          setequal(x$donor_id,d$donor_id), !anyNA(x$cluster),
          setequal(x$cluster,c("C1","C2")), !anyNA(d$condition),
          setequal(d$condition,c("control","treated")))
d$n_cells <- as.integer(table(factor(x$donor_id, levels=d$donor_id)))
d$n_C1 <- as.integer(table(factor(x$donor_id[x$cluster=="C1"],levels=d$donor_id)))
stopifnot(all(d$n_cells>0))
d$C1_percent <- 100*d$n_C1/d$n_cells
summary <- do.call(rbind,lapply(sort(unique(d$condition)),function(g){
 z <- d[d$condition==g,]
 data.frame(condition=g,n_donors=nrow(z),n_cells=sum(z$n_cells),n_C1=sum(z$n_C1),
            pooled_C1_percent=100*sum(z$n_C1)/sum(z$n_cells),
            mean_donor_C1_percent=mean(z$C1_percent))
}))
stopifnot(all(m$cluster=="C1"), valid_ids(m$marker), is.numeric(m$percent_detected),
          all(is.finite(m$percent_detected)), all(m$percent_detected>=0 & m$percent_detected<=100),
          valid_ids(ref$candidate_type), !anyNA(ref$expected_high_markers))
high <- m$marker[m$percent_detected>=50]
ref$n_expected_high <- vapply(strsplit(ref$expected_high_markers,","),length,integer(1))
ref$n_observed_high <- vapply(strsplit(ref$expected_high_markers,","),function(z){
 stopifnot(all(z %in% m$marker)); sum(z %in% high)
},integer(1))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(z,n) write.table(z,file.path(out_dir,n),sep="\t",row.names=FALSE,quote=FALSE)
write_tsv(d,"donor-summary.tsv");write_tsv(summary,"condition-summary.tsv")
write_tsv(ref,"teaching-marker-check.tsv")
capture.output(sessionInfo(),file=file.path(out_dir,"session-info.txt"))
files <- list.files(in_dir,full.names=TRUE)
write_tsv(data.frame(file=files,md5=unname(tools::md5sum(files))),"input-checksums.tsv")
print(d);print(summary);print(ref)
message("Descriptive checks only; no annotation validation or significance test performed.")
