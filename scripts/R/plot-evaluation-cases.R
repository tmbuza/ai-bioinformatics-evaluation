# Original CDI synthetic evaluation cases. Run from the repository root:
# Rscript scripts/R/plot-evaluation-cases.R
# Base R only. Regenerates five PNGs; no statistical tests are performed.
read_tsv <- function(p) read.delim(p, check.names = FALSE)
out <- "results/figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
teal <- "#036281"; gold <- "#f7c546"; ink <- "#0f172a"
figure <- function(name, draw) {
  png(file.path(out, name), width = 1800, height = 1100, res = 180)
  on.exit(dev.off())
  par(mar = c(5, 5, 4, 2), col.axis = ink, col.lab = ink, col.main = ink,
      fg = ink, family = "sans", las = 1)
  draw()
}
# 04: match by identifiers, never by row position.
m <- read_tsv("cases/microbiome/data/sample-metadata.tsv")
f <- read_tsv("cases/microbiome/data/feature-counts.tsv")
stopifnot(!anyDuplicated(m$sample_id), !anyDuplicated(f$feature_id),
          setequal(m$sample_id, names(f)[-1]), "F1" %in% f$feature_id)
x <- as.matrix(f[, m$sample_id]); stopifnot(is.numeric(x), all(x >= 0), all(colSums(x) > 0))
m$percent <- 100 * x[f$feature_id == "F1", ] / colSums(x)
ids <- sort(unique(m$participant_id)); colors <- c(teal, "#577d34", "#98622b", "#7463a0")
figure("04-microbiome-paired.png", function() {
  plot(NA, xlim = c(0.9, 2.5), ylim = c(0, 50), xaxt = "n", xlab = "Visit",
       ylab = "F1 relative abundance (%)", main = "Microbiome: follow each participant")
  axis(1, 1:2, c("Before", "After"))
  for (i in seq_along(ids)) {
    z <- m[m$participant_id == ids[i], ]; stopifnot(nrow(z) == 2, setequal(z$visit, c("before", "after")))
    y <- z$percent[match(c("before", "after"), z$visit)]
    lines(1:2, y, type = "b", pch = 19, lwd = 2, col = colors[(i-1) %% 4+1])
    text(2.06, y[2], ids[i], pos = 4, col = colors[(i-1) %% 4+1], cex = 0.85)
  }
  mtext("Synthetic counts; relative abundance does not measure absolute load", side = 3, cex = 0.8)
})
# 05: a design plot, not differential-expression evidence.
m <- read_tsv("cases/rnaseq/data/sample-metadata.tsv")
stopifnot(!anyDuplicated(m$sample_id), setequal(m$condition, c("control", "treated")), setequal(m$batch, c("B1", "B2")))
tab <- table(factor(m$condition, levels = c("control", "treated")), factor(m$batch, levels = c("B1", "B2")))
figure("05-rnaseq-design.png", function() {
  plot(NA, xlim = c(0.5, 2.5), ylim = c(0.5, 2.5), xaxt = "n", yaxt = "n",
       xlab = "Batch", ylab = "", main = "RNA-seq: condition and batch are confounded")
  axis(1, 1:2, c("B1", "B2")); axis(2, 1:2, c("Control", "Treated"))
  for (i in 1:2) for (j in 1:2) {
    rect(j-0.48, i-0.48, j+0.48, i+0.48, col = if (tab[i,j] > 0) teal else "#eef2f6", border = "white")
    text(j, i, paste(tab[i,j], "samples"), col = if (tab[i,j] > 0) "white" else ink, cex = 1.3)
  }
  mtext("No condition comparison is available within either batch", side = 3, cex = 0.8)
})
# 06: distinguish pooled cells from donor-level summaries.
d <- read_tsv("cases/singlecell/data/donors.tsv"); c <- read_tsv("cases/singlecell/data/cells.tsv")
stopifnot(!anyDuplicated(d$donor_id), !anyDuplicated(c$cell_id), all(c$donor_id %in% d$donor_id))
d$n <- vapply(d$donor_id, function(id) sum(c$donor_id == id), numeric(1))
d$p <- vapply(d$donor_id, function(id) mean(c$cluster[c$donor_id == id] == "C1") * 100, numeric(1))
stopifnot(all(d$n > 0), setequal(d$condition, c("control", "treated")))
figure("06-singlecell-weighting.png", function() {
  par(mfrow = c(1,2), mar = c(5,5,4,1))
  bars <- barplot(d$p, names.arg = d$donor_id, ylim = c(0,110), col = ifelse(d$condition == "control", teal, gold),
                  ylab = "C1 cells (%)", main = "Each donor", border = NA)
  text(bars, d$p+6, paste0("n=", d$n), cex = 0.8)
  legend("top", c("Control", "Treated"), fill = c(teal,gold), bty = "n", horiz = TRUE, cex = 0.8)
  groups <- c("control", "treated")
  pooled <- sapply(groups, function(g) weighted.mean(d$p[d$condition == g], d$n[d$condition == g]))
  means <- sapply(groups, function(g) mean(d$p[d$condition == g]))
  b <- barplot(rbind(pooled, means), beside = TRUE, names.arg = c("Control", "Treated"), ylim = c(0,110),
               col = c(teal,gold), border = NA, main = "Choice of summary", ylab = "C1 cells (%)")
  text(b, rbind(pooled,means)+5, as.character(rbind(pooled,means)), cex = 0.85)
  legend("top", c("Pooled cells", "Donor mean"), fill = c(teal,gold), bty = "n", cex = 0.8)
})
# 07: display missingness and descriptive sensitivity, not significance.
m <- read_tsv("cases/proteomics/data/sample-metadata.tsv")
p <- read_tsv("cases/proteomics/data/log2-intensities.tsv")
stopifnot(!anyDuplicated(m$sample_id), !anyDuplicated(p$protein_id), setequal(m$sample_id, names(p)[-1]))
x <- as.matrix(p[, m$sample_id]); stopifnot(is.numeric(x), "P1" %in% p$protein_id)
v <- x[p$protein_id == "P1", ]
diffs <- sapply(c(NA,18,22), function(fill) {
 z <- v; if (!is.na(fill)) z[is.na(z)] <- fill
 mean(z[m$condition == "treated"], na.rm = TRUE) - mean(z[m$condition == "control"], na.rm = TRUE)
})
figure("07-proteomics-missingness.png", function() {
  par(mfrow = c(1,2), mar = c(6,5,4,1))
  plot(NA, xlim = c(0.5,ncol(x)+0.5), ylim = c(0.5,nrow(x)+0.5), xaxt = "n", yaxt = "n",
       xlab = "Sample", ylab = "", main = "Observed or missing?")
  axis(1, seq_len(ncol(x)), colnames(x)); axis(2, seq_len(nrow(x)), p$protein_id)
  for (i in seq_len(nrow(x))) for (j in seq_len(ncol(x))) {
    rect(j-0.48,i-0.48,j+0.48,i+0.48,col = if(is.na(x[i,j])) gold else teal,border="white")
    text(j,i,if(is.na(x[i,j])) "NA" else format(x[i,j]),col=if(is.na(x[i,j])) ink else "white",cex=0.8)
  }
  mtext("Numbers are log2 intensities", side = 3, cex = 0.7)
  plot(1:3, diffs, pch = 19, col = teal, cex = 1.4, ylim = c(-4.5,3), xaxt = "n",
       xlab = "Missing-value assumption", ylab = "Treated - control (mean log2)", main = "P1: sensitivity of direction")
  axis(1,1:3,c("Observed\nonly","Fill 18","Fill 22")); abline(h=0,lty=2,col="gray50")
  text(1:3,diffs+0.45,format(round(diffs,3)),cex=0.9)
})
# 08: odds ratios on a log axis; point estimates only.
g <- read_tsv("cases/gwas/data/carrier-counts.tsv")
stopifnot(!anyNA(g), all(g$carriers > 0), all(g$noncarriers > 0),
          nrow(g) == 4, !anyDuplicated(paste(g$stratum,g$status)))
labels <- c("A","B","pooled")
ors <- sapply(labels,function(s) {
 z <- if(s=="pooled") g else g[g$stratum==s,]
 a <- colSums(z[z$status=="case",c("carriers","noncarriers"),drop=FALSE])
 b <- colSums(z[z$status=="control",c("carriers","noncarriers"),drop=FALSE])
 unname(a[1]*b[2]/(a[2]*b[1]))
})
figure("08-gwas-stratification.png", function() {
  par(mar = c(5,7,4,2))
  plot(ors,3:1,log="x",xlim=c(0.5,8),ylim=c(0.5,3.5),yaxt="n",xaxt="n",pch=19,cex=1.5,
       col=c(teal,teal,"#aa7100"),xlab="Carrier odds ratio (log scale)",ylab="",main="GWAS: pooling changes the association")
  axis(1,c(0.5,1,2,4,8),c("0.5","1","2","4","8"));axis(2,3:1,c("Stratum A","Stratum B","Pooled"))
  abline(v=1,lty=2,col="gray50");text(ors,3:1+0.25,format(round(ors,3)))
  mtext("Descriptive point estimates; no confidence intervals or significance tests",side=3,cex=0.75)
})
capture.output(sessionInfo(), file = file.path(out, "plot-session-info.txt"))
message("Saved five case figures to ", out)
