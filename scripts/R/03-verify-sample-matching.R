# ABE-003: original synthetic teaching example; base R only.
# Run from the repository root. No biological inference is performed.
# Intentionally mismatched identifiers demonstrate a failed data check.
feature_ids <- c("S01", "S02", "S03", "S04")
metadata_ids <- c("S03", "S01", "S02", "S05")

out_dir <- file.path("results", "ch03")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

valid_ids <- function(x) {
  !anyNA(x) && all(nzchar(trimws(x))) && !anyDuplicated(x)
}
inputs_valid <- valid_ids(feature_ids) && valid_ids(metadata_ids)
same_ids <- inputs_valid && setequal(feature_ids, metadata_ids)
same_order <- inputs_valid && identical(feature_ids, metadata_ids)
checks <- data.frame(
  check = c("identifiers_present_and_unique", "same_identifier_set", "same_order"),
  passed = c(inputs_valid, same_ids, same_order)
)
all_ids <- union(feature_ids, metadata_ids)
membership <- data.frame(
  sample_id = all_ids,
  in_feature_table = all_ids %in% feature_ids,
  in_metadata = all_ids %in% metadata_ids
)
write.table(checks, file.path(out_dir, "checks.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
write.table(membership, file.path(out_dir, "sample-membership.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "session-info.txt"))
write.table(
  data.frame(source = c("feature_table", "metadata"),
             sample_ids = c(paste(feature_ids, collapse = ","),
                            paste(metadata_ids, collapse = ","))),
  file.path(out_dir, "input-identifiers.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE
)
print(checks, row.names = FALSE)
if (!same_ids) {
  message("REVIEW REQUIRED: sample identifiers do not match. No alignment performed.")
} else {
  alignment <- match(feature_ids, metadata_ids)
  stopifnot(identical(metadata_ids[alignment], feature_ids))
  message("Identifier sets match; metadata can be reordered by identifier.")
}
message("Evidence saved to ", out_dir)
# Normal completion means the checks ran, not that the data passed.
