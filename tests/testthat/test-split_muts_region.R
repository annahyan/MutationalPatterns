context("test-split_muts_region")

# Read in genomic regions
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
  package = "MutationalPatterns"
))
promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
  package = "MutationalPatterns"
))
flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
  package = "MutationalPatterns"
))

# Combine the regions into a single GRangesList
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
seqlevelsStyle(regions) <- "UCSC"

# Read in some variants.
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

# Only use two samples to reduce runtime
grl <- grl[1:2]


# Run function
output <- split_muts_region(grl, regions)
output_single_gr <- split_muts_region(grl[[1]], regions)
output_single_region <- split_muts_region(grl, regions[[1]])
output_noother <- split_muts_region(grl, regions, include_other = FALSE)


test_that("Output has correct class", {
  expect_true(inherits(output, "CompressedGRangesList"))
  expect_true(inherits(output_single_gr, "CompressedGRangesList"))
  expect_true(inherits(output_single_region, "CompressedGRangesList"))
  expect_true(inherits(output_noother, "CompressedGRangesList"))
})

expected_length <- function(grl, regions) {
  exp_length <- (length(regions) + 1) * length(grl) # nr. samples * nr. regions. +1 is for the variants in 'other'
  return(exp_length)
}

test_that("Output GRangesList has correct length", {
  expect_equal(length(output), expected_length(grl, regions))
  expect_equal(length(output_single_gr), expected_length(grl[1], regions))
  expect_equal(length(output_single_region), expected_length(grl, regions[1]))
  expect_equal(length(output_noother), length(regions) * length(grl))
})

expected <- readRDS(system.file("states/grl_split_region.rds",
  package = "MutationalPatterns"
))[1:8]
test_that("Output transforms correctly", {
  expect_equal(output, expected)
})
