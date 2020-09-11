context("test-genomic_distribution")

# Read vcfs
vcfs <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))

# Load genomic regions
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
  package = "MutationalPatterns"
))

promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
  package = "MutationalPatterns"
))

flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
  package = "MutationalPatterns"
))

# Combine regions and set seqlevelstyle
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
seqlevelsStyle(regions) <- "UCSC"

# Get the callable regions
surveyed_file <- system.file("extdata/callableloci-sample.bed",
  package = "MutationalPatterns"
)

library(rtracklayer)
surveyed <- rtracklayer::import(surveyed_file)
seqlevelsStyle(surveyed) <- "UCSC"

# Use the same callable loci for all samples.
surveyed_list <- rep(list(surveyed), 9)

## Calculate the number of observed and expected number of mutations in
## each genomic regions for each sample.
output <- genomic_distribution(vcfs, surveyed_list, regions)


test_that("Output has correct class", {
  expect_true(inherits(output, c("data.frame")))
})

test_that("Output has correct size", {
  expect_equal(dim(output), c(27, 8))
})

# Test that an error is given when the vcf_list and surveyed list are not the same size
test_that("An error is given when input sizes don't match", {
  expect_error(
    {
      genomic_distribution(vcfs, surveyed_list[1:8], regions)
    },
    "vcf_list and surveyed_list must have the same length"
  )
})

# Test that an error is given when regions_list names are not set.
regions_noname <- regions
names(regions_noname) <- NULL
test_that("An error is given when regions_list names are not set", {
  expect_error(
    {
      genomic_distribution(vcfs, surveyed_list, regions_noname)
    },
    "Please set the names of region_list using"
  )
})
