context("test-context_potential_damage_analysis")


# Get contexts
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))

contexts <- rownames(mut_mat)[1:6]

# Load the corresponding reference genome.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Load transcription database
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Set gene ids
# TP53
gene_ids <- c(7157)

# Run the function
output <- context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids)

# Run the function with verbosity
output_verbose <- context_potential_damage_analysis(contexts, txdb, ref_genome, gene_ids, verbose = TRUE)

# Run unit tests
test_that("Output has correct class", {
  expect_true(inherits(output, "tbl_df"))
  expect_true(inherits(output_verbose, "tbl_df"))
})

test_that("Output has correct size", {
  expect_equal(dim(output), c(24, 5))
  expect_equal(dim(output_verbose), c(24, 5))
})

# Expected
expected <- readRDS(system.file("states/context_mismatches.rds",
  package = "MutationalPatterns"
))

test_that("Output is equal to expected", {
  expect_equal(output, expected)
})
