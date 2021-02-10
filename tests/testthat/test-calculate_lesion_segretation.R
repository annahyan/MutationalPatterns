context("test-calculate_lesion_segregation")

# To test mut_matrix, we need to load the reference genome first.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

# Load GRangesList
grl <- readRDS(system.file("states/read_vcfs_as_granges_output.rds",
  package = "MutationalPatterns"
))
sample_names <- c(
  "colon1", "colon2", "colon3",
  "intestine1", "intestine2", "intestine3",
  "liver1", "liver2", "liver3"
)


## Only look at some samples, to reduce the runtime
grl <- grl[1:2]
sample_names <- sample_names[1:2]

# Perform lesion segregation calculations
output <- calculate_lesion_segregation(grl, sample_names)
output_per_type <- calculate_lesion_segregation(grl, sample_names,
  split_by_type = TRUE, ref_genome = ref_genome
)
output_wald <- calculate_lesion_segregation(grl,
  sample_names,
  test = "wald-wolfowitz"
)

## Calculate lesion segregation using the rl20.
chromosomes <- paste0("chr", c(1:22, "X"))
output_rl20 <- calculate_lesion_segregation(grl,
  sample_names,
  test = "rl20",
  ref_genome = ref_genome,
  chromosomes = chromosomes
)

# Run tests
test_that("Output has correct class", {
  expect_true(inherits(output, c("tbl_df")))
  expect_true(inherits(output_per_type, c("tbl_df")))
  expect_true(inherits(output_wald, c("tbl_df")))
  expect_true(inherits(output_rl20, c("tbl_df")))
})

test_that("Output has correct dimensions", {
  expect_equal(dim(output), c(2, 8))
  expect_equal(dim(output_per_type), c(2, 8))
  expect_equal(dim(output_wald), c(2, 5))
  expect_equal(dim(output_rl20), c(2, 5))
})

expected <- readRDS(system.file("states/lesion_segregation.rds",
  package = "MutationalPatterns"
))
test_that("transforms correctly", {
  expect_equal(output, expected)
})

# Test that An error is thrown when the arguments are incorrectly combined.
test_that("An error is thrown when the arguments are incorrectly combined", {
  expect_error(
    {
      calculate_lesion_segregation(grl, sample_names[1])
    },
    "The vcf_list and the sample_names should be equally long."
  )
  expect_error(
    {
      calculate_lesion_segregation(grl, sample_names, test = "wald-wolfowitz", split_by_type = TRUE)
    },
    "The 'split_by_type' argument can only be used with the binomial test"
  )
  expect_error(
    {
      calculate_lesion_segregation(grl, sample_names, split_by_type = TRUE)
    },
    "The ref_genome needs to be set when"
  )
  expect_error(
    {
      calculate_lesion_segregation(grl, sample_names, test = "rl20", chromosomes = "chr1")
    },
    "The ref_genome needs to be set when"
  )
  expect_error(
    {
      calculate_lesion_segregation(grl, sample_names, test = "rl20", ref_genome = ref_genome)
    },
    "The chromosomes need to be set when"
  )
})
