context("test-extract_signatures")

# Load mutation matrix
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
  package = "MutationalPatterns"
))
mut_mat <- mut_mat[1:10,]

# extract signatures
output <- extract_signatures(mut_mat, rank = 2, nrun = 1)

# Check it also works with indels
indel_counts <- readRDS(system.file("states/blood_indel_counts.rds",
  package = "MutationalPatterns"
))
indel_counts <- indel_counts[1:8,]

output_indel <- extract_signatures(indel_counts, rank = 2, nrun = 1)

# Check that the variational bayes method works
output_bayes <- extract_signatures(mut_mat, rank = 2, nrun = 1, nmf_type = "variational_bayes")

# Tests
test_that("Outputs a list", {
  expect_true(inherits(output, c("list")))
  expect_true(inherits(output_indel, c("list")))
  expect_true(inherits(output_bayes, c("list")))
})

test_that("Output contains signatures, contribution and reconstructed", {
  expect_identical(names(output), c("signatures", "contribution", "reconstructed"))
  expect_identical(names(output_indel), c("signatures", "contribution", "reconstructed"))
  expect_identical(names(output_bayes), c("signatures", "contribution", "reconstructed"))
})

test_that("Output elements are matrixes", {
  expect_true(inherits(output$signatures, "matrix"))
  expect_true(inherits(output$contribution, "matrix"))
  expect_true(inherits(output$reconstructed, "matrix"))
  expect_true(inherits(output_indel$signatures, "matrix"))
  expect_true(inherits(output_indel$contribution, "matrix"))
  expect_true(inherits(output_indel$reconstructed, "matrix"))
  expect_true(inherits(output_bayes$signatures, "matrix"))
  expect_true(inherits(output_bayes$contribution, "matrix"))
  expect_true(inherits(output_bayes$reconstructed, "matrix"))
})

test_that("Output elements have correct size", {
  expect_equal(dim(output$signatures), c(10, 2))
  expect_equal(dim(output$contribution), c(2, 9))
  expect_equal(dim(output$reconstructed), c(10, 9))
  expect_equal(dim(output_indel$signatures), c(8, 2))
  expect_equal(dim(output_indel$contribution), c(2, 3))
  expect_equal(dim(output_indel$reconstructed), c(8, 3))
  expect_equal(dim(output_bayes$signatures), c(10, 2))
  expect_equal(dim(output_bayes$contribution), c(2, 9))
  expect_equal(dim(output_bayes$reconstructed), c(10, 9))
})

# Test that an error is given when using a incorrect rank
test_that("An error is given when used with incorrect rank", {
  expect_error(
    {
      extract_signatures(mut_mat, rank = 2.5, nrun = 1)
    },
    "Rank should be a positive integer"
  )
})

# Test that an error is given when using a rank higher than the number of columns in the mut_mat
test_that("An error is given when rank is higher than the number of mut_mat columns", {
  expect_error(
    {
      extract_signatures(mut_mat, rank = 30, nrun = 1)
    },
    paste0(
      "The rank should be smaller than the number of ",
      "samples in the input matrix."
    )
  )
})
