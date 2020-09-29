context("test-plot_profile_region")

# Read the long mutation matrix information:
input <- readRDS(system.file("states/mut_mat_longregions.rds",
  package = "MutationalPatterns"
))

## Plot the 96-profile of three samples
output <- plot_profile_region(input)
output_relative_sample_feature <- plot_profile_region(input, mode = "relative_sample_feature")
output_absolute <- plot_profile_region(input, mode = "absolute")
output_condensed <- plot_profile_region(input, condensed = TRUE)

test_that("Output has correct class", {
  expect_true(inherits(output, c("gg")))
  expect_true(inherits(output_relative_sample_feature, c("gg")))
  expect_true(inherits(output_absolute, c("gg")))
  expect_true(inherits(output_condensed, c("gg")))
})
