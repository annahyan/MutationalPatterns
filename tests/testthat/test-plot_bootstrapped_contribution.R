context("test-plot_bootstrapped_contribution")

## contri_boots
contri_boots <- readRDS(system.file("states/bootstrapped_snv_refit.rds",
  package = "MutationalPatterns"
))

## Plot contexts

test_that("Output has correct class", {
  output <- plot_bootstrapped_contribution(contri_boots)
  expect_true(inherits(output, c("gg")))

  output <- plot_bootstrapped_contribution(contri_boots, mode = "relative")
  expect_true(inherits(output, c("gg")))

  output <- plot_bootstrapped_contribution(contri_boots, plot_type = "barplot")
  expect_true(inherits(output, c("gg")))
})
