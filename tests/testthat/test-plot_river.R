context("test-plot_river")


# Get input data
mut_mat <- readRDS(system.file("states/mut_mat_data.rds",
                               package = "MutationalPatterns"
))

mut_mat_extended <- readRDS(system.file("states/mut_mat_data_extended.rds",
                                        package = "MutationalPatterns"
))


## Create riverplot of profile
output <- plot_river(mut_mat)

## Create condensed riverplot of extended profile
output_extended <- plot_river(mut_mat_extended, condensed = TRUE)


test_that("Output has correct class", {
    expect_true(inherits(output, c("gg")))
    expect_true(inherits(output_extended, c("gg")))
})
