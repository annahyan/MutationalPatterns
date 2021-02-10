context("test-get_mut_type")

# Get a grl with variants.
grl <- readRDS(system.file("states/blood_grl.rds",
  package = "MutationalPatterns"
))

# Only use two samples to reduce runtime
grl <- grl[1:2]

## Get a specific mutation type.
snv_grl <- get_mut_type(grl, "snv")
indel_grl <- get_mut_type(grl, "indel")
dbs_grl <- get_mut_type(grl, "dbs")
mbs_grl <- get_mut_type(grl, "mbs")
gr_singlesample <- get_mut_type(grl[[1]], type = "dbs")
empty_gr <- get_mut_type(grl[[1]][0], type = "dbs")
gr_nodbs <- get_mut_type(grl[[1]][1:20], type = "dbs")

# Change names of indel_grl, to make them prettier.
remove_names_gr <- function(gr) {
  names(gr) <- seq_along(gr)
  return(gr)
}
indel_grl <- purrr::map(as.list(indel_grl), remove_names_gr) %>%
  GRangesList()

expected_indel_grl <- readRDS(system.file("states/blood_grl_indel.rds",
  package = "MutationalPatterns"
))[1:2]


test_that("Output has correct class", {
  expect_true(inherits(snv_grl, c("GRanges", "CompressedGRangesList")))
  expect_true(inherits(indel_grl, c("GRanges", "CompressedGRangesList")))
  expect_true(inherits(dbs_grl, c("GRanges", "CompressedGRangesList")))
  expect_true(inherits(mbs_grl, c("GRanges", "CompressedGRangesList")))
  expect_true(inherits(gr_singlesample, c("GRanges")))
  expect_true(inherits(empty_gr, c("GRanges")))
  expect_true(inherits(gr_nodbs, c("GRanges")))
})

test_that("Output is equal to expected", {
  expect_equal(indel_grl, expected_indel_grl)
})

test_that("Empty gr is returned when a mut type is not present", {
  expect_equal(length(empty_gr), 0)
})

test_that("Empty gr as input results in a empty output gr", {
  expect_equal(length(gr_nodbs), 0)
})
