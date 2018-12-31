library(MetaGxPancreas)

context("Checking loadPancreasDatasets")


test_that("ensure datasets and duplicates are properly loaded from the hub and package", {
  dataAndDuplicates = MetaGxPancreas::loadPancreasDatasets()
  seData = dataAndDuplicates$summarizedExperiments
  expect_equal(class(seData[[1]])[1], "RangedSummarizedExperiment")
})
