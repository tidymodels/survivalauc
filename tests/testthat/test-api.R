library(testthat)

test_that("compute_auc_components doesn't crash", {
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = sample_data$status,
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    NA
  )
})
