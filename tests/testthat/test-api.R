library(testthat)

is_increasing <- function(x) {
  all(diff(x) >= 0)
}

is_decreasing <- function(x) {
  all(diff(x) <= 0)
}

test_that("compute_auc_components works as expected", {
  expect_error(
    res <- compute_auc_components(
      time = sample_data$time,
      status = sample_data$status,
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    NA
  )

  expect_named(res, c("unique_times", "sensitivity", "specificity", "auc"))

  observed_times <- sort(unique(sample_data$time[sample_data$status == 1]))

  expect_equal(
    res$unique_times,
    observed_times
  )

  expect_equal(
    dim(res$sensitivity),
    c(length(observed_times), length(sample_data$threshold))
  )

  expect_equal(
    dim(res$specificity),
    c(length(observed_times), length(sample_data$threshold))
  )

  expect_equal(
    length(res$auc),
    length(observed_times)
  )

  expect_equal(is_increasing(res$auc), TRUE)

  expect_equal(
    apply(res$sensitivity, 1, is_decreasing),
    rep(TRUE, length(observed_times))
  )

  expect_equal(
    apply(res$sensitivity, 2, is_decreasing),
    rep(TRUE, length(sample_data$threshold))
  )

  expect_equal(
    apply(res$specificity, 1, is_increasing),
    rep(TRUE, length(observed_times))
  )

  expect_equal(
    apply(res$specificity, 2, is_increasing),
    rep(TRUE, length(sample_data$threshold))
  )
})

test_that("compute_auc_components errors with bad input", {

  expect_error(
    compute_auc_components(
      time = rep("A", 100),
      status = sample_data$status,
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    "`time` must be a double vector."
  )
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = rep("A", 100),
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    "`status` must only take the values 0 and 1."
  )
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = sample_data$status,
      score = rep("A", 100),
      threshold = sample_data$threshold
    ),
    "`score` must be a double vector."
  )
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = sample_data$status,
      score = sample_data$score,
      threshold = "A"
    ),
    "`threshold` must be a double vector."
  )
  expect_error(
    compute_auc_components(
      time = rep(1, 50),
      status = sample_data$status,
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    "`time`, `status`, and `score` must be the same length."
  )

  expect_error(
    compute_auc_components(
      time = c(sample_data$time[-1], NA),
      status = sample_data$status,
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    "`time` must not have any missing values."
  )
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = c(sample_data$status[-1], NA),
      score = sample_data$score,
      threshold = sample_data$threshold
    ),
    "`status` must not have any missing values."
  )
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = sample_data$status,
      score = c(sample_data$score[-1], NA),
      threshold = sample_data$threshold
    ),
    "`score` must not have any missing values."
  )
  expect_error(
    compute_auc_components(
      time = sample_data$time,
      status = sample_data$status,
      score = sample_data$score,
      threshold = c(sample_data$threshold[-1], NA)
    ),
    "`threshold` must not have any missing values."
  )
})
