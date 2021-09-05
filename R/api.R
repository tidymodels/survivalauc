#' Compute Survival AUC
#'
#' Implementing the alternative estimators of the time-dependent area
#' under the ROC curve, time-dependent sensitivity and specificity,
#' and time-dependent ROC curves.
#'
#' @param time Numeric, survival/censoring time.
#' @param status Numeric, censoring indicator taking value 1 if survival time is
#'   observed and 0 otherwise.
#' @param score Numeric, risk score obtained from a risk prediction model.
#' @param threshold Numeric, thresholds used in estimating sensitivity and
#' specificity.
#'
#' @details
#' `time`, `status`, and `score` must have the same length.
#'
#' The corresponding estimators are presented in equations (9)-(11) in
#' Chambless and Diao (Statistics in Medicine, 2006; 25: 3474-3486).
#'
#' @return NULL
#' @export
#'
#' @examples
#' compute_auc_components(
#'   time = sample_data$time,
#'   status = sample_data$status,
#'   score = sample_data$score,
#'   threshold = sample_data$threshold
#' )
compute_auc_components <- function(time, status, score, threshold) {

  unique_time <- sort(unique(time))

  .Call(survivalauc_compute_auc_components, time, status, score, threshold,
        unique_time)
}
