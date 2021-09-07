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
#' `time`, `status`, and `score` must have the same length and be numeric.
#'
#' @return  A list with 4 elements:
#'
#' \describe{
#'   \item{unique_times}{distinct observed failure times in ascending order.}
#'   \item{sensitivity}{`m` by `k` matrix, the `(i,j)`th element corresponds
#'                      to the estimated sensitivity at time point
#'                      `unique_times[i]` with threshold `thresold[j]`}
#'   \item{specificity}{`m` by `k` matrix, the `(i,j)`th element corresponds
#'                      to the estimated specificity at time point
#'                      `unique_times[i]` with threshold `thresold[j]`}
#'   \item{auc}{the estimated AUC at time point `unique_times`}
#' }
#'
#' @references
#' Chambless, L. E. and G. Diao (2006). Estimation of time-dependent area under
#' the ROC curve for long-term risk prediction. Statistics in Medicine 25,
#' 3474--3486.
#'
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
