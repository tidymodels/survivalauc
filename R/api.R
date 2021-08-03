compute_auc_components <- function(time, status, score, threshold) {
  .Call(survivalauc_compute_auc_components, time, status, score, threshold)
}
