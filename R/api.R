compute_auc_components <- function(time, status, score, threshold) {

  unique_time <- sort(unique(time))

  .Call(survivalauc_compute_auc_components, time, status, score, threshold, unique_time)
}
