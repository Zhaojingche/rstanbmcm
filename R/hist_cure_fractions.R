
#' hist_cure_fractions
#'
#' @param stan_out Nested list
#'
#' @return
#' @export
#'
hist_cure_fractions <- function(stan_out) {

  plot_dat <- NULL

  for (i in names(stan_out)) {
    for (j in names(stan_out[[1]])) {

      plot_dat <-
        rbind(plot_dat,
              data.frame(i, j, stan_out[[i]]))
    }
  }

  names(plot_dat) <- c("event_type", "Tx", "curefraction")

  ggplot(plot_dat, aes(x = curefraction)) +
    facet_grid(event_type ~ Tx) +
    geom_histogram()
}

