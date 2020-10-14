
#' plot_S_event_type
#'
#' Plot results of running stan
#' relative survival mixture cure model.
#'
#' @param stan_out Nested list of stan output
#'
#' @importFrom purrr map
#' @importFrom reshape2 melt
#' @importFrom rstan extract
#' @importFrom dplyr mutate
#'
#' @examples
#' load("data/stan_out.RData")
#'
plot_S_event_type <- function(stan_out) {

  ##TODO:
  # text cure fractions

  fit_stan <- list()
  S_stats <- list()
  S_pred <- NULL

  event_types <- names(stan_out)
  tx_names <- names(stan_out[[1]])

  for (i in event_types) {

    fit_stan[[i]] <- list()
    S_stats[[i]] <- list()

    for (j in tx_names) {

      fit_stan[[i]][[j]] <-
        rstan::extract(stan_out[[i]][[j]])

      # rearrange to time as rows
      S_dat <-
        list(
          t(fit_stan[[i]][[j]]$S_pred) %>%
            as_tibble() %>%
            mutate(month = 1:n(),
                   type = "S_pred"),
          t(fit_stan[[i]][[j]]$S0) %>%
            as_tibble() %>%
            mutate(month = 1:n(),
                   type = "S0"),
          t(fit_stan[[i]][[j]]$S_bg) %>%
            as_tibble() %>%
            mutate(month = 1:n(),
                   type = "S_bg"))

      # means and credible intervals
      S_stats[[i]][[j]] <-
        S_dat %>%
        do.call(rbind, .) %>%
        melt(id.vars = c("month", "type")) %>%
        group_by(month, type) %>%
        summarise(mean = mean(value),
                  lower = quantile(value, probs = 0.025),
                  upper = quantile(value, probs = 0.975))
    }
  }

  # unnest
  plot_dat <-
    S_stats %>%
    map(bind_rows, .id = "Tx") %>%
    bind_rows(.id = "event_type") %>%
    mutate(scenario = paste(event_type, Tx, sep = "_"))

  ggplot(plot_dat, aes(month, mean, group = type, colour = type)) +
    geom_line() +
    # facet_grid(. ~ scenario)
    facet_grid(event_type ~ Tx) +
    ylab("Survival") +
    geom_ribbon(aes(x = month, ymin = lower, ymax = upper, fill = type),
                linetype = 0,
                alpha = 0.2)
}

