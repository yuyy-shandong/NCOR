#' A Negative Control Outcome Regression for
#' Eliminating Unobserved Confounding in Time-series Studies
#'
#' @param data an optional data frame containing the variables in the model.
#' @param coef_y3_x the coef of exposure on post-exposure outcome
#' @param coef_y1_x the coef of exposure on pre-exposure outcome
#' @param se the standard error of coef_y3_x
#' @param boots_no the number of bootstrap for estimation

#'
#' @return  \code{causal} the casual effect estimation
#' @return  \code{lag} the lag causal effect estimation



ncor_summary <- function(data = data, coef_y3_x = "coef1", coef_y1_x = "coef2",
    se = "se", boot_no = 1000) {

    new_data <- data[, c(coef_y3_x, coef_y1_x, se)]
    colnames(new_data) <- c("coef1", "coef2", "se1")
    new_data <- data.frame(new_data)

    weight <- NULL
    for (ppp in 1:length(new_data$coef1)) {
        for (qqq in 1:length(new_data$coef1)) {
            if (ppp == qqq) {
                weight[ppp] <- (1 / new_data$se1[ppp]) * (1 / new_data$se1[qqq])
            }

        }
    }

    new_data_coef <- data.frame(new_data$coef1, new_data$coef2)

    RE_model_egg <- lm(coef1 ~ coef2, data = new_data_coef, weights = weight)
    result_egg <- summary(RE_model_egg)$coefficients
    rownames(result_egg) <- c("egg_causal", "egg_diret+1")
    causal_egg <- summary(RE_model_egg)$coefficients[1, 1]
    direct_egg <- summary(RE_model_egg)$coefficients[2, 1] - 1

    casual_boot <- NULL
    direct_boot <- NULL
    NN <- length(new_data$coef1)
    for (boot in 1:boot_no) {
        num_sample <- sample(1:NN, NN, replace = T)
        weigth_www <- weight[num_sample]
        boot_sample <- new_data_coef[num_sample, ]
        boot_model_egg <- lm(coef1 ~ coef2, data = boot_sample,
            weights = weigth_www)
        casual_boot[boot] <- summary(boot_model_egg)$coefficients[1,
            1]
        direct_boot[boot] <- summary(boot_model_egg)$coefficients[2,
            1] - 1
    }
    se_causal_egg <- sd(casual_boot, na.rm = T)
    se_direct_egg <- sd(direct_boot, na.rm = T)

    p1_cau_egg <- pnorm(causal_egg, mean = 0, sd = se_causal_egg,
        lower.tail = TRUE, log.p = FALSE)
    p_cau_egg <- 2 * min(p1_cau_egg, 1 - p1_cau_egg)

    p1_dir_egg <- pnorm(direct_egg, mean = 0, sd = se_direct_egg,
        lower.tail = TRUE, log.p = FALSE)
    p_dir_egg <- 2 * min(p1_dir_egg, 1 - p1_dir_egg)

    result_egg_causal <- c(causal_egg, se_causal_egg, p_cau_egg)
    names(result_egg_causal) <- c("Est_causal", "std.Err.causal",
        "P_causal")
    result_egg_direct <- c(direct_egg, se_direct_egg, p_dir_egg)
    names(result_egg_direct) <- c("Est_lag", "std.Err.lag", "P_lag")

    result <- list()
    result$causal <- result_egg_causal
    result$lag <- result_egg_direct

    return(result)

}
