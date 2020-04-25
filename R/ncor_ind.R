#' A Negative Control Outcome Regression for
#' Eliminating Unobserved Confounding in Time-series Studies
#'
#' @param data an optional data frame containing the variables in the model.
#' @param pre_outc_name the name of pre-exposure outcome
#' @param post_outc_name the name of post-exposure outcome
#' @param expo_namem the name of exposure
#' @param boots_no the number of bootstrap for IVW estimation
#' @param centre the number of reacher centre or time series fragments
#' @param method method of estimation
#'
#' @return  \code{causal} the casual effect estimation
#' @return  \code{lag} the lag causal effect estimation

#' @examples
#'
#' NN <- 10
#' data_total <- NULL
#' for(i in 1:NN){
#' c0 <- 0.5
#' dia <- 0.5
#' diff_c1 <- 0
#' r <- 0.5
#' c1 <- rnorm(1,0.5,1)
#' c2 <- rnorm(1,0.5,1)
#' N <- 100
#' Sigma <- matrix(c(1,r,r^2,r,1,r,r^2,r,1),3,3)
#' u <- MASS::mvrnorm(n=N, mu=rep(0,3), Sigma=Sigma)
#' colnames(u) <- c('u1','u2','u3')
#' u <- data.frame(u)
#' u1 <- u$u1
#' u2 <- u$u2
#' u3 <- u$u3
#'
#' eps_y1 <- rnorm(N,0,0.1)
#' eps_y3 <- rnorm(N,0,0.1)
#' eps_x2 <-  rnorm(N,0,0.1)
#' y1 <- c1*u1 + eps_y1
#' x2 <- c2*u2 + eps_x2
#' y3 <- c0*x2 + c1*u3+diff_c1*u3 + dia*y1 + eps_y3
#' data_cenre <- data.frame(x2,u1,u2,u3,y1,y3,i)
#' data_total <- rbind(data_total,data_cenre)
#'
#' }
#'
#' model <- ncor_ind (data=data_total, pre_outc_name='y1',expo_namem = 'x2',
#'                   post_outc_name='y3',centre='i',method='IVW',boot_no=1000)
#'
#' model
#'
#'

ncor_ind <- function(data = data, pre_outc_name = "y1", expo_namem = "x",
    post_outc_name = "y3", centre = "centre", method = c("regression",
        "IVW"), boot_no = NULL) {

    new_data <- data[, c(pre_outc_name, expo_namem, post_outc_name,
        centre)]
    colnames(new_data) <- c("y1", "x", "y3", "centre")
    centre_name <- unique(new_data$centre)

    coef1 <- NULL
    coef2 <- NULL
    se1 <- NULL


    for (i in 1:length(centre_name)) {

        data_s <- new_data[new_data$centre == centre_name[i],
            ]
        model1 <- lm(y3 ~ x, data = data_s)
        coef1[i] <- summary(model1)$coefficients[2, 1]
        se1[i] <- summary(model1)$coefficients[2, 2]

        model2 <- lm(y1 ~ x, data = data_s)
        coef2[i] <- summary(model2)$coefficients[2, 1]

    }

    if (method == "regression") {
        data_coef <- data.frame(coef1, coef2)
        RE_model <- lm(coef1 ~ coef2, data = data_coef)
        causal_CR <- summary(RE_model)$coefficients[1, 1]
        se_causal_CR <- summary(RE_model)$coefficients[1, 2]
        direct_CR <- summary(RE_model)$coefficients[2, 1] - 1

        P_causal_CR <- summary(RE_model)$coefficients[1, 4]
        se_direct_CR <- summary(RE_model)$coefficients[2, 2]
        p1_dir_CR <- pnorm(direct_CR, mean = 0, sd = se_direct_CR,
            lower.tail = TRUE, log.p = FALSE)
        p_dir_CR <- 2 * min(p1_dir_CR, 1 - p1_dir_CR)

        causal_est <- c(causal_CR, se_causal_CR, P_causal_CR)
        names(causal_est) <- c("Est_causal", "std.Err.causal",
            "P_causal")

        se_est <- c(direct_CR, se_direct_CR, p_dir_CR)
        names(se_est) <- c("Est_lag", "std.Err.lag", "P_lag")

        result <- list()

        result$causal <- causal_est
        result$lag <- se_est

        return(result)

    } else {

        weight <- NULL
        for (ppp in 1:length(centre_name)) {
            for (qqq in 1:length(centre_name)) {
                if (ppp == qqq) {
                  weight[ppp] <- (1 / se1[ppp]) * (1 / se1[qqq])
                }

            }
        }

        new_data_coef <- data.frame(coef1, coef2)

        RE_model_egg <- lm(coef1 ~ coef2, data = new_data_coef,
            weights = weight)
        result_egg <- summary(RE_model_egg)$coefficients
        rownames(result_egg) <- c("egg_causal", "egg_diret+1")
        causal_egg <- summary(RE_model_egg)$coefficients[1, 1]
        direct_egg <- summary(RE_model_egg)$coefficients[2, 1] -
            1

        casual_boot <- NULL
        direct_boot <- NULL
        NN <- length(centre_name)
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
        names(result_egg_direct) <- c("Est_lag", "std.Err.lag",
            "P_lag")

        result <- list()
        result$causal <- result_egg_causal
        result$lag <- result_egg_direct

        return(result)
    }
}
