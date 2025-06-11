#' Extract information from input data in preparation for imputation.
#'
#' @description
#' Extracts information from a dataset for use in a TGIFA model.
#'
#'
#' @param data The data to be processed. If not a matrix, will be transformed in to a matrix.
#' @param scaled `TRUE`/`FALSE.` Should the input data be scaled? Defaults to `FALSE.`
#' @param coding The coding used to indicate a data entry is missing. Defaults to `NA.`
#' @param pareto_scale `TRUE`/`FALSE.` Should the input data be Pareto scaled? Defaults to `FALSE.`
#' @param log_trans `TRUE`/`FALSE.` Should the input data be log transformed? Defaults to `FALSE.`
#'
#' @return A list containing three entries.
#'
#' * data: Contains the processed data in matrix format with missing entries coded as `NA`.
#'
#' * missingness_pattern: Contains the missingness pattern of the data for use in the TGIFA model.
#'
#' * divisors: Contains a vector of the divisors used to scale the input data, if applicable. If unscaled, the value will be `NULL`.
#'
#' @importFrom stats sd
#'
#' @noRd
TGIFA_prep <- function(data, scaled = FALSE, coding = NA, pareto_scale = FALSE,
                       log_trans = FALSE) {

  # convert to matrix as required
  if (!(is.matrix(data))){
    data <- as.matrix(data)
  }

  # how missing data is encoded
  data[data == coding] <- NA

  # centre or scale as required
  if (scaled){
    data <- scale(data, center = FALSE, scale = scaled)
  }

  divisors <- NULL

  if (pareto_scale) {

    divisors <- sqrt(apply(data, 2, sd, na.rm = TRUE))

    data <- sweep(data, MARGIN = 2, STATS = divisors, FUN = "/")
  }

  if (log_trans) {
    data <- log(data)
  }

  # generate missingness pattern
  missingness_pattern <- 1 - is.na(data)

  return(list("data" = data, "missingness_pattern" = missingness_pattern, "divisors" = divisors))
}

#' Mode
#'
#' @param x A vector.
#'
#' @return The mode of x. If x is multi-modal, the first mode (by index in x) is returned.
#'
#' @noRd
Mode <- function(x) {

  options <- unique(x)

  return(options[which.max(tabulate(match(x, options)))])

}

#' Results formatting for TGIFA results
#'
#' @param TGIFA_res Draws from the MCMC chain of the TGIFA model
#' @param burn The number of MCMC iterations to be discarded in an initial burn.
#' @param thin The level of thinning to take place in the MCMC chain. E.g., `thin = 5` will retain every fifth draw, post-burn.
#' @param cred_level The desired credible interval to report for imputed values. E.g., `cred_level = 0.95` will return a 95 percent credible interval.
#'
#' @return Returns a list containing two entries.
#'
#' * imputed_dataset: The returned dataset in matrix format with the imputed values inserted.
#'
#' * imputation_info: A dataframe containing information on imputed values and associated uncertainty for each missing entry in the original dataset.
#'
#' @importFrom stats quantile
#'
#' @noRd
TGIFA_res_format <- function(TGIFA_res, burn, thin, cred_level = 0.95) {

  # extract input data from results object
  original_data <- TGIFA_res$input_data

  n <- dim(original_data)[1]
  p <- dim(original_data)[2]

  n.iters <- length((TGIFA_res$store_data[-(1:(burn/thin))]))

  # extract missingness information from results object
  missingness <- TGIFA_res$input_missingness

  miss_vec <- which(missingness == 0)
  miss_ind <- which(missingness == 0, arr.ind = TRUE)

  if (dim(miss_ind)[1] == 0) {

    print("There are no missing values in this dataset. Did you input the correct missingness coding?")

    return(list("imputed_dataset" = NA,
                "imputation_info" = NA))

  } else {

    # Z designation
    store_Z <- sapply(TGIFA_res$store_Z[-(1:(burn/thin))], FUN = matrix)

    # imputed data
    data_mat <- array(NA, dim = c(length((TGIFA_res$store_data[-(1:(burn/thin))])[[1]]),
                                  n.iters))

    for (i in 1:n.iters) {

      data_mat[ , i] <- (TGIFA_res$store_data[-(1:(burn/thin))])[[i]]

    }

    # imputed values
    just_imputed <- data_mat

    # modal designations
    mode_Z <- apply(store_Z, 1, Mode)

    MNAR_designations <- which(mode_Z == 1)
    MAR_designations <- which(mode_Z == 0)

    designations_vec <- rep("MNAR", dim(store_Z)[1])
    designations_vec[MAR_designations] <- "MAR"

    MNAR_post_mean <- c()
    MNAR_cred_upper <- c()
    MNAR_cred_lower <- c()

    cred_probs <- c((1 - cred_level) / 2,
                    1 - (1 - cred_level) / 2)

    for (mnar_ind in MNAR_designations) {

      relevant_draws <- just_imputed[mnar_ind, which(store_Z[mnar_ind, ] == mode_Z[mnar_ind])]

      post_mean <- median(relevant_draws)
      cred_int <- quantile(relevant_draws, probs = cred_probs)
      MNAR_post_mean <- c(MNAR_post_mean, post_mean)
      MNAR_cred_lower <- c(MNAR_cred_lower, cred_int[1])
      MNAR_cred_upper <- c(MNAR_cred_upper, cred_int[2])

    }

    MAR_post_mean <- c()
    MAR_cred_upper <- c()
    MAR_cred_lower <- c()

    for (mar_ind in MAR_designations) {

      relevant_draws <- just_imputed[mar_ind, which(store_Z[mar_ind, ] == mode_Z[mar_ind])]

      post_mean <- median(relevant_draws)
      cred_int <- quantile(relevant_draws, probs = cred_probs)
      MAR_post_mean <- c(MAR_post_mean, post_mean)
      MAR_cred_lower <- c(MAR_cred_lower, cred_int[1])
      MAR_cred_upper <- c(MAR_cred_upper, cred_int[2])

    }


    TGIFA_imp <- original_data
    TGIFA_imp[miss_vec][MNAR_designations] <- MNAR_post_mean
    TGIFA_imp[miss_vec][MAR_designations] <- MAR_post_mean

    cred_int_upper <- matrix(NA, nrow = n, ncol = p)
    cred_int_upper[miss_vec][MNAR_designations] <- MNAR_cred_upper
    cred_int_upper[miss_vec][MAR_designations] <- MAR_cred_upper

    cred_int_lower <- matrix(NA, nrow = n, ncol = p)
    cred_int_lower[miss_vec][MNAR_designations] <- MNAR_cred_lower
    cred_int_lower[miss_vec][MAR_designations] <- MAR_cred_lower

    # get prop of each point's Zs that are MNAR
    z_props <- apply(store_Z, 1, sum) / dim(store_Z)[2]

    # flip designated MAR to 1 - MNAR certainty
    z_props[MAR_designations] <- 1 - z_props[MAR_designations]

    # for ease of indexing
    z_uncertainty_TGIFA <- matrix(NA, nrow = n, ncol = p)
    z_uncertainty_TGIFA[miss_vec][MNAR_designations] <- (1 - z_props[MNAR_designations])
    z_uncertainty_TGIFA[miss_vec][MAR_designations] <- (1 - z_props[MAR_designations])

    # create dataframe with uncertainty information
    imputation_info <- data.frame(entry_row = miss_ind[ , 1],
                                  entry_col = miss_ind[ , 2],
                                  imputed_val = TGIFA_imp[miss_vec],
                                  cred_int_upper = cred_int_upper[miss_vec],
                                  cred_int_lower = cred_int_lower[miss_vec],
                                  miss_mech = designations_vec,
                                  miss_mech_unc = z_uncertainty_TGIFA[miss_vec]
    )


    return(list("imputed_dataset" = TGIFA_imp,
                "imputation_info" = imputation_info))

  }

}

# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~
# functions used
# ~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~

# evaluate p tilde
d_p_tilde <- function(Y, given_mu, lambda, eta, Sigma_inv, logged = TRUE) {

  d_val <- 0

  for (i in 1:nrow(Y)) {
    val <- matrix(Y[i, ]) - given_mu - mat.mult(lambda, matrix(eta[i, ]))
    d_val <- d_val + (-0.5 * mat.mult(mat.mult(t(val), Sigma_inv), val))
  }

  if (!logged) {
    d_val <- exp(d_val)
  }

  return(d_val)
}

# evaluate p tilde by col
d_p_tilde_j <- function(y_j, mu_j, lambda_j, eta_i, sigma_inv_j, logged = TRUE) {

  d_val <- 0

  for (i in 1:nrow(y_j)) {

    d_val <- d_val + (dnorm(sqrt(sigma_inv_j) * (y_j[i, ] -
                                                   mu_j -
                                                   mat.mult(Rfast::transpose(lambda_j), eta_i))))

  }

  if (!logged) {

    d_val <- exp(d_val)

  }

  return(d_val)

}

d_p_numerator_eta <- function(eta_i, y_i, mu, lambda, Sigma_inv, logged = TRUE) {

  val <- y_i - (mu + mat.mult(lambda, eta_i))

  d_p <- -0.5 * (mat.mult(transpose(eta_i), eta_i) + (mat.mult(mat.mult(transpose(val), Sigma_inv), val)))

  if (!logged) {

    d_p <- exp(d_p)

  }

  return(d_p)

}
