#' Truncated Gaussian Infinite Factor Analysis Imputation Model
#'
#' @description
#' Applies the TGIFA imputation method to the input dataset.
#'
#' @param input_data The dataset with missing value requiring imputation in matrix format.
#' @param coding The coding used to indicate a data entry is missing. Defaults to `NA`.
#' @param n.iters The number of MCMC iterations desired. Defaults to 10000.
#' @param k.star Parameter \eqn{k^{*}} in the TGIFA model. Defaults to 5.
#' @param verbose Is a readout desired? Defaults to `TRUE`.
#' @param burn The number of MCMC iterations to be discarded in an initial burn. Defaults to 5000.
#' @param thin The level of thinning to take place in the MCMC chain. Defaults to 5, meaning only every fifth draw will be kept.
#' @param kappa_1 Parameter \eqn{\kappa_1} in the TGIFA model. Defaults to 3.
#' @param kappa_2 Parameter \eqn{\kappa_2} in the TGIFA model. Defaults to 2.
#' @param a_sigma Parameter \eqn{a_{\sigma}} in the TGIFA model. Defaults to 1.
#' @param b_sigma Parameter \eqn{b_{\sigma}} in the TGIFA model. Defaults to 0.25.
#' @param a_1 Parameter \eqn{a_1} in the TGIFA model. Defaults to 2.1.
#' @param a_2 Parameter \eqn{a_2} in the TGIFA model. Defaults to 3.1.
#' @param prop_mult A list containing scalar multipliers for the covariance of parameters with Metropolis Hastings or Exchange algorithm update steps.
#' @param return_chain Should resulting draws from the MCMC chain be returned? Defaults to `FALSE`.
#' @param cred_level The desired credible interval to report for imputed values. E.g., `cred_level = 0.95` will return a 95 percent credible interval.
#' @param var_mult The desired proportion of the variance of the input data to be assigned to initial values of the idiosyncratic errors.
#'
#' @return Returns a list containing entries:
#'
#' * imputed_dataset: The imputed dataset in matrix format.
#'
#' * imputation_info: A dataframe containing information on imputed values and associated uncertainty for each missing entry in the original dataset.
#'
#' * chain: Returned only if return_chain = `TRUE`. Provides all draws from the MCMC chain and further information.
#'
#'     * store_data: MCMC draws of the imputed data entries.
#'
#'     * store_alpha: MCMC draws of the \eqn{\alpha} parameter.
#'
#'     * store_sigma_inv: MCMC draws of the \eqn{\sigma^{-2}} entries.
#'
#'     * store_lambda: MCMC draws of the loadings matrix.
#'
#'     * store_eta: MCMC draws of the latent factor scores.
#'
#'     * store_mu: MCMC draws of the mean vector.
#'
#'     * store_k_t: MCMC draws of the effective factor number (unchanging).
#'
#'     * store_phi: MCMC draws of the local shrinkage parameters.
#'
#'     * store_delta: MCMC draws of the multiplicative gamma process parameters.
#'
#'     * store_b_sigma: MCMC draws of the \eqn{b_{\sigma}} parameter (unchanging).
#'
#'     * store_Z: MCMC draws of the missingness designation indicator for each imputed entry.
#'
#'     * input_data: The dataset input into the tIFA model.
#'
#'     * input_missingness: The missingness pattern input into the tIFA model.
#'
#'     * store_MH_mu: Exchange algorithm acceptance for the mean parameter.
#'
#'     * store_MH_lambda: Exchange algorithm acceptance for the loadings matrix..
#'
#'     * store_MH_eta: Metropolis-Hastings acceptance for the latent factor scores
#'
#'     * store_MH_sigma_inv: Exchange algorithm acceptance for the \eqn{\sigma^{-2}} entries.
#'
#' @importFrom Rfast transpose mat.mult colsums
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom truncdist rtrunc
#' @importFrom stats rgamma rbeta rbinom runif prcomp cov
#' @importFrom truncnorm ptruncnorm rtruncnorm
#' @importFrom Rcpp cppFunction
#' @importFrom softImpute softImpute
#'
#' @export
#'
#' @examples
#' # generate example data
#' set.seed(1)
#' example_data <- matrix(abs(rnorm(200)), nrow = 5)
#'
#' # add missingness to example data coded as 0.001
#' example_data[4, 2] <- 0.001
#' example_data[2, 18] <- 0.001
#'
#' # run TGIFA model
#' # short chain for example
#' res <- TGIFA_model(input_data = example_data, coding = 0.001, n.iters = 100, k.star = 3,
#'                   burn = 50, thin = 5)
#'
#' # checkout imputed dataset
#' res$imputed_dataset
#'
#' # checkout further information on imputed entries
#' res$imputation_info
#'


TGIFA_model <- function(input_data, coding = NA, n.iters = 10000, k.star = 5,
                        verbose = TRUE, burn = 5000, thin = 5,
                        kappa_1 = 3L, kappa_2 = 2L,
                        a_sigma = 1L, b_sigma = 0.25, a_1 = 2.1, a_2 = 3.1,
                        prop_mult = list("mu_prop_mult" = 0.05,
                                         "eta_prop_mult" = 3L,
                                         "lambda_prop_mult" = 3L,
                                         "sigma_inv_prop_mult" = 1.5),
                        return_chain = FALSE,
                        cred_level = 0.95,
                        var_mult = 1.0) {

  Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends = "RcppArmadillo")

  # code data as expected
  # extract missingness information
  prepped_data <- TGIFA_prep(input_data, coding = coding)

  missingness <- prepped_data$missingness_pattern

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # initialisation
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # data dimensions
  n <- nrow(prepped_data$data)
  p <- ncol(prepped_data$data)

  # seeding function for rcpp functions
  get_seed <- function() {

    sample.int(.Machine$integer.max, 1)

  }

  # truncation point
  trunc.point <- min(prepped_data$data, na.rm = TRUE)

  # indices of missing data
  # col 1 contains row indices in data
  # col 2 contains col indices in data
  # so miss_ind[1, ] are the indices of the first missing point in the data
  miss_ind <- which(missingness == 0, arr.ind = TRUE)
  miss_vec <- which(missingness == 0)

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # def starting values
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # initialise imputation with svd imputation
  start_svd <- softImpute(prepped_data$data, rank.max = k.star, type = "svd")

  start_imp <- start_svd$u %*% diag(start_svd$d) %*% t(start_svd$v)

  initial_dataset <- prepped_data$data
  initial_dataset[is.na(initial_dataset)] <- start_imp[is.na(initial_dataset)]

  # data for TGIFA
  data_working <- initial_dataset

  # initialise remaining parameters

  # use pca to initialise factors and loadings
  pca_res <- prcomp(data_working, scale. = FALSE, center = FALSE)

  lambda <- pca_res$rotation[ , 1:k.star]  # p x k.star
  rownames(lambda) <- NULL
  colnames(lambda) <- NULL

  eta <- matrix(NA, nrow = n, ncol = k.star)

  for (i in 1:n) {

    eta[i, ] <- Rfast::rmvnorm(n = 1,
                               mu = rep(0, k.star),
                               sigma = diag(k.star),
                               seed = get_seed())

  }  # n x k

  mu <- matrix(apply(prepped_data$data, 2, mean, na.rm = TRUE)) - apply(Rfast::mat.mult(lambda, t(eta)), 1, mean)
  mu_tilde <- mu - 1
  mu_varphi <- 1 / (0.05 * apply(data_working, 2, mean, na.rm = TRUE))
  mu_varphi[(1:p)[-unique(miss_ind[,2])]] <- 1L

  alpha <- runif(1, 0, length(miss_vec) / (n * p))

  sigma_inv <- matrix(1 / (apply(prepped_data$data, 2, var, na.rm = TRUE) * var_mult))

  Sigma_inv <- diag(sigma_inv[ , 1])

  phi <- matrix(rgamma(n = p * k.star, shape = kappa_1, rate = kappa_2),  # p x k.star
                nrow = p)

  delta <- matrix(rep(0, k.star), nrow = k.star)  # k.star x 1
  delta[1, 1] <- rgamma(n = 1, shape = a_1, rate = 1)
  delta[-1, 1] <- truncdist::rtrunc(spec = "gamma",
                                    n = k.star - 1,
                                    shape = a_2,
                                    rate = 1,
                                    a = 1)

  tau <- matrix(cumprod(delta))  # k.star x 1

  rm(prepped_data)

  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~
  # set up storage
  # ~~~~~~~~~~~~~~~~~~~~
  # ~~~~~~~~~~~~~~~~~~~~

  # imputed dataset
  store_data <- list()

  # alpha - prob of MAR
  store_alpha <- list()

  # sigma
  store_sigma_inv <- list()

  # lambda - factor loadings
  store_lambda <- list()

  # eta
  store_eta <- list()

  # mu - mean
  store_mu <- list()

  # effective factors
  store_k_t <- list()

  # phi
  store_phi <- list()

  # delta
  store_delta <- list()

  # Z
  store_Z <- list()

  # MH iters
  store_MH_mu <- list()
  store_MH_lambda <- list()
  store_MH_eta <- list()
  store_MH_sigma_inv <- list()

  for (m in seq_len(n.iters)) {

    if (verbose) {

      # create a readout
      statement <- paste("TGIFA process running. Now on iteration ", m, " of ", n.iters)
      print(statement)

    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # lambda update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # objects for proposal distribution
    eta_t_eta <- mat.mult(Rfast::transpose(eta) , eta)  # k.star x k.star

    MH_lambda <- matrix(NA, nrow = p, ncol = 1)

    for (j in 1:p) {

      if (k.star == 1) {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau), nrow = 1))  # k.star x k.star

      } else {

        lambda_var <- armaInv(sigma_inv[j] * eta_t_eta + diag(c(phi[j, ] * tau)))  # k.star x k.star

      }

      lambda_mean <- mat.mult(lambda_var, mat.mult(
        Rfast::transpose(eta), as.matrix(data_working[ , j] - mu[j])) * sigma_inv[j])  # k.star x 1


      lambda_j_prime <- Rfast::rmvnorm(n = 1, mu = as.vector(lambda_mean),
                                       sigma = prop_mult$lambda_prop_mult * lambda_var,
                                       seed = get_seed())

      lambda_j_prime <- matrix(lambda_j_prime)

      # generate auxiliary y_jprime given lambda_j_prime
      y_j_prime <- matrix(NA, nrow = n, ncol = 1)

      for (i in 1:n) {

        y_j_prime[i, ] <- truncnorm::rtruncnorm(n = 1, a = 0, b = Inf,
                                                mean = mu[j, 1] + mat.mult(Rfast::transpose(lambda_j_prime),
                                                                           matrix(eta[i, ])),
                                                sd = sigma_inv[j, 1])

      }

      # calc acceptance prob
      # use logs for numerical stability
      A_lambda_j <- Rfast::dmvnorm(lambda[j, ], mu = lambda_mean,
                                   sigma = prop_mult$lambda_prop_mult * lambda_var, logged = TRUE) -
        Rfast::dmvnorm(Rfast::transpose(lambda_j_prime), mu = lambda_mean,
                       sigma = prop_mult$lambda_prop_mult * lambda_var, logged = TRUE) +
        Rfast::dmvnorm(Rfast::transpose(lambda_j_prime), mu = rep(0, k.star),
                       sigma = diag(c(phi[j, ] * tau)), logged = TRUE) -
        Rfast::dmvnorm(Rfast::transpose(matrix(lambda[j, ])), mu = rep(0, k.star),
                       sigma = diag(c(phi[j, ] * tau)), logged = TRUE) +
        d_p_tilde_j(y_j = matrix(data_working[ , j]), mu_j = mu[j, 1],
                    lambda_j = lambda_j_prime, eta_i = matrix(eta[i, ]),
                    sigma_inv_j = sigma_inv[j, 1], logged = TRUE) -
        d_p_tilde_j(matrix(data_working[ , j]), mu[j, 1], matrix(lambda[j, ]), matrix(eta[i, ]),
                    sigma_inv[j, 1], logged = TRUE) +
        d_p_tilde_j(y_j_prime, mu[j, 1], matrix(lambda[j, ]), matrix(eta[i, ]), sigma_inv[j, 1], logged = TRUE) -
        d_p_tilde_j(y_j_prime, mu[j, 1], lambda_j_prime, matrix(eta[i, ]), sigma_inv[j, 1], logged = TRUE)

      # accept-reject step
      if (runif(1) < exp(A_lambda_j)) {
        # update mu value
        lambda[j, ] <- lambda_j_prime
        MH_lambda[j, 1] <- 1  # update MH acceptance count
      } else {
        MH_lambda[j, 1] <- 0  # update MH acceptance count
      }

    }

    store_MH_lambda[[m]] <- MH_lambda

    # lambda is a p x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # eta update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    MH_eta <- matrix(NA, nrow = n, ncol = 1)

    # objects for and drawing proposed new eta
    lambda_t_Sigma <- mat.mult(Rfast::transpose(lambda), Sigma_inv)

    eta_var <- armaInv(diag(k.star) + mat.mult(lambda_t_Sigma, lambda))  # k.star x k.star

    for (i in 1:n) {

      eta_mean <- mat.mult(eta_var, mat.mult(lambda_t_Sigma, as.matrix(data_working[i, ] - mu)))  # k.star x 1

      eta_i_prime <-  matrix(Rfast::rmvnorm(n = 1, mu = as.vector(eta_mean),
                                            sigma = prop_mult$eta_prop_mult * eta_var,
                                            seed = get_seed()))

      A_eta_i <- Rfast::dmvnorm(eta[i, ], mu = eta_mean,
                                sigma = prop_mult$eta_prop_mult * eta_var, logged = TRUE) -
        Rfast::dmvnorm(Rfast::transpose(eta_i_prime), mu = eta_mean,
                       sigma = prop_mult$eta_prop_mult * eta_var, logged = TRUE) +
        d_p_numerator_eta(eta_i_prime, matrix(data_working[i, ]), mu, lambda, Sigma_inv, logged = TRUE) -
        d_p_numerator_eta(matrix(eta[i, ]), matrix(data_working[i, ]), mu, lambda, Sigma_inv, logged = TRUE)

      # accept-reject step
      if (runif(1) < exp(A_eta_i)) {
        # update mu value
        eta[i, ] <- eta_i_prime
        MH_eta[i, 1] <- 1  # update MH acceptance count
      } else {
        MH_eta[i, 1] <- 0  # update MH acceptance count
      }


    }

    store_MH_eta[[m]] <- MH_eta

    # eta is a n x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # sigma_inv update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # object for calculations
    lambda_eta_t <- mat.mult(lambda, Rfast::transpose(eta))  # p x n

    sigma_rate <- b_sigma + 0.5 * colsums((sweep(data_working, 2, mu, FUN = "-") -
                                             Rfast::transpose(lambda_eta_t))^2)

    MH_sigma_inv <- matrix(NA, nrow = p, ncol = 1)

    # propose new sigma_inv value
    sigma_inv_prime <- matrix(rgamma(p, shape = a_sigma + n/2,
                                     rate = prop_mult$sigma_inv_prop_mult * sigma_rate))

    for (j in 1:p) {

      # generate auxiliary y_j_prime given sigma_inv_j_prime
      y_j_prime <- matrix(NA, nrow = n, ncol = 1)

      for (i in 1:n) {

        y_j_prime[i, ] <- truncnorm::rtruncnorm(n = 1, a = 0, b = Inf,
                                                mean = mu[j, 1] + mat.mult(Rfast::transpose(matrix(lambda[j, ])),
                                                                           matrix(eta[i, ])),
                                                sd = sigma_inv_prime[j, 1])

      }

      A_sigma_j <- dgamma(sigma_inv[j, ], shape = a_sigma + n/2,
                          rate = prop_mult$sigma_inv_prop_mult * sigma_rate[j], log = TRUE) -
        dgamma(sigma_inv_prime[j, ], shape = a_sigma + n/2,
               rate = prop_mult$sigma_inv_prop_mult * sigma_rate[j], log = TRUE) +
        dgamma(sigma_inv_prime[j, ], shape = a_sigma,
               rate = b_sigma, log = TRUE) -
        dgamma(sigma_inv[j, ], shape = a_sigma,
               rate = b_sigma, log = TRUE) +
        d_p_tilde_j(y_j = matrix(data_working[ , j]), mu_j = mu[j, 1],
                    lambda_j = matrix(lambda[j, ]), eta_i = matrix(eta[i, ]),
                    sigma_inv_j = sigma_inv_prime[j, ], logged = TRUE) -
        d_p_tilde_j(y_j = matrix(data_working[ , j]), mu_j = mu[j, 1],
                    lambda_j = matrix(lambda[j, ]), eta_i = matrix(eta[i, ]),
                    sigma_inv_j = sigma_inv[j, ], logged = TRUE) +
        d_p_tilde_j(y_j = y_j_prime, mu_j = mu[j, 1],
                    lambda_j = matrix(lambda[j, ]), eta_i = matrix(eta[i, ]),
                    sigma_inv_j = sigma_inv[j, ], logged = TRUE) -
        d_p_tilde_j(y_j = y_j_prime, mu_j = mu[j, 1],
                    lambda_j = matrix(lambda[j, ]), eta_i = matrix(eta[i, ]),
                    sigma_inv_j = sigma_inv_prime[j, ], logged = TRUE)


      # accept-reject step
      if (runif(1) < exp(A_sigma_j)) {
        # update sigma_inv value
        sigma_inv[j, ] <- sigma_inv_prime[j, ]
        MH_sigma_inv[j, 1] <- 1  # update MH acceptance count
      } else {
        MH_sigma_inv[j, 1] <- 0  # update MH acceptance count
      }

    }

    store_MH_sigma_inv[[m]] <- MH_sigma_inv

    Sigma_inv <- diag(sigma_inv[ , 1])

    # sigma_inv is a 100 x 1 matrix
    # Sigma_inv is a 100 x 100 diag matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # mu update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # objects for proposal dist
    mu_temp1 <- mu_varphi * diag(p)

    mu_var <- diag(1 / (mu_varphi + n * sigma_inv[ , 1]))

    mu_temp2 <- mat.mult(Sigma_inv , matrix(colsums(data_working - mat.mult(eta, t(lambda)))))

    mu_mean <- mat.mult(mu_var, mu_temp2 + mat.mult(mu_temp1, matrix(mu_tilde)))

    # generate proposed mu_prime
    mu_prime <- Rfast::rmvnorm(n = 1, mu = as.vector(mu_mean),
                               sigma = prop_mult$mu_prop_mult * mu_var,
                               seed = get_seed())

    mu_prime <- matrix(mu_prime, nrow = p) # p x 1

    # tidy up
    rm(mu_temp1)
    rm(mu_temp2)

    # generate auxiliary Y_prime given mu_prime
    Y_prime <- matrix(NA, nrow = n, ncol = p)

    for (j in 1:p) {

      for (i in 1:n) {

        Y_prime[i, j] <- truncnorm::rtruncnorm(n = 1, a = 0, b = Inf,
                                               mean = mu_prime[j, 1] + mat.mult(Rfast::transpose(matrix(lambda[j, ])),
                                                                                matrix(eta[i, ])),
                                               sd = sigma_inv[j, 1])

      }
    }

    A_mu <- d_p_tilde(data_working, mu_prime, lambda, eta, Sigma_inv, logged = TRUE) -
      d_p_tilde(data_working, mu, lambda, eta, Sigma_inv, logged = TRUE) +
      mvtnorm::dmvnorm(as.vector(mu_prime), mean = mu_tilde, sigma = (1/mu_varphi) * diag(p), log = TRUE) -
      mvtnorm::dmvnorm(as.vector(mu), mean = mu_tilde, sigma = (1/mu_varphi) * diag(p), log = TRUE) +
      d_p_tilde(Y_prime, mu, lambda, eta, Sigma_inv, logged = TRUE) -
      d_p_tilde(Y_prime, mu_prime, lambda, eta, Sigma_inv, logged = TRUE)

    # accept-reject step
    if (runif(1) < exp(A_mu)) {
      # update mu value
      mu <- mu_prime
      store_MH_mu[m] <- 1  # update MH acceptance count
    } else {
      store_MH_mu[m] <- 0  # update MH acceptance count
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # phi update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for (h in 1:k.star) {

      phi[ , h] <- rgamma(p, shape = rep(0.5 + kappa_1, p), rate = kappa_2 + 0.5 * tau[h] * lambda[ , h]^2)

    }

    # phi is a p x k.star matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # delta and tau update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    delta[1, 1] <- rgamma(1,
                          shape = a_1 + 0.5 * k.star * p,
                          rate = 1 + 0.5 * sum((tau / delta[1, 1]) * colsums(lambda^2 * phi)))

    tau <- matrix(cumprod(delta))

    if (k.star > 1) {

      for (h in 2:k.star) {

        delta_shape <- a_2 + 0.5 * p * (k.star - h + 1)
        delta_rate <- 1 + 0.5 * sum(((tau / delta[h, 1]) * colsums(lambda^2 * phi))[h:k.star])

        delta[h, 1] <- truncdist::rtrunc(spec = "gamma",
                                         n = 1,
                                         shape = delta_shape,
                                         rate = delta_rate,
                                         a = 1)

        tau <- matrix(cumprod(delta))

      }

    }


    # delta is a k.star x 1 matrix
    # tau is a k.star x 1 matrix

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # alpha update
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    alpha_shape1 <- sum((1 - missingness) * (data_working >= trunc.point)) + 1

    alpha_shape2 <- sum((missingness) * (data_working >= trunc.point)) + 1

    alpha <- rbeta(1, shape1 = alpha_shape1, shape2 = alpha_shape2)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # imputation
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Z <- rep(NA, length(miss_vec))

    if (dim(miss_ind)[1] != 0) {

      impute_mean <- sweep(mat.mult(lambda, Rfast::transpose(eta)), 1, mu, "+")  # p x n
      impute_var <- 1 / sigma_inv  # p x p

      for (point_ind in 1:dim(miss_ind)[1]) {

        # ptruncnorm calcs lower.tail
        lower_prob <- truncnorm::ptruncnorm(trunc.point,
                                            a = 0,
                                            mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                            sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        # as truncation of dist is 0, Inf
        # and lower_prob is 0, trunc.point
        higher_prob <- 1 - lower_prob

        zij_prob <- (lower_prob) / (lower_prob + alpha * higher_prob)

        zij <- rbinom(1, 1, zij_prob)

        Z[point_ind] <- zij

        if (zij == 0) {

          # impute as MAR
          data_working[miss_ind[point_ind, ][1], miss_ind[point_ind, ][2]] <-
            truncnorm::rtruncnorm(1, a = trunc.point,
                                  b = Inf,
                                  mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                  sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        } else {

          # impute as MNAR
          data_working[miss_ind[point_ind, ][1], miss_ind[point_ind, ][2]] <-
            truncnorm::rtruncnorm(1, a = 0, b = trunc.point,
                                  mean = impute_mean[miss_ind[point_ind, ][2], miss_ind[point_ind, ][1]],
                                  sd = sqrt(impute_var[miss_ind[point_ind, ][2], 1]))

        }

      }

    }

    # ~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~
    # update storage
    # ~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~

    if (m %% thin == 0) {

      # imputed dataset
      store_data[[m / thin]] <- data_working[miss_vec]

      # alpha - prob of MAR
      store_alpha[[m / thin]] <- alpha

      # sigma
      store_sigma_inv[[m / thin]] <- sigma_inv

      # lambda - factor loadings
      store_lambda[[m / thin]] <- lambda

      # eta
      store_eta[[m / thin]] <- eta

      # mu - mean
      store_mu[[m / thin]] <- mu

      # effective factors
      store_k_t[[m / thin]] <- k.star

      # phi
      store_phi[[m / thin]] <- phi

      # delta
      store_delta[[m / thin]] <- delta

      # Z
      store_Z[[m / thin]] <- Z

    }

  }

  TGIFA_res <- list("store_data" = store_data,
                    "store_alpha" = store_alpha,
                    "store_sigma_inv" = store_sigma_inv,
                    "store_lambda" = store_lambda,
                    "store_eta" = store_eta,
                    "store_mu" = store_mu,
                    "store_k_t" = store_k_t,
                    "store_phi" = store_phi,
                    "store_delta" = store_delta,
                    "store_b_sigma" = b_sigma,
                    "store_Z" = store_Z,
                    "store_MH_mu" = store_MH_mu,
                    "store_MH_lambda" = store_MH_lambda,
                    "store_MH_eta" = store_MH_eta,
                    "store_MH_sigma_inv" = store_MH_sigma_inv,
                    "input_data" = input_data,
                    "input_missingness" = missingness,
                    "prop_mult" = prop_mult,
                    "var_mult" = var_mult
  )

  formatted_res <- TGIFA_res_format(TGIFA_res, burn = burn, thin = thin,
                                    cred_level = cred_level)

  output <- list("imputed_dataset" = formatted_res$imputed_dataset,
                 "imputation_info" = formatted_res$imputation_info)

  if (return_chain) {

    output$chain <- TGIFA_res

  }

  return(output)
}
