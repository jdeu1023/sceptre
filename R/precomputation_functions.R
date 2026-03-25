#' Run response precomputation
#'
#' Runs the precomputation for a given response vector, which consists of regressing the response expressions onto the covariate matrix.
#'
#' We (i) perform an initial Poisson regression, (ii) estimate the size parameter using the Poisson residuals, and then (iii) perform an NB regression, treating the size parameter as estimated in step (ii) as fixed and known.
#'
#' @param expressions the numeric vector of response expressions
#' @param covariate_matrix the covariate matrix on which to regress the expressions (NOTE: the matrix should contain an interecept term)
#'
#' @return a list containing the following elements: (i) "fitted_coefs": a vector of fitted coefficients; (ii) "theta": the fitted theta.
#' @noRd
perform_response_precomputation <- function(expressions, covariate_matrix, n_nonzero_thresh = Inf) {
  if (sum(expressions != 0L) < n_nonzero_thresh) {
    gp_fit <- glmGamPoi::glm_gp(data = matrix(expressions, nrow = 1), design = covariate_matrix, size_factors = FALSE,
                                 overdispersion = FALSE, overdispersion_shrinkage = FALSE, do_cox_reid_adjustment = FALSE)
    fitted_coefs <- gp_fit$Beta[1, ]
    fitted_values <- as.vector(gp_fit$Mu)
  } else {
    glm_fit <- stats::glm.fit(y = expressions, x = covariate_matrix, family = stats::poisson())
    fitted_coefs <- glm_fit$coefficients
    fitted_values <- glm_fit$fitted.values
  }
  dfr <- length(expressions) - ncol(covariate_matrix)
  response_theta_list <- estimate_theta(
    y = expressions, mu = fitted_values, dfr = dfr,
    limit = 50, eps = (.Machine$double.eps)^(1 / 4)
  )
  # check that NAs are absent from the fitted coefficient vector
  if (any(is.na(fitted_coefs))) {
    problem_covariates <- paste0(names(which(is.na(fitted_coefs))))
    stop("The coefficients corresponding to the following covariates cannot be estimated in the regression model: ",
         paste0(problem_covariates, collapse = ", "), ". Consider removing these covariates from the model (by updating `formula_object` in `set_analysis_parameters()`).")
  }
  theta <- max(min(response_theta_list[[1]], 1000), 0.01)
  result <- list(fitted_coefs = fitted_coefs, theta = theta)
  return(result)
}


perform_grna_precomputation <- function(trt_idxs, covariate_matrix, return_fitted_values) {
  indicator <- integer(length = nrow(covariate_matrix))
  indicator[trt_idxs] <- 1L
  logistic_fit <- stats::glm.fit(y = indicator, x = covariate_matrix, family = stats::binomial())
  if (return_fitted_values) {
    out <- logistic_fit$fitted.values
  } else {
    out <- logistic_fit$coefficients
  }
  return(out)
}


compute_D_matrix <- function(Zt_wZ, wZ) {
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1 / sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  return(D)
}


compute_precomputation_pieces <- function(expression_vector, covariate_matrix, fitted_coefs, theta, full_test_stat) {
  mu <- exp(as.numeric(covariate_matrix %*% fitted_coefs))
  if (full_test_stat) {
    denom <- 1 + mu / theta
    w <- mu / denom
    a <- (expression_vector - mu) / denom
    wZ <- w * covariate_matrix
    Zt_wZ <- t(covariate_matrix) %*% wZ
    D <- compute_D_matrix(Zt_wZ, wZ)
    out <- list(mu = mu, w = w, a = a, D = D)
  } else {
    a <- expression_vector - (expression_vector * mu + theta * mu) / (theta + mu)
    b <- (theta * mu) / (theta + mu)
    out <- list(mu = mu, a = a, b = b)
  }
  return(out)
}
