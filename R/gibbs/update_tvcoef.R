# R port of MCT-NYFED/functions/update_tvcoef.m
#
# Carter-Kohn (1994) forward-filter backward-sample for a time-varying
# parameter (TVP) regression:
#
#   y_t          = alpha_t * x_t + sigma_t .* eps_t,    eps ~ N(0, I_N)
#   vec(alpha_t) = vec(alpha_{t-1}) + lambda .* ups_t,  ups ~ N(0, I_{N*K})
#   alpha_0      ~ N(0, var_init)
#
# Used by NY Fed to update the time-varying loadings alpha_tau and
# alpha_eps. Our Variant Ac uses CONSTANT loadings (lambda is an
# N-vector, not a path), so this is needed only if we later switch
# to TVP loadings — port now to keep the option open.
#
# Output: T x (N*K) matrix where row t is vec(alpha_t).

#' Draw alpha_t (TVP regression coefficients) via Carter-Kohn FFBS.
#'
#' @param y        T x N matrix of regressands (NA for missing).
#' @param x        T x K matrix of regressors.
#' @param sigma    T x N matrix of obs-noise standard deviations.
#' @param lambda   Length-N vector of RW step sizes on the coefficients.
#' @param var_init (N*K) x (N*K) prior covariance for alpha_0. Default
#'                 mimics MATLAB: 0.1*I + 100 * kron(I_K, J_N).
#' @return T x (N*K) matrix; row t is vec(alpha_t).
update_tvcoef <- function(y, x, sigma, lambda, var_init = NULL) {
  if (is.null(dim(y))) y <- matrix(y, ncol = 1L)
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  if (is.null(dim(sigma))) sigma <- matrix(sigma, ncol = ncol(y))

  T_     <- nrow(y)
  N      <- ncol(y)
  K      <- ncol(x)
  n_coef <- N * K
  stopifnot(length(lambda) == N)

  if (is.null(var_init)) {
    var_init <- 0.1 * diag(n_coef) +
                100 * kronecker(diag(K), matrix(1, N, N))
  }

  Sigma_eta <- diag(lambda^2, nrow = N)
  # Sigma_eta needs to act on vec(alpha_t) of length n_coef.
  # vec(alpha_t) stacks columns; with K coefficients each Nx1, that's
  # K blocks of size N. The N RW step sizes apply within each block,
  # so the full Sigma_eta is kron(I_K, diag(lambda^2)).
  Sigma_eta <- kronecker(diag(K), diag(lambda^2, nrow = N))

  # Symmetric Cholesky inverse (matches MATLAB linsolve SYM+POSDEF)
  sym_inv <- function(A) {
    chol2inv(chol(A))
  }

  # Storage for forward pass
  X1_KF <- matrix(0, n_coef, T_)
  P1_KF <- array(0, c(n_coef, n_coef, T_))
  X2_KF <- matrix(0, n_coef, T_)
  P2_KF <- array(0, c(n_coef, n_coef, T_))

  X1 <- rep(0, n_coef)
  P1 <- var_init

  for (t in seq_len(T_)) {
    y_t <- y[t, ]
    # H = kron(x[t, ], diag(N))  shape: N x (N*K)
    H_full <- kronecker(matrix(x[t, ], 1L, K), diag(N))
    Sigma_eps_full <- diag(sigma[t, ]^2, nrow = N)

    # Predict
    X2 <- X1
    P2 <- P1 + Sigma_eta

    miss <- is.na(y_t)
    if (all(miss)) {
      X1 <- X2
      P1 <- P2
    } else {
      kept   <- !miss
      H      <- H_full[kept, , drop = FALSE]
      Sigeps <- Sigma_eps_full[kept, kept, drop = FALSE]
      e      <- y_t[kept] - as.numeric(H %*% X2)
      S      <- H %*% P2 %*% t(H) + Sigeps
      S_inv  <- tryCatch(sym_inv(.symmetrize(S)),
                         error = function(e) MASS::ginv(.symmetrize(S)))
      K_mat  <- P2 %*% t(H) %*% S_inv
      X1     <- as.numeric(X2 + K_mat %*% e)
      P1     <- (diag(n_coef) - K_mat %*% H) %*% P2
    }

    X1_KF[, t]    <- X1
    P1_KF[, , t]  <- P1
    X2_KF[, t]    <- X2
    P2_KF[, , t]  <- P2
  }

  # Backward sample (Carter-Kohn)
  alpha <- matrix(NA_real_, T_, n_coef)
  X3    <- X1
  P3    <- .symmetrize(P1)
  Xdraw <- as.numeric(MASS::mvrnorm(1, X3, P3))
  alpha[T_, ] <- Xdraw

  for (t in (T_ - 1L):1L) {
    X1 <- X1_KF[, t]
    P1 <- P1_KF[, , t]
    X2 <- X2_KF[, t + 1L]
    P2 <- P2_KF[, , t + 1L]
    # P3_denom = P2^{-1} P1 (with pinv fallback)
    P3_denom <- tryCatch(sym_inv(.symmetrize(P2)) %*% P1,
                         error = function(e) MASS::ginv(P2) %*% P1)
    X3 <- as.numeric(X1 + t(P3_denom) %*% (Xdraw - X2))
    P3 <- .symmetrize(P1 - t(P3_denom) %*% P1)
    # Numerical guard: P3 occasionally drifts slightly non-PSD
    P3_eig <- eigen(P3, symmetric = TRUE, only.values = TRUE)$values
    if (min(P3_eig) < 0) {
      P3 <- P3 + (-min(P3_eig) + 1e-10) * diag(n_coef)
    }
    Xdraw <- as.numeric(MASS::mvrnorm(1, X3, P3))
    alpha[t, ] <- Xdraw
  }

  alpha
}
