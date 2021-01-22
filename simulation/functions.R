# rho - correlation of the design (equicorrelated)
generate_SRRR_data <- function(n, p, Sigma_X, q, sparsity, rank, C_noise, Sigma_E, SNR, n_test = NULL, seed = NULL, C_true = NULL, missing_rate = NULL) {
  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  if (length(Sigma_X) > 1) {
    X_train <- mvrnorm(n, rep(0, p), Sigma_X)
  } else { # rho = Sigma_X, equicorrelated normal dist
    X_train_Y <- rnorm(n)
    X_train_Z <- matrix(rnorm(n*p), n, p)
    X_train <- sqrt(Sigma_X) * X_train_Y + sqrt(1-Sigma_X) * X_train_Z
  }

  if (is.null(C_true)) {
    C_true <- matrix(0, p, q)
    supp_true <- sample(p, sparsity)
    C_true[supp_true, ] <- tcrossprod(matrix(rnorm(sparsity*rank), sparsity, rank), matrix(rnorm(q*rank), q, rank))

    if (C_noise != 0) {
      C_true[supp_true, ] <- C_true[supp_true, ] + matrix(rnorm(sparsity*q, 0, C_noise), sparsity, q)
    }
  }

  if (length(Sigma_X) > 1) {
    Var_Mu <- sum(diag(crossprod(C_true, Sigma_X %*% C_true)))
  } else {
    Var_Mu <- Sigma_X * sum(apply(C_true, 2, sum)^2) + (1-Sigma_X) * sum(diag(crossprod(C_true)))
  }
  Var_E <- sum(diag(Sigma_E))
  sigma2_e <- Var_Mu/Var_E/SNR
  E_train <- mvrnorm(n, rep(0, q), sigma2_e * Sigma_E)

  Mu_train <- X_train %*% C_true
  Y_train <- Mu_train + E_train

  if (!is.null(missing_rate) && missing_rate > 0) {
    miss_idx <- matrix(sample(c(TRUE, FALSE), size = n*q, prob = c(missing_rate, 1-missing_rate), replace = T), n, q)
    Y_train[miss_idx] <- NA
  }

  out <- list(X_train = X_train, Mu_train = Mu_train, Y_train = Y_train, sigma2_e = sigma2_e, call = call, C_true = C_true,
              supp_true = supp_true)

  if (!is.null(n_test)) {
    if (length(Sigma_X) > 1) {
      X_test <- mvrnorm(n_test, rep(0, p), Sigma_X)
    } else { # rho = Sigma_X, equicorrelated normal dist
      X_test_Y <- rnorm(n_test)
      X_test_Z <- matrix(rnorm(n_test*p), n_test, p)
      X_test <- sqrt(Sigma_X) * X_test_Y + sqrt(1-Sigma_X) * X_test_Z
    }
    Mu_test <- X_test %*% C_true
    E_test <- mvrnorm(n_test, rep(0, q), sigma2_e * Sigma_E)
    out$X_test <- X_test
    out$Mu_test <- Mu_test
    out$Y_test <- Mu_test + E_test
  }

  out
}


norm2 <- function(x) {
  out <- sqrt(sum(x^2))
  out
}

which_row_active <- function(X) {
  out <- which(apply(X, 1, function(x) any(x != 0)))
  out
}

which_row_active_glmnet <- function(fit, ilam) {
  out <- which_row_active(sapply(fit$beta, function(mm) mm[, ilam]))
  out
}

SRRR_relaxed_path <- function(X, Y, r, SRRR_fit, max.iter = 10, thresh = 1E-7) {
  n <- nrow(X)
  p <- ncol(X)
  q <- nrow(Y)
  nlambda <- length(SRRR_fit)
  fit_list <- vector("list", nlambda)
  for (ilam in 1:nlambda) {
    ind_map <- SRRR_fit[[ilam]]$ind_map
    B_active <- which_row_active(SRRR_fit[[ilam]]$B)
    active_set <- ind_map[B_active]
    if (length(active_set) > n) {
      fit_list <- fit_list[1:(ilam-1)]
      break
    }
    if (length(active_set) != 0) {
      fit_list[[ilam]] <- RRR_iterative(X[, active_set, drop = F], Y, r, B0 = SRRR_fit[[ilam]]$B[B_active, , drop=F], niter = max.iter, thresh = thresh)
      fit_list[[ilam]]$ind_map <- active_set
    } else {
      fit_list[[ilam]] <- NA
    }
  }
  fit_list
}


SRRR_path_glmnet <- function(X, Y, r, glmnet.fit, nlambda = 50, batch.size = 100, lambda.min.ratio = 0.01, max.iter = 10, thresh = 1e-7, seed = 1223,
                      X_val = NULL, Y_val = NULL, early_stopping = FALSE, lag = 2, is.warm.start = TRUE, is.A.converge = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  lambda <- glmnet.fit$lambda
  max.nlambda <- length(lambda)
  nlambda <- min(nlambda, max.nlambda)
  
  # loop over lambda
  fit_list <- vector("list", nlambda)
  niter_list <- vector("list", nlambda)
  
  num_violations <- rep(NA, nlambda)
  B_init <- NULL
  
  for (ilam in 1:nlambda) {
    cat(paste0(ilam, " "))
    num_violations[ilam] <- 0
    lam <- lambda[ilam]
    active_glmnet <- which_row_active_glmnet(glmnet.fit, ilam)
    if (length(active_glmnet) <= 1) {
      fit_list[[ilam]] <- NA
      next
    }
    features_train <- X[, active_glmnet, drop = F]
    B_init <- matrix(sapply(glmnet.fit$beta, function(mm) mm[active_glmnet, ilam]), ncol = q)
    B_svd <- svd(B_init, r, r)
    B_init <- B_svd$u %*% diag(B_svd$d, nrow = r)
    
    strong_set <- active_glmnet
    ind_map <- active_glmnet
    check <- FALSE
    while (!check) {
      fit <- SRRR_iterative(features_train, Y, lam, r, max.iter, B_init, thresh, lambda_seq = lambda, is.warm.start = is.warm.start)
      niter_list[[ilam]] <- c(niter_list[[ilam]], fit$niter)
      prod_resid <- crossprod(X, Y %*% fit$A - features_train %*% fit$B)
      norm_prod <- apply(prod_resid, 1, norm2)
      norm_prod[strong_set] <- NA
      check <- all(norm_prod <= max(n*lam), na.rm = T)
      if (!check) {
        var_violate <- which(norm_prod > max(n*lam))
        num_violations[ilam] <- num_violations[ilam] + 1
        print(check)
        print(length(var_violate))
        features_train <- cbind(features_train, X[, var_violate, drop = F])
        B_init <- rbind(B_init, matrix(0, nrow = length(var_violate), ncol = r))
        ind_map <- c(ind_map, var_violate)
        strong_set <- c(strong_set, var_violate)
      }
    }
    cat("\n")
    active <- which_row_active(fit$B)
    fit$ind_map <- ind_map
    fit$active <- ind_map[active]
    fit_list[[ilam]] <- fit
    if (early_stopping) {
      err_val[ilam] <- mean((Y_val - X_val[, ind_map, drop = F] %*% fit$C)^2)
      if (min(tail(err_val, lag)) > min(head(err_val, max(ilam-lag, 0)))) {
        cat(paste0("Early stopping at ", ilam, ". Trailing Err: ", paste(tail(err_val, 2*lag), collapse = " "), "\n"))
        break
      }
    }
  }
  fit_list <- fit_list[1:ilam]
  names(fit_list) <- paste0("s", 0:(ilam-1))
  niter_list <- niter_list[1:ilam]
  out <- list(fit = fit_list, niter = niter_list)
  out
}




SRRR_path <- function(X, Y, r, nlambda = 50, batch.size = 100, lambda.min.ratio = 0.01, max.nlambda = 100, max.iter = 10, thresh = 1e-7, seed = 1223,
                      X_val = NULL, Y_val = NULL, early_stopping = FALSE, lag = 2, is.warm.start = TRUE, is.A.converge = TRUE) {
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  missing_Y <- NULL
  if (any(is.na(Y))) {
    missing_Y <- is.na(Y)
    Y_zero <- Y
    Y_zero[missing_Y] <- 0
    prod_init <- crossprod(X, Y_zero)
    count_no_missing_Y <- apply(1-missing_Y, 2, sum)
    prod_init <- t(t(prod_init) / count_no_missing_Y) * n
    norm_prod_init <- apply(prod_init, 1, norm2)
    for (j in 1:ncol(Y)) {
      Y[is.na(Y[, j]), j] <- mean(Y[, j], na.rm = T)
    }
  } else {
    prod_init <- crossprod(X, Y)
    norm_prod_init <- apply(prod_init, 1, norm2)
  }

  max_lambda <- max(norm_prod_init) / n
  print(max_lambda)
  lambda <- exp(seq(from = log(max_lambda), to = log(lambda.min.ratio*max_lambda), length.out = max.nlambda))

  var0 <- sum(Y^2)/(2*n)
  
  # loop over lambda
  fit_list <- vector("list", nlambda)
  niter_list <- vector("list", nlambda)
  if (!is.null(missing_Y)) A_niter_list <- vector("list", nlambda)

  num_violations <- rep(NA, nlambda)
  features_train <- NULL
  ind_map <- c()
  active <- c()
  if (early_stopping) err_val <- c()
  B_init <- NULL
  for (ilam in 1:nlambda) {
    cat(paste0(ilam, " "))
    num_violations[ilam] <- 0
    lam <- lambda[ilam]
    check <- FALSE
    if (length(active) == 0) {
      norm_prod <- norm_prod_init
    }
    if (!is.null(features_train)) {
      features_train <- features_train[, active, drop = FALSE]
      B_init <- B_init[active, , drop = FALSE]
    }
    ind_map <- ind_map[active]
    var_to_ignore <- ind_map
    idx_strong <- c()
    var_violate <- c()
    while (!check) {
      var_to_ignore <- c(var_to_ignore, idx_strong)
      norm_prod[var_to_ignore] <- NA
      order_norm <- order(norm_prod, decreasing = TRUE, na.last = NA)
      idx_strong <- c(head(order_norm, batch.size), setdiff(var_violate, order_norm))
      features_train <- cbind(features_train, X[, idx_strong, drop = FALSE])
      ind_map <- c(ind_map, idx_strong)
      if (is.null(B_init) || nrow(B_init) == 0) {
        B_init <- matrix(rnorm(r*length(idx_strong)), length(idx_strong), r)
      } else {
        B_init <- rbind(B_init, matrix(0, length(idx_strong), r))
      }
      if (is.null(missing_Y)) {
        fit <- SRRR_iterative(features_train, Y, lam, r, max.iter, B_init, thresh*var0, lambda_seq = lambda, is.warm.start = is.warm.start)
        niter_list[[ilam]] <- c(niter_list[[ilam]], fit$niter)
      } else {
        fit <- SRRR_iterative_missing(features_train, Y, missing_Y, lam, r, max.iter, B_init, thresh*var0, lambda_seq = lambda, is.warm.start = is.warm.start, is.A.converge = is.A.converge)
        Y <- fit$Y
        niter_list[[ilam]] <- c(niter_list[[ilam]], fit$niter)
        A_niter_list[[ilam]] <- c(A_niter_list[[ilam]], fit$A_niter)
      }
      # check KKT
      prod_resid <- crossprod(X, Y %*% fit$A - features_train %*% fit$B)
      norm_prod <- apply(prod_resid, 1, norm2)
      max_norm_strong <- max(norm_prod[c(var_to_ignore, idx_strong)])
      norm_prod_with_NA <- norm_prod
      norm_prod_with_NA[c(var_to_ignore, idx_strong)] <- NA
      check <- all(norm_prod_with_NA <= max(n*lam, max_norm_strong), na.rm = T)  # loose check
      if (!check) {
        var_violate <- which(norm_prod_with_NA > max(n*lam, max_norm_strong))
        num_violations[ilam] <- num_violations[ilam] + 1
        print(check)
      }
    }
    cat("\n")
    active <- which_row_active(fit$B)
    B_init <- fit$B
    fit$ind_map <- ind_map
    fit$active <- ind_map[active]
    fit_list[[ilam]] <- fit
    if (early_stopping) {
      err_val[ilam] <- mean((Y_val - X_val[, ind_map, drop = F] %*% fit$C)^2)
      if (min(tail(err_val, lag)) > min(head(err_val, max(ilam-lag, 0)))) {
        cat(paste0("Early stopping at ", ilam, ". Trailing Err: ", paste(tail(err_val, 2*lag), collapse = " "), "\n"))
        break
      }
    }
  }
  fit_list <- fit_list[1:ilam]
  names(fit_list) <- paste0("s", 0:(ilam-1))
  niter_list <- niter_list[1:ilam]
  out <- list(fit = fit_list, niter = niter_list)
  if (!is.null(missing_Y)) out$A_niter <- A_niter_list
  out
}


# solve reduced-rank regression using Chen 2012
SRRR_iterative_missing <- function(X, Y, Y_missing, lambda, r, niter, B0, thresh = 1e-7, lambda_seq = NULL, is.warm.start = FALSE, is.A.converge = TRUE) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  B <- B0

  obj_values <- rep(NA, niter)
  A_niter <- rep(NA, niter)
  message <- "Terminated"
  
  count <- 0

  for (k in 1:niter) {
    count <- count + 1
    A_niter[k] <- 0
    # fix B, solve A
    score <- X %*% B
    impute_iter_count <- 0
    while (TRUE) {
      A_niter[k] <- A_niter[k] + 1
      impute_iter_count <- impute_iter_count + 1
      crossmat <- crossprod(Y, score)
      svd_cross <- svd(crossmat)
      A <- tcrossprod(svd_cross$u, svd_cross$v)
      Y_new <- tcrossprod(score, A)
      delta <- mean((Y[Y_missing] - Y_new[Y_missing])^2)
      # browser()
      Y[Y_missing] <- Y_new[Y_missing]
      if (delta < thresh || (!is.A.converge && A_niter[k] > 0)) {
        print(impute_iter_count)
        break
      }
    }
    # fix A, solve B
    YA <- Y %*% A
    # keep a strong set, update only when KKT is violated ...
    if (!is.null(lambda_seq)) {
      lambda_seq <- sort(unique(c(lambda_seq, lambda)), decreasing = T)
      lambda_seq <- lambda_seq[lambda_seq >= lambda]
    }
    if (k == 1 && !is.warm.start) {
      mfit <- glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda_seq)
    } else {
      mfit <- glmnetPlus::glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda, beta0 = B)
    }
    beta_single <- coef(mfit, s = lambda, x = X, y = YA)
    B <- do.call(cbind, beta_single)[-1, ]
    C <- tcrossprod(B, A)
    obj_values[k] <- 1/(2*n) * sum((Y - X %*% C)^2) +
      lambda * sum(apply(B, 1, function(x) sqrt(sum(x^2))))
    if (k > 1 && abs(obj_values[k] - obj_values[k-1])/obj_values[1] < thresh) {
      message <- "Converged"
      obj_values <- obj_values[1:k]
      A_niter <- A_niter[1:k]
      print(k)
      break
    }
  }
  if (message == "Terminated") print(message)

  # browser()
  out <- list(B = B, A = A, C = C, object = obj_values, message = message, niter = count, Y = Y, A_niter = A_niter)
  out
}



# solve reduced-rank regression using Chen 2012
SRRR_iterative <- function(X, Y, lambda, r, niter, B0, thresh = 1e-7, lambda_seq = NULL, is.warm.start = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)

  B <- B0

  obj_values <- rep(NA, niter)
  message <- "Terminated"
  count <- 0

  for (k in 1:niter) {
    count <- count + 1
    # fix B, solve A
    crossmat <- crossprod(crossprod(X, Y), B)
    svd_cross <- svd(crossmat)
    A <- tcrossprod(svd_cross$u, svd_cross$v)
    # fix A, solve B
    YA <- Y %*% A
    # keep a strong set, update only when KKT is violated ...
    if (!is.null(lambda_seq)) {
      lambda_seq <- sort(unique(c(lambda_seq, lambda)), decreasing = T)
      lambda_seq <- lambda_seq[lambda_seq >= lambda]
    }
    if (k == 1 && !is.warm.start) {
      mfit <- glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda_seq)
    } else {
      mfit <- glmnetPlus::glmnet(x = X, y = YA, family = "mgaussian", standardize = F, intercept = F, lambda = lambda, beta0 = B)
    }
    beta_single <- coef(mfit, s = lambda, x = X, y = YA)
    B <- do.call(cbind, beta_single)[-1, ]
    C <- tcrossprod(B, A)
    obj_values[k] <- 1/(2*n) * sum((Y - X %*% C)^2) +
      lambda * sum(apply(B, 1, function(x) sqrt(sum(x^2))))
    if (k > 1 && abs(obj_values[k] - obj_values[k-1]) < thresh) {
      message <- "Converged"
      obj_values <- obj_values[1:k]
      print(k)
      break
    }
  }
  if (message == "Terminated") print(message)

  out <- list(B = B, A = A, C = C, object = obj_values, message = message, niter = count)
  out
}


RRR_iterative <- function(X, Y, r, niter, B0, thresh = 1e-7) {
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  
  B <- B0
  
  obj_values <- rep(NA, niter)
  message <- "Terminated"
  count <- 0
  
  for (k in 1:niter) {
    count <- count + 1
    # fix B, solve A
    crossmat <- crossprod(crossprod(X, Y), B)
    svd_cross <- svd(crossmat)
    A <- tcrossprod(svd_cross$u, svd_cross$v)
    # fix A, solve B
    YA <- Y %*% A
    # keep a strong set, update only when KKT is violated ...
    B <- solve(crossprod(X), crossprod(X, YA))
    C <- tcrossprod(B, A)
    obj_values[k] <- 1/(2*n) * sum((Y - X %*% C)^2)
    if (k > 1 && abs(obj_values[k] - obj_values[k-1])/obj_values[1] < thresh) {
      message <- "Converged"
      obj_values <- obj_values[1:k]
      print(k)
      break
    }
  }
  if (message == "Terminated") print(message)
  
  out <- list(B = B, A = A, C = C, object = obj_values, message = message, niter = count)
  out
}


find_R2_MSE_SRRR <- function(fit, X, Y, C_true = NULL) {
  ind_not_NA <- which(sapply(fit, is.list))
  pred_SRRR <- simplify2array(lapply(fit[ind_not_NA], function(ff) as.matrix(X[, ff$ind_map, drop = F] %*% tcrossprod(ff$B, ff$A))))
  if (!is.null(C_true)) {
    pred_true <- X %*% C_true
    MSE_SRRR <- rep(NA, length(fit))
    MSE_scaled_SRRR <- rep(NA, length(fit))
    MSE_SRRR[ind_not_NA] <- apply(pred_SRRR, 3, function(z) mean((pred_true-z)^2))
    MSE_scaled_SRRR[ind_not_NA] <- MSE_SRRR[ind_not_NA] / mean(scale(pred_true, scale = F)^2)
  }
  R2_SRRR <- rep(NA, length(fit))
  R2_SRRR[ind_not_NA] <- 1 - apply(pred_SRRR, 3, function(z) mean((Y-z)^2)) / mean(scale(Y, scale = F)^2)
  out <- list(R2 = R2_SRRR)
  if (!is.null(C_true)) {
    out <- list(out, list(MSE = MSE_SRRR, MSE_scaled = MSE_scaled_SRRR))
  }
  out
}


find_selection_accuracy <- function(fit, true_coef) {
  p <- nrow(true_coef)
  truth <- rep("inactive", p)
  which_true_active <- which_row_active(true_coef)
  truth[which_true_active] <- "active"
  truth <- factor(truth, levels = c("inactive", "active"))
  selection_list <- vector("list", length(fit))
  for (k in 1:length(fit)) {
    which_est_active <- fit[[k]]$ind_map[which_row_active(fit[[k]]$B)]
    estimate <- rep("inactive", p)
    estimate[which_est_active] <- "active"
    estimate <- factor(estimate, levels = c("inactive", "active"))
    selection_list[[k]] <- table(estimate = estimate, truth = truth)
  }
  selection_list
}
