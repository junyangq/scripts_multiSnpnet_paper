library(latex2exp)
library(glmnet)
library(ggplot2)
library(MASS)

args <- commandArgs(trailingOnly = TRUE)

source("functions.R")

scenario <- 10

if (length(args) == 0) {
  stop("At least one argument must be supplied for the scenario.", call. = FALSE)
} else {
  scenario <- as.integer(args[1])
}

source("configs.R")

glmnet.control(fdev = 0, devmax = 1)

set.seed(1223)

nlambda <- 100
B <- 30
rank_choices <- c(2, 3, 4, 5)

mydata_list <- vector(mode = "list", length = B)

for (b in 1:B) {
  mydata_list[[b]] <- generate_SRRR_data(
    n = n,
    p = p,
    Sigma_X = Sigma_X,
    q = q,
    sparsity = sparsity,
    rank = rank,
    C_noise = C_noise,
    Sigma_E = Sigma_E,
    SNR = SNR,
    n_test = n_test
  )
}

# saveRDS(mydata_list, paste0("../results/data_scenario_", scenario, ".rds"))

results <- vector(mode = "list", length = B)

for (b in 1:B) {
  cat(paste0("B = ", b, "\n"))
  mydata <- mydata_list[[b]]

  fit_glmnet <- glmnet(mydata$X_train, mydata$Y_train, family = "mgaussian", standardize = F, intercept = F)
  pred_glmnet_test <- predict(fit_glmnet, newx = mydata$X_test, s = fit_glmnet$lambda[1:nlambda])
  R2_glmnet_test <- 1 - apply(pred_glmnet_test, 3, function(z) mean((mydata$Y_test-z)^2)) / mean(scale(mydata$Y_test, scale = F)^2)
  fit_glmnet$R2_test <- R2_glmnet_test
  coef_array <- simplify2array(lapply(fit_glmnet$beta, function(cc) as.matrix(cc)))
  fit_glmnet$active <- apply(coef_array, 2, which_row_active)
  fit_glmnet$best_R2_test <- max(R2_glmnet_test)
  fit_glmnet$best_active_test <- fit_glmnet$active[[which.max(R2_glmnet_test)]]

  X_oracle_train <- mydata$X_train[, mydata$supp_true, drop = F]
  SXX_oracle_train <- crossprod(X_oracle_train) / nrow(X_oracle_train)
  SXY_oracle_train <- crossprod(X_oracle_train, mydata$Y_train) / nrow(X_oracle_train)
  inv_mult <- solve(SXX_oracle_train, SXY_oracle_train)
  mult_train <- crossprod(SXY_oracle_train, inv_mult)
  eigen_train <- eigen(mult_train, symmetric = TRUE)
  A_oracle <- eigen_train$vectors[, 1:rank]
  B_oracle <- inv_mult %*% A_oracle
  pred_SRRR_oracle <- mydata$X_test[, mydata$supp_true, drop = F] %*% tcrossprod(B_oracle, A_oracle)
  fit_oracle_SRRR <- list()
  fit_oracle_SRRR$B <- B_oracle
  fit_oracle_SRRR$A <- A_oracle
  fit_oracle_SRRR$R2_test <- 1 - mean((mydata$Y_test-pred_SRRR_oracle)^2) / mean(scale(mydata$Y_test, scale = F)^2)

  fit_SRRR_path <- vector(mode = "list", length = length(rank_choices))
  for (ell in 1:length(rank_choices)) {
    cat(paste0("rank: ", rank_choices[ell], "\n"))
    fit_SRRR_path[[ell]] <- list(fit = SRRR_path(mydata$X_train, mydata$Y_train, r = rank_choices[ell], nlambda = nlambda,
                                      batch.size = 20, lambda.min.ratio = ifelse(n < p, 0.01, 0.0001), max.nlambda = 100, max.iter = 20,
                                      thresh = 1e-7))
    R2_MSE_SRRR_test <- find_R2_MSE_SRRR(fit_SRRR_path[[ell]]$fit$fit, mydata$X_test, mydata$Y_test, mydata$C_true)
    fit_SRRR_path[[ell]]$MSE_test <- R2_MSE_SRRR_test$MSE
    fit_SRRR_path[[ell]]$R2_test <- R2_MSE_SRRR_test$R2
    fit_SRRR_path[[ell]]$best_R2_test <- max(R2_MSE_SRRR_test$R2)
    fit_SRRR_path[[ell]]$active <- lapply(fit_SRRR_path[[ell]]$fit, function(ff) ff$active)
    fit_SRRR_path[[ell]]$best_active_test <- fit_SRRR_path[[ell]]$fit$active[[which.max(R2_MSE_SRRR_test$R2)]]
  }
  results[[b]] <- list(glmnet = fit_glmnet, SRRR = fit_SRRR_path, SRRR_rank_choices = rank_choices, oracle = fit_oracle_SRRR)
}

# saveRDS(results, paste0("../results/results_scenario_", scenario, ".rds"))

glmnet.control(factory = TRUE)
