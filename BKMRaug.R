library(bkmr)
library(ggplot2)
library(dplyr)

# 1. Setup Simulation Parameters
set.seed(99)
n <- 300

p <- 5
n_sims <- 1
rho_corr <- 0 # Correlation between chemicals

offset <- 0
lod_quantile <- 0.2
mcmc_iter <- 1000
varsel <- TRUE

# Define RMSE functions before the loop
rmse <- function(true, pred) sqrt(mean((true - apply(pred, 2, mean))^2))

rmse_by_lod_count <- function(h_true, pred, n_below_lod, p_val = NULL) {
  if(is.null(p_val)) p_val <- max(n_below_lod)
  results <- data.frame(n_below_lod = integer(), n_obs = integer(), rmse = numeric())
  for(k in 0:p_val) {
    idx <- which(n_below_lod == k)
    if(length(idx) > 0) {
      rmse_val <- rmse(h_true[idx], pred[, idx, drop = FALSE])
      results <- rbind(results, data.frame(n_below_lod = k, n_obs = length(idx), rmse = rmse_val))
    }
  }
  total_rmse <- rmse(h_true, pred)
  results <- rbind(results, data.frame(n_below_lod = NA, n_obs = length(h_true), rmse = total_rmse))
  return(results)
}

# Initialize accumulators
all_h_true          <- c()
all_n_below_lod     <- c()
all_pred_uncensored <- NULL
all_pred_impute     <- NULL
all_pred_augmented  <- NULL

if (varsel) {
  all_pips_uncensored <- vector("list", n_sims)
  all_pips_impute <- vector("list", n_sims)
  all_pips_augmented <- vector("list", n_sims)
}

for(sim in 1:n_sims) {
  cat(sprintf("\n=== Simulation %d/%d ===\n", sim, n_sims))
  
  # Generate true underlying exposures (log-scale)
  Sigma <- matrix(rho_corr, nrow = p, ncol = p)
  diag(Sigma) <- 1
  #Z_true <- MASS::mvrnorm(n*2, mu = rep(0, p), Sigma = Sigma)
  
 
  # Z_true <- matrix(runif(2*n*p,0,4),ncol = p) + offset
  # h_true <- 4 * plogis(1/4 * (Z_true[,1] + Z_true[,2] - 2*offset + 1/2 * (Z_true[,1]-offset) * (Z_true[,2]-offset)), location = 0, scale = 0.5)
  # h_true <- (Z_true[, 1] + Z_true[, 2]+1.5)^2 - 2.25

  # Z_true <- matrix(rlnorm(2*n*p),ncol = p)

  Z_true <- matrix(rgamma(2*n*p, 1, 1), ncol = p)
  h_true <- 4 * plogis(1/4 * (Z_true[,1] + Z_true[,2] - 2*offset + 1/2 * (Z_true[,1]-offset) * (Z_true[,2]-offset)), location = 1, scale = 0.3)
  # h_true <- (.75 * (Z_true[, 1] + Z_true[, 2]+1.5))^2 - 2.25
  
  y <- h_true + rnorm(2*n, sd = 1) 
  
  # 2. Impose LoD
  lod <- apply(Z_true, 2, quantile, probs = lod_quantile)
  
  # Create observed data
  Z_obs <- Z_true
  for(j in 1:p) {
    Z_obs[Z_true[,j] < lod[j], j] <- NA
  }
  
  # 3. Prepare Three Datasets
  # A. Uncensored (Oracle/Gold Standard)
  
  Z_uncensored <- scale(Z_true)
  # Z_uncensored <- Z_true
  
  # B. Imputation (LoD / sqrt(2))
  Z_impute <- Z_obs
  for(j in 1:p) {
    Z_impute[is.na(Z_impute[,j]), j] <- lod[j] / sqrt(2)
  }
  Z_impute <- scale(Z_impute)
  
  # C. Augmented (Indicator + Continuous)
  Z_aug_cont <- Z_obs
  Z_aug_ind  <- matrix(1, nrow = n*2, ncol = p)
  for(j in 1:p) {
    Z_aug_cont[is.na(Z_aug_cont[,j]), j] <- lod[j]
    Z_aug_cont[,j] <- Z_aug_cont[,j] - lod[j]
    Z_aug_ind[is.na(Z_obs[,j]), j] <- 0
  }
  Z_aug <- cbind(scale(Z_aug_cont), Z_aug_ind)
  # Z_aug <- cbind(Z_aug_cont, Z_aug_ind)
  groups_aug <- c(1:p,1:p)
  
  
  
  #3.5 Split data
  train_split <- c(1:n)
  test_split <- c((n+1):(2*n))
  
  y_train = y[train_split]
  y_test = y[test_split]
  
  h_train = h_true[train_split]
  h_test = h_true[test_split]
  
  Z_unc_train = Z_uncensored[train_split,]
  Z_unc_test = Z_uncensored[test_split,]
  
  Z_imp_train = Z_impute[train_split,]
  Z_imp_test = Z_impute[test_split,]
  
  Z_aug_train = Z_aug[train_split,]
  Z_aug_test = Z_aug[test_split,]
  
  # 4. Run BKMR Models
  
  # m_uncensored <- kmbayes(y = y_train, Z = Z_unc_train, iter = mcmc_iter, varsel = varsel)
  # m_impute     <- kmbayes(y = y_train, Z = Z_imp_train, iter = mcmc_iter, varsel = varsel)
  m_augmented  <- kmbayes(y = y_train, Z = Z_aug_train, iter = mcmc_iter, varsel = varsel,groups = groups_aug)

  if (varsel) {
    # all_pips_uncensored[[sim]] <- ExtractPIPs(m_uncensored)
    # all_pips_impute[[sim]] <- ExtractPIPs(m_impute)
    all_pips_augmented[[sim]] <- ExtractPIPs(m_augmented)
  }

  # 5. Predict h for each model
  # pred_uncensored <- SamplePred(m_uncensored, Znew = Z_unc_test)
  # pred_impute     <- SamplePred(m_impute, Znew = Z_imp_test)
  pred_augmented  <- SamplePred(m_augmented, Znew = Z_aug_test)
  
  # Determine how many values per observation are below LoD
  n_below_lod <- rowSums(is.na(Z_obs[test_split,]))
  
  # Accumulate across simulations
  all_h_true          <- c(all_h_true, h_test)
  all_n_below_lod     <- c(all_n_below_lod, n_below_lod)
  # all_pred_uncensored <- cbind(all_pred_uncensored, pred_uncensored)
  # all_pred_impute     <- cbind(all_pred_impute, pred_impute)
  all_pred_augmented  <- cbind(all_pred_augmented, pred_augmented)
  
} # end simulation loop

# Calculate combined RMSE stratified by number of censored values
cat(sprintf("\nRMSE across %d simulations (n=%d total observations):\n", n_sims, n_sims * n))

# results_uncens  <- rmse_by_lod_count(all_h_true, all_pred_uncensored, all_n_below_lod, p)
# results_impute  <- rmse_by_lod_count(all_h_true, all_pred_impute, all_n_below_lod, p)
results_augment <- rmse_by_lod_count(all_h_true, all_pred_augmented, all_n_below_lod, p)

# cat("Uncensored:\n")
# print(results_uncens)
# cat("\nImputed:\n")
# print(results_impute)
cat("\nAugmented:\n")
print(results_augment)

if (varsel) {
  pips_by_sim <- tibble::tibble(
    sim = seq_len(n_sims),
    # uncensored = all_pips_uncensored,
    # imputed = all_pips_impute,
    augmented = all_pips_augmented
  )

  cat("\nPIPs by simulation:\n")
  print(pips_by_sim)
}



# save.image(file = "Offset2.RData")




# set.seed(9)
# n <- 30000
# p <- 3
# offset <- 0

# Z_true <- matrix(rgamma(2*n*p, 1, 1), ncol = p)
# # h_true <- 4 * plogis(1/4 * (Z_true[,1] + Z_true[,2] - 2*offset + 1/2 * (Z_true[,1]-offset) * (Z_true[,2]-offset)), location = 1, scale = 0.5)
# # h_true <- (.75 * (Z_true[, 1] + Z_true[, 2]+1.5))^2 - 2.25


# h_true <- 4 * plogis(1/4 * (Z_true[,1] + Z_true[,2] + 1/2 * (Z_true[,1]) * (Z_true[,2])), location = 2, scale = 0.3)

# library(ggplot2)
# df <- data.frame(Z1 = Z_impute[,1], Z2 = Z_impute[,2], h = h_true)
# plt <- ggplot(df, aes(x = Z1, y = Z2, color = h)) +
#   geom_point(alpha = 0.8, size = 1.2) +
#   scale_color_gradient(low = "blue", high = "red") +
#   labs(x = "Z_true[,1]", y = "Z_true[,2]", color = "h_true", title = "Z_true[,1] vs Z_true[,2] colored by h_true") +
#   theme_minimal()
# plt


debugonce(bkmr:::rdelta.group.update)
m_augmented <- kmbayes(
  y = y_train,
  Z = Z_aug_train,
  iter = mcmc_iter,
  varsel = varsel,
  groups = groups_aug
)

trace(rdelta.group.update, browser)
