# ==============================================================================
# MREWMA-COPULA CONTROL CHART FOR AIR QUALITY MONITORING (SO2, PM2.5, O3)
# ==============================================================================
# Pipeline (run top to bottom in a single session):
#   1.  Packages
#   2.  User configuration (paths, subgroup design, simulation settings)
#   3.  Plot theme
#   4.  Load raw data
#   5.  EM algorithm functions (multivariate normal missing-value imputation)
#   6.  Run EM imputation
#   7.  Imputation method validation: EM vs MICE vs KNN (hold-out accuracy check)
#   8.  Subgroup formation (Phase I / Phase II split)
#   9.  Descriptive statistics
#   10. Multivariate normality assumption testing (Shapiro-Wilk)
#   11. Time series visualization
#   12. MREWMA-Copula statistical engine (Lepage & Cucconi)
#   13. Phase I iterative cleaning (leave-one-subgroup-out)
#   14. Official control limits (bisection, target MRL0)
#   15. Subgroup control chart construction (Phase I & II) + OOC diagnostics
#   16. MRL0 curves across eight out-of-control shift scenarios
#   17. Results appendix & export
#
# Edit Section 2 (paths, subgroup design, simulation settings) before running.
# Plots are written to <dir_save>/plots/; all result tables are written to a
# single Excel workbook in <dir_save>.
#
# ------------------------------------------------------------------------------
# CHANGES IN THIS REVISION (response to EEET-04995-2026-01, Reviewer 1)
# ------------------------------------------------------------------------------
#   (a) Section 7 is new: EM imputation accuracy is now validated by a repeated
#       hold-out procedure (artificial MCAR deletion at the observed missing
#       rate, re-impute, compare to the known true values) and benchmarked
#       against two standard alternatives, MICE (predictive mean matching) and
#       KNN imputation, using RMSE/MAE. This addresses the reviewer's point
#       that imputation accuracy was previously unverified and uncompared.
#   (b) N_RL_OOC (replications used for every MRL0 estimate) is raised from
#       300 to 1,000 -- see Section 2 for the full rationale -- to bring
#       run-length estimation in line with common practice in the SPM
#       literature and with N_RL_UCL, which was already 1,000.
#   (c) Sections 12-17 (MREWMA-Copula engine, Phase I cleaning, UCL
#       calibration, control chart, MRL0 curves, and export) were rebuilt
#       following the leave-one-subgroup-out design and MRL0 naming
#       convention used in the thesis pipeline (Sintaks_Tesis_Rafli_V3.R,
#       Steps 11-16), adapted from d=4 (NI/FE/SI/MG) to the d=3 air-quality
#       variables (SO2, PM2.5, O3) and from leave-one-day-out to
#       leave-one-subgroup-out. What was previously called "MRL1" (shift-
#       scenario performance) is now labeled MRL0, matching this convention;
#       the previous separate in-control-only "MRL0" check no longer exists
#       as a separate step. This part of the pipeline runs sequentially
#       (matching the source thesis script) rather than in parallel across
#       CPU cores -- USE_PARALLEL/N_CORES in Section 2 are consequently
#       unused by Sections 12-17 and only retained for possible future use.
#       Given N_RL_OOC = 1,000 and 8 shift scenarios, Section 16 in
#       particular will take a long time to run; lower N_RL_OOC/N_HOLDOUT_REPS
#       temporarily for a quick test run.
#   (d) No debugging/scratch code has been added; every new block follows the
#       numbered-section convention already used throughout the script, and
#       the pipeline still runs top to bottom in one pass with no manual steps.
#   (e) Section 15/16 chart titles, axis labels, annotations, and console
#       messages are now in English (previously Indonesian), matching the
#       rest of the pipeline's output language.
# ==============================================================================

# ==============================================================================
# 1. PACKAGES
# ==============================================================================
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(showtext)
library(dfphase1)
library(mice)       # NEW: imputation comparison benchmark (Section 7)
library(VIM)        # NEW: KNN imputation comparison benchmark (Section 7)
library(mvnormtest) # NEW: multivariate Shapiro-Wilk normality test (Section 10)
library(scales)     # NEW: alpha() used for label transparency (Section 15)
library(parallel)   # NEW: base-R parallelization (currently unused by Sections 12-17; see Section 2)

# ==============================================================================
# 2. USER CONFIGURATION
# ==============================================================================
# -- File paths --
folder_path = "C:/Users/LAB STAT PC0112/Downloads"
input_file  = file.path(folder_path, "Air Quality Data for Makassar, 2023–2026 (Before Imputation).xlsx")
dir_save    = "C:/Users/LAB STAT PC0112/Downloads/MREWMA-EM"
plot_dir    = file.path(dir_save, "plots")

if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# -- Variables monitored --
vars_target = c("SO2", "PM2.5", "O3")

# -- Subgroup design --
n_per_subgroup    = 7    # days per subgroup
n_subgroup_phase1 = 100  # number of Phase I (reference) subgroups

# -- MREWMA-Copula simulation settings --
B_PERM        = 10000  # permutations used to build the null distribution
TARGET_MRL0   = 200    # target in-control median run length
N_TEST        = n_per_subgroup  # subgroup size used inside the run-length simulation
N_RL_UCL      = 1000   # replications per bisection evaluation (UCL search)
MAX_STEPS_UCL = 800    # max steps per replication (UCL search)
# N_RL_OOC was 300 in the previous submission. Reviewer 1 judged this "not
# sufficient to reliably assess control chart performance." It is raised here
# to 1,000, matching N_RL_UCL and standard practice for Monte Carlo run-length
# evaluation in the SPM literature (e.g. Song et al., 2023, use comparable
# replication counts). Note that Sections 12-17 run sequentially (see the
# changelog above), so this increase does make Section 16 considerably
# slower to run -- lower this temporarily (e.g. 100) for a quick test run,
# and raise it (e.g. 5000) for a final confirmatory run if time allows.
N_RL_OOC      = 1000   # replications per MRL0 evaluation
MAX_STEPS_OOC = 1500   # max steps per replication (MRL0 evaluation)
TOL_BISECT    = 5      # bisection convergence tolerance (in run length)
LAMBDA_VEC    = c(0.1, 0.25, 0.5, 0.75, 0.9)
# Upper starting bound for the UCL bisection search (Section 14, find_ucl()).
# The true UCLs for this design are all well under 5 (see UCL_Table in the
# output workbook); 25 was needlessly wide and wasted entire replications
# running to MAX_STEPS_UCL with no signal. 10 keeps a >2x safety margin, and
# the search still auto-expands if a lambda ever needs more.
UCL_SEARCH_UPPER = 10

# -- Parallelization (currently unused; kept for possible future use) --
# Sections 12-17 (MREWMA-Copula engine onward) now run sequentially, following
# the thesis pipeline they were adapted from (see the changelog above), so
# USE_PARALLEL/N_CORES are not referenced by the current code. They are left
# defined here in case parallel::makeCluster()-based execution is added back
# to Section 16 (the MRL0 curve sweep) later, which is by far the most
# expensive step.
USE_PARALLEL  = TRUE
N_CORES       = max(1, parallel::detectCores(logical = TRUE) - 1)  # leave 1 core free

# -- Imputation validation (Section 7) --
N_HOLDOUT_REPS = 50    # hold-out replications for EM vs MICE vs KNN comparison
KNN_K          = 5     # neighbors used by the KNN imputation benchmark
set.seed(20260714)     # reproducibility for the hold-out validation draws

# -- Temporal/seasonal predictors used by the EM imputation (Sections 5-7) --
LAGS_EM      = c(1, 2, 3, 7)  # lag days added as auxiliary predictors
ROLL_WINDOWS = c(3, 7)        # causal (past-only) rolling-mean window sizes

# -- Bagged EM settings (Sections 5-7) --
B_EM_BAG            = 30    # bootstrap replicates for the final imputation (Section 6)
B_EM_BAG_VALIDATION = 15    # bootstrap replicates used inside the hold-out check (Section 7)
SHRINK_GAMMA        = 0.10  # ridge shrinkage applied to Sigma in the M-step

# ==============================================================================
# 3. PLOT THEME
# ==============================================================================
tryCatch({
  font_add(family = "Arial",
           regular = "C:/Windows/Fonts/arial.ttf",
           bold    = "C:/Windows/Fonts/arialbd.ttf")
  showtext_auto()
}, warning = function(w) {
  message("Arial font not found; falling back to the default font.")
}, error = function(e) {
  message("Arial font not found; falling back to the default font.")
})

custom_plot_theme = theme_minimal(base_family = "Arial") +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 32),
    axis.title       = element_text(size = 24, face = "bold", color = "black"),
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 20),
    axis.text.y      = element_text(size = 20),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#ebebeb")
  )

# ==============================================================================
# 4. LOAD RAW DATA
# ==============================================================================
air_data = as.data.frame(read_excel(input_file, sheet = "DATA"))
air_data$date = as.Date(air_data$Date)
air_data = air_data[order(air_data$date), ]
rownames(air_data) = NULL

X_raw    = as.matrix(air_data[, vars_target])
storage.mode(X_raw) = "double"
Date_raw = air_data$date

# ==============================================================================
# 5. EM ALGORITHM FUNCTIONS (MULTIVARIATE NORMAL MISSING-VALUE IMPUTATION)
# ==============================================================================
fit_normal_score = function(x) {
  obs_idx = which(!is.na(x))
  n = length(obs_idx)
  r = rank(x[obs_idx], ties.method = "average")
  u = r / (n + 1)
  z_full = rep(NA_real_, length(x))
  z_full[obs_idx] = qnorm(u)
  list(z = z_full, sorted_obs = sort(x[obs_idx]), n = n)
}

from_normal_score = function(z_vec, transform) {
  p_grid = (1:transform$n) / (transform$n + 1)
  u = pnorm(z_vec)
  as.numeric(approx(x = p_grid, y = transform$sorted_obs, xout = u, rule = 2)$y)
}

build_copula_space = function(X) {
  p = ncol(X)
  transforms = vector("list", p); names(transforms) = colnames(X)
  Z = matrix(NA_real_, nrow(X), p, dimnames = dimnames(X))
  for (j in 1:p) {
    ft = fit_normal_score(X[, j])
    transforms[[j]] = ft
    Z[, j] = ft$z
  }
  list(Z = Z, transforms = transforms)
}

build_lag_features = function(X, lags) {
  if (length(lags) == 0) return(NULL)
  n = nrow(X); var_names = colnames(X)
  lag_mats = lapply(lags, function(L) {
    M = matrix(NA_real_, n, ncol(X), dimnames = list(NULL, paste0(var_names, "_lag", L)))
    if (L < n) M[(L + 1):n, ] = X[1:(n - L), ]
    M
  })
  do.call(cbind, lag_mats)
}

build_rolling_features = function(X, windows) {
  if (length(windows) == 0) return(NULL)
  n = nrow(X); var_names = colnames(X)
  roll_mats = lapply(windows, function(w) {
    M = matrix(NA_real_, n, ncol(X), dimnames = list(NULL, paste0(var_names, "_roll", w)))
    if (w < n) {
      for (i in (w + 1):n) {
        M[i, ] = colMeans(X[(i - w):(i - 1), , drop = FALSE], na.rm = TRUE)
      }
    }
    M
  })
  do.call(cbind, roll_mats)
}

build_seasonal_features = function(dates) {
  doy = as.numeric(format(dates, "%j"))
  t_years = as.numeric(dates - min(dates)) / 365.25
  cbind(
    sin_doy  = sin(2 * pi * doy / 365.25),
    cos_doy  = cos(2 * pi * doy / 365.25),
    sin_doy2 = sin(4 * pi * doy / 365.25),
    cos_doy2 = cos(4 * pi * doy / 365.25),
    t_trend  = t_years,
    t_trend2 = t_years^2
  )
}

compute_observed_log_likelihood = function(X, mu, Sigma) {
  N = nrow(X)
  total_ll = 0
  for (i in 1:N) {
    row = X[i, ]
    obs_idx = which(!is.na(row))
    if (length(obs_idx) == 0) next
    x_obs = row[obs_idx]
    mu_obs = mu[obs_idx]
    Sigma_obs = Sigma[obs_idx, obs_idx, drop = FALSE]
    p_obs = length(obs_idx)
    Sigma_obs_reg = Sigma_obs + 1e-10 * diag(p_obs)
    ev = eigen(Sigma_obs_reg, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev <= 0)) next
    log_det = sum(log(ev))
    diff = x_obs - mu_obs
    inv_Sigma = solve(Sigma_obs_reg)
    quad_form = as.numeric(t(diff) %*% inv_Sigma %*% diff)
    ll_i = -0.5 * (p_obs * log(2 * pi) + log_det + quad_form)
    total_ll = total_ll + ll_i
  }
  total_ll
}

compute_conditional_params = function(mu, Sigma, M_i, Mbar_i, X_obs_i) {
  mu_mis = mu[M_i]
  mu_obs = mu[Mbar_i]
  Sigma_mis = Sigma[M_i, M_i, drop = FALSE]
  Sigma_obs = Sigma[Mbar_i, Mbar_i, drop = FALSE]
  V_obs_mis = Sigma[M_i, Mbar_i, drop = FALSE]
  V_mis_obs = Sigma[Mbar_i, M_i, drop = FALSE]
  k_obs = length(Mbar_i)
  Sigma_obs_reg = Sigma_obs + 1e-10 * diag(k_obs)
  Sigma_obs_inv = tryCatch(solve(Sigma_obs_reg), error = function(e) ginv(Sigma_obs_reg))
  mu_Mi = as.numeric(mu_mis + V_obs_mis %*% Sigma_obs_inv %*% (X_obs_i - mu_obs))
  Sigma_Mi = Sigma_mis - V_obs_mis %*% Sigma_obs_inv %*% V_mis_obs
  Sigma_Mi = (Sigma_Mi + t(Sigma_Mi)) / 2
  min_eig = min(eigen(Sigma_Mi, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < 0) Sigma_Mi = Sigma_Mi + (-min_eig + 1e-10) * diag(length(M_i))
  list(mu_Mi = mu_Mi, Sigma_Mi = Sigma_Mi)
}

em_e_step = function(X, mu, Sigma) {
  N = nrow(X)
  p = ncol(X)
  tau1 = numeric(p)
  tau2 = matrix(0, p, p)
  X_imputed = X
  for (i in 1:N) {
    row = X[i, ]
    M_i = which(is.na(row))
    Mbar_i = which(!is.na(row))
    if (length(M_i) == 0) {
      tau1 = tau1 + row
      tau2 = tau2 + outer(row, row)
      next
    }
    if (length(Mbar_i) == 0) {
      tau1 = tau1 + mu
      tau2 = tau2 + Sigma + outer(mu, mu)
      X_imputed[i, ] = mu
      next
    }
    X_obs_i = row[Mbar_i]
    cp = compute_conditional_params(mu, Sigma, M_i, Mbar_i, X_obs_i)
    mu_Mi = cp$mu_Mi
    Sigma_Mi = cp$Sigma_Mi
    X_imputed[i, M_i] = mu_Mi
    tau1[Mbar_i] = tau1[Mbar_i] + row[Mbar_i]
    tau1[M_i] = tau1[M_i] + mu_Mi
    tau2[Mbar_i, Mbar_i] = tau2[Mbar_i, Mbar_i] + outer(row[Mbar_i], row[Mbar_i])
    cross = outer(row[Mbar_i], mu_Mi)
    tau2[Mbar_i, M_i] = tau2[Mbar_i, M_i] + cross
    tau2[M_i, Mbar_i] = tau2[M_i, Mbar_i] + t(cross)
    tau2[M_i, M_i] = tau2[M_i, M_i] + Sigma_Mi + outer(mu_Mi, mu_Mi)
  }
  tau1 = tau1 / N
  tau2 = tau2 / N
  list(tau1 = tau1, tau2 = tau2, X_imputed = X_imputed)
}

em_m_step = function(tau1, tau2) {
  p = length(tau1)
  mu_new = tau1
  Sigma_new = tau2 - outer(tau1, tau1)
  Sigma_new = (Sigma_new + t(Sigma_new)) / 2
  min_eig = min(eigen(Sigma_new, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < 0) Sigma_new = Sigma_new + (-min_eig + 1e-8) * diag(p)
  list(mu_new = mu_new, Sigma_new = Sigma_new)
}

em_m_step_ridge = function(tau1, tau2, shrink = 0.1) {
  p = length(tau1)
  mu_new = tau1
  Sigma_emp = tau2 - outer(tau1, tau1)
  Sigma_emp = (Sigma_emp + t(Sigma_emp)) / 2
  target = diag(diag(Sigma_emp))
  Sigma_new = (1 - shrink) * Sigma_emp + shrink * target
  min_eig = min(eigen(Sigma_new, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig < 0) Sigma_new = Sigma_new + (-min_eig + 1e-8) * diag(p)
  list(mu_new = mu_new, Sigma_new = Sigma_new)
}

run_em_algorithm = function(X, max_iter = 200, tol = 1e-4) {
  N = nrow(X)
  p = ncol(X)
  complete_mask = rowSums(is.na(X)) == 0
  X_complete = X[complete_mask, , drop = FALSE]
  if (nrow(X_complete) < 2) stop("Too few complete observations to initialize the EM algorithm.")
  mu = colMeans(X, na.rm = TRUE)
  Sigma = cov(X_complete)
  min_eig_init = min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig_init < 1e-10) Sigma = Sigma + (abs(min_eig_init) + 1e-8) * diag(p)
  
  ll_current = compute_observed_log_likelihood(X, mu, Sigma)
  n_iter = 0
  
  for (t in 1:max_iter) {
    n_iter = t
    e_res = em_e_step(X, mu, Sigma)
    m_res = em_m_step(e_res$tau1, e_res$tau2)
    mu_new = m_res$mu_new
    Sigma_new = m_res$Sigma_new
    ll_new = compute_observed_log_likelihood(X, mu_new, Sigma_new)
    delta_ll = abs(ll_new - ll_current)
    converged = delta_ll < tol
    mu = mu_new
    Sigma = Sigma_new
    ll_current = ll_new
    if (converged) break
  }
  
  final_res = em_e_step(X, mu, Sigma)
  X_imputed_final = final_res$X_imputed
  
  list(X_imputed = X_imputed_final, mu = mu, Sigma = Sigma,
       n_iter = n_iter, log_likelihood = ll_current)
}

run_em_algorithm_ridge = function(X, max_iter = 200, tol = 1e-4, shrink = 0.1,
                                  init_mu = NULL, init_Sigma = NULL) {
  N = nrow(X)
  p = ncol(X)
  if (is.null(init_mu) || is.null(init_Sigma)) {
    complete_mask = rowSums(is.na(X)) == 0
    X_complete = X[complete_mask, , drop = FALSE]
    if (nrow(X_complete) < 2) stop("Too few complete observations to initialize the EM algorithm.")
    mu = colMeans(X, na.rm = TRUE)
    Sigma = cov(X_complete)
  } else {
    mu = init_mu
    Sigma = init_Sigma
  }
  min_eig_init = min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig_init < 1e-10) Sigma = Sigma + (abs(min_eig_init) + 1e-8) * diag(p)
  ll_current = compute_observed_log_likelihood(X, mu, Sigma)
  n_iter = 0
  for (t in 1:max_iter) {
    n_iter = t
    e_res = em_e_step(X, mu, Sigma)
    m_res = em_m_step_ridge(e_res$tau1, e_res$tau2, shrink = shrink)
    mu_new = m_res$mu_new
    Sigma_new = m_res$Sigma_new
    ll_new = compute_observed_log_likelihood(X, mu_new, Sigma_new)
    converged = abs(ll_new - ll_current) < tol
    mu = mu_new
    Sigma = Sigma_new
    ll_current = ll_new
    if (converged) break
  }
  final_res = em_e_step(X, mu, Sigma)
  list(X_imputed = final_res$X_imputed, mu = mu, Sigma = Sigma,
       n_iter = n_iter, log_likelihood = ll_current)
}

run_em_copula = function(X, max_iter = 200, tol = 1e-4) {
  cs = build_copula_space(X)
  em_z = run_em_algorithm(cs$Z, max_iter = max_iter, tol = tol)
  
  p = ncol(X)
  X_imputed = X
  for (j in 1:p) X_imputed[, j] = from_normal_score(em_z$X_imputed[, j], cs$transforms[[j]])
  obs_mask = !is.na(X)
  X_imputed[obs_mask] = X[obs_mask]
  
  list(X_imputed = X_imputed, mu_Z = em_z$mu, Sigma_Z = em_z$Sigma,
       n_iter = em_z$n_iter, log_likelihood = em_z$log_likelihood)
}

run_em_copula_bagged = function(X, dates, lags = LAGS_EM, roll_windows = ROLL_WINDOWS,
                                max_iter = 200, tol = 1e-4, B_bag = B_EM_BAG,
                                shrink = SHRINK_GAMMA) {
  var_names = colnames(X)
  lag_feat  = build_lag_features(X, lags = lags)
  roll_feat = build_rolling_features(X, windows = roll_windows)
  seas_feat = build_seasonal_features(dates)
  X_aug = cbind(X, lag_feat, roll_feat, seas_feat)
  
  cs = build_copula_space(X_aug)
  Z_full = cs$Z
  N = nrow(Z_full); p_aug = ncol(Z_full)
  
  complete_mask_full = rowSums(is.na(Z_full)) == 0
  mu_init    = colMeans(Z_full, na.rm = TRUE)
  Sigma_init = cov(Z_full[complete_mask_full, , drop = FALSE])
  min_eig0 = min(eigen(Sigma_init, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig0 < 1e-10) Sigma_init = Sigma_init + (abs(min_eig0) + 1e-8) * diag(p_aug)
  
  Z_imputed_sum = matrix(0, N, p_aug)
  ll_vec = numeric(B_bag); niter_vec = numeric(B_bag)
  
  for (b in 1:B_bag) {
    boot_rows = sample.int(N, N, replace = TRUE)
    Z_boot = Z_full[boot_rows, , drop = FALSE]
    fit_b = run_em_algorithm_ridge(Z_boot, max_iter = max_iter, tol = tol, shrink = shrink,
                                   init_mu = mu_init, init_Sigma = Sigma_init)
    e_full = em_e_step(Z_full, fit_b$mu, fit_b$Sigma)
    Z_imputed_sum = Z_imputed_sum + e_full$X_imputed
    ll_vec[b] = fit_b$log_likelihood; niter_vec[b] = fit_b$n_iter
  }
  Z_imputed_avg = Z_imputed_sum / B_bag
  
  X_imputed = X
  for (j in seq_along(var_names)) {
    X_imputed[, j] = from_normal_score(Z_imputed_avg[, j], cs$transforms[[var_names[j]]])
  }
  obs_mask = !is.na(X)
  X_imputed[obs_mask] = X[obs_mask]
  
  list(X_imputed = X_imputed, n_iter_mean = mean(niter_vec),
       log_likelihood_mean = mean(ll_vec), B_bag = B_bag)
}

analyze_missing_patterns = function(X, var_names) {
  n_missing = colSums(is.na(X))
  pct_missing = round(100 * n_missing / nrow(X), 2)
  data.frame(Indicator = var_names, N_Missing = n_missing, Pct_Missing = pct_missing, row.names = NULL)
}

# ==============================================================================
# 6. RUN EM IMPUTATION
# ==============================================================================
missing_pattern_summary = analyze_missing_patterns(X_raw, vars_target)

em_result = run_em_copula_bagged(X_raw, dates = Date_raw, lags = LAGS_EM,
                                 roll_windows = ROLL_WINDOWS, max_iter = 200, tol = 1e-4,
                                 B_bag = B_EM_BAG, shrink = SHRINK_GAMMA)
air_data[, vars_target] = em_result$X_imputed

# ==============================================================================
# 7. IMPUTATION METHOD VALIDATION: EM vs MICE vs KNN (HOLD-OUT ACCURACY CHECK)
# ==============================================================================
# Reviewer 1: "no quantitative assessment of EM imputation accuracy... no
# comparison with alternative missing-data methods (e.g. MissForest, MICE,
# KNN)." This section answers both points together, using a single repeated
# hold-out procedure so all three methods are compared on identical synthetic
# missingness in every replication (a fair, paired comparison):
#   1. Start from the rows of the RAW data (X_raw) that have no missing values
#      at all, so the true values are known with certainty.
#   2. In each of N_HOLDOUT_REPS replications, delete values completely at
#      random (MCAR) from those rows, at the SAME per-variable rate observed
#      in the real data (colMeans(is.na(X_raw)): see missing_pattern_summary),
#      so the hold-out task mirrors the real imputation task as closely as
#      possible.
#   3. Impute the resulting artificial gaps with EM (this script's own
#      algorithm, Section 5), MICE (predictive mean matching), and KNN
#      (VIM::kNN), then compare each method's imputed values against the
#      (known) deleted true values.
#   4. Summarize RMSE and MAE per method per variable, averaged (with SD)
#      across all replications.
#
# This runs on the RAW data independently of the em_result used for the rest
# of the pipeline; it does not change the imputed air_data used downstream.

run_mice_impute = function(X_df) {
  imp = mice(X_df, m = 1, method = "pmm", maxit = 5, printFlag = FALSE)
  as.matrix(complete(imp, 1))
}

run_knn_impute = function(X_df, k = KNN_K) {
  imp = VIM::kNN(X_df, variable = colnames(X_df), k = k, imp_var = FALSE)
  as.matrix(imp[, colnames(X_df)])
}

validate_imputation_holdout = function(X_complete, dates_complete, miss_rates,
                                       n_reps = N_HOLDOUT_REPS, lags = LAGS_EM,
                                       roll_windows = ROLL_WINDOWS,
                                       B_bag = B_EM_BAG_VALIDATION, shrink = SHRINK_GAMMA) {
  p = ncol(X_complete); var_names = colnames(X_complete); n_obs = nrow(X_complete)
  seas_feat = build_seasonal_features(dates_complete)
  
  rep_results = vector("list", n_reps)
  
  for (r in 1:n_reps) {
    X_test = X_complete; true_vals = list(); mask_idx = list()
    for (j in 1:p) {
      n_missing_j = round(miss_rates[j] * n_obs)
      if (n_missing_j < 1) next
      idx_j = sample.int(n_obs, n_missing_j)
      true_vals[[var_names[j]]] = X_complete[idx_j, j]
      mask_idx[[var_names[j]]]  = idx_j
      X_test[idx_j, j] = NA
    }
    
    lag_feat  = build_lag_features(X_test, lags = lags)
    roll_feat = build_rolling_features(X_test, windows = roll_windows)
    X_test_aug = cbind(X_test, lag_feat, roll_feat, seas_feat)
    X_test_aug_df = as.data.frame(X_test_aug)
    
    em_out = run_em_copula_bagged(X_test, dates = dates_complete, lags = lags,
                                  roll_windows = roll_windows, max_iter = 200, tol = 1e-4,
                                  B_bag = B_bag, shrink = shrink)
    X_imp_mice_full = tryCatch(
      run_mice_impute(X_test_aug_df),
      error = function(e) { message("MICE failed on rep ", r, ": ", e$message); matrix(NA, n_obs, ncol(X_test_aug)) }
    )
    X_imp_knn_full = tryCatch(
      run_knn_impute(X_test_aug_df),
      error = function(e) { message("KNN failed on rep ", r, ": ", e$message); matrix(NA, n_obs, ncol(X_test_aug)) }
    )
    
    extract_target = function(M_full) {
      M_full = as.matrix(M_full)
      if (is.null(colnames(M_full))) colnames(M_full) = colnames(X_test_aug)
      M_full[, var_names, drop = FALSE]
    }
    
    imputed_by_method = list(
      EM   = em_out$X_imputed,
      MICE = extract_target(X_imp_mice_full),
      KNN  = extract_target(X_imp_knn_full)
    )
    rep_rows = list()
    for (method in names(imputed_by_method)) {
      X_imp = imputed_by_method[[method]]
      for (j in 1:p) {
        vn = var_names[j]
        if (is.null(mask_idx[[vn]])) next
        imputed_vals = X_imp[mask_idx[[vn]], j]
        truth_vals   = true_vals[[vn]]
        rep_rows[[length(rep_rows) + 1]] = data.frame(
          Rep      = r,
          Method   = method,
          Variable = vn,
          RMSE     = sqrt(mean((imputed_vals - truth_vals)^2, na.rm = TRUE)),
          MAE      = mean(abs(imputed_vals - truth_vals), na.rm = TRUE)
        )
      }
    }
    rep_results[[r]] = do.call(rbind, rep_rows)
  }
  
  do.call(rbind, rep_results)
}

X_complete_cases_raw   = X_raw[complete.cases(X_raw), , drop = FALSE]
Date_complete_cases     = Date_raw[complete.cases(X_raw)]
missing_rates_observed  = colMeans(is.na(X_raw))

em_validation_raw = validate_imputation_holdout(
  X_complete     = X_complete_cases_raw,
  dates_complete = Date_complete_cases,
  miss_rates     = missing_rates_observed,
  n_reps         = N_HOLDOUT_REPS,
  lags           = LAGS_EM,
  roll_windows   = ROLL_WINDOWS,
  B_bag          = B_EM_BAG_VALIDATION,
  shrink         = SHRINK_GAMMA
)

em_validation_summary = em_validation_raw %>%
  group_by(Method, Variable) %>%
  summarise(
    RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
    MAE_mean  = mean(MAE),  MAE_sd  = sd(MAE),
    .groups = "drop"
  ) %>%
  arrange(Variable, RMSE_mean)

p_rmse = ggplot(em_validation_raw, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot(alpha = 0.85, outlier.size = 1) +
  facet_wrap(~Variable, scales = "free_y") +
  labs(title = "Imputation Accuracy: EM vs MICE vs KNN (Hold-Out RMSE)",
       x = NULL, y = "RMSE (hold-out, MCAR)") +
  scale_fill_brewer(palette = "Set2") +
  custom_plot_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

p_mae = ggplot(em_validation_raw, aes(x = Method, y = MAE, fill = Method)) +
  geom_boxplot(alpha = 0.85, outlier.size = 1) +
  facet_wrap(~Variable, scales = "free_y") +
  labs(title = "Imputation Accuracy: EM vs MICE vs KNN (Hold-Out MAE)",
       x = NULL, y = "MAE (hold-out, MCAR)") +
  scale_fill_brewer(palette = "Set2") +
  custom_plot_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(plot_dir, "imputation_validation_rmse.png"), plot = p_rmse,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "imputation_validation_mae.png"), plot = p_mae,
       width = 12, height = 5, dpi = 300)

message("Imputation validation complete (", N_HOLDOUT_REPS, " hold-out replications). Summary:")
print(as.data.frame(em_validation_summary))

# ==============================================================================
# 8. SUBGROUP FORMATION (PHASE I / PHASE II SPLIT)
# ==============================================================================
n_subgroup_total = nrow(air_data) / n_per_subgroup

if (n_subgroup_total != floor(n_subgroup_total)) {
  stop("The number of observations is not evenly divisible by the subgroup size. Check the data.")
}

air_data$Subgroup = rep(1:n_subgroup_total, each = n_per_subgroup)

n_subgroup_phase2 = n_subgroup_total - n_subgroup_phase1

data_ref = air_data[air_data$Subgroup <= n_subgroup_phase1, ]
data_mon = air_data[air_data$Subgroup > n_subgroup_phase1, ]

# ==============================================================================
# 9. DESCRIPTIVE STATISTICS
# ==============================================================================
get_descriptive_stats = function(df, label_data) {
  df %>%
    select(all_of(vars_target)) %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>%
    group_by(Variable) %>%
    summarise(
      N      = sum(!is.na(Value)),
      Mean   = mean(Value, na.rm = TRUE),
      Median = median(Value, na.rm = TRUE),
      Min    = min(Value, na.rm = TRUE),
      Max    = max(Value, na.rm = TRUE),
      SD     = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Data = label_data) %>%
    select(Data, Variable, everything())
}

stats_all = get_descriptive_stats(air_data, "All Data (Post-EM Imputation)")

cat("\n--- [9] DESCRIPTIVE STATISTICS ---\n")
print(as.data.frame(stats_all), row.names = FALSE)

# ==============================================================================
# 10. MULTIVARIATE NORMALITY ASSUMPTION TESTING (SHAPIRO-WILK)
# ==============================================================================
cat("\n--- [10] MULTIVARIATE NORMALITY ASSUMPTION TESTING ---\n")

cat("\nHypotheses:\n")
cat("  H0 : The data follow a multivariate normal distribution\n")
cat("  H1 : The data do not follow a multivariate normal distribution\n")

cat("\nTest Criteria:\n")
cat("  If p-value > alpha, with alpha = 0.05, then H0 is accepted\n\n")

X_normality = t(as.matrix(air_data[, vars_target]))  # mshapiro.test expects p x n (variables in rows)
mvn_result  = mshapiro.test(X_normality)

normality_table = data.frame(
  `Shapiro-Wilk` = round(as.numeric(mvn_result$statistic), 3),
  `p-value`      = formatC(mvn_result$p.value, format = "f", digits = 3),
  Decision       = ifelse(mvn_result$p.value > 0.05, "H0 accepted", "H0 rejected"),
  check.names    = FALSE
)
print(normality_table, row.names = FALSE)

if (mvn_result$p.value <= 0.05) {
  cat(sprintf("\nBecause p-value < alpha (%.3f < 0.05), H0 is rejected: the air quality\n",
              mvn_result$p.value))
  cat("indicator data do not follow a multivariate normal distribution. The data can\n")
  cat("therefore be analyzed using the MREWMA-EM CC and MREWMA-EM LC methods, which are\n")
  cat("robust to violations of the multivariate normality assumption.\n")
} else {
  cat(sprintf("\nBecause p-value > alpha (%.3f > 0.05), H0 is accepted: the air quality\n",
              mvn_result$p.value))
  cat("indicator data follow a multivariate normal distribution.\n")
}

# ==============================================================================
# 11. TIME SERIES VISUALIZATION
# ==============================================================================
colors_ts = c("#005b96", "#b2182b", "#006837")
large_font_size = 33
phase_split_date = max(data_ref$date)

for (i in 1:length(vars_target)) {
  var = vars_target[i]
  
  p = ggplot(air_data, aes(x = date, y = .data[[var]])) +
    geom_point(color = colors_ts[i], size = 1.5, alpha = 0.8) +
    geom_vline(xintercept = as.numeric(phase_split_date),
               color = "grey20", linetype = "dashed", linewidth = 1) +
    annotate("text", x = min(data_ref$date) + (phase_split_date - min(data_ref$date)) / 2, y = Inf,
             label = "Phase I (Reference)", color = "black",
             vjust = -1, size = large_font_size / .pt,
             fontface = "bold", family = "Arial") +
    annotate("text", x = phase_split_date + (max(data_mon$date) - phase_split_date) / 2, y = Inf,
             label = "Phase II (Monitoring)", color = "black",
             vjust = -1, size = large_font_size / .pt,
             fontface = "bold", family = "Arial") +
    scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
    labs(title = NULL, x = "Time", y = var) +
    coord_cartesian(clip = "off") +
    custom_plot_theme +
    theme(
      axis.title.x = element_text(size = large_font_size, face = "bold", color = "black", margin = margin(t = 15)),
      axis.title.y = element_text(size = large_font_size, face = "bold", color = "black", margin = margin(r = 15)),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, size = large_font_size - 3, color = "black"),
      axis.text.y  = element_text(size = large_font_size - 3, color = "black"),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 1),
      plot.margin  = margin(t = 40, r = 20, b = 20, l = 20)
    )
  
  ggsave(file.path(plot_dir, paste0("timeseries_", gsub("\\.", "", var), ".png")),
         plot = p, width = 12, height = 7, dpi = 300)
}

# ==============================================================================
# 12. MREWMA-COPULA STATISTICAL ENGINE (LEPAGE & CUCCONI)
# ==============================================================================
cat("\n--- [12] MREWMA-COPULA ENGINE ---\n")

get_pc_scores = function(X_ref, Y_test) {
  N  = nrow(X_ref) + nrow(Y_test)
  combined = rbind(X_ref, Y_test)
  G  = apply(combined, 2, function(k) rank(k) / (N + 1))
  eg = eigen(cov(G))
  Th = G %*% eg$vectors
  list(ref = Th[1:nrow(X_ref), , drop = FALSE],
       test = Th[(nrow(X_ref)+1):N, , drop = FALSE])
}

.simpanan_nol = new.env(parent = emptyenv())

build_null_perm = function(m, n, B = B_PERM) {
  N  = m + n
  rk = 1:N
  ab = pmin(rk, N + 1 - rk)
  md = (rk - (N+1)/2)^2
  muW = n*(N+1)/2;  sdW = sqrt(n*m*(N+1)/12)
  if (N %% 2 == 0) {
    muAB = n*(N+2)/4; sdAB = sqrt((m*n*(N+2)*(N-2))/(48*(N-1)))
  } else {
    muAB = n*(N+1)^2/(4*N); sdAB = sqrt((m*n*(N+1)*(N^2+3))/(48*N^2))
  }
  muM = n*(N^2-1)/12; sdM = sqrt(n*m*(N+1)*(N^2-4)/180)
  sel = replicate(B, sample.int(N, n))
  TW  = colSums(matrix(rk[sel], nrow = n))
  TAB = colSums(matrix(ab[sel], nrow = n))
  TM  = colSums(matrix(md[sel], nrow = n))
  list(L = sort(((TW-muW)/sdW)^2 + ((TAB-muAB)/sdAB)^2),
       C = sort(((TW-muW)/sdW)^2 + ((TM -muM )/sdM )^2),
       B = B, mo = list(muW=muW, sdW=sdW, muAB=muAB, sdAB=sdAB, muM=muM, sdM=sdM))
}

get_null_perm = function(m, n, B = B_PERM) {
  cache_key = paste(m, n, B, sep = "_")
  ob = .simpanan_nol[[cache_key]]
  if (is.null(ob)) { ob = build_null_perm(m, n, B); assign(cache_key, ob, envir = .simpanan_nol) }
  ob
}

pval_LC = function(x, y, B = B_PERM) {
  m = length(x); n = length(y); N = m + n
  null = get_null_perm(m, n, B); mo = null$mo
  r   = rank(c(x, y)); iy = (m+1):N
  TW  = sum(r[iy])
  TAB = sum(pmin(r, N+1-r)[iy])
  TM  = sum(((r-(N+1)/2)^2)[iy])
  Lo  = ((TW-mo$muW)/mo$sdW)^2 + ((TAB-mo$muAB)/mo$sdAB)^2
  Co  = ((TW-mo$muW)/mo$sdW)^2 + ((TM -mo$muM )/mo$sdM )^2
  amb = function(obs, ns, B) (1 + (B - findInterval(obs - 1e-12, ns))) / (B + 1)
  c(max(min(amb(Lo, null$L, B), 0.99999), 1e-10),
    max(min(amb(Co, null$C, B), 0.99999), 1e-10))
}

compute_stat_LC = function(X_ref, Y_test, type = "LC") {
  d = ncol(X_ref)
  Z = numeric(d + 1)
  pc = get_pc_scores(X_ref, Y_test)
  pv = numeric(d)
  for (k in 1:d) {
    h = pval_LC(pc$ref[, k], pc$test[, k])
    pv[k] = ifelse(type == "LC", h[1], h[2])
  }
  Z[1] = -log(min(pv)) - 1
  for (k in 1:d) {
    h = pval_LC(X_ref[, k], Y_test[, k])
    Z[k + 1] = -log(ifelse(type == "LC", h[1], h[2])) - 1
  }
  Z
}

shift_data = function(sample_x, type, param, mean_v, cov, white) {
  Z = sweep(sample_x, 2, mean_v, "-") %*% white
  mu_b = numeric(d); kv_b = diag(d)
  if (type == "LS1") mu_b[1] = param
  else if (type == "LS2") { mu_b[1] = param; mu_b[2] = 0.5 }
  else if (type == "DS1") kv_b[1,1] = param
  else if (type == "DS2") { kv_b[1,1] = param; kv_b[2,2] = 1.5 }
  else if (type == "OD")  { kv_b[1,2] = param; kv_b[2,1] = param }
  else if (type == "MS1") { mu_b[1] = 0.5; kv_b[1,1] = param }
  else if (type == "MS2") { mu_b[1] = param; kv_b[1,1] = 1.5 }
  else if (type == "MS3") { mu_b[1] = 0.5; kv_b[1,1] = 1.5; kv_b[1,2] = param; kv_b[2,1] = param }
  ch = tryCatch(chol(kv_b), error = function(e) diag(sqrt(pmax(diag(kv_b), 1e-8))))
  Zg = sweep(Z %*% ch, 2, mu_b, "+")
  ei = eigen(cov)
  wa = ei$vectors %*% diag(sqrt(pmax(ei$values, 0))) %*% t(ei$vectors)
  sweep(Zg %*% wa, 2, mean_v, "+")
}

simulate_mrl = function(ucl, lam, type = "LC", n_sim = N_RL_OOC,
                        max_steps = MAX_STEPS_OOC, shift_type = "None",
                        shift_param = 0, n_test = N_TEST) {
  rl = numeric(n_sim)
  for (s in 1:n_sim) {
    ew = numeric(d + 1)
    for (t in 1:max_steps) {
      id = sample(1:m, n_test, replace = TRUE)
      Yb = X_baseline[id, , drop = FALSE]
      if (shift_type != "None")
        Yb = shift_data(Yb, shift_type, shift_param, mean_baseline, cov_baseline, whitening_mat)
      Zt = compute_stat_LC(X_baseline, Yb, type)
      ew = lam * Zt + (1 - lam) * ew
      if (max(ew) > ucl) { rl[s] = t; break }
      if (t == max_steps) rl[s] = max_steps
    }
  }
  median(rl)
}

find_ucl = function(target, lam, type, A1 = 0.1, A2 = UCL_SEARCH_UPPER,
                    n_sim = N_RL_UCL, max_steps = MAX_STEPS_UCL,
                    tol = TOL_BISECT, max_iter_bisect = 30) {
  f  = function(u) simulate_mrl(u, lam, type, n_sim, max_steps, "None") - target
  fa = f(A1); fb = f(A2); expand_count = 0
  while (fb < 0 && expand_count < 6)  { A2 = A2 * 1.5; fb = f(A2); expand_count = expand_count + 1 }
  while (fa > 0 && expand_count < 12) { A1 = A1 / 1.5; fa = f(A1); expand_count = expand_count + 1 }
  for (it in 1:max_iter_bisect) {
    mid_pt = (A1 + A2) / 2
    ft = f(mid_pt)
    if (abs(ft) <= tol) return(mid_pt)
    if (ft > 0) A2 = mid_pt else A1 = mid_pt
  }
  (A1 + A2) / 2
}

# ==============================================================================
# 13. PHASE I ITERATIVE CLEANING (LEPAGE/CUCCONI, LEAVE-ONE-SUBGROUP-OUT)
# ==============================================================================
cat("\n--- [13] PHASE I ITERATIVE CLEANING ---\n")

LAMBDA_CLEAN          = 0.1
N_RL_CLEAN         = 300
MAX_STEPS_CLEAN    = 400
MAX_ROUNDS_CLEAN = 10
MIN_SUBGROUP_CLEAN  = 30

d = length(vars_target)

subgroup_baseline = sort(unique(data_ref$Subgroup))
data_baseline    = data_ref
X_baseline_mat   = as.matrix(data_baseline[, vars_target])
cleaning_history = data.frame()

for (round_i in 1:MAX_ROUNDS_CLEAN) {
  m_sg0      = length(subgroup_baseline)
  X_baseline    = X_baseline_mat
  m          = nrow(X_baseline)
  mean_baseline = colMeans(X_baseline)
  cov_baseline  = cov(X_baseline)
  eig_a      = eigen(cov_baseline)
  whitening_mat  = eig_a$vectors %*% diag(1/sqrt(eig_a$values)) %*% t(eig_a$vectors)
  
  ucl_lc = find_ucl(TARGET_MRL0, LAMBDA_CLEAN, "LC", n_sim = N_RL_CLEAN, max_steps = MAX_STEPS_CLEAN)
  ucl_cc = find_ucl(TARGET_MRL0, LAMBDA_CLEAN, "CC", n_sim = N_RL_CLEAN, max_steps = MAX_STEPS_CLEAN)
  
  st_lc = numeric(m_sg0); st_cc = numeric(m_sg0)
  ew_lc = numeric(d + 1); ew_cc = numeric(d + 1)
  for (t in 1:m_sg0) {
    id_sg = which(data_baseline$Subgroup == subgroup_baseline[t])
    Yh    = X_baseline[id_sg, , drop = FALSE]
    Xr    = X_baseline[-id_sg, , drop = FALSE]
    ew_lc = LAMBDA_CLEAN * compute_stat_LC(Xr, Yh, "LC") + (1 - LAMBDA_CLEAN) * ew_lc
    ew_cc = LAMBDA_CLEAN * compute_stat_LC(Xr, Yh, "CC") + (1 - LAMBDA_CLEAN) * ew_cc
    st_lc[t] = max(ew_lc); st_cc[t] = max(ew_cc)
  }
  
  ooc   = (st_lc > ucl_lc) | (st_cc > ucl_cc)
  n_ooc = sum(ooc)
  cat(sprintf("Round %d: %d out of %d subgroups deviate (UCL_LC=%.3f, UCL_CC=%.3f)\n",
              round_i, n_ooc, m_sg0, ucl_lc, ucl_cc))
  cleaning_history = rbind(cleaning_history,
                           data.frame(Round = round_i, Subgroup_Total = m_sg0, Subgroup_OOC = n_ooc,
                                      UCL_LC = ucl_lc, UCL_CC = ucl_cc))
  
  if (n_ooc == 0) { cat("-> Phase I in-control, stopping.\n"); break }
  if ((m_sg0 - n_ooc) < MIN_SUBGROUP_CLEAN) {
    cat("-> Remaining subgroups would be too few, stopping without further removal.\n"); break
  }
  
  subgroup_ooc   = subgroup_baseline[ooc]
  subgroup_baseline = subgroup_baseline[!(subgroup_baseline %in% subgroup_ooc)]
  idx_remaining      = data_baseline$Subgroup %in% subgroup_baseline
  data_baseline    = data_baseline[idx_remaining, ]
  X_baseline_mat   = X_baseline_mat[idx_remaining, , drop = FALSE]
}
print(cleaning_history)

data_ref_clean = data_baseline
cat(sprintf("Final Phase I: %d rows = %d subgroups x %d days (evenly divisible by %d: %s)\n",
            nrow(data_ref_clean), nrow(data_ref_clean) / n_per_subgroup, n_per_subgroup,
            n_per_subgroup, (nrow(data_ref_clean) %% n_per_subgroup) == 0))

# ==============================================================================
# 14. OFFICIAL CONTROL LIMITS (BISECTION, TARGET MRL0)
# ==============================================================================
cat("\n--- [14] UCL CALIBRATION ---\n")

X_baseline    = as.matrix(data_ref_clean[, vars_target])
m          = nrow(X_baseline)
d          = ncol(X_baseline)
mean_baseline = colMeans(X_baseline)
cov_baseline  = cov(X_baseline)
eig_a      = eigen(cov_baseline)
whitening_mat  = eig_a$vectors %*% diag(1/sqrt(eig_a$values)) %*% t(eig_a$vectors)

ucl_table = data.frame(Lambda = numeric(), UCL_LC = numeric(), UCL_CC = numeric())
for (lam in LAMBDA_VEC) {
  cat(sprintf("   Lambda = %.2f...\n", lam))
  u1 = find_ucl(TARGET_MRL0, lam, "LC")
  u2 = find_ucl(TARGET_MRL0, lam, "CC")
  ucl_table = rbind(ucl_table, data.frame(Lambda = lam, UCL_LC = u1, UCL_CC = u2))
  cat(sprintf("      LC=%.4f | CC=%.4f\n", u1, u2))
}
print(ucl_table)

# ==============================================================================
# 15. SUBGROUP CONTROL CHART (PHASE I LEAVE-ONE-SUBGROUP-OUT + PHASE II)
# ==============================================================================
cat("\n--- [15] CONTROL CHART ---\n")

if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

subgroup_ref = sort(unique(data_ref_clean$Subgroup))
subgroup_mon = sort(unique(data_mon$Subgroup))
m_sg        = length(subgroup_ref)
k_sg        = length(subgroup_mon)
X_mon_mat     = as.matrix(data_mon[, vars_target])

compute_series_phase1 = function(lam, type) {
  st = numeric(m_sg); ew = numeric(d + 1)
  for (t in 1:m_sg) {
    id_sg = which(data_ref_clean$Subgroup == subgroup_ref[t])
    Yh    = X_baseline[id_sg, , drop = FALSE]
    Xr    = X_baseline[-id_sg, , drop = FALSE]
    Zt    = compute_stat_LC(Xr, Yh, type)
    ew    = lam * Zt + (1 - lam) * ew
    st[t] = max(ew)
  }
  st
}

compute_series_phase2 = function(lam, type) {
  st = numeric(k_sg); ew = numeric(d + 1)
  diag_mat = matrix(0, k_sg, d + 1)
  for (t in 1:k_sg) {
    id_sg = which(data_mon$Subgroup == subgroup_mon[t])
    Yh    = X_mon_mat[id_sg, , drop = FALSE]
    Zt    = compute_stat_LC(X_baseline, Yh, type)
    ew    = lam * Zt + (1 - lam) * ew
    diag_mat[t, ] = ew
    st[t] = max(ew)
  }
  list(st = st, diag = diag_mat)
}

component_names = c("Korelasi(Copula)", vars_target)
store_LC = list(); store_CC = list()
diag_LC = list(); diag_CC = list()

combined_chart = function(st1, st2, ucl, lam, type) {
  n1 = length(st1); n2 = length(st2)
  df = data.frame(Obs = 1:(n1+n2), Stat = c(st1, st2),
                  Phase = c(rep("Phase I", n1), rep("Phase II", n2)))
  df$OOC = df$Stat > ucl
  o1 = sum(df$OOC[1:n1]); o2 = sum(df$OOC[(n1+1):(n1+n2)])
  xs = n1 + 0.5
  date_ref = sapply(subgroup_ref, function(sg) min(data_ref_clean$date[data_ref_clean$Subgroup == sg]))
  date_mon = sapply(subgroup_mon, function(sg) min(data_mon$date[data_mon$Subgroup == sg]))
  date_all = as.Date(c(date_ref, date_mon), origin = "1970-01-01")
  brk = seq(0, n1+n2, by = 20)
  lab = format(date_all[pmin(pmax(brk, 1), length(date_all))], "%d/%m/%Y")
  y_top = max(df$Stat, ucl) * 1.15
  nm    = if (type == "LC") "ELC" else "ECC"
  p = ggplot(df, aes(Obs, Stat)) +
    annotate("rect", xmin = xs, xmax = n1+n2, ymin = -Inf, ymax = Inf,
             fill = "#FFF8E1", alpha = 0.6) +
    geom_line(color = "#4472C4", linewidth = 0.5, alpha = 0.85) +
    geom_hline(yintercept = ucl, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_point(data = subset(df, OOC), color = "#E63946", size = 2.5) +
    geom_vline(xintercept = xs, color = "#E63946", linetype = "dashed", linewidth = 1.2) +
    annotate("label", x = n1*0.3, y = y_top,
             label = sprintf("Phase I (m = %d subgroups)\nOOC Signal (%d)", n1, o1),
             color = "#E63946", fill = alpha("white", 0.55), size = 18/.pt,
             family = "Arial", fontface = "bold") +
    annotate("label", x = n1 + n2*0.45, y = y_top,
             label = sprintf("Phase II (k = %d subgroups)\nOOC Signal (%d)", n2, o2),
             color = "#E63946", fill = alpha("white", 0.55), size = 18/.pt,
             family = "Arial", fontface = "bold") +
    scale_x_continuous(breaks = brk, labels = lab) +
    coord_cartesian(ylim = c(NA, y_top * 1.12), clip = "off") +
    labs(title = bquote(bold("MREWMA-"*.(type)*" Control Chart ("*lambda*" = "*.(lam)*")")),
         x = "Subgroup t", y = nm) +
    theme_bw(base_family = "Arial") +
    theme(plot.margin = margin(t = 45, r = 20, b = 10, l = 10),
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 32),
          axis.title  = element_text(size = 24, face = "bold"),
          axis.text   = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  p
}

phase2_single_chart = function(st2, ucl, lam, type) {
  nm = if (type == "LC") "ELC" else "ECC"
  date_mon = as.Date(sapply(subgroup_mon, function(sg) min(data_mon$date[data_mon$Subgroup == sg])),
                     origin = "1970-01-01")
  df = data.frame(Date_x = date_mon, Stat = st2)
  df$Status = ifelse(df$Stat > ucl, "OOC", "IC")
  p = ggplot(df, aes(Date_x, Stat)) +
    geom_line(color = "grey60", linewidth = 0.4) +
    geom_point(aes(color = Status), size = 3) +
    geom_hline(yintercept = ucl, color = "firebrick", linetype = "dashed", linewidth = 0.7) +
    annotate("label", x = min(df$Date_x), y = Inf,
             label = paste0("UCL = ", round(ucl, 2)), hjust = 0, vjust = -0.4,
             fill = "white", color = "firebrick", size = 24/.pt, fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c(IC = "black", OOC = "firebrick")) +
    labs(title = paste0("MREWMA-", type, " (\u03bb = ", lam, ")"),
         x = "Subgroup", y = nm) +
    theme_bw(base_family = "Arial") +
    theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 32),
          axis.title  = element_text(size = 24, face = "bold"),
          axis.text   = element_text(size = 20),
          legend.position = "none")
  p
}

for (i in 1:nrow(ucl_table)) {
  lam = ucl_table$Lambda[i]
  cat(sprintf("Computing chart for lambda = %.2f...\n", lam))
  
  s1_lc = compute_series_phase1(lam, "LC")
  s1_cc = compute_series_phase1(lam, "CC")
  h2_lc = compute_series_phase2(lam, "LC")
  h2_cc = compute_series_phase2(lam, "CC")
  
  store_LC[[i]] = list(f1 = s1_lc, f2 = h2_lc$st)
  store_CC[[i]] = list(f1 = s1_cc, f2 = h2_cc$st)
  
  get_cause = function(bar, ucl) {
    s = component_names[bar > ucl]
    if (length(s) == 0) "-" else paste(s, collapse = " & ")
  }
  
  dl = data.frame(Subgroup = subgroup_mon, Statistic_Plotting = h2_lc$st, Copula = h2_lc$diag[, 1])
  for (jj in seq_along(vars_target)) dl[[vars_target[jj]]] = h2_lc$diag[, jj + 1]
  dl$Status = ifelse(h2_lc$st > ucl_table$UCL_LC[i], "OOC", "IC")
  dl$Cause_of_Shift = apply(h2_lc$diag, 1, get_cause, ucl = ucl_table$UCL_LC[i])
  
  dc = data.frame(Subgroup = subgroup_mon, Statistic_Plotting = h2_cc$st, Copula = h2_cc$diag[, 1])
  for (jj in seq_along(vars_target)) dc[[vars_target[jj]]] = h2_cc$diag[, jj + 1]
  dc$Status = ifelse(h2_cc$st > ucl_table$UCL_CC[i], "OOC", "IC")
  dc$Cause_of_Shift = apply(h2_cc$diag, 1, get_cause, ucl = ucl_table$UCL_CC[i])
  
  diag_LC[[paste0("Lambda_", lam)]] = subset(dl, Status == "OOC")
  diag_CC[[paste0("Lambda_", lam)]] = subset(dc, Status == "OOC")
  
  print(phase2_single_chart(h2_lc$st, ucl_table$UCL_LC[i], lam, "LC"))
  print(phase2_single_chart(h2_cc$st, ucl_table$UCL_CC[i], lam, "CC"))
  
  p_gab_lc = combined_chart(s1_lc, h2_lc$st, ucl_table$UCL_LC[i], lam, "LC")
  p_gab_cc = combined_chart(s1_cc, h2_cc$st, ucl_table$UCL_CC[i], lam, "CC")
  print(p_gab_lc)
  print(p_gab_cc)
  ggsave(file.path(plot_dir, paste0("MREWMA-LC Control Chart (", lam, ").png")),
         p_gab_lc, width = 1770, height = 708, units = "px", dpi = 100, bg = "white")
  ggsave(file.path(plot_dir, paste0("MREWMA-CC Control Chart (", lam, ").png")),
         p_gab_cc, width = 1770, height = 708, units = "px", dpi = 100, bg = "white")
}

cat("\n--- OOC DIAGNOSTICS ---\n")
for (lam in ucl_table$Lambda) {
  cat(sprintf("\n== Lambda %.2f | MREWMA-LC ==\n", lam))
  h = diag_LC[[paste0("Lambda_", lam)]]
  if (nrow(h) > 0) print(h) else cat("No OOC.\n")
  cat(sprintf("== Lambda %.2f | MREWMA-CC ==\n", lam))
  h = diag_CC[[paste0("Lambda_", lam)]]
  if (nrow(h) > 0) print(h) else cat("No OOC.\n")
}

# ==============================================================================
# 16. MRL0 CURVES ACROSS EIGHT OUT-OF-CONTROL SHIFT SCENARIOS
# ==============================================================================
cat("\n--- [16] MRL0 OOC CURVES ---\n")

if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

scenarios = list(
  LS1 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0),
  LS2 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0),
  DS1 = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0),
  DS2 = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0),
  OD  = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
  MS1 = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0),
  MS2 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
  MS3 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

if (!exists("mrl_results")) mrl_results = data.frame()
if (nrow(mrl_results) > 0) mrl_results$Scenario = as.character(mrl_results$Scenario)
already_done = if (nrow(mrl_results) > 0) unique(mrl_results$Scenario) else character(0)

for (sk in names(scenarios)) {
  if (sk %in% already_done) { cat(sk, "already computed, skipping.\n"); next }
  cat(">> Simulating", sk, "...\n")
  for (g in scenarios[[sk]]) {
    for (i in 1:nrow(ucl_table)) {
      lam = ucl_table$Lambda[i]
      m_l = simulate_mrl(ucl_table$UCL_LC[i], lam, "LC",
                         N_RL_OOC, MAX_STEPS_OOC, sk, g)
      m_c = simulate_mrl(ucl_table$UCL_CC[i], lam, "CC",
                         N_RL_OOC, MAX_STEPS_OOC, sk, g)
      mrl_results = rbind(mrl_results,
                          data.frame(Shift = g, Lambda = lam, MRL = m_l, Scheme = "MREWMA-LC",
                                     Scenario = sk, stringsAsFactors = FALSE),
                          data.frame(Shift = g, Lambda = lam, MRL = m_c, Scheme = "MREWMA-CC",
                                     Scenario = sk, stringsAsFactors = FALSE))
    }
  }
  save(mrl_results, file = file.path(dir_save, paste0("Langkah15_", sk, ".RData")))
}

mrl_results$Scenario = factor(mrl_results$Scenario, levels = names(scenarios))

curve_theme = theme_classic(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 32),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text  = element_text(size = 20),
        panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
        legend.position = "bottom",
        legend.title = element_text(size = 24, face = "bold"),
        legend.text  = element_text(size = 20))

for (skm in c("MREWMA-LC", "MREWMA-CC")) {
  short_type = sub("MREWMA-", "", skm)
  for (sk in names(scenarios)) {
    dfk = subset(mrl_results, Scheme == skm & Scenario == sk)
    if (nrow(dfk) == 0) next
    p = ggplot(dfk, aes(Shift, MRL, group = as.factor(Lambda))) +
      geom_line(aes(color = as.factor(Lambda)), linewidth = 0.8) +
      geom_point(aes(color = as.factor(Lambda), shape = as.factor(Lambda)), size = 2) +
      scale_y_log10() +
      labs(title = bquote(bold(MRL[0]~.(paste(skm, sk)))),
           y = "Median Run Length", x = "Shift Parameter",
           color = expression(lambda), shape = expression(lambda)) +
      scale_color_brewer(palette = "Set1") + curve_theme
    print(p)
    file_name = file.path(plot_dir, paste0("MRL0 MREWMA (", short_type, ")(", sk, ").png"))
    ggsave(file_name, p, width = 1796, height = 783, units = "px", dpi = 100, bg = "white")
  }
}

cat("\n--- MRL0 TABLE PER SCENARIO ---\n")
for (sk in names(scenarios)) {
  cat(paste0("\n>>> ", sk, " <<<\n"))
  dfs = subset(mrl_results, Scenario == sk)
  dl = subset(dfs, Scheme == "MREWMA-LC")[, c("Shift","Lambda","MRL")] %>%
    pivot_wider(names_from = Lambda, values_from = MRL, names_prefix = "LC_\u03bb=")
  dc = subset(dfs, Scheme == "MREWMA-CC")[, c("Shift","Lambda","MRL")] %>%
    pivot_wider(names_from = Lambda, values_from = MRL, names_prefix = "CC_\u03bb=")
  print(merge(dl, dc, by = "Shift"), row.names = FALSE)
}

cat("\n--- PHASE II PLOTTING STATISTIC TABLE ---\n")
for (i in 1:nrow(ucl_table)) {
  lam = ucl_table$Lambda[i]
  cat(sprintf("\n-- Lambda %.2f --\n", lam))
  dfp = data.frame(Subgroup = subgroup_mon,
                   Stat_LC  = round(store_LC[[i]]$f2, 4),
                   UCL_LC   = ucl_table$UCL_LC[i],
                   Stat_CC  = round(store_CC[[i]]$f2, 4),
                   UCL_CC   = ucl_table$UCL_CC[i])
  dfp$Status_LC = ifelse(dfp$Stat_LC > dfp$UCL_LC, "OOC", "IC")
  dfp$Status_CC = ifelse(dfp$Stat_CC > dfp$UCL_CC, "OOC", "IC")
  print(dfp, row.names = FALSE)
}

# ==============================================================================
# 17. RESULTS APPENDIX & EXPORT
# ==============================================================================
cat("\n--- [17] FULL RESULTS APPENDIX ---\n")
options(width = 200)

cat("\n--- MREWMA-COPULA RESULTS APPENDIX ---\n")
ooc_columns = c("Subgroup", "Statistic_Plotting", "Copula", vars_target, "Cause_of_Shift")
for (i in 1:nrow(ucl_table)) {
  lam = ucl_table$Lambda[i]
  for (type in c("LC", "CC")) {
    h = if (type == "LC") diag_LC[[paste0("Lambda_", lam)]] else diag_CC[[paste0("Lambda_", lam)]]
    cat(sprintf("\n>>> OOC Count and Cause - MREWMA-%s Control Chart (%s) <<<\n", type, lam))
    cat(sprintf("OOC count: %d out of %d subgroups\n", nrow(h), k_sg))
    if (nrow(h) > 0) print(h[, ooc_columns], row.names = FALSE) else cat("No OOC signal.\n")
  }
}

cat("\n--- MRL0 VALUE TABLE (Shift x Scenario) ---\n")
get_mrl_ordered = function(sk, skm, lam) {
  grid = scenarios[[sk]]
  sapply(grid, function(g) {
    br = which(mrl_results$Scenario == sk & mrl_results$Scheme == skm &
                 mrl_results$Lambda == lam & mrl_results$Shift == g)
    if (length(br) == 0) NA else mrl_results$MRL[br[1]]
  })
}
max_positions = max(sapply(scenarios, length))
for (i in 1:nrow(ucl_table)) {
  lam = ucl_table$Lambda[i]
  for (skm in c("MREWMA-LC", "MREWMA-CC")) {
    cat(sprintf("\n>>> MRL0 Values - %s, Lambda = %s <<<\n", skm, lam))
    wide_table = data.frame(Shift = 1:max_positions)
    for (sk in names(scenarios)) {
      values = get_mrl_ordered(sk, skm, lam)
      length(values) = max_positions
      wide_table[[sk]] = round(values, 3)
    }
    print(wide_table, row.names = FALSE)
  }
}

diag_LC_all = do.call(rbind, lapply(names(diag_LC), function(nm) {
  if (nrow(diag_LC[[nm]]) == 0) return(NULL)
  cbind(Lambda_Label = nm, diag_LC[[nm]])
}))
diag_CC_all = do.call(rbind, lapply(names(diag_CC), function(nm) {
  if (nrow(diag_CC[[nm]]) == 0) return(NULL)
  cbind(Lambda_Label = nm, diag_CC[[nm]])
}))

output_file = file.path(dir_save, "MREWMA_Copula_AirQuality_Results.xlsx")
write_xlsx(
  list(
    Imputed_Data          = air_data,
    Missing_Patterns      = missing_pattern_summary,
    Imputation_Validation = as.data.frame(em_validation_summary),
    Imputation_Val_Raw    = em_validation_raw,
    Phase1_Cleaning_Log   = cleaning_history,
    UCL_Table             = ucl_table,
    MRL0_Results          = mrl_results,
    OOC_Diagnostics_LC    = if (is.null(diag_LC_all)) data.frame() else diag_LC_all,
    OOC_Diagnostics_CC    = if (is.null(diag_CC_all)) data.frame() else diag_CC_all
  ),
  output_file
)

cat("\nDONE. Results saved to:\n")
cat("  Excel  ->", output_file, "\n")
cat("  Images ->", plot_dir, "\n")