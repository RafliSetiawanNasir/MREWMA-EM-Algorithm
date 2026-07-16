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
#       Langkah 11-16), adapted from d=4 (NI/FE/SI/MG) to the d=3 air-quality
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
#   (e) Section 10 is new: descriptive statistics (Section 9) now prints its
#       summary table, and multivariate normality is formally tested via the
#       Shapiro-Wilk test (mvnormtest::mshapiro.test) on the three imputed
#       pollutant variables jointly, with explicit hypotheses, testing
#       criterion, and decision reported. Sections 11-17 are renumbered
#       accordingly (previously 10-16).
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
dir_save    = folder_path
plot_dir    = file.path(dir_save, "plots")

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
# Upper starting bound for the UCL bisection search (Section 14, cari_ucl()).
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

if (!dir.exists(dir_save)) dir.create(dir_save, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

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

cat("\n--- [9] STATISTIK DESKRIPTIF ---\n")
print(as.data.frame(stats_all), row.names = FALSE)

# ==============================================================================
# 10. MULTIVARIATE NORMALITY ASSUMPTION TESTING (SHAPIRO-WILK)
# ==============================================================================
cat("\n--- [10] PENGUJIAN ASUMSI NORMALITAS MULTIVARIAT ---\n")

cat("\nHipotesis:\n")
cat("  H0 : Data berdistribusi normal multivariat\n")
cat("  H1 : Data tidak berdistribusi normal multivariat\n")

cat("\nKriteria Pengujian:\n")
cat("  Jika nilai p-value > alpha, dengan alpha = 0.05, maka H0 diterima\n\n")

X_normality = t(as.matrix(air_data[, vars_target]))  # mshapiro.test expects p x n (variables in rows)
mvn_result  = mshapiro.test(X_normality)

normality_table = data.frame(
  `Shapiro-Wilk` = round(as.numeric(mvn_result$statistic), 3),
  `p-value`      = formatC(mvn_result$p.value, format = "f", digits = 3),
  Keputusan      = ifelse(mvn_result$p.value > 0.05, "H0 diterima", "H0 ditolak"),
  check.names    = FALSE
)
print(normality_table, row.names = FALSE)

if (mvn_result$p.value <= 0.05) {
  cat(sprintf("\nKarena nilai p-value < alpha (%.3f < 0.05), maka H0 ditolak: data indikator\n",
              mvn_result$p.value))
  cat("pencemaran udara tidak berdistribusi normal multivariat. Selanjutnya, data dapat\n")
  cat("dianalisis menggunakan metode MREWMA-EM CC dan MREWMA-EM LC yang tahan terhadap\n")
  cat("pelanggaran asumsi distribusi normal pada data.\n")
} else {
  cat(sprintf("\nKarena nilai p-value > alpha (%.3f > 0.05), maka H0 diterima: data indikator\n",
              mvn_result$p.value))
  cat("pencemaran udara berdistribusi normal multivariat.\n")
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
cat("\n--- [11] MESIN MREWMA-COPULA ---\n")

ambil_pcs = function(X_ref, Y_uji) {
  N  = nrow(X_ref) + nrow(Y_uji)
  gb = rbind(X_ref, Y_uji)
  G  = apply(gb, 2, function(k) rank(k) / (N + 1))
  eg = eigen(cov(G))
  Th = G %*% eg$vectors
  list(ref = Th[1:nrow(X_ref), , drop = FALSE],
       uji = Th[(nrow(X_ref)+1):N, , drop = FALSE])
}

.simpanan_nol = new.env(parent = emptyenv())

bangun_nol_perm = function(m, n, B = B_PERM) {
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

ambil_nol_perm = function(m, n, B = B_PERM) {
  kunci = paste(m, n, B, sep = "_")
  ob = .simpanan_nol[[kunci]]
  if (is.null(ob)) { ob = bangun_nol_perm(m, n, B); assign(kunci, ob, envir = .simpanan_nol) }
  ob
}

pval_LC = function(x, y, B = B_PERM) {
  m = length(x); n = length(y); N = m + n
  nol = ambil_nol_perm(m, n, B); mo = nol$mo
  r   = rank(c(x, y)); iy = (m+1):N
  TW  = sum(r[iy])
  TAB = sum(pmin(r, N+1-r)[iy])
  TM  = sum(((r-(N+1)/2)^2)[iy])
  Lo  = ((TW-mo$muW)/mo$sdW)^2 + ((TAB-mo$muAB)/mo$sdAB)^2
  Co  = ((TW-mo$muW)/mo$sdW)^2 + ((TM -mo$muM )/mo$sdM )^2
  amb = function(obs, ns, B) (1 + (B - findInterval(obs - 1e-12, ns))) / (B + 1)
  c(max(min(amb(Lo, nol$L, B), 0.99999), 1e-10),
    max(min(amb(Co, nol$C, B), 0.99999), 1e-10))
}

statistik_LC = function(X_ref, Y_uji, tipe = "LC") {
  d = ncol(X_ref)
  Z = numeric(d + 1)
  pc = ambil_pcs(X_ref, Y_uji)
  pv = numeric(d)
  for (k in 1:d) {
    h = pval_LC(pc$ref[, k], pc$uji[, k])
    pv[k] = ifelse(tipe == "LC", h[1], h[2])
  }
  Z[1] = -log(min(pv)) - 1
  for (k in 1:d) {
    h = pval_LC(X_ref[, k], Y_uji[, k])
    Z[k + 1] = -log(ifelse(tipe == "LC", h[1], h[2])) - 1
  }
  Z
}

geser_data = function(sampel, tipe, prm, rata, kov, putih) {
  Z = sweep(sampel, 2, rata, "-") %*% putih
  mu_b = numeric(d); kv_b = diag(d)
  if (tipe == "LS1") mu_b[1] = prm
  else if (tipe == "LS2") { mu_b[1] = prm; mu_b[2] = 0.5 }
  else if (tipe == "DS1") kv_b[1,1] = prm
  else if (tipe == "DS2") { kv_b[1,1] = prm; kv_b[2,2] = 1.5 }
  else if (tipe == "OD")  { kv_b[1,2] = prm; kv_b[2,1] = prm }
  else if (tipe == "MS1") { mu_b[1] = 0.5; kv_b[1,1] = prm }
  else if (tipe == "MS2") { mu_b[1] = prm; kv_b[1,1] = 1.5 }
  else if (tipe == "MS3") { mu_b[1] = 0.5; kv_b[1,1] = 1.5; kv_b[1,2] = prm; kv_b[2,1] = prm }
  ch = tryCatch(chol(kv_b), error = function(e) diag(sqrt(pmax(diag(kv_b), 1e-8))))
  Zg = sweep(Z %*% ch, 2, mu_b, "+")
  ei = eigen(kov)
  wa = ei$vectors %*% diag(sqrt(pmax(ei$values, 0))) %*% t(ei$vectors)
  sweep(Zg %*% wa, 2, rata, "+")
}

mrl_simulasi = function(ucl, lam, tipe = "LC", n_sim = N_RL_OOC,
                        maks = MAX_STEPS_OOC, jenis_geser = "None",
                        prm_geser = 0, n_uji = N_TEST) {
  rl = numeric(n_sim)
  for (s in 1:n_sim) {
    ew = numeric(d + 1)
    for (t in 1:maks) {
      id = sample(1:m, n_uji, replace = TRUE)
      Yb = X_acuan[id, , drop = FALSE]
      if (jenis_geser != "None")
        Yb = geser_data(Yb, jenis_geser, prm_geser, rata_acuan, kov_acuan, mat_putih)
      Zt = statistik_LC(X_acuan, Yb, tipe)
      ew = lam * Zt + (1 - lam) * ew
      if (max(ew) > ucl) { rl[s] = t; break }
      if (t == maks) rl[s] = maks
    }
  }
  median(rl)
}

cari_ucl = function(target, lam, tipe, A1 = 0.1, A2 = UCL_SEARCH_UPPER,
                    n_sim = N_RL_UCL, maks = MAX_STEPS_UCL,
                    tol = TOL_BISECT, maks_iter = 30) {
  f  = function(u) mrl_simulasi(u, lam, tipe, n_sim, maks, "None") - target
  fa = f(A1); fb = f(A2); ek = 0
  while (fb < 0 && ek < 6)  { A2 = A2 * 1.5; fb = f(A2); ek = ek + 1 }
  while (fa > 0 && ek < 12) { A1 = A1 / 1.5; fa = f(A1); ek = ek + 1 }
  for (it in 1:maks_iter) {
    tg = (A1 + A2) / 2
    ft = f(tg)
    if (abs(ft) <= tol) return(tg)
    if (ft > 0) A2 = tg else A1 = tg
  }
  (A1 + A2) / 2
}

# ==============================================================================
# 13. PHASE I ITERATIVE CLEANING (LEPAGE/CUCCONI, LEAVE-ONE-SUBGROUP-OUT)
# ==============================================================================
cat("\n--- [12] PEMBERSIHAN FASE I ITERATIF ---\n")

LAM_BERSIH          = 0.1
N_RL_BERSIH         = 300
MAX_STEPS_BERSIH    = 400
MAKS_PUTARAN_BERSIH = 10
MIN_SUBGRUP_BERSIH  = 30

d = length(vars_target)

subgrup_acuan = sort(unique(data_ref$Subgroup))
data_acuan    = data_ref
X_acuan_mat   = as.matrix(data_acuan[, vars_target])
riwayat_bersih = data.frame()

for (put in 1:MAKS_PUTARAN_BERSIH) {
  m_sg0      = length(subgrup_acuan)
  X_acuan    = X_acuan_mat
  m          = nrow(X_acuan)
  rata_acuan = colMeans(X_acuan)
  kov_acuan  = cov(X_acuan)
  eig_a      = eigen(kov_acuan)
  mat_putih  = eig_a$vectors %*% diag(1/sqrt(eig_a$values)) %*% t(eig_a$vectors)

  ucl_lc = cari_ucl(TARGET_MRL0, LAM_BERSIH, "LC", n_sim = N_RL_BERSIH, maks = MAX_STEPS_BERSIH)
  ucl_cc = cari_ucl(TARGET_MRL0, LAM_BERSIH, "CC", n_sim = N_RL_BERSIH, maks = MAX_STEPS_BERSIH)

  st_lc = numeric(m_sg0); st_cc = numeric(m_sg0)
  ew_lc = numeric(d + 1); ew_cc = numeric(d + 1)
  for (t in 1:m_sg0) {
    id_sg = which(data_acuan$Subgroup == subgrup_acuan[t])
    Yh    = X_acuan[id_sg, , drop = FALSE]
    Xr    = X_acuan[-id_sg, , drop = FALSE]
    ew_lc = LAM_BERSIH * statistik_LC(Xr, Yh, "LC") + (1 - LAM_BERSIH) * ew_lc
    ew_cc = LAM_BERSIH * statistik_LC(Xr, Yh, "CC") + (1 - LAM_BERSIH) * ew_cc
    st_lc[t] = max(ew_lc); st_cc[t] = max(ew_cc)
  }

  ooc   = (st_lc > ucl_lc) | (st_cc > ucl_cc)
  n_ooc = sum(ooc)
  cat(sprintf("Putaran %d: %d dari %d subgrup menyimpang (UCL_LC=%.3f, UCL_CC=%.3f)\n",
              put, n_ooc, m_sg0, ucl_lc, ucl_cc))
  riwayat_bersih = rbind(riwayat_bersih,
                         data.frame(Putaran = put, Subgrup_Total = m_sg0, Subgrup_OOC = n_ooc,
                                    UCL_LC = ucl_lc, UCL_CC = ucl_cc))

  if (n_ooc == 0) { cat("-> Fase I in-control, berhenti.\n"); break }
  if ((m_sg0 - n_ooc) < MIN_SUBGRUP_BERSIH) {
    cat("-> Sisa subgrup akan terlalu sedikit, berhenti tanpa membuang.\n"); break
  }

  subgrup_ooc   = subgrup_acuan[ooc]
  subgrup_acuan = subgrup_acuan[!(subgrup_acuan %in% subgrup_ooc)]
  idx_sisa      = data_acuan$Subgroup %in% subgrup_acuan
  data_acuan    = data_acuan[idx_sisa, ]
  X_acuan_mat   = X_acuan_mat[idx_sisa, , drop = FALSE]
}
print(riwayat_bersih)

data_ref_clean = data_acuan
cat(sprintf("Fase I final: %d baris = %d subgrup x %d hari (habis dibagi %d: %s)\n",
            nrow(data_ref_clean), nrow(data_ref_clean) / n_per_subgroup, n_per_subgroup,
            n_per_subgroup, (nrow(data_ref_clean) %% n_per_subgroup) == 0))

# ==============================================================================
# 14. OFFICIAL CONTROL LIMITS (BISECTION, TARGET MRL0)
# ==============================================================================
cat("\n--- [13] KALIBRASI UCL ---\n")

X_acuan    = as.matrix(data_ref_clean[, vars_target])
m          = nrow(X_acuan)
d          = ncol(X_acuan)
rata_acuan = colMeans(X_acuan)
kov_acuan  = cov(X_acuan)
eig_a      = eigen(kov_acuan)
mat_putih  = eig_a$vectors %*% diag(1/sqrt(eig_a$values)) %*% t(eig_a$vectors)

tabel_ucl = data.frame(Lambda = numeric(), UCL_LC = numeric(), UCL_CC = numeric())
for (lam in LAMBDA_VEC) {
  cat(sprintf("   Lambda = %.2f...\n", lam))
  u1 = cari_ucl(TARGET_MRL0, lam, "LC")
  u2 = cari_ucl(TARGET_MRL0, lam, "CC")
  tabel_ucl = rbind(tabel_ucl, data.frame(Lambda = lam, UCL_LC = u1, UCL_CC = u2))
  cat(sprintf("      LC=%.4f | CC=%.4f\n", u1, u2))
}
print(tabel_ucl)

# ==============================================================================
# 15. SUBGROUP CONTROL CHART (PHASE I LEAVE-ONE-SUBGROUP-OUT + PHASE II)
# ==============================================================================
cat("\n--- [14] BAGAN KENDALI ---\n")

subgrup_ref = sort(unique(data_ref_clean$Subgroup))
subgrup_mon = sort(unique(data_mon$Subgroup))
m_sg        = length(subgrup_ref)
k_sg        = length(subgrup_mon)
X_mon_m     = as.matrix(data_mon[, vars_target])

hitung_deret_fase1 = function(lam, tipe) {
  st = numeric(m_sg); ew = numeric(d + 1)
  for (t in 1:m_sg) {
    id_sg = which(data_ref_clean$Subgroup == subgrup_ref[t])
    Yh    = X_acuan[id_sg, , drop = FALSE]
    Xr    = X_acuan[-id_sg, , drop = FALSE]
    Zt    = statistik_LC(Xr, Yh, tipe)
    ew    = lam * Zt + (1 - lam) * ew
    st[t] = max(ew)
  }
  st
}

hitung_deret_fase2 = function(lam, tipe) {
  st = numeric(k_sg); ew = numeric(d + 1)
  diag_mat = matrix(0, k_sg, d + 1)
  for (t in 1:k_sg) {
    id_sg = which(data_mon$Subgroup == subgrup_mon[t])
    Yh    = X_mon_m[id_sg, , drop = FALSE]
    Zt    = statistik_LC(X_acuan, Yh, tipe)
    ew    = lam * Zt + (1 - lam) * ew
    diag_mat[t, ] = ew
    st[t] = max(ew)
  }
  list(st = st, diag = diag_mat)
}

nama_komp = c("Korelasi(Copula)", vars_target)
simpan_LC = list(); simpan_CC = list()
diag_LC = list(); diag_CC = list()

bagan_gabungan = function(st1, st2, ucl, lam, tipe) {
  n1 = length(st1); n2 = length(st2)
  df = data.frame(Obs = 1:(n1+n2), Stat = c(st1, st2),
                  Fase = c(rep("Phase I", n1), rep("Phase II", n2)))
  df$OOC = df$Stat > ucl
  o1 = sum(df$OOC[1:n1]); o2 = sum(df$OOC[(n1+1):(n1+n2)])
  xs = n1 + 0.5
  tgl_ref = sapply(subgrup_ref, function(sg) min(data_ref_clean$date[data_ref_clean$Subgroup == sg]))
  tgl_mon = sapply(subgrup_mon, function(sg) min(data_mon$date[data_mon$Subgroup == sg]))
  tgl_semua = as.Date(c(tgl_ref, tgl_mon), origin = "1970-01-01")
  brk = seq(0, n1+n2, by = 20)
  lab = format(tgl_semua[pmin(pmax(brk, 1), length(tgl_semua))], "%d/%m/%Y")
  y_top = max(df$Stat, ucl) * 1.15
  nm    = if (tipe == "LC") "ELC" else "ECC"
  p = ggplot(df, aes(Obs, Stat)) +
    annotate("rect", xmin = xs, xmax = n1+n2, ymin = -Inf, ymax = Inf,
             fill = "#FFF8E1", alpha = 0.6) +
    geom_line(color = "#4472C4", linewidth = 0.5, alpha = 0.85) +
    geom_hline(yintercept = ucl, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_point(data = subset(df, OOC), color = "#E63946", size = 2.5) +
    geom_vline(xintercept = xs, color = "#E63946", linetype = "dashed", linewidth = 1.2) +
    annotate("label", x = n1*0.3, y = y_top,
             label = sprintf("Phase I (m = %d subgrup)\nOOC Signal (%d)", n1, o1),
             color = "#E63946", fill = alpha("white", 0.55), size = 18/.pt,
             family = "Arial", fontface = "bold") +
    annotate("label", x = n1 + n2*0.45, y = y_top,
             label = sprintf("Phase II (k = %d subgrup)\nOOC Signal (%d)", n2, o2),
             color = "#E63946", fill = alpha("white", 0.55), size = 18/.pt,
             family = "Arial", fontface = "bold") +
    scale_x_continuous(breaks = brk, labels = lab) +
    coord_cartesian(ylim = c(NA, y_top * 1.12), clip = "off") +
    labs(title = bquote(bold("Bagan Kendali MREWMA-"*.(tipe)*
                               " ("*lambda*" = "*.(lam)*")")),
         x = "Subgrup t", y = nm) +
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

bagan_fase2_tunggal = function(st2, ucl, lam, tipe) {
  nm = if (tipe == "LC") "ELC" else "ECC"
  tgl_mon = as.Date(sapply(subgrup_mon, function(sg) min(data_mon$date[data_mon$Subgroup == sg])),
                    origin = "1970-01-01")
  df = data.frame(Tanggal = tgl_mon, Stat = st2)
  df$Status = ifelse(df$Stat > ucl, "OOC", "IC")
  p = ggplot(df, aes(Tanggal, Stat)) +
    geom_line(color = "grey60", linewidth = 0.4) +
    geom_point(aes(color = Status), size = 3) +
    geom_hline(yintercept = ucl, color = "firebrick", linetype = "dashed", linewidth = 0.7) +
    annotate("label", x = min(df$Tanggal), y = Inf,
             label = paste0("UCL = ", round(ucl, 2)), hjust = 0, vjust = -0.4,
             fill = "white", color = "firebrick", size = 24/.pt, fontface = "bold") +
    coord_cartesian(clip = "off") +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    scale_color_manual(values = c(IC = "black", OOC = "firebrick")) +
    labs(title = paste0("MREWMA-", tipe, " (\u03bb = ", lam, ")"),
         x = "Subgrup", y = nm) +
    theme_bw(base_family = "Arial") +
    theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 32),
          axis.title  = element_text(size = 24, face = "bold"),
          axis.text   = element_text(size = 20),
          legend.position = "none")
  p
}

for (i in 1:nrow(tabel_ucl)) {
  lam = tabel_ucl$Lambda[i]
  cat(sprintf("Menghitung bagan lambda = %.2f...\n", lam))

  s1_lc = hitung_deret_fase1(lam, "LC")
  s1_cc = hitung_deret_fase1(lam, "CC")
  h2_lc = hitung_deret_fase2(lam, "LC")
  h2_cc = hitung_deret_fase2(lam, "CC")

  simpan_LC[[i]] = list(f1 = s1_lc, f2 = h2_lc$st)
  simpan_CC[[i]] = list(f1 = s1_cc, f2 = h2_cc$st)

  amb_sebab = function(bar, ucl) {
    s = nama_komp[bar > ucl]
    if (length(s) == 0) "-" else paste(s, collapse = " & ")
  }

  dl = data.frame(Subgroup = subgrup_mon, Statistic_Plotting = h2_lc$st, Copula = h2_lc$diag[, 1])
  for (jj in seq_along(vars_target)) dl[[vars_target[jj]]] = h2_lc$diag[, jj + 1]
  dl$Status = ifelse(h2_lc$st > tabel_ucl$UCL_LC[i], "OOC", "IC")
  dl$Cause_of_Shift = apply(h2_lc$diag, 1, amb_sebab, ucl = tabel_ucl$UCL_LC[i])

  dc = data.frame(Subgroup = subgrup_mon, Statistic_Plotting = h2_cc$st, Copula = h2_cc$diag[, 1])
  for (jj in seq_along(vars_target)) dc[[vars_target[jj]]] = h2_cc$diag[, jj + 1]
  dc$Status = ifelse(h2_cc$st > tabel_ucl$UCL_CC[i], "OOC", "IC")
  dc$Cause_of_Shift = apply(h2_cc$diag, 1, amb_sebab, ucl = tabel_ucl$UCL_CC[i])

  diag_LC[[paste0("Lambda_", lam)]] = subset(dl, Status == "OOC")
  diag_CC[[paste0("Lambda_", lam)]] = subset(dc, Status == "OOC")

  print(bagan_fase2_tunggal(h2_lc$st, tabel_ucl$UCL_LC[i], lam, "LC"))
  print(bagan_fase2_tunggal(h2_cc$st, tabel_ucl$UCL_CC[i], lam, "CC"))

  p_gab_lc = bagan_gabungan(s1_lc, h2_lc$st, tabel_ucl$UCL_LC[i], lam, "LC")
  p_gab_cc = bagan_gabungan(s1_cc, h2_cc$st, tabel_ucl$UCL_CC[i], lam, "CC")
  print(p_gab_lc)
  print(p_gab_cc)
  ggsave(file.path(plot_dir, paste0("Bagan Kendali MREWMA-LC (", lam, ").png")),
         p_gab_lc, width = 1770, height = 708, units = "px", dpi = 100, bg = "white")
  ggsave(file.path(plot_dir, paste0("Bagan Kendali MREWMA-CC (", lam, ").png")),
         p_gab_cc, width = 1770, height = 708, units = "px", dpi = 100, bg = "white")
}

cat("\n--- DIAGNOSTIK OOC ---\n")
for (lam in tabel_ucl$Lambda) {
  cat(sprintf("\n== Lambda %.2f | MREWMA-LC ==\n", lam))
  h = diag_LC[[paste0("Lambda_", lam)]]
  if (nrow(h) > 0) print(h) else cat("Tidak ada OOC.\n")
  cat(sprintf("== Lambda %.2f | MREWMA-CC ==\n", lam))
  h = diag_CC[[paste0("Lambda_", lam)]]
  if (nrow(h) > 0) print(h) else cat("Tidak ada OOC.\n")
}

# ==============================================================================
# 16. MRL0 CURVES ACROSS EIGHT OUT-OF-CONTROL SHIFT SCENARIOS
# ==============================================================================
cat("\n--- [15] KURVA MRL0 OOC ---\n")

skenario = list(
  LS1 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0),
  LS2 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0),
  DS1 = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0),
  DS2 = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0),
  OD  = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
  MS1 = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0),
  MS2 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
  MS3 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

if (!exists("hasil_mrl")) hasil_mrl = data.frame()
if (nrow(hasil_mrl) > 0) hasil_mrl$Skenario = as.character(hasil_mrl$Skenario)
sudah = if (nrow(hasil_mrl) > 0) unique(hasil_mrl$Skenario) else character(0)

for (sk in names(skenario)) {
  if (sk %in% sudah) { cat(sk, "sudah ada, lewati.\n"); next }
  cat(">> Simulasi", sk, "...\n")
  for (g in skenario[[sk]]) {
    for (i in 1:nrow(tabel_ucl)) {
      lam = tabel_ucl$Lambda[i]
      m_l = mrl_simulasi(tabel_ucl$UCL_LC[i], lam, "LC",
                         N_RL_OOC, MAX_STEPS_OOC, sk, g)
      m_c = mrl_simulasi(tabel_ucl$UCL_CC[i], lam, "CC",
                         N_RL_OOC, MAX_STEPS_OOC, sk, g)
      hasil_mrl = rbind(hasil_mrl,
                        data.frame(Shift = g, Lambda = lam, MRL = m_l, Skema = "MREWMA-LC",
                                   Skenario = sk, stringsAsFactors = FALSE),
                        data.frame(Shift = g, Lambda = lam, MRL = m_c, Skema = "MREWMA-CC",
                                   Skenario = sk, stringsAsFactors = FALSE))
    }
  }
  save(hasil_mrl, file = file.path(dir_save, paste0("Langkah15_", sk, ".RData")))
}

hasil_mrl$Skenario = factor(hasil_mrl$Skenario, levels = names(skenario))

tema_kurva = theme_classic(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 32),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text  = element_text(size = 20),
        panel.grid.major = element_line(color = "grey90", linetype = "dotted"),
        legend.position = "bottom",
        legend.title = element_text(size = 24, face = "bold"),
        legend.text  = element_text(size = 20))

for (skm in c("MREWMA-LC", "MREWMA-CC")) {
  tipe_singkat = sub("MREWMA-", "", skm)
  for (sk in names(skenario)) {
    dfk = subset(hasil_mrl, Skema == skm & Skenario == sk)
    if (nrow(dfk) == 0) next
    p = ggplot(dfk, aes(Shift, MRL, group = as.factor(Lambda))) +
      geom_line(aes(color = as.factor(Lambda)), linewidth = 0.8) +
      geom_point(aes(color = as.factor(Lambda), shape = as.factor(Lambda)), size = 2) +
      scale_y_log10() +
      labs(title = bquote(bold(MRL[0]~.(paste(skm, sk)))),
           y = "Median Run Length", x = "Shift Parameter",
           color = expression(lambda), shape = expression(lambda)) +
      scale_color_brewer(palette = "Set1") + tema_kurva
    print(p)
    nama_file = file.path(plot_dir, paste0("MRL0 MREWMA (", tipe_singkat, ")(", sk, ").png"))
    ggsave(nama_file, p, width = 1796, height = 783, units = "px", dpi = 100, bg = "white")
  }
}

cat("\n--- TABEL MRL0 PER SKENARIO ---\n")
for (sk in names(skenario)) {
  cat(paste0("\n>>> ", sk, " <<<\n"))
  dfs = subset(hasil_mrl, Skenario == sk)
  dl = subset(dfs, Skema == "MREWMA-LC")[, c("Shift","Lambda","MRL")] %>%
    pivot_wider(names_from = Lambda, values_from = MRL, names_prefix = "LC_\u03bb=")
  dc = subset(dfs, Skema == "MREWMA-CC")[, c("Shift","Lambda","MRL")] %>%
    pivot_wider(names_from = Lambda, values_from = MRL, names_prefix = "CC_\u03bb=")
  print(merge(dl, dc, by = "Shift"), row.names = FALSE)
}

cat("\n--- TABEL STATISTIK PLOTTING FASE II ---\n")
for (i in 1:nrow(tabel_ucl)) {
  lam = tabel_ucl$Lambda[i]
  cat(sprintf("\n-- Lambda %.2f --\n", lam))
  dfp = data.frame(Subgroup = subgrup_mon,
                   Stat_LC  = round(simpan_LC[[i]]$f2, 4),
                   UCL_LC   = tabel_ucl$UCL_LC[i],
                   Stat_CC  = round(simpan_CC[[i]]$f2, 4),
                   UCL_CC   = tabel_ucl$UCL_CC[i])
  dfp$Status_LC = ifelse(dfp$Stat_LC > dfp$UCL_LC, "OOC", "IC")
  dfp$Status_CC = ifelse(dfp$Stat_CC > dfp$UCL_CC, "OOC", "IC")
  print(dfp, row.names = FALSE)
}

# ==============================================================================
# 17. RESULTS APPENDIX & EXPORT
# ==============================================================================
cat("\n--- [16] LAMPIRAN HASIL LENGKAP ---\n")
options(width = 200)

cat("\n--- LAMPIRAN HASIL MREWMA-COPULA ---\n")
kolom_ooc = c("Subgroup", "Statistic_Plotting", "Copula", vars_target, "Cause_of_Shift")
for (i in 1:nrow(tabel_ucl)) {
  lam = tabel_ucl$Lambda[i]
  for (tipe in c("LC", "CC")) {
    h = if (tipe == "LC") diag_LC[[paste0("Lambda_", lam)]] else diag_CC[[paste0("Lambda_", lam)]]
    cat(sprintf("\n>>> Jumlah dan Penyebab OOC Bagan Kendali MREWMA-%s (%s) <<<\n", tipe, lam))
    cat(sprintf("Jumlah OOC: %d dari %d subgrup\n", nrow(h), k_sg))
    if (nrow(h) > 0) print(h[, kolom_ooc], row.names = FALSE) else cat("Tidak ada sinyal OOC.\n")
  }
}

cat("\n--- TABEL NILAI MRL0 (Shift x Skenario) ---\n")
ambil_mrl_urut = function(sk, skm, lam) {
  grid = skenario[[sk]]
  sapply(grid, function(g) {
    br = which(hasil_mrl$Skenario == sk & hasil_mrl$Skema == skm &
                 hasil_mrl$Lambda == lam & hasil_mrl$Shift == g)
    if (length(br) == 0) NA else hasil_mrl$MRL[br[1]]
  })
}
maks_posisi = max(sapply(skenario, length))
for (i in 1:nrow(tabel_ucl)) {
  lam = tabel_ucl$Lambda[i]
  for (skm in c("MREWMA-LC", "MREWMA-CC")) {
    cat(sprintf("\n>>> Nilai MRL0 - %s, Lambda = %s <<<\n", skm, lam))
    tabel_lebar = data.frame(Shift = 1:maks_posisi)
    for (sk in names(skenario)) {
      nilai = ambil_mrl_urut(sk, skm, lam)
      length(nilai) = maks_posisi
      tabel_lebar[[sk]] = round(nilai, 3)
    }
    print(tabel_lebar, row.names = FALSE)
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
    Descriptive_Statistics = as.data.frame(stats_all),
    Normality_Test        = normality_table,
    Phase1_Cleaning_Log   = riwayat_bersih,
    UCL_Table             = tabel_ucl,
    MRL0_Results          = hasil_mrl,
    OOC_Diagnostics_LC    = if (is.null(diag_LC_all)) data.frame() else diag_LC_all,
    OOC_Diagnostics_CC    = if (is.null(diag_CC_all)) data.frame() else diag_CC_all
  ),
  output_file
)

cat("\nSELESAI. Hasil disimpan di:\n")
cat("  Excel  ->", output_file, "\n")
cat("  Gambar ->", plot_dir, "\n")
