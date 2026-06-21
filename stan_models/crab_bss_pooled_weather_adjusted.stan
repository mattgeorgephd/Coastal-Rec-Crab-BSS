// =============================================================================
// Pooled CPUE Crab Creel Model — WEATHER-ADJUSTED VARIANT
// =============================================================================
//
// This is a surgical extension of crab_bss_pooled.stan. It adds daily tide and
// weather covariates to the log-linear predictors for effort (lambda_E) and
// CPUE (lambda_C), and adds pointwise log-likelihood generation (log_lik) so
// that PSIS-LOO can be used to compare baseline vs weather-adjusted fits.
//
// All structural features of the baseline pooled model are preserved verbatim:
//   - Adaptive-resolution AR(1) over P_n periods
//   - Cholesky correlation structure on AR innovations (Lcorr_E, Lcorr_C)
//   - B1_C weekend CPUE effect
//   - L_effective non-centered lognormal (shore) or fixed (boat)
//   - Lognormal R_G prior (gear-per-crabber, can exceed 1)
//   - Beta R_T prior (trailer-per-group)
//   - Sparse per-observation gamma overdispersion eps_E_H_obs
//   - I/E direct effort integration via lognormal anchor
//   - Negative-binomial CPUE likelihood with dispersion r_C
//
// Collapses to the baseline pooled model when K_E = K_C = 0 (zero-column
// covariate matrices). Posteriors in that regime should be statistically
// identical to crab_bss_pooled.stan up to MCMC noise. Verify before trusting
// covariate-augmented inference.
//
// -----------------------------------------------------------------------------
// CHANGES FROM BASELINE (search "WEATHER-ADJ" to locate each)
// -----------------------------------------------------------------------------
//   1. data block:     K_E, K_C, X_E, X_C, prior_sd_gamma
//   2. parameters:     gamma_E[K_E], gamma_C[K_C]
//   3. transformed:    cov_contrib_E and cov_contrib_C added inside exp() of
//                      lambda_E_S[s][d,g] and lambda_C_S[s][d,g]
//   4. model block:    gamma_E, gamma_C priors with K>0 guards
//   5. generated:      log_lik vector for PSIS-LOO
//
// -----------------------------------------------------------------------------
// NOTE ON log_lik AND OVERDISPERSION
// -----------------------------------------------------------------------------
// As of B1.5 (see version note below), effort counts (Gear_I, T_I) use the
// marginal Negative Binomial 2 form, neg_binomial_2(lambda * R, r_E), directly
// in the MODEL block. Earlier versions used the conditional Poisson with per-
// observation random effects eps_E_H_obs[i] ~ Gamma(r_E, r_E); the gamma-
// Poisson mixture integrates out to exactly this neg_binomial_2, so the change
// is inference-preserving and removed a centered funnel. log_lik has always
// used this marginal form (correct for the predictive density of a NEW
// observation whose eps_E_H_obs would be unknown; Vehtari et al. 2017), so the
// model block and log_lik are now identical in form for effort counts. The
// ELPD difference between baseline and weather-adjusted models is valid as long
// as both fits use the same log_lik construction, which they do here.
//
// Interview catches use neg_binomial_2(..., r_C) in both the model block and
// log_lik. Gear expansion, trailer expansion, and I/E anchors have no latent
// random effects; log_lik for those uses their native likelihood directly.
//
// v6.6 (B1.3): the AR(1) initial states omega_E_0 / omega_C_0 are non-centered
//   (omega_*_0 = stationary SD x raw), matching crab_bss_pooled.stan, to remove
//   the centered funnel behind the boat divergences. Inference is unchanged.
//
// v6.7 (B1.5): effort-count overdispersion marginalized to neg_binomial_2 in
//   the model block, matching crab_bss_pooled.stan. The eps_E_H_obs parameters
//   and their Gamma(r_E, r_E) prior are removed; the gamma-Poisson mixture
//   equals neg_binomial_2(lambda * R, r_E) exactly, so inference is unchanged
//   while a centered high-dimensional funnel is removed. log_lik was already
//   on this marginal form, so the model block and log_lik now match.
//
// v6.8 (B1.6): sigma_IE given a proper prior unconditionally (matching
//   crab_bss_pooled.stan). Fixes the improper flat sigma_IE direction when
//   IE_n = 0 (the boat), which drifted to ~1e307 and drove divergences.
//   Inference-preserving for the reported quantities.
//
// v6.9 (B1.7): single-cell hierarchical scale collapsed (matching
//   crab_bss_pooled.stan). sigma_mu_E / eps_mu_E and the C counterparts removed;
//   mu_E = mu_mu_E ~ normal(prior). With G = S = 1 the scale was an unidentified
//   single-cell funnel and, after B1.6, the boat's last divergence source. NOT
//   inference-preserving (removes the HalfCauchy variance layer, tightening the
//   still-wide prior tails on mu_E); restore the hierarchy for any G>1 / S>1 use.
// =============================================================================

data {
  int<lower=1> D;                        // Number of days in sub-season
  int<lower=1> G;                        // Number of gear groups (1 in pooled)
  int<lower=1> S;                        // Number of sections (1)
  int<lower=1> P_n;                      // Number of AR periods
  int<lower=1> period[D];               // Day-to-period mapping
  vector<lower=0,upper=1>[D] w;          // Weekend indicator
  vector<lower=0,upper=1>[D] holiday;    // Holiday indicator
  real<lower=0> O[D,S,G];               // Open/closed

  // --- Day length ---
  vector<lower=0>[D] L_data;
  int<lower=0,upper=1> estimate_L;
  vector<lower=0>[D] L_prior_sigma;

  // --- Effort observations (sparse overdispersion) ---
  int<lower=0> n_effort_obs;

  int<lower=0> Gear_n;
  int<lower=1> day_Gear[Gear_n];
  int<lower=1> section_Gear[Gear_n];
  int<lower=0> Gear_I[Gear_n];

  int<lower=0> T_n;
  int<lower=1> day_T[T_n];
  int<lower=1> section_T[T_n];
  int<lower=0> T_I[T_n];

  int<lower=0> Crab_n;
  int<lower=1> day_Crab[Crab_n];
  int<lower=1> section_Crab[Crab_n];
  int<lower=0> Crab_I[Crab_n];
  real<lower=0,upper=1> p_I_crab;

  // --- Interview CPUE ---
  int<lower=0> IntC;
  int<lower=1> day_IntC[IntC];
  int<lower=1> gear_IntC[IntC];
  int<lower=1> section_IntC[IntC];
  int<lower=0> c[IntC];
  vector<lower=0>[IntC] h;

  // --- Gear-per-crabber interviews ---
  int<lower=0> IntA_gear;
  int<lower=0> Gear_A[IntA_gear];
  int<lower=1> A_A_gear[IntA_gear];

  // --- Trailer-per-group interviews ---
  int<lower=0> IntA_trailer;
  int<lower=0> T_A_int[IntA_trailer];
  int<lower=1> A_A_trailer[IntA_trailer];

  // --- I/E direct effort observations ---
  int<lower=0> IE_n;
  int<lower=1> day_IE[IE_n];
  int<lower=1> section_IE[IE_n];
  vector<lower=0>[IE_n] IE_crabber_hours;

  // --- Hyperparameters ---
  real value_cauchyDF_sigma_eps_E;
  real value_cauchyDF_sigma_eps_C;
  real value_cauchyDF_sigma_r_E;
  real value_cauchyDF_sigma_r_C;
  real value_betashape_phi_E_scaled;
  real value_betashape_phi_C_scaled;
  real value_normal_sigma_B1;
  real value_normal_sigma_B2;
  real value_normal_sigma_B1_C;
  real value_normal_mu_mu_C;
  real value_normal_sigma_mu_C;
  real value_normal_mu_mu_E;
  real value_normal_sigma_mu_E;
  real value_cauchyDF_sigma_mu_C;
  real value_cauchyDF_sigma_mu_E;

  // --- Data-driven priors ---
  real<lower=0> R_G_prior_mu;
  real<lower=0> R_G_prior_sigma;
  real<lower=0> R_T_alpha;
  real<lower=0> R_T_beta;

  // === WEATHER-ADJ: Covariate inputs ========================================
  int<lower=0> K_E;                      // Number of effort covariates (0 = baseline)
  int<lower=0> K_C;                      // Number of CPUE covariates (0 = baseline)
  matrix[D, K_E] X_E;                    // Standardized daily effort covariates
  matrix[D, K_C] X_C;                    // Standardized daily CPUE covariates
  real<lower=0> prior_sd_gamma;          // Prior SD on gamma_E and gamma_C (typ 0.35)
  // ==========================================================================
}

transformed data {
  // === WEATHER-ADJ: pre-compute total observation count for log_lik =========
  // Declared here (not in generated quantities) so the vector[N_obs] log_lik
  // declaration has a resolvable dimension at Stan compile time.
  int N_obs = Gear_n + T_n + IntC + IntA_gear + IntA_trailer + IE_n;
  // ==========================================================================
}

parameters {
  real B1;
  real B2;
  real B1_C;

  real<lower=0> sigma_eps_E;
  cholesky_factor_corr[G*S] Lcorr_E;
  real<lower=0> sigma_r_E;
  real<lower=0,upper=1> phi_E_scaled;
  matrix[P_n-1, G*S] eps_E;
  matrix[G,S] omega_E_0_raw;   // B1.3: non-centered AR(1) initial state (raw)
  real mu_mu_E[G];
  // B1.7: sigma_mu_E and eps_mu_E removed. With G=1, S=1 (this pooled model is
  //       always single-cell) the hierarchical scale is a single-cell funnel:
  //       mu_E = mu_mu_E + eps_mu_E*sigma_mu_E is 3 parameters for 1 identified
  //       quantity, sigma_mu_E prior-dominated and unidentified. mu_E now equals
  //       mu_mu_E directly (see transformed parameters).

  // B1.5: eps_E_H_obs[n_effort_obs] removed; effort-count overdispersion is now
  //       marginalized into neg_binomial_2 in the model block. n_effort_obs is
  //       retained in the data block only for R-interface compatibility.

  real<lower=0> R_G;
  real<lower=0,upper=1> R_T;

  real<lower=0> sigma_IE;

  real<lower=0> sigma_eps_C;
  cholesky_factor_corr[G*S] Lcorr_C;
  real<lower=0,upper=1> phi_C_scaled;
  real<lower=0> sigma_r_C;
  matrix[P_n-1, G*S] eps_C;
  matrix[G,S] omega_C_0_raw;   // B1.3: non-centered AR(1) initial state (raw)
  real mu_mu_C[G];
  // B1.7: sigma_mu_C and eps_mu_C removed (same single-cell funnel; mu_C = mu_mu_C).

  vector[D * estimate_L] L_raw;

  // === WEATHER-ADJ: Covariate coefficients ==================================
  vector[K_E] gamma_E;                   // Log-scale effect of each X_E column on lambda_E
  vector[K_C] gamma_C;                   // Log-scale effect of each X_C column on lambda_C
  // ==========================================================================
}

transformed parameters {
  matrix[G,S] mu_E;
  real<lower=-1,upper=1> phi_E;
  matrix[P_n, G*S] omega_E;
  matrix[G,S] omega_E_0;
  matrix<lower=0>[D,G] lambda_E_S[S];
  real<lower=0> r_E;

  matrix[G,S] mu_C;
  real<lower=-1,upper=1> phi_C;
  real<lower=0> r_C;
  matrix[P_n, G*S] omega_C;
  matrix[G,S] omega_C_0;
  matrix<lower=0>[D,G] lambda_C_S[S];

  vector<lower=0>[D] L;

  // --- Compute L ---
  if (estimate_L == 1) {
    for (d in 1:D)
      L[d] = L_data[d] * exp(L_prior_sigma[d] * L_raw[d]);
  } else {
    L = L_data;
  }

  r_E = 1 / square(sigma_r_E);
  r_C = 1 / square(sigma_r_C);
  phi_E = (phi_E_scaled * 2) - 1;
  phi_C = (phi_C_scaled * 2) - 1;

  // --- B1.3: non-centered AR(1) initial state (see crab_bss_pooled.stan). ---
  omega_E_0 = sqrt(square(sigma_eps_E) / (1 - square(phi_E))) * omega_E_0_raw;
  omega_C_0 = sqrt(square(sigma_eps_C) / (1 - square(phi_C))) * omega_C_0_raw;

  // --- AR(1) over P_n periods ---
  omega_E[1,] = to_row_vector(omega_E_0);
  omega_C[1,] = to_row_vector(omega_C_0);
  for (p in 2:P_n) {
    omega_E[p,] = to_row_vector(phi_E * to_vector(omega_E[p-1,]) +
      diag_pre_multiply(rep_vector(sigma_eps_E, G*S), Lcorr_E) * to_vector(eps_E[p-1,]));
    omega_C[p,] = to_row_vector(phi_C * to_vector(omega_C[p-1,]) +
      diag_pre_multiply(rep_vector(sigma_eps_C, G*S), Lcorr_C) * to_vector(eps_C[p-1,]));
  }

  for (g in 1:G) {
    for (s in 1:S) {
      mu_E[g,s] = mu_mu_E[g];   // B1.7: single-cell collapse (S=1, G=1)
      mu_C[g,s] = mu_mu_C[g];   // B1.7: single-cell collapse
    }
    for (d in 1:D) {
      for (s in 1:S) {
        // === WEATHER-ADJ: covariate contribution inlined into exp() ========
        // When K_E = 0 or K_C = 0, the row-vector-by-vector product of empty
        // objects is 0.0 in Stan; the ternary is a safety belt for versions
        // that warn on zero-dimensional operations. This line is equivalent
        // to the baseline pooled model when K_E = K_C = 0.
        // ===================================================================

        // Effort: AR deviation (at period resolution) + weekend + holiday + covariates
        lambda_E_S[s][d,g] = exp(mu_E[g,s] +
          to_matrix(omega_E[period[d],], G, S)[g,s] + B1 * w[d] + B2 * holiday[d] +
          ((K_E > 0) ? (X_E[d,] * gamma_E) : 0.0)) * O[d,s,g];
        // CPUE: AR deviation + weekend CPUE effect + covariates
        lambda_C_S[s][d,g] = exp(mu_C[g,s] +
          to_matrix(omega_C[period[d],], G, S)[g,s] + B1_C * w[d] +
          ((K_C > 0) ? (X_C[d,] * gamma_C) : 0.0)) * O[d,s,g];
      }
    }
  }
}

model {
  sigma_eps_E ~ cauchy(0, value_cauchyDF_sigma_eps_E);
  sigma_eps_C ~ cauchy(0, value_cauchyDF_sigma_eps_C);
  Lcorr_E ~ lkj_corr_cholesky(1);
  Lcorr_C ~ lkj_corr_cholesky(1);
  phi_E_scaled ~ beta(value_betashape_phi_E_scaled, value_betashape_phi_E_scaled);
  phi_C_scaled ~ beta(value_betashape_phi_C_scaled, value_betashape_phi_C_scaled);
  sigma_r_E ~ cauchy(0, value_cauchyDF_sigma_r_E);
  sigma_r_C ~ cauchy(0, value_cauchyDF_sigma_r_C);
  // B1.7: sigma_mu_E / sigma_mu_C priors removed (parameters no longer exist).
  //       value_cauchyDF_sigma_mu_* remain in the data block for R-interface
  //       compatibility and are now unused.
  B1 ~ normal(0, value_normal_sigma_B1);
  B2 ~ normal(0, value_normal_sigma_B2);
  B1_C ~ normal(0, value_normal_sigma_B1_C);

  to_vector(eps_E) ~ std_normal();
  to_vector(eps_C) ~ std_normal();

  R_G ~ lognormal(log(R_G_prior_mu), R_G_prior_sigma);

  if (T_n > 0 || IntA_trailer > 0) {
    R_T ~ beta(R_T_alpha, R_T_beta);
  }

  if (estimate_L == 1) {
    L_raw ~ std_normal();
  }

  for (g in 1:G) {
    mu_mu_E[g] ~ normal(value_normal_mu_mu_E, value_normal_sigma_mu_E);
    mu_mu_C[g] ~ normal(value_normal_mu_mu_C, value_normal_sigma_mu_C);
    for (s in 1:S) {
      omega_E_0_raw[g,s] ~ std_normal();   // B1.3: prior on raw; omega_*_0 scaled in TP
      omega_C_0_raw[g,s] ~ std_normal();
    }
  }

  // B1.5: effort-count overdispersion marginalized to neg_binomial_2 (see the
  //       gear/trailer loops below and the header note). eps_E_H_obs and its
  //       Gamma(r_E, r_E) prior are removed; r_E = 1 / sigma_r_E^2 is unchanged.

  // === WEATHER-ADJ: covariate priors ========================================
  if (K_E > 0) gamma_E ~ normal(0, prior_sd_gamma);
  if (K_C > 0) gamma_C ~ normal(0, prior_sd_gamma);
  // ==========================================================================

  for (i in 1:Gear_n) {
    Gear_I[i] ~ neg_binomial_2(
      lambda_E_S[section_Gear[i]][day_Gear[i], 1] * R_G, r_E
    );
  }

  for (i in 1:T_n) {
    T_I[i] ~ neg_binomial_2(
      lambda_E_S[section_T[i]][day_T[i], G] * R_T, r_E
    );
  }

  for (a in 1:IntC) {
    c[a] ~ neg_binomial_2(
      lambda_C_S[section_IntC[a]][day_IntC[a], gear_IntC[a]] * h[a], r_C
    );
  }

  for (a in 1:IntA_gear) {
    Gear_A[a] ~ poisson(A_A_gear[a] * R_G);
  }

  for (a in 1:IntA_trailer) {
    T_A_int[a] ~ bernoulli(R_T);
  }

  // B1.6: sigma_IE gets a proper prior unconditionally. When IE_n = 0 (the boat
  //       has no I/E observations) the old code left sigma_IE with no prior and
  //       no likelihood: an improper flat direction that drifted to ~1e307 and
  //       was the boat's dominant divergence source (the divergence diagnostic
  //       found sigma_IE at the floating-point ceiling for the boat). sigma_IE
  //       is decoupled from effort and catch, so this is inference-preserving
  //       for E and C; it only makes the posterior proper and the sampler sane.
  sigma_IE ~ exponential(5);
  if (IE_n > 0) {
    for (i in 1:IE_n) {
      IE_crabber_hours[i] ~ lognormal(
        log(lambda_E_S[section_IE[i]][day_IE[i], 1] * L[day_IE[i]]),
        sigma_IE
      );
    }
  }
}

generated quantities {
  matrix[G*S, G*S] Omega_C;
  matrix[G*S, G*S] Omega_E;

  matrix<lower=0>[D,G] lambda_Ctot_S[S];
  matrix<lower=0>[D,G] C_expected[S];
  real<lower=0> C_expected_sum;

  matrix<lower=0>[D,G] C[S];
  real<lower=0> C_sum;

  matrix<lower=0>[D,G] E[S];
  real<lower=0> E_sum;

  real R_G_out;
  real sigma_IE_out;
  real B1_C_out;
  vector[D] L_out;

  // === WEATHER-ADJ: pointwise log likelihood for PSIS-LOO ===================
  // N_obs is computed in transformed data. Ordering (must match R-side
  // mask_stan_data() in BSS-GH-weather-tide-covariates.Rmd):
  //   1. Gear counts           1                              .. Gear_n
  //   2. Trailer counts        Gear_n+1                       .. Gear_n+T_n
  //   3. Interview catches     Gear_n+T_n+1                   .. Gear_n+T_n+IntC
  //   4. Gear expansions       Gear_n+T_n+IntC+1              .. ... +IntA_gear
  //   5. Trailer expansions    ...+IntA_gear+1                .. ... +IntA_trailer
  //   6. I/E anchor obs        ...+IntA_trailer+1             .. ... +IE_n
  vector[N_obs] log_lik;
  // ==========================================================================

  Omega_C = multiply_lower_tri_self_transpose(Lcorr_C);
  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  sigma_IE_out = sigma_IE;
  B1_C_out = B1_C;
  L_out = L;

  C_sum = 0;
  C_expected_sum = 0;
  E_sum = 0;

  for (g in 1:G) {
    for (d in 1:D) {
      for (s in 1:S) {
        lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * L[d] * lambda_C_S[s][d,g];
        C_expected[s][d,g] = lambda_Ctot_S[s][d,g];
        C_expected_sum = C_expected_sum + C_expected[s][d,g];

        if (lambda_Ctot_S[s][d,g] < 1e9) {
          C[s][d,g] = poisson_rng(lambda_Ctot_S[s][d,g]);
        } else {
          C[s][d,g] = lambda_Ctot_S[s][d,g];
        }
        C_sum = C_sum + C[s][d,g];

        E[s][d,g] = lambda_E_S[s][d,g] * L[d];
        E_sum = E_sum + E[s][d,g];
      }
    }
  }

  // === WEATHER-ADJ: pointwise log likelihood construction ===================
  // CRITICAL: these expressions must produce a valid predictive density for
  // each observation. For effort counts (Gear_I, T_I) we use the MARGINAL
  // neg_binomial_2 form rather than the conditional Poisson used in the model
  // block, because eps_E_H_obs is observation-specific and unavailable for a
  // held-out observation. The Poisson(lambda * eps_E_H * R) with
  // eps_E_H ~ Gamma(r_E, r_E) integrates to neg_binomial_2(lambda * R, r_E).
  //
  // Other observations have no latent random effects; log_lik matches the
  // model block directly.
  {
    int idx = 0;

    // 1. Gear counts (marginal neg_binomial_2 over eps_E_H)
    for (i in 1:Gear_n) {
      idx = idx + 1;
      log_lik[idx] = neg_binomial_2_lpmf(Gear_I[i] |
        lambda_E_S[section_Gear[i]][day_Gear[i], 1] * R_G,
        r_E);
    }

    // 2. Trailer counts (marginal neg_binomial_2 over eps_E_H)
    for (i in 1:T_n) {
      idx = idx + 1;
      log_lik[idx] = neg_binomial_2_lpmf(T_I[i] |
        lambda_E_S[section_T[i]][day_T[i], G] * R_T,
        r_E);
    }

    // 3. Interview catches (neg_binomial_2 direct)
    for (a in 1:IntC) {
      idx = idx + 1;
      log_lik[idx] = neg_binomial_2_lpmf(c[a] |
        lambda_C_S[section_IntC[a]][day_IntC[a], gear_IntC[a]] * h[a],
        r_C);
    }

    // 4. Gear-per-crabber expansion (Poisson direct)
    for (a in 1:IntA_gear) {
      idx = idx + 1;
      log_lik[idx] = poisson_lpmf(Gear_A[a] | A_A_gear[a] * R_G);
    }

    // 5. Trailer-per-group expansion (Bernoulli direct)
    for (a in 1:IntA_trailer) {
      idx = idx + 1;
      log_lik[idx] = bernoulli_lpmf(T_A_int[a] | R_T);
    }

    // 6. I/E anchor (lognormal direct)
    if (IE_n > 0) {
      for (i in 1:IE_n) {
        idx = idx + 1;
        log_lik[idx] = lognormal_lpdf(IE_crabber_hours[i] |
          log(lambda_E_S[section_IE[i]][day_IE[i], 1] * L[day_IE[i]]),
          sigma_IE);
      }
    }
  }
}
