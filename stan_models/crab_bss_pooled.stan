// =============================================================================
// Pooled CPUE Crab Creel Model
// Adaptive-resolution AR(1) + CPUE Day-Type Effect + L_effective Uncertainty
//
// The AR(1) processes for effort and CPUE evolve over P_n temporal periods.
// The R preprocessing selects the resolution (daily, weekly, or monthly)
// based on data density for each population × sub-season, then sets:
//   P_n = number of AR periods (= D for daily, ~D/7 for weekly, etc.)
//   period[d] = mapping from each day to its AR period index
//
// When P_n = D and period[d] = d, this is equivalent to a daily AR(1).
// When P_n = n_months and period maps days to months, it behaves as
// a monthly AR(1). The Stan code is identical in both cases.
//
// Additional features:
//   - B1_C: weekend CPUE effect
//   - L_effective as parameter with lognormal prior (shore only)
//   - Data-driven R_G prior from interview data
//   - Informative R_T prior Beta(alpha, beta)
//   - Sparse overdispersion (observation-indexed)
//   - I/E direct effort integration
//   - Dual reporting: expected catch + predictive draws
//
// v6.6 (B1.3): the AR(1) initial states omega_E_0 / omega_C_0 are non-centered
//   (omega_*_0 = stationary SD x raw) to remove the centered funnel that drove
//   the boat divergences (~98% of iterations). The implied prior is unchanged,
//   so the posterior is the same; only the sampling geometry improves.
//
// v6.7 (B1.5): the per-observation effort overdispersion is marginalized. The
//   gamma-Poisson form Poisson(lambda * eps_E_H_obs * R) with
//   eps_E_H_obs ~ Gamma(r_E, r_E) is replaced by its exact marginal,
//   neg_binomial_2(lambda * R, r_E), for both gear and trailer counts. This
//   removes the n_effort_obs latent eps_E_H_obs parameters, a centered
//   high-dimensional funnel (neck at large r_E / small sigma_r) that survived
//   the B1.3 non-centering and produced the residual shore divergences. The
//   marginalization is exact (gamma-Poisson == negative binomial), so the
//   posterior over every reported quantity is unchanged; only nuisance
//   parameters and the funnel are removed. The interview catch likelihood
//   already used neg_binomial_2(.., r_C); the effort counts now match it.
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
  real<lower=0> sigma_mu_E;
  matrix[G,S] eps_mu_E;

  // B1.5: eps_E_H_obs[n_effort_obs] removed. The effort-count overdispersion is
  //       now marginalized into neg_binomial_2 in the model block, so the per-
  //       observation latent effects (and their funnel) no longer exist as
  //       parameters. n_effort_obs is retained in the data block only for
  //       R-interface compatibility; it is no longer referenced.

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
  real<lower=0> sigma_mu_C;
  matrix[G,S] eps_mu_C;

  vector[D * estimate_L] L_raw;
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

  // --- B1.3: non-centered AR(1) initial state. omega_*_0 = stationary SD x raw,
  //     which reproduces the original normal(0, sqrt(sigma_eps^2/(1-phi^2)))
  //     prior but moves the sigma_eps/phi-dependent scale into a deterministic
  //     transform. This removes the centered funnel that produced near-total
  //     boat divergences (treedepth 0; adapt_delta 0.99 could not fix it). ---
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
      mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E;
      mu_C[g,s] = mu_mu_C[g] + eps_mu_C[g,s] * sigma_mu_C;
    }
    for (d in 1:D) {
      for (s in 1:S) {
        // Effort: AR deviation (at period resolution) + weekend + holiday
        lambda_E_S[s][d,g] = exp(mu_E[g,s] +
          to_matrix(omega_E[period[d],], G, S)[g,s] + B1 * w[d] + B2 * holiday[d]) * O[d,s,g];
        // CPUE: AR deviation + weekend CPUE effect
        lambda_C_S[s][d,g] = exp(mu_C[g,s] +
          to_matrix(omega_C[period[d],], G, S)[g,s] + B1_C * w[d]) * O[d,s,g];
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
  sigma_mu_E ~ cauchy(0, value_cauchyDF_sigma_mu_E);
  sigma_mu_C ~ cauchy(0, value_cauchyDF_sigma_mu_C);
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
      eps_mu_E[g,s] ~ std_normal();
      eps_mu_C[g,s] ~ std_normal();
    }
  }

  // B1.5: effort-count overdispersion marginalized to neg_binomial_2. The
  //       previous form was Gear_I ~ Poisson(lambda * eps_E_H_obs * R_G) with
  //       eps_E_H_obs ~ Gamma(r_E, r_E); integrating out eps_E_H_obs gives
  //       neg_binomial_2(lambda * R_G, r_E) exactly (mean lambda*R_G, variance
  //       lambda*R_G + (lambda*R_G)^2 / r_E). r_E = 1 / sigma_r_E^2 is
  //       unchanged, so the overdispersion is identical; only the latent per-
  //       observation eps parameters (and their centered funnel) are removed.
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

  if (IE_n > 0) {
    sigma_IE ~ exponential(5);
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
}
