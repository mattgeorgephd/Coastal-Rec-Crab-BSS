// =============================================================================
// Pooled CPUE Crab Creel Model — Sparse Overdispersion + I/E Integration
// (crab_bss_pooled.stan)
//
// RESTRUCTURED eps_E_H: Within-day overdispersion parameters are allocated
// ONLY for actual effort observations, not for all D*H*S*G combinations.
//
// INGRESS/EGRESS INTEGRATION (Option 2):
// On days with I/E surveys, observed crabber-hours enter as a direct
// lognormal observation of lambda_E * L, bypassing R_G and day-length
// assumptions. This calibrates the gear-count pathway and provides
// high-confidence anchor points for the effort trajectory.
//
// Single CPUE process shared across all gear types.
// Holiday effect B2 separates holiday effort from regular weekends.
// Effort: log(lambda_E) = mu + omega + B1*weekend + B2*holiday
// =============================================================================

data {
  int<lower=1> D;
  int<lower=1> G;
  int<lower=1> S;
  int<lower=1> P_n;
  int<lower=1> period[D];
  vector<lower=0,upper=1>[D] w;
  vector<lower=0,upper=1>[D] holiday;
  vector<lower=0>[D] L;
  real<lower=0> O[D,S,G];

  // Total effort observations (Gear_n + T_n). Each gets one eps_E_H_obs.
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

  int<lower=0> IntC;
  int<lower=1> day_IntC[IntC];
  int<lower=1> gear_IntC[IntC];
  int<lower=1> section_IntC[IntC];
  int<lower=0> c[IntC];
  vector<lower=0>[IntC] h;

  int<lower=0> IntA_gear;
  int<lower=0> Gear_A[IntA_gear];
  int<lower=1> A_A_gear[IntA_gear];

  int<lower=0> IntA_trailer;
  int<lower=0> T_A_int[IntA_trailer];
  int<lower=1> A_A_trailer[IntA_trailer];

  // Ingress/egress direct crabber-hours observations
  int<lower=0> IE_n;
  int<lower=1> day_IE[IE_n];
  int<lower=1> section_IE[IE_n];
  vector<lower=0>[IE_n] IE_crabber_hours;

  real value_cauchyDF_sigma_eps_E;
  real value_cauchyDF_sigma_eps_C;
  real value_cauchyDF_sigma_r_E;
  real value_cauchyDF_sigma_r_C;
  real value_betashape_phi_E_scaled;
  real value_betashape_phi_C_scaled;
  real value_normal_sigma_B1;
  real value_normal_sigma_B2;
  real value_normal_mu_mu_C;
  real value_normal_sigma_mu_C;
  real value_normal_mu_mu_E;
  real value_normal_sigma_mu_E;
  real value_cauchyDF_sigma_mu_C;
  real value_cauchyDF_sigma_mu_E;
}

parameters {
  real B1;
  real B2;
  real<lower=0> sigma_eps_E;
  cholesky_factor_corr[G*S] Lcorr_E;
  real<lower=0> sigma_r_E;
  real<lower=0,upper=1> phi_E_scaled;
  matrix[P_n-1, G*S] eps_E;
  matrix[G,S] omega_E_0;
  real mu_mu_E[G];
  real<lower=0> sigma_mu_E;
  matrix[G,S] eps_mu_E;

  // Sparse overdispersion: one per actual effort observation
  vector<lower=0>[n_effort_obs] eps_E_H_obs;

  real<lower=0> R_G;
  real<lower=0,upper=1> R_T;

  // I/E measurement error (log scale)
  real<lower=0> sigma_IE;

  real<lower=0> sigma_eps_C;
  cholesky_factor_corr[G*S] Lcorr_C;
  real<lower=0,upper=1> phi_C_scaled;
  real<lower=0> sigma_r_C;
  matrix[P_n-1, G*S] eps_C;
  matrix[G,S] omega_C_0;
  real mu_mu_C[G];
  real<lower=0> sigma_mu_C;
  matrix[G,S] eps_mu_C;
}

transformed parameters {
  matrix[G,S] mu_E;
  real<lower=-1,upper=1> phi_E;
  matrix[P_n, G*S] omega_E;
  matrix<lower=0>[D,G] lambda_E_S[S];
  real<lower=0> r_E;

  matrix[G,S] mu_C;
  real<lower=-1,upper=1> phi_C;
  real<lower=0> r_C;
  matrix[P_n, G*S] omega_C;
  matrix<lower=0>[D,G] lambda_C_S[S];

  r_E = 1 / square(sigma_r_E);
  r_C = 1 / square(sigma_r_C);
  phi_E = (phi_E_scaled * 2) - 1;
  phi_C = (phi_C_scaled * 2) - 1;

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
        lambda_E_S[s][d,g] = exp(mu_E[g,s] +
          to_matrix(omega_E[period[d],], G, S)[g,s] + B1 * w[d] + B2 * holiday[d]) * O[d,s,g];
        lambda_C_S[s][d,g] = exp(mu_C[g,s] +
          to_matrix(omega_C[period[d],], G, S)[g,s]) * O[d,s,g];
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

  to_vector(eps_E) ~ std_normal();
  to_vector(eps_C) ~ std_normal();

  R_G ~ lognormal(log(1.3), 0.3);
  if (T_n > 0 || IntA_trailer > 0) {
    R_T ~ beta(0.5, 0.5);
  }

  for (g in 1:G) {
    mu_mu_E[g] ~ normal(value_normal_mu_mu_E, value_normal_sigma_mu_E);
    mu_mu_C[g] ~ normal(value_normal_mu_mu_C, value_normal_sigma_mu_C);
    for (s in 1:S) {
      omega_E_0[g,s] ~ normal(0, sqrt(square(sigma_eps_E) / (1 - square(phi_E))));
      omega_C_0[g,s] ~ normal(0, sqrt(square(sigma_eps_C) / (1 - square(phi_C))));
      eps_mu_E[g,s] ~ std_normal();
      eps_mu_C[g,s] ~ std_normal();
    }
  }

  // Sparse overdispersion prior: vectorized over actual observations only
  eps_E_H_obs ~ gamma(r_E, r_E);

  // --- Gear counts: obs indices 1..Gear_n ---
  for (i in 1:Gear_n) {
    Gear_I[i] ~ poisson(
      lambda_E_S[section_Gear[i]][day_Gear[i], 1] * eps_E_H_obs[i] * R_G
    );
  }

  // --- Trailer counts: obs indices Gear_n+1..Gear_n+T_n ---
  for (i in 1:T_n) {
    T_I[i] ~ poisson(
      lambda_E_S[section_T[i]][day_T[i], G] * eps_E_H_obs[Gear_n + i] * R_T
    );
  }

  // --- Interview CPUE ---
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

  // --- I/E direct effort observations ---
  // I/E crabber-hours are a direct measurement of lambda_E * L (daily effort)
  // with lognormal measurement error. This bypasses R_G and day-length assumptions.
  if (IE_n > 0) {
    sigma_IE ~ exponential(5);  // prior: ~0.2 on log scale (±20% measurement error)
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
  matrix<lower=0>[D,G] C[S];
  matrix<lower=0>[D,G] E[S];
  real<lower=0> C_sum;
  real<lower=0> E_sum;
  real R_G_out;
  real sigma_IE_out;

  Omega_C = multiply_lower_tri_self_transpose(Lcorr_C);
  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  sigma_IE_out = sigma_IE;
  C_sum = 0;
  E_sum = 0;

  for (g in 1:G) {
    for (d in 1:D) {
      for (s in 1:S) {
        lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * L[d] * lambda_C_S[s][d,g];
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
