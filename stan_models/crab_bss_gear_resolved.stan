// =============================================================================
// Gear-Resolved CPUE Crab Creel Model (crab_bss_gear_resolved.stan)
//
// Adapted from WDFW FW CreelEstimates BSS model
// Per-gear-type CPUE processes with independent AR(1) dynamics.
// Holiday effect B2 separates holiday effort from regular weekends.
// Per-gear-type NegBin overdispersion (r_C_gear) for each gear type.
//
// ARCHITECTURE: Shared effort × separate CPUE per gear type
//   - Single latent effort process lambda_E (total crabbers present)
//   - G_gear independent CPUE AR(1) processes lambda_C_gear[d, g_gear]
//   - Time-varying gear-type proportions pi_gear[period] (Dirichlet)
//   - Daily catch per gear type:
//     C_gear[d,g] = lambda_E[d] × L[d] × pi_gear[period[d],g] × lambda_C_gear[d,g]
//
// EFFORT SIDE: Unchanged from v2.0
//   - Gear counts (docks) with R_G expansion
//   - Trailer counts (boat launch) with R_T expansion
//   - Weekend/holiday effect B1
//   - Within-day gamma overdispersion eps_E_H
//
// CPUE SIDE: Replaced from v2.0
//   - G_gear separate AR(1) processes (shared phi_C, sigma_eps_C)
//   - Each interview's catch depends on its gear type
//   - Categorical likelihood for gear-type assignment informs pi_gear
//   - Negative binomial catch likelihood per gear type
//
// GENERATED QUANTITIES:
//   - C_gear[D, G_gear]: daily catch by gear type
//   - C_sum_gear[G_gear]: season total by gear type
//   - C[S][D,G]: total catch (backward compatible)
//   - E[S][D,G], E_sum: effort (unchanged)
// =============================================================================

data {
  // --- Core dimensions ---
  int<lower=1> D;              // days in sub-season
  int<lower=1> G;              // crabber types (always 1; shore or boat per fit)
  int<lower=1> S;              // sections (always 1)
  int<lower=1> H;              // max count sequences per day
  int<lower=1> P_n;            // time periods
  int<lower=1> period[D];      // period index per day
  vector<lower=0,upper=1>[D] w;       // weekend/holiday indicator (1 for both)
  vector<lower=0,upper=1>[D] holiday;  // holiday-only indicator (1 for holidays, 0 otherwise)
  vector<lower=0>[D] L;                // day length (hours, from suncalc)
  real<lower=0> O[D,S,G];      // open/closed status

  // --- Gear-type dimension ---
  int<lower=1> G_gear;         // number of gear types present in this sub-season

  // --- GEAR index counts (docks: counting crab gear in water) ---
  int<lower=0> Gear_n;
  int<lower=1> day_Gear[Gear_n];
  int<lower=1> section_Gear[Gear_n];
  int<lower=1> countnum_Gear[Gear_n];
  int<lower=0> Gear_I[Gear_n];

  // --- TRAILER index counts (boat launch) ---
  int<lower=0> T_n;
  int<lower=1> day_T[T_n];
  int<lower=1> section_T[T_n];
  int<lower=1> countnum_T[T_n];
  int<lower=0> T_I[T_n];

  // --- Direct crabber counts (reserved for future jetty counts) ---
  int<lower=0> Crab_n;
  int<lower=1> day_Crab[Crab_n];
  int<lower=1> section_Crab[Crab_n];
  int<lower=1> countnum_Crab[Crab_n];
  int<lower=0> Crab_I[Crab_n];
  real<lower=0,upper=1> p_I_crab;

  // --- Interview CPUE data (with gear type) ---
  int<lower=0> IntC;
  int<lower=1> day_IntC[IntC];
  int<lower=1> gear_IntC[IntC];                    // crabber type G-index (always 1)
  int<lower=1> section_IntC[IntC];
  int<lower=0> c[IntC];                            // crab caught
  vector<lower=0>[IntC] h;                         // hours fished
  int<lower=1,upper=G_gear> gear_type_IntC[IntC];  // gear type per interview

  // --- Interview expansion: gear per crabber ---
  int<lower=0> IntA_gear;
  int<lower=0> Gear_A[IntA_gear];
  int<lower=1> A_A_gear[IntA_gear];

  // --- Interview expansion: trailers per group ---
  int<lower=0> IntA_trailer;
  int<lower=0> T_A_int[IntA_trailer];
  int<lower=1> A_A_trailer[IntA_trailer];

  // --- Gear-type proportions prior ---
  vector<lower=0>[G_gear] pi_gear_alpha;  // Dirichlet concentration per gear type

  // --- Hyperparameters ---
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
  // --- Effort process ---
  real B1;                             // weekend/holiday effect
  real B2;                             // additional holiday effect (beyond B1)
  real<lower=0> sigma_eps_E;
  cholesky_factor_corr[G*S] Lcorr_E;
  real<lower=0> sigma_r_E;
  real<lower=0,upper=1> phi_E_scaled;
  matrix[P_n-1, G*S] eps_E;
  matrix[G,S] omega_E_0;
  matrix<lower=0>[D,G] eps_E_H[S,H];
  real mu_mu_E[G];
  real<lower=0> sigma_mu_E;
  matrix[G,S] eps_mu_E;

  // --- Expansion parameters ---
  real<lower=0> R_G;                   // gear per crabber
  real<lower=0,upper=1> R_T;          // trailers per boat group

  // --- CPUE process per gear type ---
  real mu_C_gear[G_gear];              // CPUE intercept per gear type (log scale)
  real<lower=0> sigma_eps_C_gear;      // shared CPUE process error SD
  real<lower=0,upper=1> phi_C_gear_scaled;  // shared AR(1) coefficient (0-1, rescaled)
  real<lower=0> sigma_r_C_gear[G_gear]; // NegBin overdispersion scale per gear type
  matrix[P_n-1, G_gear] eps_C_gear;   // CPUE innovations [periods-1 × gear types]
  vector[G_gear] omega_C_gear_0;       // initial CPUE state per gear type

  // --- Gear-type proportions ---
  simplex[G_gear] pi_gear[P_n];       // per-period gear-type mix
}

transformed parameters {
  // --- Effort (unchanged) ---
  matrix[G,S] mu_E;
  real<lower=-1,upper=1> phi_E;
  matrix[P_n, G*S] omega_E;
  matrix<lower=0>[D,G] lambda_E_S[S];
  matrix<lower=0>[D,G] lambda_E_S_I[S,H];
  real<lower=0> r_E;

  // --- CPUE by gear type ---
  real<lower=-1,upper=1> phi_C_gear;
  vector<lower=0>[G_gear] r_C_gear;
  matrix[P_n, G_gear] omega_C_gear;
  matrix<lower=0>[D, G_gear] lambda_C_gear;

  r_E = 1 / square(sigma_r_E);
  for (gg in 1:G_gear) {
    r_C_gear[gg] = 1 / square(sigma_r_C_gear[gg]);
  }
  phi_E = (phi_E_scaled * 2) - 1;
  phi_C_gear = (phi_C_gear_scaled * 2) - 1;

  // --- Effort AR(1) process (unchanged) ---
  omega_E[1,] = to_row_vector(omega_E_0);
  for (p in 2:P_n) {
    omega_E[p,] = to_row_vector(phi_E * to_vector(omega_E[p-1,]) +
      diag_pre_multiply(rep_vector(sigma_eps_E, G*S), Lcorr_E) * to_vector(eps_E[p-1,]));
  }

  // --- Gear-type CPUE AR(1) processes ---
  for (gg in 1:G_gear) {
    omega_C_gear[1, gg] = omega_C_gear_0[gg];
  }
  for (p in 2:P_n) {
    for (gg in 1:G_gear) {
      omega_C_gear[p, gg] = phi_C_gear * omega_C_gear[p-1, gg] +
        sigma_eps_C_gear * eps_C_gear[p-1, gg];
    }
  }

  // --- Daily effort rates (unchanged) ---
  for (g in 1:G) {
    for (s in 1:S) {
      mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E;
    }
    for (d in 1:D) {
      for (s in 1:S) {
        lambda_E_S[s][d,g] = exp(mu_E[g,s] +
          to_matrix(omega_E[period[d],], G, S)[g,s] + B1 * w[d] + B2 * holiday[d]) * O[d,s,g];
        for (i in 1:H) {
          lambda_E_S_I[s,i][d,g] = lambda_E_S[s][d,g] * eps_E_H[s,i][d,g];
        }
      }
    }
  }

  // --- Daily CPUE rates by gear type ---
  for (d in 1:D) {
    for (gg in 1:G_gear) {
      lambda_C_gear[d, gg] = exp(mu_C_gear[gg] + omega_C_gear[period[d], gg]);
    }
  }
}

model {
  // === EFFORT PRIORS (unchanged) ===
  sigma_eps_E ~ cauchy(0, value_cauchyDF_sigma_eps_E);
  Lcorr_E ~ lkj_corr_cholesky(1);
  phi_E_scaled ~ beta(value_betashape_phi_E_scaled, value_betashape_phi_E_scaled);
  sigma_r_E ~ cauchy(0, value_cauchyDF_sigma_r_E);
  sigma_mu_E ~ cauchy(0, value_cauchyDF_sigma_mu_E);
  B1 ~ normal(0, value_normal_sigma_B1);
  B2 ~ normal(0, value_normal_sigma_B2);  // holiday effect: additional boost beyond B1
  to_vector(eps_E) ~ std_normal();

  // Gear per crabber
  R_G ~ lognormal(log(1.3), 0.3);
  // Trailers per group (guarded)
  if (T_n > 0 || IntA_trailer > 0) {
    R_T ~ beta(0.5, 0.5);
  }

  for (g in 1:G) {
    mu_mu_E[g] ~ normal(value_normal_mu_mu_E, value_normal_sigma_mu_E);
    for (d in 1:D) {
      for (s in 1:S) {
        for (i in 1:H) {
          eps_E_H[s,i][d,g] ~ gamma(r_E, r_E);
        }
      }
    }
    for (s in 1:S) {
      omega_E_0[g,s] ~ normal(0, sqrt(square(sigma_eps_E) / (1 - square(phi_E))));
      eps_mu_E[g,s] ~ std_normal();
    }
  }

  // === CPUE PRIORS (gear-type) ===
  sigma_eps_C_gear ~ cauchy(0, value_cauchyDF_sigma_eps_C);
  phi_C_gear_scaled ~ beta(value_betashape_phi_C_scaled, value_betashape_phi_C_scaled);
  for (gg in 1:G_gear) {
    sigma_r_C_gear[gg] ~ cauchy(0, value_cauchyDF_sigma_r_C);
  }
  to_vector(eps_C_gear) ~ std_normal();

  for (gg in 1:G_gear) {
    mu_C_gear[gg] ~ normal(value_normal_mu_mu_C, value_normal_sigma_mu_C);
    omega_C_gear_0[gg] ~ normal(0, sqrt(square(sigma_eps_C_gear) / (1 - square(phi_C_gear))));
  }

  // Gear-type proportions per period
  for (p in 1:P_n) {
    pi_gear[p] ~ dirichlet(pi_gear_alpha);
  }

  // === EFFORT LIKELIHOODS (unchanged) ===

  // Gear counts (docks)
  for (i in 1:Gear_n) {
    Gear_I[i] ~ poisson(
      lambda_E_S_I[section_Gear[i], countnum_Gear[i]][day_Gear[i], 1] * R_G
    );
  }

  // Trailer counts (boat launch)
  for (i in 1:T_n) {
    T_I[i] ~ poisson(
      lambda_E_S_I[section_T[i], countnum_T[i]][day_T[i], G] * R_T
    );
  }

  // === CPUE LIKELIHOODS (gear-type indexed) ===
  for (a in 1:IntC) {
    // Gear-type assignment informs pi_gear proportions
    gear_type_IntC[a] ~ categorical(pi_gear[period[day_IntC[a]]]);
    // Catch informs gear-specific CPUE rate
    c[a] ~ neg_binomial_2(lambda_C_gear[day_IntC[a], gear_type_IntC[a]] * h[a], r_C_gear[gear_type_IntC[a]]);
  }

  // === EXPANSION LIKELIHOODS (unchanged) ===
  for (a in 1:IntA_gear) {
    Gear_A[a] ~ poisson(A_A_gear[a] * R_G);
  }
  for (a in 1:IntA_trailer) {
    T_A_int[a] ~ bernoulli(R_T);
  }
}

generated quantities {
  // Effort (unchanged structure)
  matrix[G*S, G*S] Omega_E;
  matrix<lower=0>[D,G] E[S];
  real<lower=0> E_sum;
  real R_G_out;
  vector[G_gear] r_C_gear_out;

  // Catch by gear type
  matrix<lower=0>[D, G_gear] C_gear;
  vector<lower=0>[G_gear] C_sum_gear;
  real<lower=0> C_sum;

  // Backward-compatible total catch in C[S][D,G] format
  matrix<lower=0>[D,G] C[S];

  // Pi_gear summary (expose for monitoring)
  matrix[P_n, G_gear] pi_gear_out;

  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  r_C_gear_out = r_C_gear;
  E_sum = 0;
  C_sum = 0;
  C_sum_gear = rep_vector(0, G_gear);

  // Copy pi_gear to output matrix
  for (p in 1:P_n) {
    for (gg in 1:G_gear) {
      pi_gear_out[p, gg] = pi_gear[p][gg];
    }
  }

  // Daily effort (unchanged)
  for (g in 1:G) {
    for (d in 1:D) {
      for (s in 1:S) {
        E[s][d,g] = lambda_E_S[s][d,g] * L[d];
        E_sum = E_sum + E[s][d,g];
      }
    }
  }

  // Daily catch by gear type
  for (d in 1:D) {
    real total_catch_d = 0;
    for (gg in 1:G_gear) {
      real rate = lambda_E_S[1][d,1] * L[d] * pi_gear[period[d]][gg] * lambda_C_gear[d, gg];
      if (rate > 0 && rate < 1e9) {
        C_gear[d, gg] = poisson_rng(rate);
      } else if (rate >= 1e9) {
        C_gear[d, gg] = rate;  // overflow protection
      } else {
        C_gear[d, gg] = 0;
      }
      C_sum_gear[gg] = C_sum_gear[gg] + C_gear[d, gg];
      total_catch_d = total_catch_d + C_gear[d, gg];
    }
    C_sum = C_sum + total_catch_d;
    // Backward compat: total catch across all gear types
    for (s in 1:S) {
      for (g in 1:G) {
        C[s][d,g] = total_catch_d;
      }
    }
  }
}
