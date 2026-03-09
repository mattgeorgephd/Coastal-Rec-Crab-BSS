// =============================================================================
// BSS Crab Creel Model v1.1
// Adapted from WDFW FW CreelEstimates BSS model (2024-07-24)
// 
// v1.1 changes:
//   - Added if(T_n > 0) guard around R_T prior to prevent divergences
//     when no trailer data is present (shore-only fits)
// 
// Designed for recreational Dungeness crab fishery monitoring.
//
// KEY CHANGES FROM FW MODEL:
//   1. GEAR COUNTS (docks): counts crab gear in water → R_G (gear/crabber)
//   2. TRAILER COUNTS (boat launch): counts trailers → R_T (trailer/group)
//   3. DIRECT CRABBER COUNTS (jetty): counts people crabbing → no expansion
//   4. Census blocks REMOVED (discrete sites = full spatial coverage)
//   5. Vehicle counts and bias parameter REMOVED
//   6. Poisson_rng overflow protection in generated quantities
//
// CRABBER TYPES (G):
//   g=1: Shore crabbers (dock, jetty) — informed by gear + crabber counts
//   g=2: Boat crabbers — informed by trailer counts
//   Can run G=1 (dock-only) by setting T_n=0, Crab_n=0
// =============================================================================

data {
  // --- Core dimensions ---
  int<lower=1> D;              // days in season
  int<lower=1> G;              // crabber types
  int<lower=1> S;              // sections
  int<lower=1> H;              // max count sequences per day
  int<lower=1> P_n;            // time periods
  int<lower=1> period[D];      // period index per day
  vector<lower=0,upper=1>[D] w; // weekend indicator (0=weekday, 1=weekend/holiday)
  vector<lower=0>[D] L;        // day length (hours)
  real<lower=0> O[D,S,G];      // open/closed status

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

  // --- DIRECT CRABBER counts (jetty) ---
  int<lower=0> Crab_n;
  int<lower=1> day_Crab[Crab_n];
  int<lower=1> section_Crab[Crab_n];
  int<lower=1> countnum_Crab[Crab_n];
  int<lower=0> Crab_I[Crab_n];
  real<lower=0,upper=1> p_I_crab;  // proportion of jetty visible (1.0 if full)

  // --- Interview CPUE data ---
  int<lower=0> IntC;
  int<lower=1> day_IntC[IntC];
  int<lower=1> gear_IntC[IntC];
  int<lower=1> section_IntC[IntC];
  int<lower=0> c[IntC];
  vector<lower=0>[IntC] h;

  // --- Interview expansion: gear per crabber ---
  int<lower=0> IntA_gear;
  int<lower=0> Gear_A[IntA_gear];
  int<lower=1> A_A_gear[IntA_gear];

  // --- Interview expansion: trailers per group ---
  int<lower=0> IntA_trailer;
  int<lower=0> T_A_int[IntA_trailer];
  int<lower=1> A_A_trailer[IntA_trailer];

  // --- Hyperparameters ---
  real value_cauchyDF_sigma_eps_E;
  real value_cauchyDF_sigma_eps_C;
  real value_cauchyDF_sigma_r_E;
  real value_cauchyDF_sigma_r_C;
  real value_betashape_phi_E_scaled;
  real value_betashape_phi_C_scaled;
  real value_normal_sigma_B1;
  real value_normal_mu_mu_C;
  real value_normal_sigma_mu_C;
  real value_normal_mu_mu_E;
  real value_normal_sigma_mu_E;
  real value_cauchyDF_sigma_mu_C;
  real value_cauchyDF_sigma_mu_E;
}

parameters {
  // --- Effort process ---
  real B1;
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
  real<lower=0> R_G;                   // gear per crabber (lognormal prior, typically 1-3)
  real<lower=0,upper=1> R_T;           // trailers per boat crabber group (beta prior)

  // --- CPUE process ---
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
  matrix<lower=0>[D,G] lambda_E_S_I[S,H];
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

  // AR(1) process
  omega_E[1,] = to_row_vector(omega_E_0);
  omega_C[1,] = to_row_vector(omega_C_0);
  for (p in 2:P_n) {
    omega_E[p,] = to_row_vector(phi_E * to_vector(omega_E[p-1,]) +
      diag_pre_multiply(rep_vector(sigma_eps_E, G*S), Lcorr_E) * to_vector(eps_E[p-1,]));
    omega_C[p,] = to_row_vector(phi_C * to_vector(omega_C[p-1,]) +
      diag_pre_multiply(rep_vector(sigma_eps_C, G*S), Lcorr_C) * to_vector(eps_C[p-1,]));
  }

  // Daily effort and CPUE rates
  for (g in 1:G) {
    for (s in 1:S) {
      mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E;
      mu_C[g,s] = mu_mu_C[g] + eps_mu_C[g,s] * sigma_mu_C;
    }
    for (d in 1:D) {
      for (s in 1:S) {
        lambda_E_S[s][d,g] = exp(mu_E[g,s] +
          to_matrix(omega_E[period[d],], G, S)[g,s] + B1 * w[d]) * O[d,s,g];
        lambda_C_S[s][d,g] = exp(mu_C[g,s] +
          to_matrix(omega_C[period[d],], G, S)[g,s]) * O[d,s,g];
        for (i in 1:H) {
          lambda_E_S_I[s,i][d,g] = lambda_E_S[s][d,g] * eps_E_H[s,i][d,g];
        }
      }
    }
  }
}

model {
  // === HYPERPRIORS ===
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

  // === PRIORS ===
  to_vector(eps_E) ~ std_normal();
  to_vector(eps_C) ~ std_normal();

  // Gear per crabber: prior centered ~1.3, data from 3000+ interviews will dominate
  R_G ~ lognormal(log(1.3), 0.3);
  // Trailers per group: only evaluate prior when trailer data exists
  // (avoids divergences from U-shaped beta prior with no data)
  if (T_n > 0) {
    R_T ~ beta(0.5, 0.5);
  }

  for (g in 1:G) {
    mu_mu_E[g] ~ normal(value_normal_mu_mu_E, value_normal_sigma_mu_E);
    mu_mu_C[g] ~ normal(value_normal_mu_mu_C, value_normal_sigma_mu_C);
    for (d in 1:D) {
      for (s in 1:S) {
        for (i in 1:H) {
          eps_E_H[s,i][d,g] ~ gamma(r_E, r_E);
        }
      }
    }
    for (s in 1:S) {
      omega_E_0[g,s] ~ normal(0, sqrt(square(sigma_eps_E) / (1 - square(phi_E))));
      omega_C_0[g,s] ~ normal(0, sqrt(square(sigma_eps_C) / (1 - square(phi_C))));
      eps_mu_E[g,s] ~ std_normal();
      eps_mu_C[g,s] ~ std_normal();
    }
  }

  // === LIKELIHOODS ===

  // --- Gear counts (docks): gear_observed ~ Poisson(effort_shore * R_G) ---
  for (i in 1:Gear_n) {
    Gear_I[i] ~ poisson(
      lambda_E_S_I[section_Gear[i], countnum_Gear[i]][day_Gear[i], 1] * R_G
    );
  }

  // --- Trailer counts (boat launch): trailer_observed ~ Poisson(effort_boat * R_T) ---
  for (i in 1:T_n) {
    T_I[i] ~ poisson(
      lambda_E_S_I[section_T[i], countnum_T[i]][day_T[i], G] * R_T
    );
  }

  // --- Direct crabber counts (jetty): crabber_observed ~ Poisson(effort_shore * p_I) ---
  for (i in 1:Crab_n) {
    Crab_I[i] ~ poisson(
      lambda_E_S_I[section_Crab[i], countnum_Crab[i]][day_Crab[i], 1] * p_I_crab
    );
  }

  // --- Interview CPUE: crab_caught ~ NegBin(CPUE * hours, r_C) ---
  for (a in 1:IntC) {
    c[a] ~ neg_binomial_2(
      lambda_C_S[section_IntC[a]][day_IntC[a], gear_IntC[a]] * h[a], r_C
    );
  }

  // --- Interview expansion: gear per crabber ---
  for (a in 1:IntA_gear) {
    Gear_A[a] ~ poisson(A_A_gear[a] * R_G);
  }

  // --- Interview expansion: trailers per group ---
  for (a in 1:IntA_trailer) {
    T_A_int[a] ~ bernoulli(R_T);
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
  real R_G_out;  // expose R_G for monitoring

  Omega_C = multiply_lower_tri_self_transpose(Lcorr_C);
  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  C_sum = 0;
  E_sum = 0;

  for (g in 1:G) {
    for (d in 1:D) {
      for (s in 1:S) {
        lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * L[d] * lambda_C_S[s][d,g];
        // Overflow protection for poisson_rng
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
