// =============================================================================
// Gear-Resolved CPUE Crab Creel Model - Gear-Hours Formulation
// (crab_bss_gear_resolved.stan)
//
// CPUE STRUCTURE: one CPUE process PER GEAR TYPE, not a single pooled process.
//   mu_C is a [G, S] matrix and each gear g carries its own hierarchical mean
//   mu_mu_C[g]; the AR(1) deviations omega_C run over G*S with a Cholesky
//   correlation across the gear-by-section blocks, and the interview catch
//   likelihood is indexed by each interview's gear type (gear_IntC). Gear-type
//   catch therefore carries full posterior uncertainty, rather than being
//   apportioned after the fact from interview proportions (which is what the
//   pooled model, crab_bss_pooled.stan, does). This is the only structural
//   difference from the pooled model; everything below is shared with it.
//
// GEAR-HOURS (boat populations):
//   For boat crabbers using pots, gear is deployed continuously (24 hrs/day)
//   while the trailer remains at the boat ramp. The previous formulation used
//   crabber-hours with day_length (9-16 hrs), creating a unit mismatch that
//   systematically underestimated boat catch by ~2x.
//
//   This version:
//     - lambda_E represents GEAR IN THE WATER (shore: via R_G; boats: directly)
//     - For boats: R_G_boat (gear per group) replaces R_T (trailers per crabber)
//     - CPUE denominator h = gear-hours (boats) or crabber-hours (shore)
//     - L[d] = 24 for boats (gear fishes 24/7), day_length for shore
//     - Trailer obs: T_I ~ Poisson(lambda_E / R_G_boat * eps)
//
//   Shore model is unchanged: lambda_E = crabbers, h = crabber-hours, L = day_length.
//
// SPARSE eps_E_H: overdispersion allocated only for actual observations.
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
  int<lower=1> Gear_A_boat[IntA_trailer];  // number of gear per boat group (replaces T_A_int, A_A_trailer)

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
  real<lower=0> R_G_boat;  // gear per boat group (replaces R_T)

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
    R_G_boat ~ lognormal(log(4), 0.5);  // ~4 gear per group, with range ~2-8
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
  // lambda_E = gear in water; trailers = gear / R_G_boat = groups
  for (i in 1:T_n) {
    T_I[i] ~ poisson(
      lambda_E_S[section_T[i]][day_T[i], G] / R_G_boat * eps_E_H_obs[Gear_n + i]
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

  // --- Gear per boat group: learn R_G_boat from interview data ---
  for (a in 1:IntA_trailer) {
    Gear_A_boat[a] ~ poisson(R_G_boat);
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
  real R_G_boat_out;

  Omega_C = multiply_lower_tri_self_transpose(Lcorr_C);
  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  R_G_boat_out = R_G_boat;
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
