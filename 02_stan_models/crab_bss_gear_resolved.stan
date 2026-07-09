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
//     - Trailer obs: T_I ~ neg_binomial_2(lambda_E / R_G_boat, r_E)
//
//   Shore model is unchanged: lambda_E = crabbers, h = crabber-hours, L = day_length.
//
// Holiday effect B2 separates holiday effort from regular weekends.
// Effort: log(lambda_E) = mu + omega + B1*weekend + B2*holiday
//
// PORTED FROM THE POOLED MODEL (parity work):
//   B1.3  Non-centered AR(1) initial state. omega_*_0 = stationary_SD * raw with
//         raw ~ std_normal(), replacing the centered
//         omega_*_0 ~ normal(0, sqrt(sigma_eps^2 / (1 - phi^2))). Reproduces the
//         same prior exactly while removing the funnel. Inference-preserving.
//
//   B1.5  Effort-count overdispersion marginalized to neg_binomial_2. The former
//         Poisson-Gamma mixture (per-observation latent eps_E_H_obs ~
//         Gamma(r_E, r_E)) is integrated out analytically, giving identical mean
//         and variance with r_E = 1 / sigma_r_E^2 unchanged. Removes n_effort_obs
//         latent parameters and their funnel. n_effort_obs is retained in the
//         data block for R-interface compatibility only; it is unused here.
//
//   B1.8  C_expected_sum: deterministic expected-catch total. The scale-aware
//         convergence gate (03_R_functions/bss_convergence_gate.R) tests whether
//         divergent draws move this total, so it must not carry the Poisson
//         predictive noise that C_sum does.
//
//   log_lik_gear / log_lik_trailer / log_lik_catch: pointwise log-likelihood for
//         PSIS-LOO (loo package) and Pareto-k influence diagnostics. Empty when a
//         stream is absent (trailer for shore, gear for the boat).
// =============================================================================

data {
  int<lower=1> D;
  int<lower=1> G;
  int<lower=1> S;
  int<lower=1> P_n;
  int<lower=1> period[D];
  vector<lower=0,upper=1>[D] w;
  vector<lower=0,upper=1>[D] holiday;

  // 5b: shore day length. L_data is the point value (I/E-derived L_mu for shore,
  // 24 for boats). When estimate_L = 1, L becomes a PARAMETER with a lognormal
  // prior centred on L_data and log-scale SD L_prior_sigma, so the I/E prediction
  // uncertainty propagates into the posterior instead of being asserted as known.
  // When estimate_L = 0 (boats), L is fixed at L_data.
  vector<lower=0>[D] L_data;
  int<lower=0, upper=1> estimate_L;
  vector<lower=0>[D] L_prior_sigma;

  // 5b: ingress/egress observations, a SECOND observation stream on effort.
  // IE_crabber_hours[i] measures total crabber-hours on day_IE[i]; the model
  // predicts lambda_E * L, so these observations jointly identify L and effort.
  // IE_n = 0 for boat fits (no shore I/E), which is why sigma_IE needs an
  // unconditional prior (see B1.6 note in the model block).
  int<lower=0> IE_n;
  int<lower=1> day_IE[IE_n];
  int<lower=1> section_IE[IE_n];
  vector<lower=0>[IE_n] IE_crabber_hours;

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
  real value_normal_sigma_B1_C;
  // B1.9 parity: switch the weekend CPUE effect on (1) or off (0). When 0, B1_C
  // is still sampled from its prior but drops out of the likelihood, so
  // log_lik_catch is that of the reduced model and elpd_loo is directly
  // comparable across a use_B1_C = 1 vs 0 pair of runs.
  int<lower=0, upper=1> use_B1_C;
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
  real B1_C;   // B1.9 parity: weekend/holiday effect on CPUE (pooled model L244)

  // 5b: non-centered lognormal deviation for L. Size 0 when estimate_L = 0, so
  // the boat fits carry no extra parameters at all.
  vector[D * estimate_L] L_raw;

  // 5b: lognormal SD of the I/E crabber-hour observations.
  real<lower=0> sigma_IE;
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
  //       R-interface compatibility; it is no longer referenced here.

  real<lower=0> R_G;
  real<lower=0> R_G_boat;  // gear per boat group (replaces R_T)

  real<lower=0> sigma_eps_C;
  cholesky_factor_corr[G*S] Lcorr_C;
  real<lower=0,upper=1> phi_C_scaled;
  real<lower=0> sigma_r_C;
  matrix[P_n-1, G*S] eps_C;
  matrix[G,S] omega_C_0_raw;   // B1.3: non-centered AR(1) initial state (raw)
  real mu_mu_C[G];
  real<lower=0> sigma_mu_C;
  matrix[G,S] eps_mu_C;
}

transformed parameters {
  matrix[G,S] mu_E;
  real<lower=-1,upper=1> phi_E;
  matrix[P_n, G*S] omega_E;
  matrix[G,S] omega_E_0;        // B1.3: scaled from omega_E_0_raw below
  matrix<lower=0>[D,G] lambda_E_S[S];
  real<lower=0> r_E;

  matrix[G,S] mu_C;
  real<lower=-1,upper=1> phi_C;
  real<lower=0> r_C;
  vector<lower=0>[D] L;   // 5b: parameter when estimate_L = 1, else fixed = L_data
  matrix[P_n, G*S] omega_C;
  matrix[G,S] omega_C_0;        // B1.3: scaled from omega_C_0_raw below
  matrix<lower=0>[D,G] lambda_C_S[S];

  r_E = 1 / square(sigma_r_E);
  r_C = 1 / square(sigma_r_C);
  phi_E = (phi_E_scaled * 2) - 1;
  phi_C = (phi_C_scaled * 2) - 1;

  // 5b: L ~ lognormal(log(L_data), L_prior_sigma), non-centered. Identical in
  // distribution to the centered form, without the funnel.
  if (estimate_L == 1) {
    for (d in 1:D) {
      L[d] = L_data[d] * exp(L_prior_sigma[d] * L_raw[d]);
    }
  } else {
    L = L_data;
  }

  // B1.3: non-centered AR(1) initial state. omega_*_0 = stationary SD x raw,
  //       reproducing normal(0, sqrt(sigma_eps^2 / (1 - phi^2))) exactly while
  //       removing the centered funnel. Inference-preserving.
  omega_E_0 = sqrt(square(sigma_eps_E) / (1 - square(phi_E))) * omega_E_0_raw;
  omega_C_0 = sqrt(square(sigma_eps_C) / (1 - square(phi_C))) * omega_C_0_raw;

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
          to_matrix(omega_C[period[d],], G, S)[g,s]
          + use_B1_C * B1_C * w[d]) * O[d,s,g];
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
  // Proper prior regardless of use_B1_C. When use_B1_C = 0 this parameter does
  // not enter the likelihood, so it simply samples its prior.
  B1_C ~ normal(0, value_normal_sigma_B1_C);

  to_vector(eps_E) ~ std_normal();
  to_vector(eps_C) ~ std_normal();

  // 5b: prior on the non-centered L deviation. Vacuous when estimate_L = 0
  // (L_raw has length 0).
  if (estimate_L == 1) {
    L_raw ~ std_normal();
  }

  // 5b / B1.6: sigma_IE gets a PROPER prior UNCONDITIONALLY. In the pooled model
  // this prior originally sat inside `if (IE_n > 0)`, so a fit with no I/E data
  // (the boat) left sigma_IE with neither prior nor likelihood: an improper flat
  // direction that drifted to ~1e307 and became the boat's dominant divergence
  // source. Do not move this inside the guard below.
  sigma_IE ~ exponential(5);

  if (IE_n > 0) {
    for (i in 1:IE_n) {
      IE_crabber_hours[i] ~ lognormal(
        log(lambda_E_S[section_IE[i]][day_IE[i], 1] * L[day_IE[i]]),
        sigma_IE
      );
    }
  }

  R_G ~ lognormal(log(1.3), 0.3);
  if (T_n > 0 || IntA_trailer > 0) {
    R_G_boat ~ lognormal(log(4), 0.5);  // ~4 gear per group, with range ~2-8
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

  // B1.5: effort-count overdispersion marginalized to neg_binomial_2. The prior
  //       form was Gear_I ~ Poisson(lambda * eps_E_H_obs * R_G) with
  //       eps_E_H_obs ~ Gamma(r_E, r_E); integrating out eps_E_H_obs gives
  //       neg_binomial_2(lambda * R_G, r_E) exactly (identical mean and variance).
  //       r_E = 1 / sigma_r_E^2 is unchanged, so the overdispersion is the same;
  //       only the per-observation latent eps parameters (and their funnel) are
  //       removed. The trailer keeps the gear-resolved / R_G_boat structure.

  // --- Gear counts ---
  for (i in 1:Gear_n) {
    Gear_I[i] ~ neg_binomial_2(
      lambda_E_S[section_Gear[i]][day_Gear[i], 1] * R_G, r_E
    );
  }

  // --- Trailer counts: lambda_E = gear in water; trailers = gear / R_G_boat = groups ---
  for (i in 1:T_n) {
    T_I[i] ~ neg_binomial_2(
      lambda_E_S[section_T[i]][day_T[i], G] / R_G_boat, r_E
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
  real<lower=0> C_expected_sum;   // B1.8 gate: expected (deterministic) catch sum
  real<lower=0> E_sum;
  real R_G_out;
  real R_G_boat_out;
  real sigma_IE_out;      // 5b: exposed for diagnostics (pooled parity)
  vector<lower=0>[D] L_out;   // 5b: realized day length per day

  // Pointwise log-likelihood for PSIS-LOO (loo package), one entry per obs in
  // each stream, mirroring the model-block likelihood terms exactly. Empty when
  // a stream is absent (log_lik_trailer for shore, log_lik_gear for the boat).
  vector[Gear_n] log_lik_gear;
  vector[T_n] log_lik_trailer;
  vector[IntC] log_lik_catch;

  Omega_C = multiply_lower_tri_self_transpose(Lcorr_C);
  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  R_G_boat_out = R_G_boat;
  sigma_IE_out = sigma_IE;
  L_out = L;

  // Mirror the model-block likelihood terms exactly (Gear_I, T_I, c).
  for (i in 1:Gear_n) {
    log_lik_gear[i] = neg_binomial_2_lpmf(
      Gear_I[i] | lambda_E_S[section_Gear[i]][day_Gear[i], 1] * R_G, r_E
    );
  }
  for (i in 1:T_n) {
    log_lik_trailer[i] = neg_binomial_2_lpmf(
      T_I[i] | lambda_E_S[section_T[i]][day_T[i], G] / R_G_boat, r_E
    );
  }
  for (a in 1:IntC) {
    log_lik_catch[a] = neg_binomial_2_lpmf(
      c[a] | lambda_C_S[section_IntC[a]][day_IntC[a], gear_IntC[a]] * h[a], r_C
    );
  }

  C_sum = 0;
  C_expected_sum = 0;
  E_sum = 0;

  for (g in 1:G) {
    for (d in 1:D) {
      for (s in 1:S) {
        lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * L[d] * lambda_C_S[s][d,g];
        C_expected_sum = C_expected_sum + lambda_Ctot_S[s][d,g];
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
