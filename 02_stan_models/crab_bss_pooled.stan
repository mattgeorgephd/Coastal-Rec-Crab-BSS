// -----------------------------------------------------------------------------
// Part of Coastal-Rec-Crab-BSS: recreational Dungeness crab creel estimation
// for Grays Harbor / Westport (WDFW).
// Copyright (C) 2024-2026 Washington Department of Fish and Wildlife.
//
// Adapted from CreelEstimates, the WDFW freshwater creel estimation framework:
//   https://github.com/dfw-wa/CreelEstimates   (licensed GPL-3.0).
// Substantial portions of the methodology, structure, and R/Stan code originate
// in CreelEstimates and remain (C) their authors under GPL-3.0; changes for
// recreational crab are by WDFW.
//
// This program is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License, version 3, as published by the Free
// Software Foundation. It is distributed WITHOUT ANY WARRANTY; without even the
// implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details. You should have received a copy of
// the GNU General Public License along with this program (see the LICENSE file);
// if not, see <https://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------
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
//   - L_effective / turnover tau as a parameter with lognormal prior
//   - Data-driven R_G prior from interview data
//   - R_G_boat (gear per boat group) with a lognormal prior (POOL-1)
//   - effort_scale_gear: E carries the unit of h (P1 / POOL-3)
//   - Sparse overdispersion (observation-indexed)
//   - I/E direct effort integration
//   - Dual reporting: expected catch + predictive draws
//
// v7.6 (POOL-1 + POOL-3, 2026-07-10): the BOAT is moved onto the gear-deployment
//   scale, matching crab_bss_gear_resolved.stan. R_T is replaced by R_G_boat with
//   trailer counts T_I ~ neg_binomial_2(lambda_E / R_G_boat, r_E) and interviews
//   Gear_A_boat ~ poisson(R_G_boat), so lambda_E is GEAR (not groups). effort_scale_gear
//   / E_scale generalize E = lambda_E * E_scale * L; the driver sends h = number_of_gear
//   and L = tau_boat for the boat. This removes the POOL-1 (R_T pinned at 1) and
//   POOL-3 (flat L = 24 gear-hours) defects; the boat number moves and must be
//   confirmed by a run. Shore is unchanged (crabber-hours, E_scale = 1).
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
//
// v6.8 (B1.6): sigma_IE is given a proper prior unconditionally. Previously the
//   prior sat inside if (IE_n > 0), so for a fit with no I/E data (the boat)
//   sigma_IE had no prior and no likelihood, an improper flat direction that
//   drifted to ~1e307 and was the boat's dominant divergence source. sigma_IE
//   is decoupled from effort and catch, so this is inference-preserving for the
//   reported quantities; it only makes the posterior proper.
//
// v6.9 (B1.7) ATTEMPTED AND REVERTED (2026-06-21). Collapsing the single-cell
//   hierarchical scale (removing sigma_mu_E / eps_mu_E and the C counterparts so
//   mu_E = mu_mu_E) cleared the boat funnel in offline diagnosis, but in the
//   production run it made the shore all-gear fit hang. That fit is a daily AR
//   over 289 days with ~50% of days unobserved; removing the level term that was
//   decoupled from the AR forced the overall level to reconcile directly against
//   the high-dimensional AR process, a long thin ridge the sampler traverses with
//   maximal-length trajectories (max_treedepth saturation, ~2^14 leapfrog steps
//   per iteration). shore ring-net (76 days, denser) completed in ~8 min; shore
//   all-gear did not finish in 24 h, whereas v6.8 completed it with treedepth 0%.
//   The hierarchy is restored (this file is back to the v6.8 structure). The boat
//   stays on PE; the durable fix is a more informative boat effort series, not
//   parameter surgery. See documentation Section 14.
// =============================================================================

data {
  int<lower=1> D;                        // Number of days in sub-season
  int<lower=1> G;                        // Number of gear groups (1 in pooled)
  int<lower=1> S;                        // Number of sections (1)
  int<lower=1> P_n;                      // Number of AR periods
  int<lower=1> period[D];               // Day-to-period mapping
  vector<lower=0,upper=1>[D] w;          // Weekend indicator
  vector<lower=0,upper=1>[D] holiday;    // Holiday indicator
  // improvement 4 (2026-08-25): OTHER-FISHERY OPENER EFFORT COVARIATES.
  // Generalizes the single razor-dig indicator into a K_open-column design matrix on
  // the EFFORT process, so any opener (razor dig, MA2 halibut / salmon / bottomfish)
  // can enter either population's effort model. A razor-only run with K_open = 1 is
  // mathematically identical to the retired B3 * razor[d] term.
  //
  // ON "IDENTICAL", PRECISELY. At K_open = 0 the term drops out and no B_open parameter
  // exists, so the POSTERIOR over every reported quantity is unchanged. It is NOT
  // bit-identical to a pre-2026-08-25 pooled run: the retired `real B3` was declared
  // unconditionally with a proper prior and, when razor was all zeros, entered no
  // likelihood -- a genuine sampled dimension that is now gone. Removing it changes the
  // HMC trajectories and the RNG stream, so a fixed-seed rerun will not diff to zero
  // against Run 6. Judge a baseline-reproduction run on medians and intervals within
  // Monte Carlo error, not by comparing CSVs byte for byte. (The gear-resolved model
  // never carried B3 and only gains zero-size containers, so IT is bit-identical.)
  // The driver owns column selection (manual, or automatic on the spillover diagnostic's
  // multiplicity-adjusted p-values) and passes the column labels back out of band.
  // CPUE deliberately gets no opener covariate: every opener-vs-catch-rate test in the
  // 2024-25 diagnostic was null.
  int<lower=0> K_open;
  // Passed FLAT (column-major, length D * K_open) rather than as a matrix, so the data
  // block never contains a zero-extent container: rstan marshals a zero-length vector
  // reliably (the same pattern the empty effort streams already use), whereas a D x 0
  // matrix is an unnecessary edge case. Rebuilt as a matrix in transformed data.
  vector[D * K_open] X_open_flat;
  real<lower=0> O[D,S,G];               // Open/closed

  // POOL-4: single-cell mu-hierarchy collapse lever. 0 (default) keeps the v6.8
  // decoupled level term (mu = mu_mu + eps_mu * sigma_mu); 1 collapses it to
  // mu = mu_mu, i.e. the B1.7 experiment that cleared the boat funnel offline but
  // hung the shore all-gear fit in production (289-day daily AR; see header v6.9
  // note and documentation Section 14). Exposed as data so the collapse can be
  // toggled PER FIT from the driver (params$collapse_mu_hier) without editing this
  // file mid-investigation. The durable fix for the funnel is a more informative
  // boat effort series, not this parameter surgery; this is only a safe lever.
  int<lower=0,upper=1> collapse_mu_hier;

  // --- Day length / effort expansion factor L ---
  //   shore: effective day length in HOURS (I/E-derived L_mu).
  //   boat (POOL-3): tau, the gear-deployment TURNOVER (dimensionless), a
  //         PARAMETER (estimate_L = 1) prior-centred on L_data with log-scale SD
  //         L_prior_sigma, replacing the old flat L = 24 h.
  vector<lower=0>[D] L_data;
  int<lower=0,upper=1> estimate_L;
  vector<lower=0>[D] L_prior_sigma;

  // P1 (POOL-3): does the effort expansion count CRABBERS (lambda_E, scale 0) or
  // GEAR (lambda_E * R_G, scale 1)? Shore crabber-hours and the boat (lambda_E is
  // already gear, via R_G_boat) use 0; the shore gear-hours / gear-deployments
  // options use 1. E and the catch expansion are multiplied by E_scale so E always
  // carries the same unit as the CPUE denominator h.
  int<lower=0, upper=1> effort_scale_gear;

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

  // --- Phase 1: OSP boat-count observations (second boat effort stream) ---
  // OSP_I is the DAILY BOAT TOTAL (all private boats) on OSP-sampled days. It
  // observes the same latent lambda_E as the trailer stream, via a scale
  // coefficient kappa_OSP (~ the OSP/trailer overlap ratio). Empty (OSP_n = 0)
  // for shore fits and when use_osp_boat_counts = FALSE.
  int<lower=0> OSP_n;
  int<lower=1> day_OSP[OSP_n];
  int<lower=1> section_OSP[OSP_n];
  int<lower=0> OSP_I[OSP_n];

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

  // --- Gear-per-boat-group interviews (POOL-1) ---
  // Replaces T_A_int / A_A_trailer. The model learns R_G_boat (mean gear deployed
  // per boat group) directly from observed number_of_gear via
  // Gear_A_boat[a] ~ poisson(R_G_boat), mirroring crab_bss_gear_resolved.stan.
  int<lower=0> IntA_trailer;
  int<lower=1> Gear_A_boat[IntA_trailer];

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
  real value_normal_sigma_B2_C;   // item 6a: holiday CPUE effect (B2_C) prior scale
  real value_normal_sigma_B_open; // improvement 4: opener effort-covariate prior scale
  // item 6b: same-day-effort CPUE density-dependence (OFF by default). estimate_cpue_density
  // toggles the gamma_C term; log_E_ref centers the log-effort covariate (per fit, = the
  // effort level prior mu_E) so gamma_C does not soak up the CPUE intercept.
  int<lower=0,upper=1> estimate_cpue_density;
  real log_E_ref;
  real value_normal_mu_mu_C;
  real value_normal_sigma_mu_C;
  real value_normal_mu_mu_E;
  real value_normal_sigma_mu_E;
  real value_cauchyDF_sigma_mu_C;
  real value_cauchyDF_sigma_mu_E;

  // --- Data-driven priors ---
  real<lower=0> R_G_prior_mu;
  real<lower=0> R_G_prior_sigma;
  // POOL-1: R_T_alpha / R_T_beta removed. R_G_boat now carries a fixed lognormal
  // prior (see model block), matching crab_bss_gear_resolved.stan.

  // Phase 1: kappa_OSP (OSP-to-trailer scale) lognormal prior, centered on the
  // OSP/trailer overlap ratio (~3) measured by diagnose_osp_trailer_overlap().
  real<lower=0> osp_scale_prior_mu;
  real<lower=0> osp_scale_prior_sigma;

  // Phase 2: crabbing fraction f (boat only). apply_crab_fraction = 0 pins f = 1
  // (behavior-neutral: shore, or use_crab_fraction = FALSE). When estimate = 1, f is a
  // Beta(alpha0, beta0) parameter (prior centered on the set value) optionally updated by
  // the I/E crab-vs-total classification (n_crab of n_total); when estimate = 0, f is
  // pinned to crab_fraction_value (a hard set value). f enters ONLY generated quantities,
  // scaling boat effort and catch, so it never perturbs the effort/CPUE sampling.
  int<lower=0,upper=1> apply_crab_fraction;
  int<lower=0,upper=1> crab_fraction_estimate;
  // Phase 3: f is per STRATUM (n_f_strata; = 1 reproduces the Phase 2 scalar).
  // f_stratum[d] maps each day to its stratum (month / day_type / none).
  int<lower=1> n_f_strata;
  int<lower=1,upper=n_f_strata> f_stratum[D];
  vector<lower=0,upper=1>[n_f_strata] crab_fraction_value;
  vector<lower=0>[n_f_strata] crab_fraction_alpha0;
  vector<lower=0>[n_f_strata] crab_fraction_beta0;
  int<lower=0> crab_fraction_n_total[n_f_strata];
  int<lower=0> crab_fraction_n_crab[n_f_strata];
  // improvement 8 (2026-08-25): OSP CRAB-ONLY COUNTS AS A HARD LOWER BOUND ON f.
  // OSP can report how many boats were labelled crabbing ONLY; it cannot see combo
  // trips (a boat crabbing AND fishing something else is labelled by the other
  // fishery), so the OSP crab-only share is a LOWER BOUND on f, never an estimate of
  // it. The bound is imposed by construction:
  //     f_lower[k] ~ Beta(1,1),  osp_f_n_crab[k] ~ Binomial(osp_f_n_total[k], f_lower[k])
  //     f[k]       = f_lower[k] + (1 - f_lower[k]) * theta[k]
  // theta[k] is the share of the NOT-crab-labelled boats that were also crabbing. With
  // osp_crab_lower = 0, or in a stratum with no OSP classification, f_lower[k] is pinned
  // to 0 and f[k] = theta[k], which reproduces the Phase 2/3 scalar/stratum f EXACTLY
  // (the R side then sets crab_fraction_alpha0/beta0 to the ordinary f prior rather than
  // the combo prior, so the degradation is exact, not approximate).
  int<lower=0,upper=1> osp_crab_lower;
  // Per-stratum totals: used ONLY as the has-data flag and for reporting. The
  // likelihood below is per DAY, not on these sums -- see the next paragraph.
  int<lower=0> osp_f_n_total[n_f_strata];
  int<lower=0> osp_f_n_crab[n_f_strata];
  // WHY PER-DAY AND WHY BETA-BINOMIAL. A Binomial on the stratum SUM would treat
  // every boat-day as an independent Bernoulli trial: summing ~150 operating days at
  // ~50 boats gives n ~ 7,500 and an f_lower posterior SD near 0.005, which is not
  // remotely supported by ~150 correlated daily observations. The crab share moves
  // day to day with weather, tide and which other fishery is open, so the honest
  // model is one observation PER DAY with an estimated overdispersion:
  //     osp_f_crab[i] ~ beta_binomial(osp_f_total[i], f_lower*kappa, (1-f_lower)*kappa)
  // kappa is identified by the spread of the daily shares (large kappa -> binomial,
  // small kappa -> strongly overdispersed), so the bound's precision is measured
  // rather than assumed. This matters because f multiplies the boat total linearly
  // and the bound is HARD: an over-precise f_lower would force f up with unearned
  // confidence.
  int<lower=0> OSPF_n;
  int<lower=1,upper=n_f_strata> osp_f_stratum[OSPF_n];
  int<lower=0> osp_f_total[OSPF_n];
  int<lower=0> osp_f_crab[OSPF_n];
  real<lower=0> osp_f_kappa_prior_mu;
  // Phase 3: OSP-informs-tau toggle. 0 (default) keeps the free kappa_OSP scale (Phase 1);
  // 1 makes the OSP mean use L (= tau_boat) as the turnover, so the dense OSP series
  // identifies the boat turnover. See the dev note for the kappa_OSP (~2.7) vs
  // tau_boat prior (~1.2) tension; validate by run before using.
  int<lower=0,upper=1> osp_scale_is_tau;
}

transformed data {
  matrix[D, K_open] X_open;   // improvement 4: rebuilt from the flat opener design vector
  for (k in 1:K_open)
    for (d in 1:D)
      X_open[d, k] = X_open_flat[(k - 1) * D + d];
}

parameters {
  real B1;
  real B2;
  real B1_C;
  real B2_C;      // item 6a: holiday CPUE effect (symmetric to B2 on effort)
  vector[K_open] B_open;   // improvement 4: opener effort covariates (length 0 when unused)
  real gamma_C;   // item 6b: same-day-effort CPUE density effect (inert when estimate_cpue_density = 0)

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
  real<lower=0> R_G_boat;   // POOL-1: gear per boat group (replaces R_T)
  real<lower=0> kappa_OSP;    // Phase 1: OSP-to-trailer scale (within-day boat turnover)
  real<lower=0> sigma_r_OSP;  // Phase 1: OSP observation overdispersion (own r)
  // Phase 3: per-stratum crabbing fraction. Length 0 when pinned/off (f_crab is then data).
  // improvement 8: this is now theta, the share of the NOT-crab-labelled boats that were
  // also crabbing, not f itself. With no OSP lower bound they coincide (f = theta).
  vector<lower=0,upper=1>[n_f_strata * crab_fraction_estimate] f_theta;
  // improvement 8: the OSP-observed crab-only share, the hard lower bound on f.
  vector<lower=0,upper=1>[n_f_strata * osp_crab_lower] f_lower_param;
  // Beta-binomial concentration for the DAILY OSP crab-only shares. Proper prior
  // unconditionally, like kappa_OSP / R_G_boat / sigma_IE, so it is decoupled but
  // proper when the stream is absent.
  real<lower=0> osp_f_kappa;

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
  real<lower=0> r_OSP;   // Phase 1: OSP observation overdispersion

  matrix[G,S] mu_C;
  real<lower=-1,upper=1> phi_C;
  real<lower=0> r_C;
  matrix[P_n, G*S] omega_C;
  matrix[G,S] omega_C_0;
  matrix<lower=0>[D,G] lambda_C_S[S];

  vector<lower=0>[D] L;
  real<lower=0> E_scale;   // P1/POOL-3: 1 (crabbers) or R_G (gear); see effort_scale_gear
  vector<lower=0,upper=1>[n_f_strata] f_crab;   // Phase 3: per-stratum crabbing fraction
  vector<lower=0,upper=1>[n_f_strata] f_lower;  // improvement 8: OSP crab-only lower bound

  // P1/POOL-3: convert lambda_E into the unit of the CPUE denominator h. For the
  // boat, lambda_E is already gear (via R_G_boat) so effort_scale_gear = 0 and
  // E_scale = 1; for shore crabber-hours it is also 0. E = lambda_E * E_scale * L.
  E_scale = (effort_scale_gear == 1) ? R_G : 1.0;

  // Phase 3 + improvement 8: resolve per-stratum f_crab.
  //   apply = 0                  -> f = 1 (shore, or the feature off)
  //   estimate = 1               -> f = f_lower + (1 - f_lower) * theta, so f can never
  //                                 fall below the crab-only share OSP directly observed
  //   estimate = 0               -> the pinned set values (sensitivity lever)
  // f_lower is 0 whenever the OSP bound is off or that stratum carries no OSP
  // classification, in which case f = theta exactly (Phase 2/3 behaviour).
  for (k in 1:n_f_strata) {
    if (osp_crab_lower == 1 && osp_f_n_total[k] > 0) f_lower[k] = f_lower_param[k];
    else                                             f_lower[k] = 0.0;

    if (apply_crab_fraction == 0)         f_crab[k] = 1.0;
    else if (crab_fraction_estimate == 1) f_crab[k] = f_lower[k] + (1 - f_lower[k]) * f_theta[k];
    else                                  f_crab[k] = crab_fraction_value[k];
  }

  // --- Compute L ---
  if (estimate_L == 1) {
    for (d in 1:D)
      L[d] = L_data[d] * exp(L_prior_sigma[d] * L_raw[d]);
  } else {
    L = L_data;
  }

  r_E = 1 / square(sigma_r_E);
  r_OSP = 1 / square(sigma_r_OSP);   // Phase 1
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
      // POOL-4 lever: collapse_mu_hier == 1 removes the decoupled single-cell
      // level (mu = mu_mu), the B1.7 collapse. Default 0 reproduces the v6.8
      // hierarchy EXACTLY, so the default posterior is unchanged. When collapsed,
      // eps_mu_* and sigma_mu_* keep their priors (below) but enter no likelihood,
      // so they are proper and decoupled (like sigma_IE at IE_n = 0), not a funnel.
      if (collapse_mu_hier == 1) {
        mu_E[g,s] = mu_mu_E[g];
        mu_C[g,s] = mu_mu_C[g];
      } else {
        mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E;
        mu_C[g,s] = mu_mu_C[g] + eps_mu_C[g,s] * sigma_mu_C;
      }
    }
    for (d in 1:D) {
      for (s in 1:S) {
        // Effort: AR deviation (at period resolution) + weekend + holiday + any opener
        //   covariates (improvement 4). With K_open = 0 the opener term is absent and this
        //   is byte-identical to the pre-2026-08-25 effort process.
        lambda_E_S[s][d,g] = exp(mu_E[g,s] +
          to_matrix(omega_E[period[d],], G, S)[g,s] + B1 * w[d] + B2 * holiday[d] +
          (K_open > 0 ? dot_product(X_open[d], B_open) : 0.0)) * O[d,s,g];
        // CPUE: AR deviation + weekend + holiday (B2_C, item 6a) + optional same-day-effort
        //   density-dependence (gamma_C, item 6b; only when estimate_cpue_density = 1). The
        //   density covariate is the same day's effort intensity lambda_E_S (assigned just
        //   above), centered at log_E_ref. In the pooled model O = 1, so log(lambda_E_S) is finite.
        lambda_C_S[s][d,g] = exp(mu_C[g,s] +
          to_matrix(omega_C[period[d],], G, S)[g,s] + B1_C * w[d] + B2_C * holiday[d] +
          (estimate_cpue_density == 1 ? gamma_C * (log(lambda_E_S[s][d,g]) - log_E_ref) : 0.0)) * O[d,s,g];
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
  sigma_r_OSP ~ cauchy(0, value_cauchyDF_sigma_r_E);   // Phase 1: OSP overdispersion
  sigma_r_C ~ cauchy(0, value_cauchyDF_sigma_r_C);
  sigma_mu_E ~ cauchy(0, value_cauchyDF_sigma_mu_E);
  sigma_mu_C ~ cauchy(0, value_cauchyDF_sigma_mu_C);
  B1 ~ normal(0, value_normal_sigma_B1);
  B2 ~ normal(0, value_normal_sigma_B2);
  B1_C ~ normal(0, value_normal_sigma_B1_C);
  B2_C ~ normal(0, value_normal_sigma_B2_C);   // item 6a: holiday CPUE effect
  B_open ~ normal(0, value_normal_sigma_B_open);  // improvement 4: opener effort covariates
  gamma_C ~ normal(0, 1);   // item 6b: proper prior; enters the likelihood only when estimate_cpue_density = 1

  to_vector(eps_E) ~ std_normal();
  to_vector(eps_C) ~ std_normal();

  // R_G (shore gear-per-crabber) keeps its data-driven lognormal prior
  // unconditionally. In a boat fit R_G enters no likelihood term (Gear_n = 0,
  // IntA_gear = 0, and effort_scale_gear = 0 so E_scale = 1), so it simply samples
  // its prior; the prior is proper and scale-matched, so it is harmless.
  R_G ~ lognormal(log(R_G_prior_mu), R_G_prior_sigma);

  // F1 / POOL-1: R_G_boat (gear per boat group) gets a PROPER prior UNCONDITIONALLY;
  // it replaces R_T. Because R_G_boat is declared real<lower=0>, an unguarded flat
  // prior would be IMPROPER (Stan samples exp(z) with a +z Jacobian, so the log
  // density rises without bound as z -> inf; that is the run-away bug F1 fixed for
  // gear-resolved). A lognormal is proper. In a shore fit R_G_boat enters no
  // likelihood (T_n = 0, IntA_trailer = 0) and simply samples this prior. Do NOT
  // move it inside a guard.
  R_G_boat ~ lognormal(log(4), 0.5);   // ~4 gear per group, range ~2-8

  // Phase 1: kappa_OSP proper prior UNCONDITIONALLY (proper even when OSP_n = 0,
  // like R_G_boat / sigma_IE), centered on the OSP/trailer overlap ratio.
  kappa_OSP ~ lognormal(log(osp_scale_prior_mu), osp_scale_prior_sigma);

  // Phase 2: crabbing fraction. Beta prior on f (centered on the set value) + optional
  // Binomial from the I/E crab-vs-total classification. Decoupled from effort/catch
  // sampling (f enters only generated quantities), so no identifiability interaction.
  if (crab_fraction_estimate == 1) {
    for (k in 1:n_f_strata) {
      // improvement 8: the Beta prior now sits on theta. In a stratum with no OSP bound
      // f_lower = 0 so f = theta and this IS the historical prior on f; in an
      // OSP-informed stratum the R side hands over the combo-share prior instead.
      f_theta[k] ~ beta(crab_fraction_alpha0[k], crab_fraction_beta0[k]);
      // WPT/WBL egress classification: boats seen crabbing (combo trips INCLUDED) out of
      // boats classified. This binds on f itself, so it is what pins theta once the OSP
      // bound is in play, and what carries f alone when OSP is not in port.
      if (crab_fraction_n_total[k] > 0)
        crab_fraction_n_crab[k] ~ binomial(crab_fraction_n_total[k], f_crab[k]);
    }
  }
  // improvement 8: the OSP crab-only stream. Binomial on the OBSERVED OSP daily totals,
  // so f stays entirely out of the effort/CPUE likelihoods: the boat remains exactly
  // linear in f and the model CPUE remains invariant (the validated Phase 2/3 property).
  // Proper prior unconditionally, so an unused f_lower_param is decoupled-but-proper,
  // exactly like sigma_IE at IE_n = 0.
  osp_f_kappa ~ lognormal(log(osp_f_kappa_prior_mu), 0.75);
  if (osp_crab_lower == 1) {
    for (k in 1:n_f_strata) f_lower_param[k] ~ beta(1, 1);
    // One observation per OSP day. osp_f_stratum only ever points at strata that
    // cleared the minimum, so f_lower_param there is the live parameter (never the
    // pinned 0 that f_lower carries for data-free strata).
    for (i in 1:OSPF_n)
      osp_f_crab[i] ~ beta_binomial(osp_f_total[i],
                                    f_lower_param[osp_f_stratum[i]] * osp_f_kappa,
                                    (1 - f_lower_param[osp_f_stratum[i]]) * osp_f_kappa);
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

  // POOL-1: trailer counts. lambda_E is gear in the water; trailers = groups =
  //         lambda_E / R_G_boat. (Was lambda_E * R_T, which pinned R_T at ~1 and
  //         left lambda_E in group units mismatched to gear-hours h.)
  for (i in 1:T_n) {
    T_I[i] ~ neg_binomial_2(
      lambda_E_S[section_T[i]][day_T[i], G] / R_G_boat, r_E
    );
  }

  // Phase 1: OSP boat-count stream. OSP_I is the daily total of ALL private boats,
  // = (crab boats present = lambda_E / R_G_boat) * kappa_OSP, where kappa_OSP is the
  // within-day boat turnover (OSP/trailer overlap ratio). Own overdispersion r_OSP.
  for (i in 1:OSP_n) {
    // Phase 3: osp_scale_is_tau = 1 uses L (tau_boat) as the turnover so OSP identifies
    // tau; 0 keeps the free kappa_OSP (Phase 1). See dev note (kappa_OSP ~2.7 vs tau ~1.2).
    OSP_I[i] ~ neg_binomial_2(
      (lambda_E_S[section_OSP[i]][day_OSP[i], G] / R_G_boat)
        * (osp_scale_is_tau == 1 ? L[day_OSP[i]] : kappa_OSP), r_OSP
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

  // POOL-1: learn R_G_boat (mean gear per boat group) from interview data. Was
  //         T_A_int[a] ~ bernoulli(R_T) on a vector of literal ones, which pinned
  //         R_T at 1.00 and made the trailer expansion degenerate.
  for (a in 1:IntA_trailer) {
    Gear_A_boat[a] ~ poisson(R_G_boat);
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
  real R_G_boat_out;   // POOL-1: gear per boat group (replaces R_T reporting)
  real kappa_OSP_out;  // Phase 1: OSP-to-trailer scale
  vector[n_f_strata] f_crab_out;   // Phase 3: per-stratum crabbing fraction
  vector[n_f_strata] f_lower_out;
  real osp_f_kappa_out;        // improvement 8: daily-share overdispersion of the OSP bound  // improvement 8: OSP crab-only lower bound on f
  real sigma_IE_out;
  real B1_C_out;
  real B2_C_out;       // item 6a
  vector[K_open] B_open_out;   // improvement 4 (labels travel out of band, see the driver)
  real gamma_C_out;    // item 6b
  vector[D] L_out;

  // Pointwise log-likelihood for PSIS-LOO (loo package) and Pareto-k influence
  // diagnostics, one entry per observation in each data stream. Enables the
  // project's primary model-selection tool (PSIS-LOO) on the pooled model, e.g.
  // comparing AR resolutions for the boat, and flags influential observations
  // (high Pareto-k) such as sparse-month interviews. Empty when a stream is
  // absent (e.g. log_lik_trailer for shore, log_lik_gear for the boat). Note:
  // these add a [draws x n_obs] matrix per stream to the fit; for the shore
  // all-gear fit (~3000 obs) that is a few hundred MB in memory. If that is a
  // constraint, gate these behind a data flag (available on request).
  vector[Gear_n] log_lik_gear;
  vector[T_n] log_lik_trailer;
  vector[OSP_n] log_lik_osp;   // Phase 1
  vector[IntC] log_lik_catch;

  Omega_C = multiply_lower_tri_self_transpose(Lcorr_C);
  Omega_E = multiply_lower_tri_self_transpose(Lcorr_E);
  R_G_out = R_G;
  R_G_boat_out = R_G_boat;
  kappa_OSP_out = kappa_OSP;   // Phase 1
  f_crab_out = f_crab;         // Phase 2
  f_lower_out = f_lower;
  osp_f_kappa_out = osp_f_kappa;       // improvement 8
  sigma_IE_out = sigma_IE;
  B1_C_out = B1_C;
  B2_C_out = B2_C;
  B_open_out = B_open;
  gamma_C_out = gamma_C;
  L_out = L;

  // Mirror the model-block likelihood terms exactly (lines for Gear_I, T_I, c).
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
  for (i in 1:OSP_n) {   // Phase 1
    log_lik_osp[i] = neg_binomial_2_lpmf(
      OSP_I[i] | (lambda_E_S[section_OSP[i]][day_OSP[i], G] / R_G_boat)
        * (osp_scale_is_tau == 1 ? L[day_OSP[i]] : kappa_OSP), r_OSP
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
        // P1/POOL-3: E_scale converts lambda_E to the unit of h (see effort_scale_gear).
        lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * E_scale * L[d] * lambda_C_S[s][d,g] * f_crab[f_stratum[d]];
        C_expected[s][d,g] = lambda_Ctot_S[s][d,g];
        C_expected_sum = C_expected_sum + C_expected[s][d,g];

        if (lambda_Ctot_S[s][d,g] < 1e9) {
          C[s][d,g] = poisson_rng(lambda_Ctot_S[s][d,g]);
        } else {
          C[s][d,g] = lambda_Ctot_S[s][d,g];
        }
        C_sum = C_sum + C[s][d,g];

        E[s][d,g] = lambda_E_S[s][d,g] * E_scale * L[d] * f_crab[f_stratum[d]];
        E_sum = E_sum + E[s][d,g];
      }
    }
  }
}
