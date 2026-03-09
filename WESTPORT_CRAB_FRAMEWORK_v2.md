# Westport Crab Creel Estimation Framework — Design Specification v2

## 1. Overview

This document specifies the full estimation framework for recreational Dungeness crab harvest at Westport / Grays Harbor, incorporating all crabbing modes, gear restrictions, and the commercial/charter fleet.

Season totals are computed as the **sum of three independently-estimated population components**, each with its own effort and catch model:

```
Total Westport Harvest = Shore Crabbers + Private Boat Crabbers + Commercial/Charter Boats
```

Each component produces posterior draws (BSS) or point estimates (PE) that can be summed for the port total with propagated uncertainty.

---

## 2. Season Structure

The 2024-25 season (Sep 16, 2024 – Sep 15, 2025) is split into **two gear-regime sub-seasons** for the BSS model:

| Sub-season | Dates | Gear Allowed | Key Characteristics |
|---|---|---|---|
| **Ring-net only** | Sep 16 – Nov 30 (76 days) | Ring nets, snares, traps (foldable), handlines — NO POTS | Lower CPUE, moderate effort |
| **All-gear (pot open)** | Dec 1 – Sep 15 (289 days) | All gear including pots | Higher CPUE, higher effort |

Each BSS model run covers one sub-season independently. PE estimates are computed for both sub-seasons and summed.

### Commercial/Charter Season
Commercial vessels participate recreationally primarily **before the commercial opener**. The commercial opener date is a user-specified parameter (default: Jan 15, 2025 for 2024-25 season). After this date, commercial/charter harvest is assumed negligible.

---

## 3. Population Components

### 3.1 Shore Crabbers (BSS model, G=1)

**Who:** People crabbing from docks (Float 20, Float 17-21), the Westport Jetty, and any beach locations.

**Effort data:**
- Gear counts at docks (primary: 493 obs across 186 days)
- No jetty effort counts in 2024-25 (28 jetty interviews contribute to CPUE only)
- No beach effort counts (1 interview, negligible)

**Interview data:** 3,751 interviews (Dock + Jetty + Beach modes at all GH sites)

**Effort observation model:**
```
Gear_count ~ Poisson(lambda_E_shore * eps_H * R_G)
```
Where R_G (gear per crabber, ~1.27) is estimated from interviews.

**CPUE model:**
```
crab_caught ~ NegBin(lambda_C_shore * crabber_hours, r_C)
```

**Run twice:** Once for ring-net sub-season, once for all-gear sub-season.

### 3.2 Private Boat Crabbers (BSS model, G=1 or G=2 joint with shore)

**Who:** Private recreational boats launching from Westport Boat Launch (and Ocean Shores BL).

**Effort data:**
- Trailer counts at Westport BL (177 non-zero obs across 174 days)
- Trailer counts at Ocean Shores BL (4 obs, sparse)

**Interview data:** ~405 interviews
- 301 at Westport BL with mode=Boat
- ~104 at docks with mode=Boat (private boats that pulled up to dock for interview)

**Effort observation model:**
```
Trailer_count ~ Poisson(lambda_E_boat * eps_H * R_T)
```
Where R_T (trailers per boat group, ~0.5–1.0) is estimated from interviews.

**CPUE model:**
```
crab_caught ~ NegBin(lambda_C_boat * crabber_hours, r_C)
```
Note: Boat crabbers have substantially higher CPUE than shore crabbers (~4.8 vs ~1.1 crab/trip), so separate CPUE processes are appropriate.

**Run twice:** Ring-net sub-season and all-gear sub-season.

### 3.3 Commercial/Charter Boats (Census/Tally method, NOT BSS)

**Who:** Commercial crab vessels and charter boats engaged in recreational crabbing, operating from Westport Marina.

**Effort data:**
- Daily commercial/charter/private tally from the "wes commercial tally" sheet
- Very limited gear count data from Marina (5 obs)

**Interview data:** 201 interviews at Westport Marina with mode=Boat

**Estimation method:** Direct census expansion, not BSS:
```
Daily harvest = (vessels_observed × interview_rate × mean_crab_per_vessel)
             + (vessels_missed × imputed_crab_per_vessel)
```
OR simply:
```
Daily harvest = total_vessels × mean_crab_per_interviewed_vessel
```

**Active period:** Sep 16 through commercial opener (Jan 15, 2025 this season). Negligible after.

---

## 4. Parameter Structure

```r
params <- list(
  # --- Metadata ---
  project_name      = "Coastal Recreational Crab",
  fishery_name      = "Rec Crab Grays Harbor Westport 2024-25",
  est_date_start    = "2024-09-16",
  est_date_end      = "2025-09-15",
  season_filter     = "2024-25",

  # --- Gear restriction dates ---
  pot_open_date     = "2024-12-01",   # Pots allowed starting this date
  pot_close_date    = "2025-09-16",   # Pots prohibited starting (next season)
  # Sub-season 1 (ring-net only): est_date_start to pot_open_date - 1
  # Sub-season 2 (all-gear):      pot_open_date  to est_date_end

  # --- Commercial/charter season ---
  commercial_opener = "2025-01-15",   # Commercial season opens; charter/commercial
                                      # recreational effort drops to negligible

  # --- Catch groups ---
  est_catch_groups = data.frame(
    species = c("Dungeness", "Red_Rock"),
    fate    = c("Kept", "Kept"),
    stringsAsFactors = FALSE
  ),

  # --- Day type ---
  days_wkend        = c("Friday", "Saturday", "Sunday"),
  min_fishing_time  = 0.5,

  # --- Time stratification ---
  period_pe         = "month",
  period_bss        = "month",

  # --- Sections ---
  # S=1 for now (all Westport); expandable to S=2+ for finer spatial strata
  sections          = c(1),

  # --- BSS model ---
  bss_model_file    = "BSS_crab_model_01.stan",
  bss_chains        = 4,
  bss_iter          = 2000,
  bss_warmup        = 1000,
  bss_adapt_delta   = 0.9,
  bss_max_treedepth = 10,
  bss_cores         = 4,
  bss_max_interviews = 1000,  # NULL to use all; integer to subsample
  bss_max_count_seq  = 3,

  # --- Export ---
  model_used        = "Both models",
  data_grade        = "provisional",
  export            = "local"
)
```

---

## 5. Estimation Workflow

```
┌─────────────────────────────────────┐
│  1. Load & Transform Data           │
│     fetch_crab_data(params)         │
│     → dwg (effort, interview, catch)│
└─────────────┬───────────────────────┘
              │
              ▼
┌─────────────────────────────────────┐
│  2. Prep Days & Summarize           │
│     prep_days_crab()                │
│     prep_dwg_crab_summary()         │
└─────────────┬───────────────────────┘
              │
    ┌─────────┼──────────┐
    ▼         ▼          ▼
┌────────┐ ┌────────┐ ┌──────────────┐
│ SHORE  │ │ BOAT   │ │ COMM/CHARTER │
│ BSS+PE │ │ BSS+PE │ │ Census PE    │
└───┬────┘ └───┬────┘ └──────┬───────┘
    │          │              │
    │   ┌──────┘              │
    │   │  For each:          │
    │   │  ├─ Sub-season 1    │
    │   │  │  (ring-net only) │
    │   │  └─ Sub-season 2    │
    │   │     (all-gear)      │
    │   │                     │
    ▼   ▼                     ▼
┌─────────────────────────────────────┐
│  3. Combine Estimates               │
│     Shore + Boat + Comm/Charter     │
│     Sum posterior draws for BSS     │
│     Sum PEs with propagated SE      │
│     → Port total with uncertainty   │
└─────────────┬───────────────────────┘
              │
              ▼
┌─────────────────────────────────────┐
│  4. Output                          │
│     Tables: mode × sub-season       │
│     Plots: daily, monthly, totals   │
│     CSVs: stratum-level results     │
└─────────────────────────────────────┘
```

---

## 6. Mode-to-Data Mapping

### Effort Count Assignment

| Effort Site | Count Type | Assigned To | Notes |
|---|---|---|---|
| Westport Docks Float 20 | Gear count | Shore (g=1) | Primary shore effort |
| Westport Docks Float 17-21 | Gear count | Shore (g=1) | Paired with F20 by time |
| Westport Boat Launch | Trailer count | Private boat (g=2) | Primary boat effort |
| Ocean Shores Boat Launch | Trailer count | Private boat (g=2) | Sparse but included |
| Westport Marina | Gear count | Commercial/charter | Separate census approach |
| Westport Jetty | None in 2024-25 | Shore CPUE only | No effort counts available |
| Damon Point | None | Shore CPUE only | 1 beach interview |

### Interview Assignment

| Interview Area × Mode | Assigned To | Rationale |
|---|---|---|
| Float 20/17-21 × Dock | Shore (g=1) | Dock crabbers |
| Float 20/17-21 × Jetty | Shore (g=1) | Jetty crabbers (28 int) |
| Float 20/17-21 × Boat (Private) | Private boat (g=2) | Private boats interviewed at dock |
| Westport BL × Boat | Private boat (g=2) | Primary boat crabber source |
| Westport Marina × Boat (Private) | Private boat (g=2) | Private boats moored at marina |
| Westport Marina × Boat (Comm/Charter) | Commercial/charter | Separate estimation |
| Westport Marina × Dock | Shore (g=1) | Dock crabbers at marina area |
| Damon Point × Beach | Shore (g=1) | Beach crabber |

**Note on Marina interviews:** Need to split by `boat_type` field:
- Private → Private boat population
- Commercial/Charter → Commercial/charter population
- Dock mode at Marina → Shore population

---

## 7. BSS Model Runs Required

For a complete port estimate with both CPUE methods and both catch groups:

| Run | Population | Sub-season | Catch Group | ~Runtime (4 cores) |
|---|---|---|---|---|
| 1 | Shore | Ring-net (76 days) | Dungeness | ~15 min |
| 2 | Shore | All-gear (289 days) | Dungeness | ~45 min |
| 3 | Shore | Ring-net | Red Rock | ~15 min |
| 4 | Shore | All-gear | Red Rock | ~45 min |
| 5 | Boat | Ring-net | Dungeness | ~10 min |
| 6 | Boat | All-gear | Dungeness | ~30 min |
| 7 | Boat | Ring-net | Red Rock | ~10 min |
| 8 | Boat | All-gear | Red Rock | ~30 min |

Total: ~3.5 hours for full analysis. Can be reduced by:
- Dropping Red Rock (runs 3,4,7,8 → saves ~2 hrs)
- Subsampling interviews more aggressively
- Combining sub-seasons with pot-season covariate (halves the number of fits)

---

## 8. Next Steps

1. **Immediate:** Update `prep_inputs_bss_crab_v2.R` to handle sub-season splitting and mode assignment
2. **Immediate:** Update `fetch_crab_data()` to assign interview mode correctly (especially Marina split by boat_type)
3. **Short-term:** Build commercial/charter census estimation function
4. **Short-term:** Build port-total combination function (sum posteriors)
5. **Future:** Add ingress/egress data when digitized (improves effort expansion)
6. **Future:** Add Westport Jetty effort counts (when collected)
7. **Future:** Expand to Willapa Bay and Columbia River ports
