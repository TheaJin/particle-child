# --- 0) Setup ---------------------------------------------------------------
library(tidyverse)

# You already have these; keep as-is
# read_airway(path)      -> returns airway rows with columns incl. gen, nj_loss_*
# region_of_gen(gen_int) -> returns "proximal"/"mid"/"distal"

# --- 1) Declare all runs in ONE tibble -------------------------------------
# Put every (subject, version, dp) with its exact subpath under root.
root_refine <- "/hpc/gjin250/refine_1d/particle_development/output"
root_takeout <- "/hpc/gjin250/Takeout/Thea_PhD/Paper 1_ method"

runs <- tribble(
  ~sbj,     ~Sex, ~Flow,      ~version, ~dp,    ~root,        ~subpath,
  # H653 (F) – both versions under refine_1d, note the "(copy)" for orig
  "H653",   "F",  "baseline", "orig",    0.01,  root_refine,  file.path("H653 (copy)", "0.01", "airway_result.txt"),
  "H653",   "F",  "baseline", "rev",     0.01,  root_refine,  file.path("H653",         "0.01", "airway_result.txt"),
  
  # H12816 (M) – orig under refine_1d; rev under Takeout
  "H12816", "M",  "baseline", "orig",    0.01,  root_refine,  file.path("H12816", "0.01", "airway_result.txt"),
  "H12816", "M",  "baseline", "rev",     0.01,  root_takeout, file.path("H12816", "250ml", "d1.0", "p001", "airway_result.txt")
) |>
  mutate(airway_path = file.path(root, subpath),
         exists = file.exists(airway_path))

# Fail fast (or switch to a soft warning if you prefer)
if (!all(runs$exists)) {
  warning("Missing files:\n",
          paste(runs$airway_path[!runs$exists], collapse = "\n"))
  stop("Fix the paths above and retry.")
}

# --- 2) Load airway rows, tag regions --------------------------------------
all_bronch <- runs |>
  mutate(data = map(airway_path, read_airway)) |>
  unnest(data) |>
  mutate(
    gen    = as.integer(gen),
    region = region_of_gen(gen),
    zone   = "bronchial"
  )

# --- 3) Aggregate to bronchi × mechanism × region --------------------------
bronch_mech <- all_bronch |>
  group_by(version, sbj, Sex, Flow, dp, region) |>
  summarise(
    Diffusion     = sum(nj_loss_dif, na.rm = TRUE),
    Impaction     = sum(nj_loss_imp, na.rm = TRUE),
    Sedimentation = sum(nj_loss_sed, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(c(Diffusion, Impaction, Sedimentation),
               names_to = "mech", values_to = "mass")

# Shares within the bronchi (per version/subject/region)
bronch_share <- bronch_mech |>
  group_by(version, sbj, Sex, Flow, dp, region) |>
  mutate(share = mass / sum(mass)) |>
  ungroup()

# --- 4) Pair orig vs rev & compute deltas ----------------------------------
delta_001 <- bronch_share |>
  filter(abs(dp - 0.01) < 1e-9) |>
  select(version, sbj, Sex, Flow, dp, region, mech, share, mass) |>
  pivot_wider(names_from = version, values_from = c(share, mass)) |>
  mutate(
    delta_pp      = 100 * (share_rev - share_orig),                         # percentage points
    delta_masspct = 100 * (mass_rev  - mass_orig) / pmax(1e-30, mass_orig)  # percent change
  )

# Inspect tidy deltas
delta_001 |>
  arrange(region, mech)

# --- 5) Plots ---------------------------------------------------------------
# Stacked shares (orig vs rev)
ggplot(filter(bronch_share, abs(dp - 0.01) < 1e-9),
       aes(version, share, fill = mech)) +
  geom_col(position = "fill") +
  facet_wrap(~ region, nrow = 1) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(
    title = "0.01 µm: bronchial mechanism split (orig vs revised)",
    x = NULL, y = "Fraction of bronchial deposition", fill = NULL
  ) +
  theme_classic(base_size = 12)

# Δ share in percentage points
ggplot(delta_001, aes(mech, delta_pp, fill = mech)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col() +
  facet_wrap(~ region, nrow = 1) +
  scale_fill_viridis_d(end = 0.9) +
  labs(
    title = "Change in bronchial mechanism share (revised − original) at 0.01 µm",
    x = NULL, y = "Δ share (pp)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")



delta_001_pretty <- delta_001 %>%
  mutate(
    delta_ppm = 1e6 * (share_rev - share_orig),  # parts per million of bronchial share
    delta_pp  = 100  * (share_rev - share_orig)  # percentage points (what you had)
  )

delta_001_clean <- delta_001_pretty %>%
  mutate(
    sym_pct = 100 * (mass_rev - mass_orig) / if_else(mass_rev + mass_orig > 0,
                                                     0.5*(mass_rev + mass_orig), NA_real_),
    mass_note = case_when(
      mass_orig == 0 & mass_rev  > 0 ~ "0 → tiny",
      mass_rev  == 0 & mass_orig > 0 ~ "tiny → 0",
      TRUE                           ~ ""
    )
  )
delta_001_clean

ggplot(delta_001_pretty, aes(mech, delta_ppm, fill = mech)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col() +
  facet_grid(sbj ~ region) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = "Change in bronchial mechanism share (revised − original) at 0.01 µm",
       x = NULL, y = "Δ share (ppm of bronchial deposition)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")


audit <- delta_001_pretty %>%
  group_by(sbj, region) %>%
  summarise(max_abs_ppm = max(abs(delta_ppm), na.rm = TRUE),
            .groups = "drop")

audit






