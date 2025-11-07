# --- 0) Setup ---------------------------------------------------------------
library(tidyverse)

# You already have these in your session:
#   read_airway(path)  -> tibble with gen, nj_loss_* columns
#   region_of_gen(int) -> "proximal"/"mid"/"distal"

# --- 1) Declare all runs for multiple dp -----------------------------------
root_refine  <- "/hpc/gjin250/refine_1d/particle_development/output"
root_takeout <- "/hpc/gjin250/Takeout/Thea_PhD/Paper 1_ method"

# Choose the diameters you want to audit
dp_list <- c(0.01, 0.05, 0.1, 1.0)

# Helper: map dp -> terminal leaf like p001/p005/p010/p100
dp_leaf <- function(dp) {
  if (abs(dp - 0.01) < 1e-9) return("p001")
  if (abs(dp - 0.05) < 1e-9) return("p005")
  if (abs(dp - 0.10) < 1e-9) return("p010")
  if (abs(dp - 1.00) < 1e-9) return("p100")
  sprintf("p%03d", as.integer(round(dp * 100))) # fallback
}
dp_leaf

# Build the run list by expanding over sbj, version, dp
runs <- tidyr::expand_grid(
  sbj     = c("H653", "H12816"),
  version = c("orig", "rev"),
  dp      = dp_list
) %>%
  mutate(
    Sex  = if_else(sbj == "H653", "Female", "Male"),
    Flow = "baseline",
    # root location depends on subject+version
    root = case_when(
      sbj == "H12816" & version == "orig" ~ root_takeout,
      TRUE                               ~ root_refine
    ),
    # subpath rule per subject+version
    subpath = pmap_chr(list(sbj, version, dp), function(s, v, d) {
      d2 <- d#  sprintf("%.2f", d)
      if (s == "H653"   && v == "orig") return(file.path("H653 (copy)", d2, "airway_result.txt"))
      if (s == "H653"   && v == "rev")  return(file.path("H653",         d2, "airway_result.txt"))
      if (s == "H12816" && v == "orig") return(file.path("H12816", "250ml", "d1.0", dp_leaf(d), "airway_result.txt"))
      if (s == "H12816" && v == "rev")  return(file.path("H12816",       d2, "airway_result.txt"))
      NA_character_
    }),
    airway_path = file.path(root, subpath),
    exists      = file.exists(airway_path)
  )

if (!all(runs$exists)) {
  warning("Missing files:\n", paste(runs$airway_path[!runs$exists], collapse = "\n"))
  stop("Fix the paths above and retry.")
}


fmt_dp <- function(d) sprintf("%.2f", d)   # "1.00", "0.10", "0.05", "0.01"

runs <- tidyr::expand_grid(
  sbj     = c("H653", "H12816"),
  version = c("orig", "rev"),
  dp      = dp_list
) %>%
  mutate(
    Sex  = if_else(sbj == "H653", "Female", "Male"),
    Flow = "baseline",
    
    # H12816: rev = Takeout; orig = refine_1d
    root = dplyr::case_when(
      sbj == "H12816" & version == "rev"  ~ root_takeout,
      TRUE                                 ~ root_refine
    ),
    
    subpath = pmap_chr(list(sbj, version, dp), function(s, v, d) {
      d2 <- fmt_dp(d)
      if (s == "H653"   && v == "orig") return(file.path("H653 (copy)", d2, "airway_result.txt"))
      if (s == "H653"   && v == "rev")  return(file.path("H653",         d2, "airway_result.txt"))
      if (s == "H12816" && v == "orig") return(file.path("H12816",       d2, "airway_result.txt"))
      if (s == "H12816" && v == "rev")  return(file.path("H12816", "250ml", "d1.0", dp_leaf(d), "airway_result.txt"))
      NA_character_
    }),
    
    airway_path = file.path(root, subpath),
    exists      = file.exists(airway_path)
  )

# Helpful sanity check
runs %>% filter(dp == 1.0) %>% select(sbj, version, airway_path, exists) %>% print(n = Inf)

if (!all(runs$exists)) {
  warning("Missing files:\n", paste(runs$airway_path[!runs$exists], collapse = "\n"))
  stop("Fix the paths above and retry.")
}





# --- 2) Load airway rows & tag region --------------------------------------
all_bronch <- runs %>%
  mutate(data = map(airway_path, read_airway)) %>%
  unnest(data) %>%
  mutate(
    gen    = as.integer(gen),
    region = region_of_gen(gen),
    zone   = "bronchial",
    # nice ordered factors for facets
    Sex    = factor(Sex, levels = c("Female","Male")),
    region = factor(region, levels = c("distal","mid","proximal")),
    dp_f   = factor(sprintf("%.2f", dp), levels = sprintf("%.2f", sort(unique(dp))))
  )

# --- 3) Aggregate: bronchi × mechanism × region ----------------------------
bronch_mech <- all_bronch %>%
  group_by(version, sbj, Sex, Flow, dp, dp_f, region) %>%
  summarise(
    Diffusion     = sum(nj_loss_dif, na.rm = TRUE),
    Impaction     = sum(nj_loss_imp, na.rm = TRUE),
    Sedimentation = sum(nj_loss_sed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Diffusion, Impaction, Sedimentation),
               names_to = "mech", values_to = "mass")

# Shares within bronchi (per version/subject/dp/region)
bronch_share <- bronch_mech %>%
  group_by(version, sbj, Sex, Flow, dp, dp_f, region) %>%
  mutate(share = mass / sum(mass)) %>%
  ungroup()

# --- 4) Pair orig vs rev & compute deltas (for ALL dps) --------------------
delta_all <- bronch_share %>%
  select(version, sbj, Sex, Flow, dp, dp_f, region, mech, share, mass) %>%
  pivot_wider(names_from = version, values_from = c(share, mass)) %>%
  mutate(
    delta_pp     = 100 * (share_rev - share_orig),                 # percentage points
    delta_ppm    = 1e6 * (share_rev - share_orig),                 # ppm of bronchial share
    delta_masspct = 100 * (mass_rev - mass_orig) / pmax(1e-30, mass_orig)
  )

# Quick audit: max |Δshare| (ppm) by dp/sex/region
audit <- delta_all %>%
  group_by(dp_f, Sex, region) %>%
  summarise(max_abs_ppm = max(abs(delta_ppm), na.rm = TRUE), .groups = "drop") %>%
  arrange(dp_f, Sex, region)
print(audit, n = Inf)

# --- 5) Plots ---------------------------------------------------------------
# (A) Stacked shares (orig vs rev) — facet by Sex × Region × dp
ggplot(bronch_share, aes(version, share, fill = mech)) +
  geom_col(position = "fill") +
  facet_grid(rows = vars(Sex), cols = vars(region, dp_f)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(
    title = "Bronchial mechanism split (orig vs revised) across particle sizes",
    x = NULL, y = "Fraction of bronchial deposition", fill = NULL
  ) +
  theme_classic(base_size = 12)

# (B) Δ share in ppm — facet by Sex × Region × dp
ggplot(delta_all, aes(mech, delta_ppm, fill = mech)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col() +
  facet_grid(rows = vars(Sex), cols = vars(region, dp_f)) +
  scale_fill_viridis_d(end = 0.9) +
  labs(
    title = "Change in bronchial mechanism share (revised − original)",
    subtitle = "Values in ppm of bronchial deposition; facets by Sex × Region × dₚ",
    x = NULL, y = "Δ share (ppm)"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

