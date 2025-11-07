# prints "1.0", "0.1", "0.05", "0.01"
fmt_dp <- function(d) {
  if (abs(d - 1.0) < 1e-12 || abs(d - 0.1) < 1e-12) {
    formatC(d, format = "f", digits = 1)
  } else {
    formatC(d, format = "f", digits = 2)
  }
}

root_refine  <- "/hpc/gjin250/refine_1d/particle_development/output"
root_takeout <- "/hpc/gjin250/Takeout/Thea_PhD/Paper 1_ method"

dp_list <- c(0.01, 0.05, 0.1, 1.0)

# Map dp -> terminal leaf for the Takeout layout
dp_leaf <- function(dp) {
  if (abs(dp - 0.01) < 1e-9) return("p001")
  if (abs(dp - 0.05) < 1e-9) return("p005")
  if (abs(dp - 0.10) < 1e-9) return("p010")
  if (abs(dp - 1.00) < 1e-9) return("p100")
  sprintf("p%03d", as.integer(round(dp * 100)))
}

runs <- tidyr::expand_grid(
  sbj     = c("H653", "H12816"),
  version = c("orig", "rev"),
  dp      = dp_list
) %>%
  mutate(
    Sex  = if_else(sbj == "H653", "Female", "Male"),
    Sex  = factor(Sex, levels = c("Female","Male")),   # order in facets
    Flow = "baseline",
    
    # H12816: rev lives under Takeout; others under refine_1d
    root = dplyr::case_when(
      sbj == "H12816" & version == "rev"  ~ root_takeout,
      TRUE                                 ~ root_refine
    ),
    
    subpath = pmap_chr(list(sbj, version, dp), function(s, v, d) {
      d_str <- fmt_dp(d)
      if (s == "H653"   && v == "orig") return(file.path("H653 (copy)", d_str, "airway_result.txt"))
      if (s == "H653"   && v == "rev")  return(file.path("H653",         d_str, "airway_result.txt"))
      if (s == "H12816" && v == "orig") return(file.path("H12816",       d_str, "airway_result.txt"))
      if (s == "H12816" && v == "rev")  return(file.path("H12816", "250ml", "d1.0", dp_leaf(d), "airway_result.txt"))
      NA_character_
    }),
    
    airway_path = file.path(root, subpath),
    exists      = file.exists(airway_path)
  )

# sanity check
runs %>% arrange(sbj, version, dp) %>% select(sbj, version, dp, airway_path, exists) %>% print(n = Inf)

stopifnot(all(runs$exists))


# 3) Use the same analysis, but let you pick any dp


dp_pick <- 1.0 # change to 0.05, 0.1, 1.0, etc.

all_bronch <- runs %>%
  mutate(data = purrr::map(airway_path, read_airway)) %>%
  tidyr::unnest(data) %>%
  mutate(gen = as.integer(gen),
         region = region_of_gen(gen),
         zone = "bronchial")

bronch_mech <- all_bronch %>%
  group_by(version, sbj, Sex, Flow, dp, region) %>%
  summarise(
    Diffusion     = sum(nj_loss_dif, na.rm=TRUE),
    Impaction     = sum(nj_loss_imp, na.rm=TRUE),
    Sedimentation = sum(nj_loss_sed, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(c(Diffusion, Impaction, Sedimentation),
               names_to="mech", values_to="mass")

bronch_share <- bronch_mech %>%
  group_by(version, sbj, Sex, Flow, dp, region) %>%
  mutate(share = mass / sum(mass)) %>%
  ungroup()

delta_dp <- bronch_share %>%
  filter(abs(dp - dp_pick) < 1e-12) %>%
  select(version, sbj, Sex, Flow, dp, region, mech, share, mass) %>%
  pivot_wider(names_from = version, values_from = c(share, mass)) %>%
  mutate(
    delta_pp      = 100 * (share_rev - share_orig),
    delta_ppm     = 1e6 * (share_rev - share_orig),
    delta_masspct = 100 * (mass_rev - mass_orig) / pmax(1e-30, mass_orig)
  )
delta_dp

# 4) Facet labels: use Sex instead of subject IDs

# Stacked shares (orig vs rev) at the chosen dp, faceted by Sex
ggplot(filter(bronch_share, abs(dp - dp_pick) < 1e-12),
       aes(version, share, fill = mech)) +
  geom_col(position = "fill") +
  facet_grid(Sex ~ region) +                      # <— here
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = sprintf("%.2g µm: bronchial mechanism split (orig vs revised)", dp_pick),
       x = NULL, y = "Fraction of bronchial deposition", fill = NULL) +
  theme_classic(base_size = 12)

# Δ share (ppm) faceted by Sex
ggplot(delta_dp, aes(mech, delta_ppm, fill = mech)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col() +
  facet_grid(Sex ~ region) +                      # <— here
  scale_fill_viridis_d(end = 0.9) +
  labs(title = sprintf("Change in bronchial mechanism share (rev − orig) at %.2g µm", dp_pick),
       x = NULL, y = "Δ share (ppm of bronchial deposition)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

