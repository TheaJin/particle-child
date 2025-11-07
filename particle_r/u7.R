library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(scales)


base_dir <- "/hpc/gjin250/refine_1d/particle_development/output"
adult_runs <- tribble(
  ~sbj,   ~Sex, ~Flow,      ~version, ~dp,   ~airway_path,
  "H653", "F",  "baseline", "orig",    0.01, file.path(base_dir, "H653 (copy)", "0.01", "airway_result.txt"),
  "H653", "F",  "baseline", "rev",     0.01, file.path(base_dir, "H653",        "0.01", "airway_result.txt")
) %>%
  mutate(exists = file.exists(airway_path))


base_dir2 <- "/hpc/gjin250/Takeout/Thea_PhD/Paper 1_ method/"
adult_runs <- tribble(
  ~sbj,   ~Sex, ~Flow,      ~version, ~dp,   ~airway_path,
  "H12816", "M",  "baseline", "orig",    0.01, file.path(base_dir, "H12816", "0.01", "airway_result.txt"),
  "H12816", "M",  "baseline", "rev",     0.01, file.path(base_dir2, "H12816", "250ml", "d1.0",  "p001", "airway_result.txt")
) %>%
  mutate(exists = file.exists(airway_path))

adult_runs
stopifnot(all(adult_runs$exists))   # fail fast if a path is wrong


# --- load all airway rows ---
all_bronch <- adult_runs %>%
  mutate(data = map(airway_path, read_airway)) %>%
  unnest(data) %>%
  mutate(region = region_of_gen(gen),
         zone   = "bronchial")

# --- aggregate by mechanism in bronchi ---
bronch_mech <- all_bronch %>%
  group_by(version, sbj, Sex, Flow, dp, region) %>%
  summarise(
    dif = sum(nj_loss_dif, na.rm=TRUE),
    imp = sum(nj_loss_imp, na.rm=TRUE),
    sed = sum(nj_loss_sed, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(dif:sed, names_to="mech", values_to="mass") %>%
  mutate(mech = recode(mech, dif="Diffusion", imp="Impaction", sed="Sedimentation"))

# --- normalise to shares within bronchi ---
bronch_share <- bronch_mech %>%
  group_by(version, sbj, Sex, Flow, dp, region) %>%
  mutate(total_bronch = sum(mass),
         share = ifelse(total_bronch > 0, mass/total_bronch, NA_real_)) %>%
  ungroup()

# --- pair orig vs rev and compute deltas ---
delta_001 <- bronch_share %>%
  filter(abs(dp - 0.01) < 1e-9) %>%
  select(version, sbj, Sex, Flow, dp, region, mech, share, mass) %>%
  pivot_wider(names_from = version, values_from = c(share, mass)) %>%
  mutate(
    delta_pp      = 100 * (share_rev - share_orig),                         # change in share (pp)
    delta_masspct = 100 * (mass_rev  - mass_orig) / pmax(1e-30, mass_orig)  # % mass change
  )

delta_001 %>% arrange(region, mech)


# Stacked shares (orig vs rev)
ggplot(bronch_share %>% filter(abs(dp-0.01)<1e-9),
       aes(version, share, fill=mech)) +
  geom_col(position="fill") +
  facet_wrap(~ region, nrow=1) +
  scale_y_continuous(labels=percent) +
  labs(title="0.01 µm: bronchial mechanism split (orig vs revised)",
       x=NULL, y="Fraction of bronchial deposition", fill=NULL) +
  theme_classic(base_size=12) +
  scale_fill_viridis_d(end=0.9)

# ∆ share (pp)
ggplot(delta_001, aes(mech, delta_pp, fill=mech)) +
  geom_hline(yintercept=0, linewidth=0.3) +
  geom_col() +
  facet_wrap(~ region, nrow=1) +
  labs(title="Change in bronchial mechanism share (revised − original) at 0.01 µm",
       x=NULL, y="∆ share (pp)") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  scale_fill_viridis_d(end=0.9)










# --- helpers (reuse your existing readers) ---
air_cols <- c("ne","x","y","z","distance","gen","horsf_ord","lobe","ne_length","ne_radius",
              "ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel")

read_airway <- function(path){
  df <- utils::read.table(path, header=FALSE, sep="", fill=FALSE, comment.char="")
  colnames(df) <- air_cols
  num <- c("ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel")
  df[num] <- lapply(df[num], function(x) as.numeric(gsub("D","E",x)))
  df %>% mutate(ne=as.integer(ne))
}


region_of_gen <- function(g){
  cut(as.integer(g),
      breaks=c(-Inf,6,14,Inf),
      labels=c("proximal (0–6)","mid (7–14)","distal (15–28)"),
      right=TRUE)
}

# --- load one run (bronchial only) ---
load_bronchial <- function(run_row){
  aw <- read_airway(run_row$airway_path) %>%
    mutate(
      version = run_row$version,
      sbj     = run_row$sbj,
      Sex     = run_row$Sex,
      Flow    = run_row$Flow,   # e.g. "High"/"Low"
      dp      = run_row$dp,
      region  = region_of_gen(gen),
      zone    = "bronchial"
    )
  aw
}

# runs: a data.frame with columns version/sbj/Sex/Flow/dp/airway_path
all_bronch <- map_dfr(split(runs, seq_len(nrow(runs))), ~load_bronchial(.x))

# --- 1) aggregate by mechanism within bronchi ---
bronch_mech <- all_bronch %>%
  group_by(version, sbj, Sex, Flow, dp, region) %>%
  summarise(
    dif = sum(nj_loss_dif, na.rm=TRUE),
    imp = sum(nj_loss_imp, na.rm=TRUE),
    sed = sum(nj_loss_sed, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(dif:sed, names_to="mech", values_to="mass") %>%
  mutate(mech = recode(mech, dif="Diffusion", imp="Impaction", sed="Sedimentation"))

# --- 2) normalise to shares within bronchi (per run) ---
bronch_share <- bronch_mech %>%
  group_by(version, sbj, Sex, Flow, dp, region) %>%
  mutate(total_bronch = sum(mass),
         share = ifelse(total_bronch > 0, mass/total_bronch, NA_real_)) %>%
  ungroup()

# --- 3) pair orig vs rev and compute deltas ---
delta_share <- bronch_share %>%
  select(version, sbj, Sex, Flow, dp, region, mech, share, mass) %>%
  pivot_wider(names_from=version, values_from=c(share, mass)) %>%
  mutate(
    delta_pp     = 100 * (share_rev - share_orig),               # percentage points
    delta_mass_% = 100 * (mass_rev - mass_orig) / pmax(1e-30, mass_orig)
  )

# --- 4) focus on dp = 0.01 µm (primary) ---
delta_001 <- delta_share %>% filter(abs(dp - 0.01) < 1e-9)

# Summary table (mean ± SD pp) by Sex × Flow × region × mech
sum_tbl <- delta_001 %>%
  group_by(Sex, Flow, region, mech) %>%
  summarise(
    mean_pp = mean(delta_pp, na.rm=TRUE),
    sd_pp   = sd(delta_pp,   na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(Sex, Flow, region, mech)

# --- 5) plots ---
# a) stacked shares per version (visual split), faceted by Sex/Flow
plot_share <- bronch_share %>%
  filter(abs(dp - 0.01) < 1e-9) %>%
  group_by(version, Sex, Flow, mech) %>%
  summarise(share = mean(share, na.rm=TRUE), .groups="drop") %>%  # avg over subjects
  ggplot(aes(x=version, y=share, fill=mech)) +
  geom_col(position="fill") +
  facet_grid(Sex ~ Flow) +
  scale_y_continuous(labels=percent) +
  labs(title="0.01 µm: bronchial mechanism split (orig vs revised)",
       x=NULL, y="Fraction of bronchial deposition", fill=NULL) +
  theme_classic(base_size=12) +
  scale_fill_viridis_d(end=0.9)

# b) ∆ share (pp) bars by mechanism, faceted
plot_delta <- delta_001 %>%
  group_by(Sex, Flow, region, mech) %>%
  summarise(mean_pp = mean(delta_pp, na.rm=TRUE), .groups="drop") %>%
  ggplot(aes(mech, mean_pp, fill=mech)) +
  geom_hline(yintercept=0, linewidth=0.3) +
  geom_col() +
  facet_grid(Sex ~ Flow + region) +
  labs(title="0.01 µm: change in bronchial mechanism share (revised – original)",
       x=NULL, y="∆ share (percentage points)") +
  theme_classic(base_size=12) +
  theme(legend.position="none") +
  scale_fill_viridis_d(end=0.9)

# c) (optional) absolute mass change %
plot_mass <- delta_001 %>%
  group_by(Sex, Flow, region, mech) %>%
  summarise(mean_pct = mean(delta_mass_%, na.rm=TRUE), .groups="drop") %>%
                            ggplot(aes(mech, mean_pct, fill=mech)) +
                              geom_hline(yintercept=0, linewidth=0.3) +
                              geom_col() +
                              facet_grid(Sex ~ Flow + region) +
                              labs(title="0.01 µm: % change in bronchial mechanism mass (revised vs original)",
                                   x=NULL, y="% change") +
                              theme_classic(base_size=12) +
                              theme(legend.position="none") +
                              scale_fill_viridis_d(end=0.9)
                            
                            # --- 6) quick paired test (per mechanism) ---
                            stats_tbl <- delta_001 %>%
                              group_by(Sex, Flow, region, mech) %>%
                              summarise(
                                n     = sum(!is.na(delta_pp)),
                                p_wil = tryCatch(wilcox.test(share_rev, share_orig, paired=TRUE)$p.value, error=function(e) NA_real_),
                                .groups="drop"
                              )
                            