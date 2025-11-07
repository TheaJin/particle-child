# A) Setup helpers (once per session)

# packages
library(dplyr); library(tidyr); library(ggplot2)
library(readr); library(purrr); library(stringr)
library(scales); library(viridisLite); library(conflicted)
conflicted::conflicts_prefer(dplyr::filter)  # avoid stats::filter collisions


mm     <- 1e-3        # m per mm
mu_air <- 1.8e-5      # Pa·s
rho_p  <- 1000        # kg/m^3
g      <- 9.81        # m/s^2

Cc <- function(dp_m){ 
  Kn <- 2*6.6e-8/dp_m
  1 + Kn*(1.257 + 0.4*exp(-1.1/Kn))
}
brownian_D <- function(dp_m, T=298, mu_visc=mu_air){
  kB <- 1.380649e-23
  kB*T*Cc(dp_m)/(3*pi*mu_visc*dp_m)
}
calc_Stk <- function(dp_m, U, r_m, rho=rho_p, mu_visc=mu_air){
  rho*dp_m^2*U/(18*mu_visc*r_m)
}
settling_velocity <- function(dp_m, rho_p=rho_p, rho_air=1.2, mu_visc=mu_air){
  # Stokes settling with slip correction
  (rho_p - rho_air)*g*dp_m^2*Cc(dp_m)/(18*mu_visc)
}
region_of_gen <- function(g){
  dplyr::case_when(g <= 6 ~ "proximal",
                   g <= 14 ~ "mid",
                   TRUE    ~ "distal")
}

column_names <- c(
  "ne","x","y","z","distance","gen","horsf_ord","lobe","ne_length","ne_radius",
  "ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel"
)

read_airway <- function(path) {
  df <- utils::read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
  colnames(df) <- column_names
  num_cols <- c("ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel")
  df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(gsub("D","E", x)))
  df %>% mutate(ne = as.integer(ne))
}

# B) Tell me where your adult runs live

#base_dir <- "/hpc/gjin250/refine_1d/particle_development/output"  # EDIT if needed
base_dir <- "/hpc/gjin250/Takeout/Thea_PhD/Paper 1_ method/"
adult_runs <- tibble::tribble(
  ~sbj,    ~sex, ~flow,  ~dp,   ~path,
  "H653",   "F",  "Low",  0.01, file.path(base_dir, "H653",  "1G", "100ml", "p001", "airway_result.txt"),
  "H653",   "F",  "High", 0.01, file.path(base_dir, "H653", "1G",   "500ml", "p001",  "airway_result.txt"),
  "H12816", "M",  "Low",  0.01, file.path(base_dir, "H12816", "100ml", "p001",  "airway_result.txt"),
  "H12816", "M",  "High", 0.01, file.path(base_dir, "H12816", "500ml", "p001", "airway_result.txt"),
  "H5977",  "M",  "Low",  0.01, file.path(base_dir, "H5977", "100ml", "p001", "airway_result.txt"),
  "H5977",  "M",  "High", 0.01, file.path(base_dir, "H5977",  "500ml", "p001",  "airway_result.txt")
)

adult_runs <- adult_runs |> dplyr::mutate(exists = file.exists(path))
print(adult_runs, n = Inf)



# C) Load everything and tag metadata
adult_all <- adult_runs %>%
  mutate(data = purrr::map(path, read_airway)) %>%
  tidyr::unnest(data) %>%
  mutate(
    dp_m   = dp*1e-6,
    r_m    = as.numeric(ne_radius)*mm,
    A      = pi*r_m^2,
    U      = as.numeric(ne_Vdot)/A,
    region = region_of_gen(as.integer(gen)),
    Sex    = dplyr::recode(sex, "F"="Female", "M"="Male")
  )

# D) Core deliverables for Rec 5.2.2
# D1. Bronchial vs Alveolar fractions by sex & flow (dp = 0.01 µm)
# bd_ad_by_sex <- adult_all %>%
#   dplyr::filter(abs(dp - 0.01) < 1e-9) %>%
#   group_by(Sex, flow, region) %>%
#   summarise(dif = sum(nj_loss_dif, na.rm=TRUE),
#             imp = sum(nj_loss_imp, na.rm=TRUE),
#             sed = sum(nj_loss_sed, na.rm=TRUE), .groups="drop") %>%
#   mutate(bronch = dif + imp + sed,
#          zone   = ifelse(region=="distal","alveolar","bronchial")) %>%
#   group_by(Sex, flow, zone) %>%
#   summarise(dep = sum(bronch), .groups="drop") %>%
#   group_by(Sex, flow) %>%
#   mutate(frac = dep/sum(dep)) %>%
#   ungroup()
# 
# ggplot(bd_ad_by_sex, aes(flow, frac, fill = zone)) +
#   geom_col(position = "fill", width=0.7) +
#   facet_wrap(~ Sex) +
#   scale_y_continuous(labels = scales::percent) +
#   scale_fill_viridis_d(end = 0.9, name = NULL) +
#   labs(title = "0.01 µm: Bronchial vs alveolar deposition by sex across flow",
#        y = "Fraction of total", x = "Flow setting") +
#   theme_classic(base_size = 12) +
#   theme(legend.position = "bottom")

bd_conducting <- adult_all %>%          # <- built from airway_result.txt
  dplyr::filter(dp == 0.01) %>%
  mutate(region = region_of_gen(as.integer(gen)),   # "proximal/mid/distal" for gens
         zone   = factor(region,
                         levels = c("proximal","mid","distal"),
                         labels = c("proximal (0–6)","mid (7–14)","distal (15–28)"))) %>%
  group_by(Sex, flow, zone) %>%
  summarise(dep = sum(nj_loss_dif + nj_loss_imp + nj_loss_sed, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Sex, flow) %>%
  mutate(frac = dep / sum(dep)) %>%
  ungroup()

ggplot(bd_conducting, aes(flow, frac, fill = zone)) +
  geom_col(position = "fill") +
  facet_wrap(~ Sex) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = "0.01 µm: conducting-airway deposition by sex across flow",
       y = "Fraction of conducting deposition", x = "Flow setting", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")




# D2. Mechanism contributions in the bronchial zone (dp = 0.01 µm)
mech_bronch <- adult_all %>%
  dplyr::filter(abs(dp - 0.01) < 1e-9, region != "distal") %>%
  group_by(Sex, flow) %>%
  summarise(
    Diffusion   = sum(nj_loss_dif, na.rm=TRUE),
    Impaction   = sum(nj_loss_imp, na.rm=TRUE),
    Sedimentation = sum(nj_loss_sed, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(-c(Sex,flow), names_to="mechanism", values_to="mass") %>%
  group_by(Sex, flow) %>%
  mutate(frac = mass/sum(mass)) %>%
  ungroup()

ggplot(mech_bronch, aes(flow, frac, fill = mechanism)) +
  geom_col(position="fill", width=0.7) +
  facet_wrap(~ Sex) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9, name = NULL) +
  labs(title = "0.01 µm (bronchial only): mechanism fractions by sex vs flow",
       y = "Fraction within bronchial zone", x = "Flow setting") +
  theme_classic(base_size=12) + theme(legend.position="bottom")

# D3. Prove impaction is negligible at 0.01 µm
pe_stk_tbl <- adult_all %>%
  dplyr::filter(abs(dp - 0.01) < 1e-9) %>%
  mutate(
    Stk = calc_Stk(dp_m, U, r_m),
    Pe  = U*(2*r_m)/brownian_D(dp_m)
  ) %>%
  group_by(Sex, flow) %>%
  summarise(
    Stk_p95 = quantile(Stk, 0.95, na.rm=TRUE),
    Pe_med  = median(Pe, na.rm=TRUE),
    .groups="drop"
  )

pe_stk_tbl


# D4. Explain the sex difference via Péclet distributions

library(ggridges)

pe_branch <- adult_all %>%
  dplyr::filter(abs(dp - 0.01) < 1e-9) %>%
  mutate(Pe = U*(2*r_m)/brownian_D(dp_m),
         log10Pe = log10(Pe))

ggplot(pe_branch, aes(x = log10Pe, y = interaction(Sex, flow), fill = Sex)) +
  ggridges::geom_density_ridges(alpha=0.7, scale=1.1, rel_min_height=0.01) +
  scale_fill_viridis_d(end=0.9, name=NULL) +
  labs(title="0.01 µm: distribution of log10(Pe) by sex and flow",
       x = expression(log[10]*"(Pe)"), y = "Sex × Flow") +
  theme_classic(base_size=12) +
  theme(legend.position="bottom")

# Δ bronchial fraction (Low − High), by sex, at dp = 0.01 µm
delta_bronch <- bd_ad_by_sex %>%
  dplyr::filter(zone == "bronchial") %>%
  dplyr::select(Sex, flow, frac) %>%      # <-- drop dep
  dplyr::distinct() %>%
  dplyr::mutate(flow = factor(flow, levels = c("High","Low"))) %>%
  tidyr::pivot_wider(names_from = flow, values_from = frac) %>%
  dplyr::mutate(
    delta_pp = 100 * (Low - High),
    High     = scales::percent(High, accuracy = 0.1),
    Low      = scales::percent(Low,  accuracy = 0.1),
    delta_pp = scales::number(delta_pp, accuracy = 0.1, signed = TRUE, suffix = " pp")
  )

delta_bronch





# Option B — True bronchial vs alveolar

term_colnames <- c(
  "np", "x", "y", "z", "distance", "lobe",
  "nu_vol", "nu_vdot0",
  "nu_mass", "nu_loss_dif", "nu_loss_sed", "ne_part_vel"
)

base_dir <- "/hpc/gjin250/Takeout/Thea_PhD/Paper 1_ method/"
path_to_terminal_txt <- tibble::tribble(
  ~sbj,    ~sex, ~flow,  ~dp,   ~path,
  "H653",   "F",  "Low",  0.01, file.path(base_dir, "H653",  "1G", "100ml", "p001", "terminal.txt"),
  "H653",   "F",  "High", 0.01, file.path(base_dir, "H653", "1G",   "500ml", "p001",  "terminal.txt"),
  "H12816", "M",  "Low",  0.01, file.path(base_dir, "H12816", "100ml", "p001",  "terminal.txt"),
  "H12816", "M",  "High", 0.01, file.path(base_dir, "H12816", "500ml", "p001", "terminal.txt"),
  "H5977",  "M",  "Low",  0.01, file.path(base_dir, "H5977", "100ml", "p001", "terminal.txt"),
  "H5977",  "M",  "High", 0.01, file.path(base_dir, "H5977",  "500ml", "p001",  "terminal.txt")
)

path_to_terminal_txt <- path_to_terminal_txt |> dplyr::mutate(exists = file.exists(path))
print(path_to_terminal_txt, n = Inf)


read_terminal <- function(path) {
  df <- utils::read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
  colnames(df) <- term_colnames
  num_cols <- c("nu_vol","nu_vdot0","nu_mass","nu_loss_dif","nu_loss_sed","ne_part_vel")
  df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(gsub("D","E", x)))
  df %>% dplyr::mutate(np = as.integer(np))
}

# airway (conducting) -> bronchial
bronch_sum <- adult_all %>%
  dplyr::summarise(bronch_mass = sum(nj_loss_dif + nj_loss_imp + nj_loss_sed, na.rm = TRUE))

# terminal file -> alveolar
term_df    <- read_terminal(path_to_terminal_txt)
alv_sum    <- term_df %>%
  dplyr::summarise(alv_mass = sum(nu_loss_dif + nu_loss_sed, na.rm = TRUE))

dep_zone <- dplyr::bind_cols(bronch_sum, alv_sum) %>%
  tidyr::pivot_longer(everything(), names_to = "zone", values_to = "mass") %>%
  dplyr::mutate(frac = mass / sum(mass),
                zone = dplyr::recode(zone,
                                     bronch_mass = "bronchial",
                                     alv_mass    = "alveolar"))

ggplot(dep_zone, aes(x = zone, y = frac, fill = zone)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(x = NULL, y = "Fraction of total", fill = NULL)









# -------------------------%

read_terminal <- function(path) {
  # Generic reader; adjust if your terminal file has headers
  df <- utils::read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
  # Try to coerce all that look numeric (D-notation to E)
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(gsub("D","E", x))))
  # Keep numeric cols and sum them as deposition; refine if you know exact names
  num <- dplyr::select(df, where(is.numeric))
  tibble::tibble(alv_dep = rowSums(num, na.rm = TRUE))
}

# Build one table per run
bronch_tbl <- adult_runs %>%
  dplyr::mutate(air_path = path) %>%
  rowwise() %>%
  dplyr::mutate(bronch = {
    air <- read_airway(air_path)
    sum(air$nj_loss_dif + air$nj_loss_imp + air$nj_loss_sed, na.rm = TRUE)
  }) %>%
  ungroup() %>%
  dplyr::select(sbj, sex, flow, dp, bronch)

alv_tbl <- adult_runs %>%
  dplyr::mutate(term_path = file.path(dirname(path), "terminal.txt")) %>%  # adjust name if needed
  rowwise() %>%
  dplyr::mutate(alveolar = {
    if (file.exists(term_path)) sum(read_terminal(term_path)$alv_dep, na.rm = TRUE) else NA_real_
  }) %>%
  ungroup() %>%
  dplyr::select(sbj, sex, flow, dp, alveolar)

bd_ad_by_sex <- bronch_tbl %>%
  dplyr::inner_join(alv_tbl, by = c("sbj","sex","flow","dp")) %>%
  dplyr::filter(dp == 0.01) %>%
  dplyr::mutate(Sex = dplyr::recode(sex, F = "Female", M = "Male"),
                FlowSetting = flow) %>%
  tidyr::pivot_longer(c(bronch, alveolar), names_to = "zone", values_to = "dep") %>%
  dplyr::group_by(Sex, FlowSetting) %>%
  dplyr::mutate(frac = dep / sum(dep)) %>%
  dplyr::ungroup()

ggplot(bd_ad_by_sex, aes(FlowSetting, frac, fill = zone)) +
  geom_col(position = "fill") +
  facet_wrap(~ Sex) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = "0.01 µm: Bronchial vs alveolar deposition by sex across flow",
       y = "Fraction of total", x = "Flow setting", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")









