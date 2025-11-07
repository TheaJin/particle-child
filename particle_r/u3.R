library(dplyr); library(tidyr); library(ggplot2); library(scales); library(viridis)

# --- helpers ---------------------------------------------------------------
# if units are mm for radius/length, set these
mm  <- 1e-3  # m
mu  <- 1.8e-5 # air dynamic viscosity (Pa·s)
rho <- 1000   # particle density (kg/m^3), adjust if needed

# Cunningham slip + Brownian D for spherical particle of diameter dp (m), 298 K
Cc <- function(dp_m, P=101325){ # Hutchins 1995 form
  Kn <- 2*6.6e-8 / dp_m
  1 + Kn*(1.257 + 0.4*exp(-1.1/Kn))
}
D_B <- function(dp_m, T=298, mu=mu){ # m^2/s
  kB <- 1.380649e-23
  kB*T*Cc(dp_m) / (3*pi*mu*dp_m)
}
Stk <- function(dp_m, U, r_m){       # Stokes number proxy with tube radius
  rho*dp_m^2*U/(18*mu*r_m)
}
Pe  <- function(U, r_m, dp_m){        # tube Pe using 2r as length
  U*(2*r_m)/D_B(dp_m)
}

# classify region by generation
region_of_gen <- function(g){
  dplyr::case_when(
    g <= 6             ~ "proximal",
    g >= 7 & g <= 14   ~ "mid",
    TRUE               ~ "distal"     # 15+
  )
}



### Rec 5.2
# --- 1) FLOW & VELOCITY BY LOBE (Ref vs LLL-target) -----------------------
flow_cmp <- function(all_data_one_dp){
  all_data_one_dp %>%
    mutate(area = pi*(as.numeric(ne_radius)*mm)^2,
           U    = as.numeric(ne_Vdot)/area) %>%   # m/s if ne_Vdot in m^3/s
    group_by(condition, lobe) %>%
    summarise(Vdot_lobe = sum(ne_Vdot, na.rm=TRUE),
              U_mean    = weighted.mean(U, w = area, na.rm=TRUE),
              .groups="drop") %>%
    group_by(condition) %>%
    mutate(Vdot_share = Vdot_lobe/sum(Vdot_lobe)) %>%
    ungroup()
}

# example: pick dp=1.00 files to pull the flow fields
one_dp <- all_data %>%  dplyr::filter(dp == 1.00)   # your combined dataframe
flow_tbl <- flow_cmp(one_dp)

# pretty comparison Ref vs Cluster for un-obstructed lobes
flow_tbl %>%
  dplyr::filter(lobe != "LLL", condition %in% c("VC (1.0)","Cluster 12")) %>%
  select(condition, lobe, Vdot_lobe, Vdot_share, U_mean) %>%
  arrange(lobe, condition)


# Total flow check (is global V̇ different?)
flow_tbl %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(Vdot_total = sum(Vdot_lobe),
                   U_global   = stats::weighted.mean(U_mean, w = Vdot_lobe))

# Per-lobe percent deltas (Cluster12 vs VC (1.0))
delta_flow <- flow_tbl %>%
  dplyr::filter(lobe != "LLL",
                condition %in% c("VC (1.0)", "Cluster 12")) %>%
  dplyr::group_by(lobe) %>%
  dplyr::summarise(
    Vdot_ref = Vdot_lobe[condition == "VC (1.0)"][1],
    Vdot_c12 = Vdot_lobe[condition == "Cluster 12"][[1]],
    U_ref    = U_mean    [condition == "VC (1.0)"][1],
    U_c12    = U_mean    [condition == "Cluster 12"][[1]],
    pct_Vdot = 100 * (Vdot_c12 / Vdot_ref - 1),
    pct_U    = 100 * (U_c12    / U_ref    - 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(lobe)

# pretty print
delta_flow %>%
  dplyr::mutate(
    pct_Vdot = scales::number(pct_Vdot, accuracy = 0.1, suffix = " %"),
    pct_U    = scales::number(pct_U,    accuracy = 0.1, suffix = " %")
  )



## b)
library(dplyr); library(tidyr)

## ---- constants (units only matter for ratios; keep consistent)
mm     <- 1e-3        # m per mm
mu_air <- 1.8e-5      # Pa·s (dynamic viscosity of air)
rho_p  <- 1000        # kg/m^3 (particle density) – adjust if you use another

## ---- physics helpers (renamed to avoid collisions)
Cc <- function(dp_m){                        # Cunningham slip
  Kn <- 2*6.6e-8 / dp_m
  1 + Kn*(1.257 + 0.4*exp(-1.1/Kn))
}
brownian_D <- function(dp_m, temp_k = 298, mu_visc = mu_air){
  kB <- 1.380649e-23
  kB*temp_k*Cc(dp_m) / (3*pi*mu_visc*dp_m)
}
calc_Stk <- function(dp_m, U, r_m, rho = rho_p, mu_visc = mu_air){
  rho*dp_m^2*U / (18*mu_visc*r_m)
}
calc_Pe <- function(U, r_m, dp_m){
  U*(2*r_m) / brownian_D(dp_m)
}

## ---- aggregate Stokes & Péclet by lobe/condition for selected dp
dp_list <- c(1.0, 3.0, 5.0)  # µm

phys_tbl <- all_data %>%
  dplyr::filter(dp %in% dp_list) %>%
  dplyr::mutate(
    r_m  = as.numeric(ne_radius) * mm,
    A    = pi*r_m^2,
    U    = as.numeric(ne_Vdot) / A,
    dp_m = dp * 1e-6,
    Stk  = calc_Stk(dp_m, U, r_m),
    Pe   = calc_Pe(U, r_m, dp_m)
  ) %>%
  dplyr::group_by(condition, dp, lobe) %>%
  dplyr::summarise(
    Stk_bar = stats::weighted.mean(Stk, w = A, na.rm = TRUE),
    Pe_bar  = stats::weighted.mean(Pe,  w = A, na.rm = TRUE),
    .groups = "drop"
  )

## ---- deltas: Cluster 12 vs VC (1.0) at dp = 1 µm (non-LLL lobes)

library(dplyr)
# Cluster 12 vs VC (1.0), dp = 1 µm, non-LLL lobes
delta_phys <- phys_tbl %>%
  dplyr::filter(dp == 1.00,
         lobe != "LLL",
         condition %in% c("VC (1.0)", "Cluster 12")) %>%
  group_by(lobe) %>%
  summarise(
    Stk_ref = Stk_bar[condition == "VC (1.0)"][1],
    Stk_c12 = Stk_bar[condition == "Cluster 12"][1],
    Pe_ref  = Pe_bar [condition == "VC (1.0)"][1],
    Pe_c12  = Pe_bar [condition == "Cluster 12"][1],
    pct_Stk = 100 * (Stk_c12 / Stk_ref - 1),
    pct_Pe  = 100 * (Pe_c12 / Pe_ref  - 1),
    .groups = "drop"
  ) %>%
  arrange(lobe)

# Pretty print
delta_phys %>%
  mutate(
    pct_Stk = scales::number(pct_Stk, accuracy = 0.1, signed = TRUE, suffix = " %"),
    pct_Pe  = scales::number(pct_Pe,  accuracy = 0.1, signed = TRUE, suffix = " %")
  )

library(dplyr); library(tidyr); library(knitr); library(kableExtra)

tab_delta <- delta_flow %>%
  select(lobe, pct_U) %>%
  rename(dU = pct_U) %>%
  left_join(delta_phys %>% select(lobe, pct_Stk, pct_Pe), by="lobe") %>%
  mutate(across(c(dU, pct_Stk, pct_Pe),
                ~ scales::number(.x, accuracy=0.1, signed=TRUE, suffix=" %"))) %>%
  rename(`ΔU`=dU, `ΔStk`=pct_Stk, `ΔPe`=pct_Pe)

kable(tab_delta, "latex", booktabs=TRUE, linesep="") %>%
  kable_styling(latex_options=c("hold_position","striped")) %>%
  add_header_above(c(" " = 1, "Clustered vs VC-Ref (dp = 1 µm)" = 3))

library(ggplot2)
tab_delta_long <- delta_phys %>%
  select(lobe, pct_Stk, pct_Pe) %>%
  left_join(delta_flow %>% select(lobe, pct_U), by="lobe") %>%
  pivot_longer(-lobe, names_to="metric", values_to="pct")

ggplot(tab_delta_long, aes(lobe, pct, fill=metric)) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_col(position="dodge", width=0.7) +
  scale_fill_viridis_d(end=0.9, labels=c(pct_U="ΔU", pct_Stk="ΔStk", pct_Pe="ΔPe")) +
  labs(title="Clustered vs VChild-Ref at dp = 1 µm", y="Percent change", x="Lobe", fill=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position="bottom", axis.text.x=element_text(angle=40, hjust=1))



# --- 2) WHERE DOES SED ↑ ? (conducting vs alveolar) -----------------------
sed_region <- all_data %>%
  mutate(region = region_of_gen(as.integer(gen))) %>%
  group_by(condition, dp, lobe, region) %>%
  summarise(sed = sum(nj_loss_sed, na.rm=TRUE),
            imp = sum(nj_loss_imp, na.rm=TRUE),
            dif = sum(nj_loss_dif, na.rm=TRUE),
            .groups="drop") %>%
  group_by(condition, dp, lobe) %>%
  mutate(sed_frac_region = sed/sum(sed)) %>%  # share of lobe's sedimentation by region
  ungroup()

# focus on un-obstructed lobes at dp=1.0
sed_region %>%
  dplyr::filter(dp == 1.00, lobe != "LLL", condition %in% c("VC (1.0)","Cluster 12")) %>%
  select(condition, lobe, region, sed_frac_region) %>%
  arrange(lobe, region, condition)

# quick plot

plot_df <- sed_region |>
  dplyr::filter(dp == 1.00, lobe != "LLL",
                condition %in% c("VC (1.0)", "Cluster 12")) |>
  dplyr::mutate(
    condition = factor(condition,
                       levels = c("Cluster 12", "VC (1.0)"),
                       labels = c("LLL-target", "VChild-Ref")),
    # optional: control region order
    region = factor(region, levels = c("proximal", "mid", "distal"))
  )

ggplot(plot_df, aes(region, sed, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ lobe, nrow = 1) +
  scale_y_continuous(labels = scales::label_number(scale = 1e-6, suffix = " mg", accuracy = 0.1)) +
  scale_fill_viridis_d(end = 0.9, name = NULL) +
  labs(title = "Sedimentation by region (dp = 1 µm)",
       y = "Sedimented mass", x = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")


# --- 3) STOKES & PECLET PROXIES (why imp ↓, sed ↑) -----------------------
dp_list <- c(1.0, 3.0, 5.0) # µm
phys_tbl <- all_data %>%
  dplyr::filter(dp %in% dp_list) %>%
  mutate(r_m = as.numeric(ne_radius)*mm,
         area = pi*r_m^2,
         U    = as.numeric(ne_Vdot)/area,
         dp_m = dp*1e-6,
         Stk  = Stk(dp_m, U, r_m),
         Pe   = Pe(U, r_m, dp_m)) %>%
  group_by(condition, dp, lobe) %>%
  summarise(Stk_bar = weighted.mean(Stk, w=area, na.rm=TRUE),
            Pe_bar  = weighted.mean(Pe,  w=area, na.rm=TRUE),
            .groups="drop")

phys_tbl %>%
  dplyr::filter(lobe!="LLL", condition %in% c("VC (1.0)","Cluster 12")) %>%
  arrange(dp, lobe, condition)


# --- 4) 0.01 µm + sex: why women ↑ BD, men ↓ BD when flow ↓? -------------
# Build subject-level table: pick your adult datasets & flow indicators.
# Compute bronchi vs alveoli fractions by sex across flow settings.

# Example skeleton (adapt column names if needed):

# Read in adult data
base_dir <- "/hpc/gjin250/refine_1d/particle_development/output"
#/hpc/gjin250/refine_1d/particle_development/output/H12816/0.01
# sbj = H12816, H653, dp = 0.01, 0.1,1.0,10.0

conds <- tibble::tibble(
  condition = c("VC (1.0)", "VC (0.90)", "VC (0.85)", "VC (0.80)", "Cluster 12"),
  subdir    = c("025-virtual","025-constrict-0.9","025-constrict-0.85","025-constrict","025-cluster-12")
)
dp_vals <- c(0.80, 1.00, 3.00, 5.00)
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

# summarise to lobe totals per condition & dp
lobe_summary <- all_data %>%
  group_by(condition, dp, lobe) %>%
  summarise(
    total_sed = sum(nj_loss_sed, na.rm = TRUE),
    total_imp = sum(nj_loss_imp, na.rm = TRUE),
    total_dif = sum(nj_loss_dif, na.rm = TRUE),
    total_loss_sum = sum(total_loss, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  mutate(
    sed_percent = 100 * total_sed / sum(total_loss_sum),
    imp_percent = 100 * total_imp / sum(total_loss_sum),
    dif_percent = 100 * total_dif / sum(total_loss_sum),
    lobe_percent = 100 * total_loss_sum / sum(total_loss_sum)
  ) %>%
  ungroup()

# ---------- 3) tidy for stacked bars ----------
plot_df <- lobe_summary %>%
  select(condition, dp, lobe, total_sed, total_imp, total_dif) %>%
  pivot_longer(c(total_sed, total_imp, total_dif),
               names_to = "mechanism", values_to = "mass") %>%
  mutate(
    mechanism = factor(mechanism, levels = c("total_sed","total_imp","total_dif"),
                       labels = c("Sedimentation","Impaction","Diffusion")),
    condition = factor(condition, levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)")),
    dp = factor(sprintf("%.2f", dp), levels = sprintf("%.2f", sort(dp_vals))),
    # order lobes nicely if you have the standard set:
    lobe = factor(lobe, levels = c("Trach","RUL","RML","RLL","LUL","LLL"))
  )




bd_ad_by_sex <- your_adult_dataframe %>%   # bind_rows of male/female runs
  filter(dp == 0.01) %>%
  mutate(region = region_of_gen(as.integer(gen))) %>%
  group_by(Sex, FlowSetting, region) %>%   # FlowSetting: however you encode low/med/high
  summarise(dif = sum(nj_loss_dif, na.rm=TRUE),
            imp = sum(nj_loss_imp, na.rm=TRUE),
            sed = sum(nj_loss_sed, na.rm=TRUE),
            .groups="drop") %>%
  mutate(bronch = dif+imp+sed,
         zone   = ifelse(region=="distal","alveolar","bronchial")) %>%
  group_by(Sex, FlowSetting, zone) %>%
  summarise(dep = sum(bronch), .groups="drop") %>%
  group_by(Sex, FlowSetting) %>%
  mutate(frac = dep/sum(dep)) %>%
  ungroup()

# Plot BD vs AD by sex and flow
ggplot(bd_ad_by_sex, aes(FlowSetting, frac, fill=zone)) +
  geom_col(position="fill") +
  facet_wrap(~ Sex) +
  scale_y_continuous(labels=percent) +
  scale_fill_viridis_d(end=0.9) +
  labs(title="0.01 µm: bronchial vs alveolar deposition by sex across flow",
       y="Fraction of total", x="Flow setting", fill=NULL) +
  theme_classic(base_size=12) + theme(legend.position="bottom")


