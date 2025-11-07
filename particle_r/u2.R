library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)

# --- constants ---
base_dir <- "/hpc/gjin250/refine_1d/particle_development/output"
conds <- tibble::tibble(
  condition = c("VC (1.0)", "VC (0.90)", "VC (0.85)", "VC (0.80)", "Cluster 12"),
  subdir    = c("025-virtual","025-constrict-0.9","025-constrict-0.85","025-constrict","025-cluster-12")
)
dp_vals <- c(0.80, 1.00, 3.00, 5.00)
column_names <- c(
  "ne","x","y","z","distance","gen","horsf_ord","lobe","ne_length","ne_radius",
  "ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel"
)

# --- helper: read a result file safely ---
read_airway <- function(path) {
  df <- utils::read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
  colnames(df) <- column_names
  num_cols <- c("ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel")
  df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(gsub("D","E", x)))
  df %>% mutate(ne = as.integer(ne))
}

# --- 1) build lobe map from your suggested file ---
lobe_map_path <- file.path(base_dir, "025-constrict", "5.00", "airway_result.txt")
stopifnot(file.exists(lobe_map_path))

lobe_map_df <- read_airway(lobe_map_path) %>%
  transmute(ne, lobe_ref = as.character(lobe))

# quick sanity check on uniqueness
stopifnot(!any(duplicated(lobe_map_df$ne)))


# --- 2) loader that applies the lobe map by joining on `ne` ---

load_one <- function(path, condition, dp, lobe_map_df) {
  read_airway(path) %>%
    dplyr::left_join(lobe_map_df, by = "ne") %>%
    dplyr::mutate(
      # prefer mapped lobe; only fall back to file lobe if map missing
      lobe = dplyr::coalesce(lobe_ref, as.character(lobe)),
      # now clean up any unknowns
      lobe = dplyr::case_when(
        lobe %in% c("Unkno","Unknown","UNK", NA_character_) ~ "Trach",
        TRUE ~ lobe
      ),
      mapped_ne = !is.na(lobe_ref),
      condition = condition,
      dp = dp,
      total_loss = nj_loss_dif + nj_loss_imp + nj_loss_sed
    ) %>%
    dplyr::select(-lobe_ref)
}



# --- 3) read ALL files with the lobe map applied ---
all_data <- purrr::map_dfr(dp_vals, function(dp) {
  dp_str <- sprintf("%.2f", dp)
  purrr::pmap_dfr(conds, function(condition, subdir) {
    path <- file.path(base_dir, subdir, dp_str, "airway_result.txt")
    if (!file.exists(path)) stop("Missing file: ", path)
    load_one(path, condition, dp, lobe_map_df)
  })
})

# after building all_data:
coverage <- all_data %>%
  dplyr::group_by(condition, dp) %>%
  dplyr::summarise(n = dplyr::n(), matched = sum(mapped_ne), pct = 100*matched/n, .groups="drop")
print(coverage)


all_data %>% dplyr::filter(ne == 2380) %>% dplyr::count(condition, dp, lobe)


# all_data %>% dplyr::filter(is.na(lobe)) %>% dplyr::count()  # should be 0
# 
# # totals per panel (mass should equal sum of mechanisms)
# panel_totals <- all_data %>%
#   dplyr::group_by(condition, dp) %>%
#   dplyr::summarise(
#     total = sum(nj_loss_dif + nj_loss_imp + nj_loss_sed, na.rm=TRUE),
#     dif   = sum(nj_loss_dif, na.rm=TRUE),
#     imp   = sum(nj_loss_imp, na.rm=TRUE),
#     sed   = sum(nj_loss_sed, na.rm=TRUE),
#     .groups="drop"
#   )
# print(panel_totals, n=Inf)

ggplot(plot_df, aes(x = lobe, y = mass, fill = mechanism)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(dp), cols = vars(condition), switch = "y",
             scales = "free_y") +  # <- remove this if you prefer a common y
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(title = "Deposition by Lobe and Condition",
       x = "Lobe (Region)", y = "Total Deposition (mass)", fill = "Mechanism") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95", colour = NA))


# --- 4) summarise by lobe/condition/dp ---
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

# --- 5) long format for stacked bars; tidy labels & order ---
plot_df <- lobe_summary %>%
  select(condition, dp, lobe, total_sed, total_imp, total_dif) %>%
  pivot_longer(c(total_sed, total_imp, total_dif),
               names_to = "mechanism", values_to = "mass") %>%
  mutate(
    mechanism = factor(mechanism,
                       levels = c("total_sed","total_imp","total_dif"),
                       labels = c("Sedimentation","Impaction","Diffusion")),
    condition = factor(condition, levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)")),
    dp = factor(sprintf("%.2f", dp), levels = sprintf("%.2f", sort(dp_vals))),
    lobe = factor(lobe, levels = c("Trach","RUL","RML","RLL","LUL","LLL"))
  )

# --- 6) plot (viridis) ---
p <- ggplot(plot_df, aes(x = lobe, y = mass, fill = mechanism)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(dp), cols = vars(condition), switch = "y", scales = "free_y") +
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  labs(
    title = "Deposition by Lobe and Condition",
    x = "Lobe (Region)", y = "Total Deposition (mass)", fill = "Mechanism"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.spacing.x = unit(8, "pt"),
    panel.spacing.y = unit(8, "pt")
  )

print(p)

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(scales)

plot_pct <- lobe_summary %>%
  select(condition, dp, lobe, sed_percent, imp_percent, dif_percent) %>%
  pivot_longer(c(sed_percent, imp_percent, dif_percent),
               names_to = "mechanism", values_to = "pct") %>%
  mutate(
    mechanism = factor(mechanism,
                       levels = c("sed_percent","imp_percent","dif_percent"),
                       labels = c("Sedimentation","Impaction","Diffusion")),
    condition = factor(condition, levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)")),
    dp = factor(sprintf("%.2f", dp), levels = sprintf("%.2f", sort(unique(dp)))),
    lobe = factor(lobe, levels = c("Trach","RUL","RML","RLL","LUL","LLL"))
  )

ggplot(plot_pct, aes(lobe, pct/100, fill = mechanism)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(dp), cols = vars(condition), switch = "y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = "Deposition Fraction by Lobe and Condition",
       x = "Lobe (Region)", y = "Share of Total Deposition",
       fill = "Mechanism") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95", colour = NA))




plot_lobe_share <- lobe_summary %>%
  mutate(lobe_share = lobe_percent/100) %>%
  select(condition, dp, lobe, lobe_share) %>%
  mutate(
    condition = factor(condition, levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)")),
    dp = factor(sprintf("%.2f", dp), levels = sprintf("%.2f", sort(unique(dp)))),
    lobe = factor(lobe, levels = c("Trach","RUL","RML","RLL","LUL","LLL"))
  )

ggplot(plot_lobe_share, aes(lobe, lobe_share, fill = condition)) +
  geom_col(position = "dodge", width = 0.75) +
  facet_grid(rows = vars(dp), switch = "y") +
  scale_fill_viridis_d(option = "D", end = 0.95) +
  scale_y_continuous(labels = label_percent()) +
  labs(title = "Lobe Share of Total Deposition",
       x = "Lobe (Region)", y = "Fraction of Total Deposition", fill = "Condition") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95", colour = NA))




ref <- lobe_summary %>%
  dplyr::filter(condition == "VC (1.0)") %>%
  dplyr::select(dp, lobe, ref_share = lobe_percent)

delta_tbl <- lobe_summary %>%
  dplyr::select(condition, dp, lobe, lobe_percent) %>%
  dplyr::left_join(ref, by = c("dp","lobe")) %>%
  dplyr::mutate(delta_pp = lobe_percent - ref_share) %>%
  dplyr::arrange(dp, lobe, condition)


# View a compact summary at dp = 1.00, for example:
delta_tbl %>%
  dplyr::filter(dp == 1.00) %>%
  dplyr::select(condition, lobe, delta_pp)


library(dplyr); library(tidyr); library(ggplot2); library(viridis); library(scales)

# use your lobe_summary already built
plot_pct <- lobe_summary %>%
  dplyr::select(condition, dp, lobe, sed_percent, imp_percent, dif_percent) %>%
  tidyr::pivot_longer(c(sed_percent, imp_percent, dif_percent),
                      names_to = "mechanism", values_to = "pct") %>%
  dplyr::mutate(
    mechanism = factor(mechanism,
                       levels = c("sed_percent","imp_percent","dif_percent"),
                       labels = c("Sedimentation","Impaction","Diffusion")),
    condition = factor(condition,
                       levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)"),
                       labels = c("LLL-target","VChild-0.80","VChild-0.85","Vchild-0.90","VChild-Ref")),
    dp = factor(sprintf("%.2f", dp), levels = sprintf("%.2f", sort(unique(dp)))),
    lobe = factor(lobe, levels = c("Trach","RUL","RML","RLL","LUL","LLL"))
  )

p_frac <- ggplot(plot_pct, aes(lobe, pct/100, fill = mechanism)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(dp), cols = vars(condition), switch = "y") +
  scale_fill_viridis_d(option = "D") +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(title = "Deposition Fraction by Lobe and Condition",
       x = "Lobe (Region)", y = "Share of Total Deposition", fill = "Mechanism") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95", colour = NA),
        legend.position = "bottom",    # put legend under the plot
        legend.box = "horizontal"      # keep items in a row
        )
p_frac




ref <- lobe_summary %>%
  dplyr::filter(condition == "VC (1.0)") %>%
  dplyr::select(dp, lobe, ref_share = lobe_percent)

delta_tbl <- lobe_summary %>%
  dplyr::select(condition, dp, lobe, lobe_percent) %>%
  dplyr::left_join(ref, by = c("dp","lobe")) %>%
  dplyr::mutate(delta_pp = lobe_percent - ref_share) %>%
  dplyr::arrange(dp, lobe, condition)

# Pretty view at dp = 1 μm
delta_tbl %>%
  dplyr::filter(dp == 1.00) %>%
  dplyr::select(condition, lobe, delta_pp) %>%
  dplyr::mutate(delta_pp = scales::number(delta_pp, accuracy = 0.1))

delta_tbl |>
  dplyr::filter(dp == 1.00, condition == "Cluster 12", lobe %in% c("LLL","LUL")) |>
  dplyr::select(lobe, delta_pp)

library(dplyr); library(tidyr); library(knitr); library(kableExtra)

delta_1um_wide <- delta_tbl %>%
  dplyr::filter(dp == 1.00) %>%
  dplyr::mutate(delta_pp = scales::number(delta_pp, accuracy = 0.1, trim = TRUE),
                delta_pp = ifelse(substr(delta_pp,1,1) == "-", delta_pp, paste0("+", delta_pp))) %>%
  dplyr::select(condition, lobe, delta_pp) %>%
  tidyr::pivot_wider(names_from = lobe, values_from = delta_pp) %>%
  dplyr::arrange(factor(condition, levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)")))

kable(delta_1um_wide, "latex", booktabs = TRUE, linesep = "") %>%
  kable_styling(latex_options = c("hold_position","striped")) %>%
  add_header_above(c(" " = 1, "Δ share vs VC (1.0), pp" = ncol(delta_1um_wide)-1)) %>%
  column_spec(2:ncol(delta_1um_wide), width = "10mm")

library(ggplot2); library(dplyr); library(viridis); library(scales)

# overall lobe fraction at dp=1 (not split by mechanism)
frac_1um <- lobe_summary %>%
  dplyr::filter(dp == 1.00) %>%
  dplyr::transmute(condition, lobe, share = lobe_percent/100)

# reference shares (VC 1.0)
ref_1um <- frac_1um %>%
  dplyr::filter(condition == "VC (1.0)") %>%
  dplyr::select(lobe, ref = share)

# deltas for labels on Cluster 12 bars
lab_1um <- frac_1um %>%
  dplyr::filter(condition == "Cluster 12") %>%
  dplyr::left_join(ref_1um, by = "lobe") %>%
  dplyr::mutate(delta = share - ref,
                label = scales::label_number(accuracy = 0.1, signed = TRUE)(delta*100),
                y = pmax(share, ref) + 0.01)   # place a bit above the taller of the two

# order for display
frac_1um <- frac_1um %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Cluster 12","VC (0.80)","VC (0.85)","VC (0.90)","VC (1.0)")),
    lobe = factor(lobe, levels = c("Trach","RUL","RML","RLL","LUL","LLL"))
  )
lab_1um$lobe <- frac_1um$lobe[match(lab_1um$lobe, frac_1um$lobe)]

p_1um <- ggplot(frac_1um, aes(lobe, share, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_viridis_d(option = "D", end = 0.95) +
  scale_y_continuous(labels = label_percent()) +
  labs(title = "Lobe share of total deposition (dp = 1 µm)",
       x = "Lobe (Region)", y = "Fraction of total deposition", fill = "Condition") +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        axis.text.x = element_text(angle = 50, hjust = 1))

# add Δ labels for Cluster 12 bars
p_1um +
  geom_text(data = lab_1um,
            aes(x = lobe, y = y, label = paste0(label, " pp")),
            color = "gray20", size = 3.5, fontface = 2, vjust = 0) +
  coord_cartesian(ylim = c(0, max(frac_1um$share) + 0.08))


# ggsave("fig_5_16_fraction.png", p_frac, width = 12, height = 10, dpi = 300)

mk_sentence <- function(dp_num = 1.00) {
  d <- delta_tbl |> dplyr::filter(dp == dp_num, condition == "Cluster 12")
  lll <- d |> dplyr::filter(lobe == "LLL") |> dplyr::pull(delta_pp)
  lul <- d |> dplyr::filter(lobe == "LUL") |> dplyr::pull(delta_pp)
  rll <- d |> dplyr::filter(lobe == "RLL") |> dplyr::pull(delta_pp)
  glue::glue("At $d_p={dp_num}\\,\\mu\\text{{m}}$, the LLL share decreases by {scales::number(lll, accuracy=0.1)} pp, while LUL and RLL increase by {scales::number(lul, accuracy=0.1)} and {scales::number(rll, accuracy=0.1)} pp, respectively.")
}

cat(mk_sentence(1.00))

delta_1um_wide <- delta_tbl %>%
  dplyr::filter(dp == 1.00) %>%
  dplyr::mutate(sign = ifelse(delta_pp >= 0, "\u2191", "\u2193"),
                delta_txt = paste0(sign, " ", scales::number(delta_pp, accuracy=0.1), " pp")) %>%
  dplyr::select(condition, lobe, delta_txt) %>%
  tidyr::pivot_wider(names_from = lobe, values_from = delta_txt)



# ggsave("deposition_lobe_condition_dp.png", p, width = 12, height = 10, dpi = 300)
