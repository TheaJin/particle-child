## deps
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)

# ---------- 0) describe your layout once ----------
base_dir <- "/hpc/gjin250/refine_1d/particle_development/output"

conds <- tibble::tibble(
  condition = c("VC (1.0)", "VC (0.90)", "VC (0.85)", "VC (0.80)", "Cluster 12"),
  subdir    = c("025-virtual","025-constrict-0.9","025-constrict-0.85","025-constrict","025-cluster-12")
)

dp_vals <- c(0.80, 1.00, 3.00, 5.00)

# these are the raw column names in airway_result.txt
column_names <- c(
  "ne","x","y","z","distance","gen","horsf_ord","lobe","ne_length","ne_radius",
  "ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel"
)

# ---------- 1) helper to load one file ----------
load_one <- function(path, condition, dp, lobe_map = NULL) {
  df <- utils::read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
  colnames(df) <- column_names
  
  # numeric columns sometimes with Fortran "D" exponents -> swap to "E"
  num_cols <- c("ne_Vdot","ne_mass","nj_loss_dif","nj_loss_imp","nj_loss_sed","nj_conc1","ne_part_vel")
  df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(gsub("D","E", x)))
  
  # apply external lobe mapping if supplied; else just fix "Unkno" -> "Trach"
  if (!is.null(lobe_map)) df$lobe <- lobe_map
  df <- df %>% mutate(lobe = ifelse(lobe %in% c("Unkno","Unknown","UNK","NA"), "Trach", as.character(lobe)))
  
  # (optional) standardise lobe labels
  df <- df %>%
    mutate(lobe = dplyr::recode(lobe,
                                "RUL"="RUL","RML"="RML","RLL"="RLL","LUL"="LUL","LLL"="LLL",
                                "Trach"="Trach", .default = lobe
    ))
  
  df %>%
    mutate(
      condition = condition,
      dp = dp,
      total_loss = nj_loss_dif + nj_loss_imp + nj_loss_sed
    )
}

# ---------- 2) build the full file list, read & summarise ----------
# if you *do* have a reliable lobe vector from one file, set lobe_map <- that vector; else leave NULL
lobe_map <- NULL

all_data <- purrr::map_dfr(dp_vals, function(dp) {
  dp_str <- sprintf("%.2f", dp)
  purrr::pmap_dfr(conds, function(condition, subdir) {
    path <- file.path(base_dir, subdir, dp_str, "airway_result.txt")
    if (!file.exists(path)) stop("Missing file: ", path)
    load_one(path, condition, dp, lobe_map)
  })
})

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

# ---------- 4) plot (viridis) ----------
ggplot(plot_df, aes(x = lobe, y = mass, fill = mechanism)) +
  geom_col(width = 0.8) +
  facet_grid(rows = vars(dp), cols = vars(condition), switch = "y") +
  scale_fill_viridis_d(option = "D", end = 0.9) +
  labs(
    title = "Deposition by Lobe and Condition",
    x = "Lobe (Region)",
    y = "Total Deposition (mass)",
    fill = "Mechanism"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = NA),
    panel.spacing.x = unit(8, "pt"),
    panel.spacing.y = unit(8, "pt")
  )
