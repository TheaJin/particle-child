library(dplyr)
library(purrr)
library(tidyr)



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

# 1) reader (as you wrote)
term_colnames <- c(
  "np","x","y","z","distance","lobe",
  "nu_vol","nu_vdot0",
  "nu_mass","nu_loss_dif","nu_loss_sed","ne_part_vel"
)
read_terminal <- function(path) {
  df <- utils::read.table(path, header = FALSE, sep = "", fill = FALSE, comment.char = "")
  colnames(df) <- term_colnames
  num_cols <- c("nu_vol","nu_vdot0","nu_mass","nu_loss_dif","nu_loss_sed","ne_part_vel")
  df[num_cols] <- lapply(df[num_cols], function(x) as.numeric(gsub("D","E", x)))
  df %>% mutate(np = as.integer(np), lobe = trimws(lobe))
}

# 2) read ALL terminal files (map row-by-row)
term_df <- pmap_dfr(path_to_terminal_txt,
                    function(sbj, sex, flow, dp, path, exists) {
                      df <- read_terminal(path)
                      df %>% mutate(sbj = sbj, Sex = sex, FlowSetting = flow, dp = dp)
                    }
)

# 3) alveolar (terminal) deposition per run
alv_sum <- term_df %>%
  group_by(sbj, Sex, FlowSetting, dp) %>%
  summarise(alv_mass = sum(nu_loss_dif + nu_loss_sed, na.rm = TRUE), .groups = "drop")

# 4) bronchial (airway_result) deposition per run
# assumes `adult_all` already contains these columns: sbj, Sex, FlowSetting, dp
bronch_sum <- adult_all %>%
  group_by(sbj, Sex, flow, dp) %>%
  summarise(bronch_mass = sum(nj_loss_dif + nj_loss_imp + nj_loss_sed, na.rm = TRUE), .groups = "drop")

# 5) combine to get fractions bronchial vs alveolar
dep_zone <- bronch_sum %>%
  dplyr::rename(FlowSetting = flow) %>%              # bronch_sum had `flow`
  dplyr::full_join(
    alv_sum,   # alv_sum had `flowsetting`
    by = c("sbj","Sex","FlowSetting","dp")
  ) %>%
  dplyr::mutate(bronch_mass = coalesce(bronch_mass, 0),
                alv_mass    = coalesce(alv_mass,    0)) %>%
  tidyr::pivot_longer(c(bronch_mass, alv_mass),
                      names_to = "zone", values_to = "mass") %>%
  dplyr::mutate(zone = dplyr::recode(zone,
                                     bronch_mass = "bronchial",
                                     alv_mass    = "alveolar")) %>%
  dplyr::group_by(sbj, Sex, FlowSetting, dp) %>%
  dplyr::mutate(frac = mass / sum(mass)) %>%
  dplyr::ungroup()

ggplot(dep_zone %>% filter(dp == 0.01),
       aes(FlowSetting, frac, fill = zone)) +
  geom_col(position = "fill") +
  facet_wrap(~ Sex) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = "0.01 µm: bronchial vs alveolar deposition by sex across flow (with terminal files)",
       y = "Fraction of total", x = "Flow setting", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")



# 1) normalize keys in each source -----------------------------
norm_sex  <- function(x) dplyr::recode(x,
                                       "F"="Female","female"="Female","f"="Female",
                                       "M"="Male","male"="Male","m"="Male", .default = x
)
norm_flow <- function(x) stringr::str_to_title(x)  # "low"->"Low", etc.

bronch_sum <- bronch_sum %>%
  dplyr::mutate(
    Sex = norm_sex(Sex),
    FlowSetting = norm_flow(flow)
  )

alv_sum <- alv_sum %>%
  dplyr::mutate(
    Sex = norm_sex(Sex),
    FlowSetting = norm_flow(FlowSetting)
  )

# 2) join and compute fractions -------------------------------
dep_zone <- bronch_sum %>%
  dplyr::full_join(alv_sum,
                   by = c("sbj","Sex","FlowSetting","dp")
  ) %>%
  dplyr::mutate(
    bronch_mass = dplyr::coalesce(bronch_mass, 0),
    alv_mass    = dplyr::coalesce(alv_mass, 0)
  ) %>%
  tidyr::pivot_longer(c(bronch_mass, alv_mass),
                      names_to = "zone", values_to = "mass") %>%
  dplyr::mutate(zone = dplyr::recode(zone,
                                     bronch_mass = "bronchial",
                                     alv_mass    = "alveolar"),
                Sex = factor(Sex, levels = c("Female","Male"))) %>%
  dplyr::group_by(sbj, Sex, FlowSetting, dp) %>%
  dplyr::mutate(frac = mass / sum(mass)) %>%
  dplyr::ungroup()

# (optional sanity check: should be exactly 2 rows per group)
dep_zone %>% dplyr::count(sbj, Sex, FlowSetting, dp)
# n should be 2 everywhere

# 3) plot (use mass + position='fill' to avoid double-normalizing)
ggplot(dep_zone, aes(FlowSetting, mass, fill = zone)) +
  geom_col(position = "fill") +
  facet_wrap(~ Sex) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_viridis_d(end = 0.9) +
  labs(title = "0.01 µm: bronchial vs alveolar deposition by sex across flow",
       y = "Fraction of total", x = "Flow setting", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

# FRACTIONS (percentage-point change)
# Fractions of total (by Sex × FlowSetting), plus Δpp = Low − High (percentage points)
frac_tbl <- dep_zone %>%
  dplyr::filter(zone %in% c("bronchial","alveolar")) %>%
  dplyr::group_by(Sex, FlowSetting, zone) %>%
  dplyr::summarise(mass = sum(mass), .groups = "drop") %>%     # aggregate mass
  dplyr::group_by(Sex, FlowSetting) %>%
  dplyr::mutate(frac = mass / sum(mass)) %>%                    # fraction within Sex×FlowSetting
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from = FlowSetting, values_from = frac) %>%
  dplyr::mutate(delta_pp = 100 * (Low - High)) %>%
  dplyr::arrange(Sex, zone)

# Pretty print
frac_tbl %>%
  dplyr::mutate(
    High    = scales::percent(High, accuracy = 0.1),
    Low     = scales::percent(Low,  accuracy = 0.1),
    delta_pp = scales::number(delta_pp, accuracy = 0.1, signed = TRUE)
  )


# ABSOLUTE MASSES (percent change High→Low)
abs_tbl <- dep_zone %>%
  group_by(Sex, FlowSetting, zone) %>%
  summarise(mass = sum(mass), .groups="drop") %>%
  tidyr::pivot_wider(names_from = FlowSetting, values_from = mass) %>%
  mutate(delta_pct = 100*(Low/High - 1)) %>%
  arrange(Sex, zone)

abs_tbl

# pool masses across subjects within Sex×Flow×zone, then make fractions
frac_tbl <- dep_zone %>%
  dplyr::filter(zone %in% c("bronchial","alveolar")) %>%
  dplyr::group_by(Sex, FlowSetting, zone) %>%
  dplyr::summarise(mass = sum(mass), .groups = "drop") %>%
  dplyr::group_by(Sex, FlowSetting) %>%
  dplyr::mutate(frac = mass / sum(mass)) %>%
  dplyr::ungroup() %>%
  dplyr::select(Sex, FlowSetting, zone, frac) %>%        # drop mass before widening
  tidyr::pivot_wider(
    id_cols   = c(Sex, zone),                            # <- important!
    names_from  = FlowSetting,
    values_from = frac
  ) %>%
  dplyr::mutate(delta_pp = 100 * (Low - High)) %>%
  dplyr::arrange(Sex, zone)

# pretty view
frac_tbl %>%
  dplyr::mutate(
    High    = scales::percent(High, accuracy = 0.1),
    Low     = scales::percent(Low,  accuracy = 0.1),
    delta_pp = scales::number(delta_pp, accuracy = 0.1, signed = TRUE)
  )


# B) Compute per-subject fractions, then average (nice if n>1 per sex)
by_sbj <- dep_zone %>%
  dplyr::filter(zone %in% c("bronchial","alveolar")) %>%
  dplyr::group_by(sbj, Sex, FlowSetting, zone) %>%
  dplyr::summarise(mass = sum(mass), .groups = "drop") %>%
  dplyr::group_by(sbj, Sex, FlowSetting) %>%
  dplyr::mutate(frac = mass / sum(mass)) %>%
  dplyr::ungroup()

frac_mean <- by_sbj %>%
  dplyr::group_by(Sex, zone, FlowSetting) %>%
  dplyr::summarise(mean_frac = mean(frac), sd = sd(frac), n = dplyr::n(),
                   se = sd/sqrt(n), .groups = "drop") %>%
  tidyr::pivot_wider(
    id_cols   = c(Sex, zone),
    names_from  = FlowSetting,
    values_from = mean_frac
  ) %>%
  dplyr::mutate(delta_pp = 100 * (Low - High))

frac_mean





