# old_adult and new_adult: bound data.frames with columns:
# sbj, Sex, FlowSetting, dp, nj_loss_dif, nj_loss_sed, nj_loss_imp, (etc.)
sum_bronch <- function(df, sbj_id, dp_sel = 0.01) {
  df |>
    dplyr::filter(sbj == sbj_id, dp == dp_sel) |>
    dplyr::group_by(Sex, FlowSetting) |>
    dplyr::summarise(
      bronch_dif = sum(nj_loss_dif, na.rm = TRUE),
      bronch_imp = sum(nj_loss_imp, na.rm = TRUE),
      bronch_sed = sum(nj_loss_sed, na.rm = TRUE),
      .groups = "drop"
    )
}
br_old <- sum_bronch(old_adult, "H653", 0.01)
br_new <- sum_bronch(new_adult, "H653", 0.01)
dplyr::left_join(br_old, br_new, by = c("Sex","FlowSetting"),
                 suffix = c("_old","_new")) |>
  dplyr::mutate(pct_delta_dif = 100*(bronch_dif_new/bronch_dif_old - 1))


# term_old and term_new: bound terminal tables read with your read_terminal()
sum_alv <- function(df, sbj_id, flow) {
  df |>
    dplyr::filter(sbj == sbj_id, flow == !!flow) |>
    dplyr::summarise(alv_mass = sum(nu_loss_dif + nu_loss_sed, na.rm = TRUE),
                     .groups = "drop")
}

mk_zone_fracs <- function(bronch_tbl, alv_tbl) {
  dplyr::left_join(bronch_tbl, alv_tbl, by = "FlowSetting") |>
    dplyr::transmute(FlowSetting, bronch = bronch_dif + bronch_sed + bronch_imp,
                     alv = alv_mass) |>
    tidyr::pivot_longer(c(bronch, alv), names_to = "zone", values_to = "mass") |>
    dplyr::group_by(FlowSetting) |>
    dplyr::mutate(frac = mass/sum(mass)) |>
    dplyr::ungroup()
}
fr_old <- mk_zone_fracs(br_old, alv_old_H653)  # make alv_old_H653 from terminal
fr_new <- mk_zone_fracs(br_new, alv_new_H653)
dplyr::left_join(fr_old, fr_new, by = c("FlowSetting","zone"),
                 suffix = c("_old","_new")) |>
  dplyr::mutate(delta_pp = 100*(frac_new - frac_old))
