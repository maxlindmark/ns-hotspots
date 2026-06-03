prepare_species_data <- function(file, threshold, depth_mean, depth_sd, ...) {
  load(paste0(home, "/data/datras_data_cod_had_whg_pok/", file, ".RData"))
  
  # Estimate length-weight parameters from CA data by year and survey
  # FIXME: I suspect lots of errors here, some estimates are wild.
  lw <- dat$CA |>
    drop_na(IndWgt, LngtClas) |>
    filter(IndWgt > 0, LngtClas > 0) |>
    slice(rep(1:n(), times = NoAtALK)) |>
    mutate(length_cm = ifelse(LngtCode == ".", LngtClas / 10, LngtClas)) |>
    group_by(Year, Survey) |>
    group_modify(~ broom::tidy(lm(log(IndWgt) ~ log(length_cm), data = .x))) |>
    ungroup() |>
    dplyr::select(Year, Survey, term, estimate) |>
    pivot_wider(names_from = term, values_from = estimate) |>
    mutate(a = exp(`(Intercept)`), b = `log(length_cm)`) |>
    dplyr::select(Year, Survey, a, b) |>
    janitor::clean_names() |>
    mutate(
      b = if_else(b < 2.8 | b > 3.2, NA, b),
      a = if_else(a < quantile(a, 0.05) | a > quantile(a, 0.95), NA, a)
    )
  
  a_fallback <- median(lw$a, na.rm = TRUE)
  b_fallback <- median(lw$b, na.rm = TRUE)
  
  by_size <- dat$HL |>
    left_join(lw, by = c("Year" = "year", "Survey" = "survey")) |>
    mutate(
      a = replace_na(a, a_fallback),
      b = replace_na(b, b_fallback),
      size_class = if_else(LngtCm < threshold, "small", "large"),
      biomass_g = Count * a * LngtCm^b
    ) |>
    summarise(
      n = sum(Count, na.rm = TRUE),
      biomass_g = sum(biomass_g, na.rm = TRUE),
      a = first(a),
      b = first(b),
      .by = c(haul.id, size_class)
    )
  
  d <- expand_grid(
    haul.id = dat$HH$haul.id,
    size_class = c("small", "large")
  ) |>
    left_join(by_size, by = c("haul.id", "size_class")) |>
    mutate(
      n = replace_na(n, 0),
      biomass_g = replace_na(biomass_g, 0),
      a = replace_na(a, a_fallback),
      b = replace_na(b, b_fallback)
    ) |>
    left_join(
      dat$HH |> dplyr::select(haul.id, Year, Survey, Gear, lon, lat, Depth,
                              TimeShotHour, Ship, HaulDur, Month),
      by = "haul.id"
    ) |>
    mutate(
      fish_size = factor(size_class, levels = c("small", "large")),
      Year = as.factor(Year),
      Ship = as.factor(Ship),
      Gear = as.factor(Gear)
    ) |>
    dplyr::select(-size_class) |>
    janitor::clean_names() |>
    add_utm_columns(ll_names = c("lon", "lat")) |>
    filter(!is.na(haul_dur)) |>
    mutate(
      month_f = as.factor(month),
      year_f = as.factor(year),
      gear = stringr::str_to_lower(gear),
      gear = fct_collapse(gear, other = c("bt4p", "bt4s", "dht", "grt")),
      gear_f = as.factor(gear),
      depth_sc = (depth - depth_mean) / depth_sd,
      depth_sc2 = depth_sc^2,
      quarter_f = as.factor(ceiling(month / 3))
    ) |> 
    filter(quarter_f %in% c("3", "4"))
  
  d
}
