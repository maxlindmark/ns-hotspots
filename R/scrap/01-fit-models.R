library(forcats)
library(stringr)
library(viridis)
library(sf)
library(sdmTMB)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
#library(ggsidekick)
#theme_set(theme_sleek())

home <- here::here()

##### ------------ Species metadata grid
## Length-weight parameters
# TODO: estimate these every year / spatial location / life stage using built in CA data
a <- 0.01
b <- 3.0

species_grid <- tibble(
  species   = c("cod", "had", "pok", "whg"),
  file      = c("datras_data_cod", "datras_data_had", "datras_data_pok", "datras_data_whg"),
  threshold = c(35, 30, 35, 27),
  cutoff    = 25
)


##### ------------ Prepare data
prepare_species_data <- function(file, threshold, ...) {
  load(paste0(home, "/data/datras_data_cod_had_whg_pok/", file, ".RData"))
  
  # Compute depth scaling params from haul data
  depth_mean <- mean(dat$HH$Depth, na.rm = TRUE)
  depth_sd   <- sd(dat$HH$Depth,   na.rm = TRUE)
  
  by_size <- dat$HL |>
    mutate(
      size_class = if_else(LngtCm < threshold, "small", "large"),
      biomass_g  = Count * a * LngtCm^b
    ) |>
    summarise(
      n         = sum(Count,     na.rm = TRUE),
      biomass_g = sum(biomass_g, na.rm = TRUE),
      .by = c(haul.id, size_class)
    )
  
  d <- expand_grid(
    haul.id    = dat$HH$haul.id,
    size_class = c("small", "large")
  ) |>
    left_join(by_size, by = c("haul.id", "size_class")) |>
    mutate(
      n         = replace_na(n, 0),
      biomass_g = replace_na(biomass_g, 0)
    ) |>
    left_join(
      dat$HH |> dplyr::select(
        haul.id, Year, Gear, lon, lat, Depth,
        TimeShotHour, Ship, HaulDur
      ),
      by = "haul.id"
    ) |>
    mutate(
      fish_size = factor(size_class, levels = c("small", "large")),
      Year      = as.factor(Year),
      Ship      = as.factor(Ship),
      Gear      = as.factor(Gear)
    ) |>
    dplyr::select(-size_class) |>
    janitor::clean_names() |>
    add_utm_columns(ll_names = c("lon", "lat")) |>
    filter(!is.na(haul_dur)) |>
    mutate(
      gear      = stringr::str_to_lower(gear),   # <-- FIX
      gear      = fct_collapse(gear, other = c("bt4p", "bt4s", "dht", "grt")),  # <-- adjusted to match lower case
      depth_sc  = (depth - depth_mean) / depth_sd,
      depth_sc2 = depth_sc^2
    )
  
  # Return data and scaling params together
  list(data = d, depth_mean = depth_mean, depth_sd = depth_sd)
}

species_grid <- species_grid |>
  mutate(prepared   = pmap(list(file = file, threshold = threshold), prepare_species_data)) |>
  mutate(
    data       = map(prepared, "data"),
    depth_mean = map_dbl(prepared, "depth_mean"),
    depth_sd   = map_dbl(prepared, "depth_sd")
  ) |>
  dplyr::select(-prepared) |>
  mutate(
    data_small = map(data, ~ filter(.x, fish_size == "small")),
    data_large = map(data, ~ filter(.x, fish_size == "large"))
  ) |>
  pivot_longer(
    cols      = c(data_small, data_large),
    names_to  = "size_class",
    values_to = "data_filtered"
  ) |>
  mutate(
    size_class    = str_remove(size_class, "data_"),
    data_filtered = map(data_filtered, ~ {
      keep_gears <- names(which(tapply(.x$biomass_g, .x$gear, sum) > 0))
      mutate(.x, gear = fct_other(gear, keep = keep_gears),
             gear = fct_recode(gear, other = "Other"))
    })
  )


##### ------------ Make prediction grid
# Use cod as the reference for the spatial domain
cod_data <- species_grid |>
  filter(species == "cod", size_class == "small") |>
  pull(data_filtered) |>
  pluck(1)

depth_model <- mgcv::gam(
  depth ~ s(X, Y, bs = "tp", k = 200),
  data = cod_data |> distinct(X, Y, depth)
)

x <- cod_data$X
y <- cod_data$Y
z <- chull(x, y)
coords <- cbind(x[z], y[z])
coords <- rbind(coords, coords[1, ])

sp_poly <- sp::SpatialPolygons(
  list(sp::Polygons(list(sp::Polygon(coords)), ID = 1))
)
sp_poly_df <- sp::SpatialPolygonsDataFrame(
  sp_poly,
  data = data.frame(ID = 1)
)

cell_width <- 5

spatial_grid <- expand.grid(
  X = seq(min(cod_data$X), max(cod_data$X), cell_width),
  Y = seq(min(cod_data$Y), max(cod_data$Y), cell_width)
)

sp::coordinates(spatial_grid) <- c("X", "Y")
inside <- !is.na(sp::over(spatial_grid, as(sp_poly_df, "SpatialPolygons")))
spatial_grid <- as.data.frame(spatial_grid[inside, ])

spatial_grid$depth <- predict(depth_model, newdata = spatial_grid)
spatial_grid <- spatial_grid |> filter(depth > 0)

make_pred_grid <- function(dat, depth_mean, depth_sd) {
  expand_grid(
    spatial_grid,
    year = sort(unique(dat$year))
  ) |>
    mutate(
      depth_sc       = (depth - depth_mean) / depth_sd,
      depth_sc2      = depth_sc^2,
      gear           = as.factor("gov"),
      ship           = factor(dat$ship[1]),
      time_shot_hour = 12,
      haul_dur       = rep(30, n())  # <- ensures same length as rows
    )
}

species_grid <- species_grid |>
  mutate(
    pred_grid = pmap(
      list(dat = data_filtered, depth_mean = depth_mean, depth_sd = depth_sd),
      make_pred_grid
    )
  )

# world <- ne_countries(scale = "medium", returnclass = "sf")
# world_utm <- world |> st_transform(crs = 32631)
# ggplot() +
#   geom_raster(
#     data = species_grid |>
#       filter(species == "cod", size_class == "small") |>
#       pull(pred_grid) |> pluck(1) |> distinct(X, Y, depth),
#     aes(x = X * 1000, y = Y * 1000, fill = depth)) +
#   geom_sf(data = world_utm) +
#   scale_fill_viridis() +
#   coord_sf(
#     xlim = range(cod_data$X) * 1000,
#     ylim = range(cod_data$Y) * 1000,
#     crs = 32631, expand = 0
#   ) +
#   labs(x = "Easting (km)", y = "Northing (km)") +
#   theme(axis.text.x = element_text(angle = 90))


##### ------------ Fit models
# species_grid$data_filtered |> 
#   purrr::map(~summary(.x$gear))
# 
# species_grid$data_filtered |> 
#   purrr::map(~summary(.x$haul_dur))
# 
# species_grid$data_filtered |>
#   purrr::map(~sum(.x$biomass_g > 0))
# 
# species_grid$data_filtered |>
#   purrr::map(~table(.x$gear))
# 
# species_grid$data_filtered |>
#   purrr::map(~table(.x$ship))
# 
# species_grid$data_filtered |>
#   purrr::map(~sum(is.na(.x)))
# 
# species_grid$data_filtered |>
#   purrr::map(~c(
#     biomass = sum(!is.finite(.x$biomass_g)),
#     depth1  = sum(!is.finite(.x$depth_sc)),
#     depth2  = sum(!is.finite(.x$depth_sc2))
#   ))



fit_model <- function(dat, pred_grid, species, size_class, cutoff, ...) {
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)
  
  fit_time <- system.time(
    m <- sdmTMB(
      biomass_g ~ year + depth_sc + depth_sc2 + gear,# + (1|ship),
      # s(time_shot_hour, bs = 'cc') +
      data = dat,
      mesh = mesh,
      family = delta_gamma(type = "poisson-link"),
      offset = log(dat$haul_dur + 5),
      spatial = "on",
      spatiotemporal = "iid",
      time = "year"
    )
  )
  
  print(sanity(m))
  
  saveRDS(m, paste0("output/m_", species, "_", size_class, ".rds"))
  
  # Residuals
  res <- simulate(m, nsim = 300, type = "mle-mvn") |>
    dharma_residuals(m, plot = FALSE)

  p_res <- res |>
    ggplot(aes(observed, expected)) +
    geom_point(color = "grey30", shape = 21, size = 0.5) +
    geom_abline(col = "tomato3", linewidth = 0.6) +
    theme(aspect.ratio = 1) +
    labs(
      x = "Observed", y = "Expected",
      title = paste("DHARMa residuals —", species, size_class)
    )
  ggsave(paste0("output/residuals_", species, "_", size_class, ".png"), p_res)
  
  # Predict
  p <- predict(
    m,
    newdata     = pred_grid,
    offset      = log(as.numeric(pred_grid$haul_dur) + 5),
    re_form_iid = NA
  ) |>
    mutate(est = est1 + est2)
  
  saveRDS(p, paste0("output/pred_", species, "_", size_class, ".rds"))
  
  list(
    model         = m,
    predictions   = p,
    residual_plot = p_res,
    fit_time      = fit_time["elapsed"]
  )
}

# Now run all
species_grid <- species_grid |>
  mutate(
    results = pmap(
      list(
        dat        = data_filtered,
        pred_grid  = pred_grid,
        species    = species,
        size_class = size_class,
        cutoff     = cutoff
      ),
      fit_model
    )
  )


## Test
# dat        <- species_grid$data_filtered[[1]]
# pred_grid  <- species_grid$pred_grid[[1]]
# species    <- species_grid$species[[1]]
# size_class <- species_grid$size_class[[1]]
# cutoff     <- species_grid$cutoff[[1]]
# mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)
# fit_time <- system.time(
#   m <- sdmTMB(
#     biomass_g ~ year + depth_sc + depth_sc2 + gear,# + (1|ship),
#     # s(time_shot_hour, bs = 'cc') +
#     data = dat,
#     mesh = mesh,
#     family = delta_gamma(type = "poisson-link"),
#     offset = log(dat$haul_dur + 5),
#     spatial = "on",
#     spatiotemporal = "off",
#     time = "year"))
# d <- sdmTMB::pcod
# mesh_pcod <- make_mesh(d, c("X", "Y"), cutoff = 10)
# m_test <- sdmTMB(density ~ 1, data = d, mesh = mesh_pcod,
#                  family = tweedie(), spatial = "on")


# Access via e.g.
# species_grid |> filter(species == "cod", size_class == "small") |> pull(results) |> pluck(1, "model")
# species_grid |> filter(species == "cod", size_class == "small") |> pull(results) |> pluck(1, "predictions")
# species_grid |> mutate(fit_time = map_dbl(results, ~ .x$fit_time))

