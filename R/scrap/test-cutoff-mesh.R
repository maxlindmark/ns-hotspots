library(concaveman)
library(marmap)
library(forcats)
library(stringr)
library(sdmTMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(ggsidekick);theme_set(theme_sleek())
library(viridis)

home <- here::here()

# Modelling notebook:
# 1. Count models are not estimated very well (delta_truncated nbinom best QQ, but some estimates explode)
# 2. Converting to biomass using survey + year specific a and b, but must be lots of errors in this data
# because some estimates are bonkers. I therefore revert to mean a and b parameters.
# 3. I filter from 1991 and onwards. I don't have catches of haddock in the first years and this doesnt 
# play well with a delta gamma (can't estimate probability of occurence if no catches)

load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_pok.RData"))

unique(dat$HH$Survey)

crs_lcc <- 3035  # ETRS89-LAEA

to_lcc <- function(df, lon = "lon", lat = "lat") {
  df |>
    sf::st_as_sf(coords = c(lon, lat), crs = 4326) |>
    sf::st_transform(crs_lcc) |>
    dplyr::mutate(
      X = sf::st_coordinates(geometry)[, 1] / 1000,
      Y = sf::st_coordinates(geometry)[, 2] / 1000
    ) |>
    sf::st_drop_geometry()
}

# Species metadata grid with settings
species_grid <- tibble(
  species = c("cod", "had", "whg"),
  file = c("datras_data_cod", "datras_data_had", "datras_data_whg"), # "datras_data_pok"),
  threshold = c(35, 30, 27),
  cutoff = 10 # for mesh
)

# Make prediction grid based on domain in cod data (but should be the same)
load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_cod.RData"))

cod_hh <- dat$HH

# Using NOAA bathy from marmap instead...
cod_hh <- cod_hh |>
  janitor::clean_names() |>
  to_lcc(lon = "lon", lat = "lat")

# This is a fiddly grid but my normal approach of convex hull doesnt play well with UK
cod_sf <- cod_hh |> 
  distinct(X, Y) |> 
  st_as_sf(coords = c("X", "Y"), crs = NA)
hull <- concaveman(cod_sf, concavity = 2)
plot(hull)
hull_buf <- st_buffer(hull, dist = 20)
plot(hull_buf)
coords <- st_coordinates(hull_buf)[, 1:2]

sp_poly <- sp::SpatialPolygons(
  list(sp::Polygons(list(sp::Polygon(coords)), ID = 1))
)
sp_poly_df <- sp::SpatialPolygonsDataFrame(
  sp_poly,
  data = data.frame(ID = 1)
)

cell_width <- 3

spatial_grid <- expand.grid(
  X = seq(min(cod_hh$X), max(cod_hh$X), cell_width),
  Y = seq(min(cod_hh$Y), max(cod_hh$Y), cell_width)
)

sp::coordinates(spatial_grid) <- c("X", "Y")
inside <- !is.na(sp::over(spatial_grid, as(sp_poly_df, "SpatialPolygons")))
spatial_grid <- as.data.frame(spatial_grid[inside, ])

pred_grid <- expand_grid(
  spatial_grid,
  year = sort(unique(cod_hh$year))
) |>
  mutate(
    lon = st_coordinates(st_transform(st_as_sf(data.frame(X = X * 1000, Y = Y * 1000), 
                                               coords = c("X", "Y"), crs = crs_lcc), 4326))[, 1],
    lat = st_coordinates(st_transform(st_as_sf(data.frame(X = X * 1000, Y = Y * 1000), 
                                               coords = c("X", "Y"), crs = crs_lcc), 4326))[, 2],
    gear = as.factor("gov"),
    ship_f = factor(cod_hh$ship[1]),
    time_shot_hour = 12,
    haul_dur = 30,
    haul_dur_sc = 0,
    year_f = as.factor(year),
    gear_f = as.factor(gear),
    month = 8,
    month_f = as.factor(month),
    quarter_f = as.factor(3)
  )

ggplot(pred_grid |> filter(year == 2002), aes(X, Y)) + 
  geom_raster() + 
  geom_point(data = cod_hh, aes(X, Y, color = "data"))

depth_box <- getNOAA.bathy(min(pred_grid$lon) - .1, max(pred_grid$lon) + .1, min(pred_grid$lat) - .1, max(pred_grid$lat) + .1)
pred_grid$depth <- get.depth(depth_box, x = pred_grid$lon, y = pred_grid$lat, locator = F)$depth
pred_grid$depth <- pred_grid$depth * (-1)

ggplot(pred_grid |>
         filter(year == 2000 & depth < 200 & depth > 1),
       aes(X, Y, fill = depth)) +
  geom_raster()

pred_grid <- pred_grid |> 
  filter(depth < 300 & depth > 1)

ggplot(pred_grid,
       aes(X, Y, fill = depth)) +
  geom_raster()

# Scale depth (use that for data also)
depth_mean <- mean(pred_grid$depth, na.rm = TRUE)
depth_sd <- sd(pred_grid$depth, na.rm = TRUE)

pred_grid <- pred_grid |> 
  mutate(depth_sc = (depth - depth_mean) / depth_sd,
         depth_sc2 = depth_sc^2)

# Test plot
world <- ne_countries(scale = "medium", returnclass = "sf")
world_lcc <- world |> st_transform(crs = crs_lcc)

ggplot() +
  geom_raster(
    data = pred_grid,
    aes(x = X * 1000, y = Y * 1000, fill = depth)
  ) +
  facet_wrap(~year, ncol = 7) +
  scale_fill_viridis() +
  geom_sf(data = world_lcc) +
  coord_sf(
    xlim = range(pred_grid$X) * 1000,
    ylim = range(pred_grid$Y) * 1000,
    crs = crs_lcc,
    expand = 0
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    position = "inside", name = "Depth"
  )) +
  labs(x = "Easting (km)", y = "Northing (km)") +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.direction = "horizontal",
    legend.position.inside = c(0.63, 0.057),
    legend.key.width = unit(1, "cm")
  )


# Prepare data
prepare_species_data <- function(file, threshold, depth_mean, depth_sd, ...) {
  load(paste0(home, "/data/datras_data_cod_had_whg_pok/", file, ".RData"))
  
  # Estimate length-weight parameters from CA data by year and survey
  lw <- dat$CA |>
    drop_na(IndWgt, LngtClas) |>
    filter(IndWgt > 0, LngtClas > 0) |>
    tidyr::uncount(NoAtALK) |>
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
      b = if_else(b < 2.9 | b > 3.1, NA, b),
      a = if_else(a < quantile(a, 0.2) | a > quantile(a, 0.8), NA, a)
    )
  
  a_fallback <- median(lw$a, na.rm = TRUE)
  b_fallback <- median(lw$b, na.rm = TRUE)
  
  by_size <- dat$HL |>
    left_join(lw, by = c("Year" = "year", "Survey" = "survey")) |>
    mutate(
      a = replace_na(a, a_fallback),
      b = replace_na(b, b_fallback),
      size_class = if_else(LngtCm < threshold, "small", "large"),
      biomass_g = Count * a * LngtCm^b,
      biomass_kg = biomass_g / 1000
    ) |>
    summarise(
      n = sum(Count, na.rm = TRUE),
      biomass_kg = sum(biomass_kg, na.rm = TRUE),
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
      biomass_kg = replace_na(biomass_kg, 0),
      a = replace_na(a, a_fallback),
      b = replace_na(b, b_fallback)
    ) |>
    left_join(
      dat$HH |> dplyr::select(
        haul.id, Year, Survey, Gear2, lon, lat, Depth,
        TimeShotHour, Ship, HaulDur, Month
      ),
      by = "haul.id"
    ) |>
    rename(Gear = Gear2) |> 
    mutate(
      fish_size = factor(size_class, levels = c("small", "large")),
      Year = as.factor(Year),
      Gear = as.factor(Gear)
    ) |>
    dplyr::select(-size_class) |>
    janitor::clean_names() |>
    to_lcc(lon = "lon", lat = "lat") |>
    filter(!is.na(haul_dur)) |>
    mutate(
      month_f = as.factor(month),
      ship_f = as.factor(ship),
      year_f = as.factor(year),
      gear = stringr::str_to_lower(gear),
      gear_f = as.factor(gear),
      haul_dur_sc = (haul_dur - mean(haul_dur)) / sd(haul_dur),
      quarter_f = as.factor(ceiling(month / 3))
    ) |>
    filter(quarter_f %in% c("3", "4"))
  
  d
}

# Load and combine all species
dat <- list()

for (i in 1:nrow(species_grid)) {
  tryCatch(
    {
      sp <- species_grid$species[i]
      file <- species_grid$file[i]
      threshold <- species_grid$threshold[i]
      cutoff <- species_grid$cutoff[i]
      
      dat[[i]] <- prepare_species_data(
        file = file,
        threshold = threshold,
        depth_mean = depth_mean,
        depth_sd = depth_sd
      ) |>
        mutate(
          species = sp,
          threshold = threshold,
          cutoff = cutoff
        )
      
      cat("Prepared:", sp, "\n")
    },
    error = function(e) {
      cat("ERROR preparing", sp, ":", conditionMessage(e), "\n")
    }
  )
}

dat <- bind_rows(dat)
dat$sp_sz <- paste(dat$species, dat$fish_size, sep = "_")

# Final preparations to data and grid
dat <- dat |>
  mutate(
    lon = st_coordinates(st_transform(st_as_sf(data.frame(X = X * 1000, Y = Y * 1000),
                                               coords = c("X", "Y"), crs = crs_lcc), 4326))[, 1],
    lat = st_coordinates(st_transform(st_as_sf(data.frame(X = X * 1000, Y = Y * 1000),
                                               coords = c("X", "Y"), crs = crs_lcc), 4326))[, 2]
  )

dat$depth <- get.depth(depth_box, x = dat$lon, y = dat$lat, locator = F)$depth
dat$depth <- dat$depth * (-1)

# Filter years from 1990+
dat <- dat |>
  mutate(year_num = as.numeric(as.character(year))) |>
  tidylog::filter(year_num >= 1991)

pred_grid <- pred_grid |>
  mutate(year_num = as.numeric(as.character(year))) |>
  tidylog::filter(year_num >= 1991)

dat <- dat |> 
  mutate(depth_sc = (depth - depth_mean) / depth_sd,
         depth_sc2 = depth_sc^2)

ggplot(dat, aes(a)) +
  geom_histogram()

ggplot(dat, aes(b)) +
  geom_histogram()

# Plot data, loop through all species
data_plot <- ggplot() +
  facet_wrap(~year, ncol = 7) +
  geom_sf(data = world_lcc) +
  coord_sf(
    xlim = range(pred_grid$X),
    ylim = range(pred_grid$Y),
    crs = crs_lcc,
    expand = 0
  ) +
  guides(color = guide_legend(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  labs(x = "Easting (km)", y = "Northing (km)") +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.direction = "horizontal",
    legend.position.inside = c(0.93, 0.08),
    legend.key.width = unit(1, "cm")
  )

for(i in unique(dat$sp_sz)){
  
  dd <- filter(dat, sp_sz == i) |> 
    mutate(Presence = ifelse(biomass_kg == 0, "0", "1"))
  
  data_plot +
    geom_point(
      data = dd,
      size = 0.75,
      aes(x = X * 1000, y = Y * 1000, color = Presence)
    ) +
    geom_sf(data = world_lcc) +
    coord_sf(
      xlim = range(pred_grid$X) * 1000,
      ylim = range(pred_grid$Y) * 1000,
      crs = crs_lcc,
      expand = 0
    )
  
  ggsave(paste0("figures/data_all_yrs_", i, ".pdf"), width = 20, height = 16)
  
}

ggplot(dat, aes(biomass_kg)) +
  geom_histogram() +
  facet_wrap(~sp_sz, scales = "free")

dat |>
  group_by(sp_sz) |>
  summarise(
    n = n(),
    min = min(biomass_kg, na.rm = TRUE),
    q01 = quantile(biomass_kg, 0.01, na.rm = TRUE),
    q05 = quantile(biomass_kg, 0.05, na.rm = TRUE),
    q25 = quantile(biomass_kg, 0.25, na.rm = TRUE),
    median = quantile(biomass_kg, 0.5, na.rm = TRUE),
    q75 = quantile(biomass_kg, 0.75, na.rm = TRUE),
    q95 = quantile(biomass_kg, 0.95, na.rm = TRUE),
    q99 = quantile(biomass_kg, 0.99, na.rm = TRUE),
    max = max(biomass_kg, na.rm = TRUE)
  )


## Fit models
dat$sp_sz

dd <- filter(dat, sp_sz == "cod_small")

mesh <- make_mesh(dd, xy_cols = c("X", "Y"), cutoff = cutoff)
plot(mesh)

m <- sdmTMB(
  biomass_kg ~ year_f + depth_sc + depth_sc2 + gear_f + quarter_f + 
    s(month, bs = "cc", k = 5) +
    (1|ship_f),
  data = dd,
  mesh = mesh,
  family = delta_gamma(type = "poisson-link"),
  offset = log(dd$haul_dur + 5),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year"
)

p <- predict(
  m,
  newdata = pred_grid,
  offset = log(as.numeric(pred_grid$haul_dur) + 5),
  re_form_iid = NA
) |>
  mutate(est = est1 + est2,
         species = species,
         size_class = size_class)





