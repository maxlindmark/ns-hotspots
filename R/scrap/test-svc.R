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
# 2. Converting to biomass using survey + year specific a and b, but must be lots of errors in this data, 
# because some estimates are bonkers. I therefore revert to mean a and b parameters if they are "off"
# 3. I filter from >=1991, because in the data, I don't have catches of haddock in the first years.
# No positive observations means I can't estimate the first model component in the delta gamma


##------ Load and set up data and prediction grid
load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_pok.RData"))

unique(dat$HH$Survey)
class(dat)

# Convert lat/lon to ETRS89-LAEA (EPSG:3035) projected coordinates to minimize distortion in the domain,
# which is important because the model estimates the range of spatial autocorrelation,
# specifically distance where correlation <0.13
# Change to  them from metres to kilometres to help estimation
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
  file = c("datras_data_cod", "datras_data_had", "datras_data_whg"), # "datras_data_pok"), # skipping saithe becase few juveniles caught
  threshold = c(35, 30, 27),
  cutoff = 15 # for mesh. Minimum distance between mesh verticies
)

# Make prediction grid based on domain in cod data (but should be the same for all species)
load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_cod.RData"))

cod_hh <- dat$HH

# Using NOAA bathy from marmap to get depth raster. Likely there are more regional shapefiles
# that have a higher resolution, but for now this will have to do
cod_hh <- cod_hh |>
  janitor::clean_names() |>
  to_lcc(lon = "lon", lat = "lat")

# This is a fiddly grid but my normal approach of a convex hull of the sampling data doesn't play well with UK(!)
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

# Size of cell width (in km)
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


##------ Iterate through all species and summarise data
# Function to process data to get standardized response variable (kg/h)
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
        TimeShotHour, Ship, HaulDur, Month, Day
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
    ) #|>
    #filter(quarter_f %in% c("3", "4"))
  
  d
}

# Loop function across all species
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







# Explore seasonal models:
dat_sub <- dat |>
  filter(sp_sz == "cod_small" & between(year_num, 2015, 2020)) |> 
  mutate(quarter2_f = ifelse(as.character(as.numeric(quarter_f)) %in% c(1, 2),
                             "1.5", "3.5"),
         quarter2_f = as.factor(quarter2_f))

mesh <- make_mesh(dat_sub, xy_cols = c("X", "Y"), cutoff = 20)

plot(mesh)

tic()
m1 <- sdmTMB(
  biomass_kg ~ 0 + year_f + depth_sc + depth_sc2 + gear_f + quarter2_f + (1|ship_f),
  data = dat_sub,
  mesh = mesh,
  family = tweedie(),
  #family = delta_gamma(type = "poisson-link"),
  offset = log(dat_sub$haul_dur + 5),
  spatial = "off",
  spatial_varying = ~0 + quarter2_f,
  spatiotemporal = "off",
  time = "year"
)
toc()

pq <- pred_grid |>
  filter(year == max(dat_sub$year_num)) |> 
  replicate_df("quarter2_f", as.factor(c(1.5, 3.5)))

p1 <- predict(
  m1,
  newdata = pq,
  offset = log(as.numeric(pq$haul_dur) + 5),
  re_form_iid = NA
)

p1 |> 
  pivot_longer(starts_with("zeta_s"), values_to = "svc", names_to = "quarter") |> 
  ggplot(aes(X, Y, fill = svc)) + 
  facet_wrap(~quarter) +
  scale_fill_gradient2() +
  geom_raster()

sanity(m1)

# Make effect of day of the year vary spatially!
library(lubridate)

dat_sub$jday <- yday(make_date(dat_sub$year_num, dat_sub$month, dat_sub$day))
dat_sub$days_in_year <- yday(make_date(dat_sub$year_num, 12, 31))
dat_sub$angle <- 2 * pi * jday / days_in_year

dat_sub$sin_Jday_2pi <- sin(angle)
dat_sub$cos_Jday_2pi <- cos(angle)

tic()
m2 <- sdmTMB(
  biomass_kg ~ 0 + year_f + depth_sc + depth_sc2 + gear_f + sin_Jday_2pi + cos_Jday_2pi + (1|ship_f),
  data = dat_sub,
  mesh = mesh,
  family = tweedie(),
  offset = log(dat_sub$haul_dur + 5),
  spatial = "off",
  spatial_varying = ~0 + sin_Jday_2pi + cos_Jday_2pi,
  spatiotemporal = "off",
  time = "year"
)
toc()

sanity(m2)

ps <- pred_grid |>
  filter(year == max(dat_sub$year_num)) |>
  replicate_df("jday", yday(make_date(2020, 1:12, 15))) |>
  mutate(
    date = as_date(jday - 1, origin = make_date(2020, 1, 1)),
    angle = 2 * pi * jday / yday(make_date(2020, 12, 31)),
    sin_Jday_2pi = sin(angle),
    cos_Jday_2pi = cos(angle),
    month2 = month(date)
  ) |> 
  mutate(quarter2_f = ifelse(month2 >= 1 & month2 <= 6,
                             "1.5", "3.5"))

p2 <- predict(
  m2,
  newdata = ps,
  offset = log(as.numeric(ps$haul_dur) + 5),
  re_form_iid = NA
)

summary(exp(p2$est))

p2 |> 
  pivot_longer(starts_with("zeta_s"), values_to = "svc", names_to = "month3") |> 
  ggplot(aes(X, Y, fill = svc)) + 
  facet_wrap(~month3) +
  scale_fill_gradient2() +
  geom_raster()

p2 |> 
  ggplot(aes(X, Y, fill = est_rf)) + 
  facet_wrap(~month2) +
  scale_fill_gradient2() +
  geom_raster()

pp2 <- p2 |> 
  mutate(est = ifelse(est > quantile(est, 0.9999), quantile(est, 0.9999), est)) |> 
  ggplot(aes(X, Y, fill = exp(est))) + 
  facet_wrap(~month2, nrow = 2) +
  scale_fill_viridis(trans = "fourth_root_power") +
  geom_raster() +
  guides(fill = "none") +
  ggtitle("Spatially varying cyclical slope of day")

pp1 <- predict(m1,
  newdata = ps,
  offset = log(as.numeric(ps$haul_dur) + 5),
  re_form_iid = NA
  ) |> 
  mutate(est = ifelse(est > quantile(est, 0.9999), quantile(est, 0.9999), est)) |> 
  ggplot(aes(X, Y, fill = exp(est))) + 
  facet_wrap(~month2, nrow = 2) +
  scale_fill_viridis(trans = "fourth_root_power") +
  geom_raster() +
  guides(fill = "none") +
  ggtitle("Spatially varying quarter")

library(patchwork)
(pp2 / pp1) & theme(legend.position = "bottom")


# Now add data as points...
dat_sub <- dat_sub |>
  mutate(
    date = as_date(jday - 1, origin = make_date(2020, 1, 1)),
    month2 = month(date)
  )

pp2 <- p2 |> 
  mutate(est = ifelse(est > quantile(est, 0.9999), quantile(est, 0.9999), est)) |> 
  ggplot(aes(X, Y, fill = exp(est))) + 
  facet_wrap(~month2, nrow = 2) +
  scale_fill_viridis(trans = "fourth_root_power") +
  geom_raster() +
  geom_point(data = dat_sub,
             aes(X, Y), inherit.aes = FALSE, size = 0.01, color = "red") +
  guides(fill = "none") +
  ggtitle("Spatially varying cyclical slope of day")

pp2

pp1 <- predict(m1,
               newdata = ps,
               offset = log(as.numeric(ps$haul_dur) + 5),
               re_form_iid = NA
) |> 
  mutate(est = ifelse(est > quantile(est, 0.9999), quantile(est, 0.9999), est)) |> 
  ggplot(aes(X, Y, fill = exp(est))) + 
  facet_wrap(~month2, nrow = 2) +
  scale_fill_viridis(trans = "fourth_root_power") +
  geom_raster() +
  geom_point(data = dat_sub,
             aes(X, Y), inherit.aes = FALSE, size = 0.01, color = "red") +
  guides(fill = "none") +
  ggtitle("Spatially varying quarter")

library(patchwork)
(pp2 / pp1) & theme(legend.position = "bottom")

AIC(m1, m2)
