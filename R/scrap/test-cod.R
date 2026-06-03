library(forcats)
library(tictoc)
library(viridis)
library(sf)
library(sdmTMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsidekick)
theme_set(theme_sleek())

home <- here::here()

##### ------------ Prepare data
## Length-weight parameters
# TODO: estimate these every year / spatial location / life stage using built in CA data
a <- 0.01
b <- 3.0

# Cod
load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_cod.RData"))

cod_by_size <- dat$HL |>
  mutate(
    size_class = if_else(LngtCm < 35, "small", "large"),
    biomass_g  = Count * a * LngtCm^b
  ) |>
  summarise(
    n = sum(Count, na.rm = TRUE),
    biomass_g = sum(biomass_g, na.rm = TRUE),
    .by = c(haul.id, size_class)
  )

cod <- expand_grid(
  haul.id    = dat$HH$haul.id,
  size_class = c("small", "large")
) |>
  left_join(cod_by_size, by = c("haul.id", "size_class")) |>
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
  mutate(n = round(n)) |>
  janitor::clean_names() |>
  add_utm_columns(ll_names = c("lon", "lat")) |>
  filter(!is.na(haul_dur)) |>
  mutate(gear = fct_collapse(
    gear,
    other = c("BT4P", "BT4S", "DHT", "GRT")
  )) |>
  filter(!n > 10000)

# Split by size class
cod_small <- cod |> filter(fish_size == "small")
cod_large <- cod |> filter(fish_size == "large")

# ! Haddock
# load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_had.RData"))
# 
# # ! Pollock
# load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_pok.RData"))
# 
# # ! Whiting
# load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_pok.RData"))



##### ------------ Make prediction grid
# Fit depth surface from observed hauls
depth_model <- mgcv::gam(
  depth ~ s(X, Y, bs = "tp", k = 300),
  data = cod |> distinct(X, Y, depth)
)

# Base spatial grid clipped to convex hull
x <- cod$X
y <- cod$Y
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

cell_width <- 10

spatial_grid <- expand.grid(
  X = seq(min(cod$X), max(cod$X), cell_width),
  Y = seq(min(cod$Y), max(cod$Y), cell_width)
)

sp::coordinates(spatial_grid) <- c("X", "Y")
inside <- !is.na(sp::over(spatial_grid, as(sp_poly_df, "SpatialPolygons")))
spatial_grid <- as.data.frame(spatial_grid[inside, ])

# Impute depth and remove land
spatial_grid$depth <- predict(depth_model, newdata = spatial_grid)
spatial_grid <- spatial_grid |> filter(depth > 0)

# Expand by year, add covariates (no fish_size — separate models)
pred_grid <- expand_grid(
  spatial_grid,
  year = sort(unique(cod$year))
) |>
  mutate(
    gear           = factor("GOV", levels = levels(cod$gear)),
    ship           = cod$ship[1],
    time_shot_hour = 12,
    haul_dur       = 30
  )


# world <- ne_countries(scale = "medium", returnclass = "sf")
# world_utm <- world |> st_transform(crs = 32631)
# ggplot() +
#   geom_raster(
#     data = pred_grid |> distinct(X, Y, depth),
#     aes(x = X * 1000, y = Y * 1000, fill = depth)) +
#   geom_sf(data = world_utm) +
#   scale_fill_viridis() +
#   coord_sf(
#     xlim = range(cod$X) * 1000,
#     ylim = range(cod$Y) * 1000,
#     crs = 32631, expand = 0
#   ) +
#   labs(x = "Easting (km)", y = "Northing (km)") +
#   theme(axis.text.x = element_text(angle = 90))


##### ------------ Fit models
tic()
mcod_small <- sdmTMB(
  biomass_g ~ year + gear,
  s(depth, k = 3) +
  # s(time_shot_hour, bs = 'cc') +
  (1|ship),
  data = cod_small,
  mesh = mesh,
  family = delta_gamma(type = "poisson-link"),
  offset = log(cod_small$haul_dur + 5),
  spatial = "on",
  spatiotemporal = "iid",
  time = "year"
)
toc()

sanity(mcod_small)

saveRDS(mcod_small, "output/mcod_small.rds")


##### ------------ Plot residuals
res <- simulate(mcod_small, nsim = 300, type = "mle-mvn") |>
  dharma_residuals(mcod_small, plot = FALSE)

res |> 
  ggplot(aes(observed, expected)) +
  geom_point(color = "grey30", shape = 21, size = 0.5) +
  geom_abline(col = "tomato3", linewidth = 0.6) +
  theme(aspect.ratio = 1) +
  labs(x = "Observed", y = "Expected")

##### ------------ Predict on grid
pcod_small <- predict(
  mcod_small,
  newdata = pred_grid,
  offset  = log(pred_grid$haul_dur + 5),
  re_form_iid = NA,
  ) |> 
  mutate(est = est1 + est2)

saveRDS(pcod_small, "output/pcod_small.rds")



### extra
world <- ne_countries(scale = "medium", returnclass = "sf")
world_utm <- world |> st_transform(crs = 32631)

## Plot haul locations coloured by depth
p_small |> 
  #filter(year == 2023) |> 
  ggplot() +
  geom_raster(
    aes(x = X * 1000, y = Y * 1000, fill = exp(est))
  ) +
  geom_sf(data = world_utm) +
  facet_wrap(~year) +
  scale_fill_viridis(trans = "sqrt") +
  scale_fill_viridis() +
  coord_sf(
    xlim = range(cod$X) * 1000,
    ylim = range(cod$Y) * 1000,
    crs = 32631, expand = 0
  ) +
  labs(x = "Easting (km)", y = "Northing (km)") +
  theme(axis.text.x = element_text(angle = 90))
