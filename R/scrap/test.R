library(concaveman)
library(forcats)
library(stringr)
library(sdmTMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(viridis)

home <- here::here()

# Modelling notebook:
# 1. Count models are not estimated very well (delta_truncated nbinom best QQ, but some estimates explode)
# 2. Converting to biomass using survey + year specific a and b, but must be lots of errors in this data
# because some estimates are bonkers. I therefore revert to mean a and b parameters.

load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_pok.RData"))

# Lambert Conformal Conic optimized for North Sea (because it spans too many UTM zones)
crs_lcc <- "+proj=lcc +lat_1=50 +lat_2=60 +lat_0=55 +lon_0=2 +datum=WGS84 +units=km +no_defs"

to_lcc <- function(df, lon = "lon", lat = "lat") {
  df |>
    sf::st_as_sf(coords = c(lon, lat), crs = 4326) |>
    sf::st_transform(crs_lcc) |>
    dplyr::mutate(
      X = sf::st_coordinates(geometry)[, 1],
      Y = sf::st_coordinates(geometry)[, 2]
    ) |>
    sf::st_drop_geometry()
}

# Species metadata grid with settings
species_grid <- tibble(
  species = c("cod", "had", "whg"),
  file = c("datras_data_cod", "datras_data_had", "datras_data_whg"), # "datras_data_pok"),
  threshold = c(35, 30, 27),
  cutoff = 25 # for mesh
)

# Make prediction grid based on domain in cod data (but should be the same)
load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_cod.RData"))

cod_hh <- dat$HH

depth_mean <- mean(cod_hh$Depth, na.rm = TRUE)
depth_sd <- sd(cod_hh$Depth, na.rm = TRUE)

cod_hh <- cod_hh |>
  janitor::clean_names() |>
  to_lcc(lon = "lon", lat = "lat")

depth_model <- mgcv::gam(
  depth ~ s(X, Y, bs = "tp", k = 200),
  data = cod_hh |> distinct(X, Y, depth)
)

# x <- cod_hh$X
# y <- cod_hh$Y
# z <- chull(x, y)
# plot(z)
# coords <- cbind(x[z], y[z])
# coords <- rbind(coords, coords[1, ])

cod_sf <- cod_hh |> 
  distinct(X, Y) |> 
  st_as_sf(coords = c("X", "Y"), crs = crs_lcc)  # already projected
hull <- concaveman(cod_sf, concavity = 2)
plot(hull)
hull_buf <- st_buffer(hull, dist = 50)
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

spatial_grid$depth <- predict(depth_model, newdata = spatial_grid)

ggplot(spatial_grid, aes(X, Y)) + 
  geom_raster() + 
  geom_point(data = cod_hh, aes(X, Y, color = depth))


pred_grid_sf <- st_as_sf(pred_grid, coords = c("X", "Y"), crs = crs_lcc)
pred_grid_ll <- st_transform(pred_grid_sf, 4326)  # WGS84 lon/lat
pred_grid$lon <- st_coordinates(pred_grid_ll)[,1]
pred_grid$lat <- st_coordinates(pred_grid_ll)[,2]



