home <- here::here()

##------ Load and set up data and prediction grid
load(paste0(home, "/data/datras_data_cod_had_whg_pok/datras_data_pok.RData"))

unique(dat$HH$Survey)

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
