library(concaveman)
library(marmap)
library(forcats)
library(stringr)
library(sdmTMB)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(readr)
library(ggsidekick);theme_set(theme_sleek())
library(viridis)

home <- here::here()


##------ Load the Norwegian data and prepare it so that I can merge it with the DATRAS data
# Read raw file lines
lines <- read_lines(
  paste0(home, "/data/NO_shrimp_survey_gadoids_DATRAS.csv")
)

# Find all header rows
header_idx <- which(str_detect(lines, "^RecordType"))

# Add artificial end position
header_end <- c(header_idx[-1] - 1, length(lines))

# Function to parse one block
parse_block <- function(start, end) {
  
  header <- lines[start]
  
  # Determine record type from first data row
  first_data <- lines[start + 1]
  rec_type <- str_split(first_data, "\t|,", simplify = TRUE)[1]
  
  # Extract lines for this block
  block_lines <- lines[start:end]
  
  # Create temporary text connection
  block_text <- paste(block_lines, collapse = "\n")
  
  # Read as dataframe
  df <- read_delim(
    I(block_text),
    delim = ",",
    show_col_types = FALSE,
    trim_ws = TRUE
  )
  
  list(
    type = rec_type,
    data = df
  )
}

# Parse all blocks
blocks <- map2(header_idx, header_end, parse_block)

# Combine blocks of same RecordType
no_HH <- blocks |>
  keep(~ .x$type == "HH") |>
  map("data") |>
  bind_rows()

no_HL <- blocks |>
  keep(~ .x$type == "HL") |>
  map("data") |>
  bind_rows()

no_CA <- blocks |>
  keep(~ .x$type == "CA") |>
  map("data") |>
  bind_rows()












##------ Iterate through all species and summarise data
# Function to process data to get standardized response variable (kg/h)
prepare_species_data <- function(file, threshold, depth_mean, depth_sd, ...) {
  load(paste0(home, "/data/datras_data_cod_had_whg_pok/", file, ".RData"))
  
  CA <- dat$CA
  
  # Estimate length-weight parameters from CA data by year and survey
  lw <- CA |>
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

ggplot() +
  scale_fill_manual(values = c("grey95", "forestgreen")) +
  geom_point(data = dat, aes(X*1000, Y*1000, color = as.factor(survey)),
             alpha = 0.25, size = 0.75) +
  geom_sf(data = world_lcc) +
  theme(legend.position = "bottom") +
  coord_sf(
    xlim = range(dat$X) * 1000,
    ylim = range(dat$Y) * 1000,
    crs = crs_lcc,
    expand = 0
  ) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.direction = "horizontal",
    #legend.key.size = unit(5, "cm"),
    legend.position.inside = c(0.96, 0.07),
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "Longitude", y = "Latitude") +
  scale_color_brewer(palette = "Dark2", direction = -1)

ggsave(paste0("figures/survey_map.png"), width = 13, height = 14.5, units = "cm")



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


##------ Fit model to all species, save predictions
# Function to fit model and predict
fit_model <- function(dat, pred_grid, species, size_class, cutoff, ...) {
  mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoff)
  
  m <- sdmTMB(
    biomass_kg ~ year_f + depth_sc + depth_sc2 + gear_f + quarter_f + 
      s(month, bs = "cc", k = 5) +
      (1|ship_f),
    data = dat,
    mesh = mesh,
    family = delta_gamma(type = "poisson-link"),
    offset = log(dat$haul_dur + 5),
    spatial = "on",
    spatiotemporal = "iid",
    time = "year"
  )
  
  # print(sanity(m))
  
  saveRDS(m, paste0("output/m_", species, "_", size_class, ".rds"))
  
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
  
  p <- predict(
    m,
    newdata = pred_grid,
    offset = log(as.numeric(pred_grid$haul_dur) + 5),
    re_form_iid = NA
  ) |>
    mutate(est = est1 + est2,
           species = species,
           size_class = size_class)
  
  saveRDS(p, paste0("output/pred_", species, "_", size_class, ".rds"))
  
  list(
    model = m,
    predictions = p,
    residual_plot = p_res
  )
}


# Fit to all species
tictoc::tic()
for (i in unique(dat$sp_sz)) {
  tryCatch(
    {
      cat("Fitting:", i, "\n")
      
      d <- filter(dat, sp_sz == i)
      
      fit_model(
        dat = d,
        pred_grid = pred_grid,
        species = unique(d$species),
        size_class = unique(d$fish_size),
        cutoff = unique(d$cutoff)
      )
      
      cat("Done:", i, "\n")
    },
    error = function(e) {
      cat("ERROR fitting", i, ":", conditionMessage(e), "\n")
    }
  )
}
tictoc::toc()

##------ Plotting is done in script 02

# Just a final check on the model convergence
# All important stuff is good. The ln_smooth_sigma is typically large,
# and might be better estimate on the full data set (all months)
list.files(path = paste0(home, "/output"),
           pattern = "^m_.*\\.rds$",
           full.names = TRUE) |>
  set_names(~ tools::file_path_sans_ext(basename(.))) |>
  map(readRDS) |>
  map(sanity)

