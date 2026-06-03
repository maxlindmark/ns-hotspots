library(ggplot2)
library(sf)
library(rnaturalearth)
library(viridis)
library(tidyr)
library(dplyr)
library(purrr)
library(forcats)
library(ggsidekick)
theme_set(theme_sleek())

home <- here::here()
crs_lcc <- 3035 # ETRS89-LAEA

# Load predictions
# see this: https://stackoverflow.com/questions/38797714/reading-multiple-rds-files
preds <- list.files(
  path = paste0(home, "/output"),
  pattern = "^pred.*\\.rds$",
  full.names = TRUE
) |>
  map(readRDS) |>
  bind_rows() |>
  mutate(
    sp_sz_code = fct_recode(
      sp_sz,
      "cod_small" = "Cod < 35 cm",
      "had_small" = "Haddock < 30 cm",
      "whg_small" = "Whiting < 27 cm",
      "cod_large" = "Cod ≥ 35 cm",
      # "had_large" = "Haddock ≥ 30 cm",
      "whg_large" = "Whiting ≥ 27 cm"
    )
  ) |>
  mutate(
    est2 = ifelse(est > quantile(est, probs = 0.9995), quantile(est, probs = 0.9995), est),
    .by = c(sp_sz, quarter_merge)
  )

# Map stuff
world <- ne_countries(scale = "medium", returnclass = "sf")
world_lcc <- world |> st_transform(crs = crs_lcc)

# Plot species
data_plot <- ggplot() +
  #facet_wrap(~year, ncol = 7) +
  geom_sf(data = world_lcc) +
  # scale_fill_viridis(trans = "fourth_root_power") +
  scale_fill_viridis(
    trans = "sqrt",
    breaks = scales::pretty_breaks(n = 3)
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5
  )) +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.8, "cm")
  ) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "kg per 30 min"
  )


### --- Plot predictions

for (i in unique(preds$sp_sz_code)) {
  dd <- filter(preds, sp_sz_code == i)

  for (j in unique(dd$quarter_merge)) {
    # All years
    data_plot +
      geom_raster(
        data = dd |> filter(quarter_merge == j),
        aes(x = X * 1000, y = Y * 1000, fill = exp(est))
      ) +
      facet_wrap(~year, ncol = 7) +
      geom_sf(data = world_lcc) +
      coord_sf(
        xlim = range(preds$X) * 1000,
        ylim = range(preds$Y) * 1000,
        crs = crs_lcc,
        expand = 0
      ) +
      ggtitle(paste0(dd$sp_sz[1], ", Q", j))

    ggsave(paste0(home, "/figures/predictions/pred_all_yrs_", i, j, ".png"), width = 22, height = 20, units = "cm")

    # Selected years
    data_plot +
      geom_raster(
        data = dd |>
          filter(year_f %in% as.factor(2019:2024)) |>
          filter(quarter_merge == j),
        aes(x = X * 1000, y = Y * 1000, fill = exp(est))
      ) +
      facet_wrap(~year, nrow = 2) +
      geom_sf(data = world_lcc) +
      coord_sf(
        xlim = range(preds$X) * 1000,
        ylim = range(preds$Y) * 1000,
        crs = crs_lcc,
        expand = 0
      ) +
      ggtitle(paste0(dd$sp_sz[1], ", Q", j))

    ggsave(paste0(home, "/figures/predictions/pred_recent_yrs_", i, j, ".png"), width = 19, height = 14, units = "cm")
  }
}


### --- Plot hotspots

# Only above quantile
for (i in unique(preds$sp_sz_code)) {
  dd <- filter(preds, sp_sz_code == i)

  for (j in unique(dd$quarter_merge)) {
    hotspot <- data_plot +
      geom_raster(
        data = dd |>
          filter(year_f %in% as.factor(2019:2024)) |>
          mutate(`Above 90% quantile` = ifelse(est > quantile(est, probs = 0.9), "1", "0")) |>
          filter(quarter_merge == j),
        aes(x = X * 1000, y = Y * 1000, fill = `Above 90% quantile`)
      ) +
      scale_fill_manual(values = c("grey95", "forestgreen")) +
      facet_wrap(~year, nrow = 2) +
      geom_sf(data = world_lcc) +
      theme(legend.position = "bottom") +
      coord_sf(
        xlim = range(preds$X) * 1000,
        ylim = range(preds$Y) * 1000,
        crs = crs_lcc,
        expand = 0
      ) +
      guides(fill = guide_colorbar(
        title.position = "top", title.hjust = 0.5,
        position = "inside"
      )) +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.direction = "horizontal",
        legend.key.width = unit(0.3, "cm"),
        legend.position.inside = c(0.95, 0.07),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5)
      ) +
      guides(fill = guide_legend(
        title.position = "top", title.hjust = 0.5,
        position = "inside"
      )) +
      labs(
        x = "Longitude", y = "Latitude",
        fill = "Hotspot",
        title = (paste0(dd$sp_sz[1], ", Q", j)),
        caption = "Hotspot calculated across year×quarter"
      )

    hotspot

    ggsave(paste0(home, "/figures/hotspots/hotspot_recent_yrs_", i, j, ".png"), width = 19, height = 14, units = "cm")

    hotspot +
      geom_raster(
        data = dd |>
          filter(year_f %in% as.factor(2019:2024)) |>
          mutate(
            `Above 90% quantile` = ifelse(est > quantile(est, probs = 0.9), "1", "0"),
            .by = c(year, quarter_merge)
          ) |>
          filter(quarter_merge == j),
        aes(x = X * 1000, y = Y * 1000, fill = `Above 90% quantile`)
      ) +
      geom_sf(data = world_lcc) +
      theme(legend.position = "bottom") +
      coord_sf(
        xlim = range(preds$X) * 1000,
        ylim = range(preds$Y) * 1000,
        crs = crs_lcc,
        expand = 0
      ) +
      labs(
        caption = "Hotspot calculated by year×quarter"
      )

    ggsave(paste0(home, "/figures/hotspots/hotspot_recent_yrs_byyear_", i, j, ".png"), width = 19, height = 14, units = "cm")
  }
}

# Same as above but all years!
# TODO


### --- CV of predictions
# Load predictions
# see this: https://stackoverflow.com/questions/38797714/reading-multiple-rds-files
cvs <- list.files(
  path = paste0(home, "/output"),
  pattern = "^cv.*\\.rds$",
  full.names = TRUE
) |>
  map(readRDS) |>
  bind_rows() |>
  mutate(
    sp_sz_code = fct_recode(
      sp_sz,
      "cod_small" = "Cod < 35 cm",
      "had_small" = "Haddock < 30 cm",
      "whg_small" = "Whiting < 27 cm",
      "cod_large" = "Cod ≥ 35 cm",
      # "had_large" = "Haddock ≥ 30 cm",
      "whg_large" = "Whiting ≥ 27 cm"
    )
  ) |> 
  mutate(quarter_merge = ifelse(quarter_merge == 1.5, "Quarter=1.5", "Quarter=2.5"))

# cvs <- list.files(
#   path = paste0(home, "/output"),
#   pattern = "^cv.*\\.rds$",
#   full.names = TRUE
# ) |>
#   map(readRDS) |>
#   bind_rows()
# 
# cvs$sp_sz <- rep(c("Cod < 35 cm", "Haddock < 30 cm", "Whiting < 27 cm", "Cod ≥ 35 cm", "Whiting ≥ 27 cm"),
#                  each = nrow(cvs)/5)
# 
# cvs <- cvs |>
#   mutate(
#     sp_sz_code = fct_recode(
#       sp_sz,
#       "cod_small" = "Cod < 35 cm",
#       "had_small" = "Haddock < 30 cm",
#       "whg_small" = "Whiting < 27 cm",
#       "cod_large" = "Cod ≥ 35 cm",
#       # "had_large" = "Haddock ≥ 30 cm",
#       "whg_large" = "Whiting ≥ 27 cm"
#     )
#   )

for (i in unique(cvs$sp_sz_code)) {
  dd <- filter(cvs, sp_sz_code == i)

  data_plot +
    geom_raster(
      data = dd,
      aes(x = X * 1000, y = Y * 1000, fill = est)
    ) +
    facet_wrap(~quarter_merge) +
    scale_fill_viridis(option = "magma") +
    labs(
      fill = "CV (prediction)",
      title = dd$sp_sz[1]
    ) +
    geom_sf(data = world_lcc) +
    coord_sf(
      xlim = range(cvs$X) * 1000,
      ylim = range(cvs$Y) * 1000,
      crs = crs_lcc,
      expand = 0
    )

  ggsave(paste0(home, "/figures/cv/cv_est_recent_yrs_", i, ".png"), width = 19, height = 13, units = "cm")
}


### --- Mean hotspot value across sims

sims <- list.files(
  path = paste0(home, "/output"),
  pattern = "^sim.*\\.rds$",
  full.names = TRUE
) |>
  map(readRDS) |>
  bind_rows() |>
  mutate(
    sp_sz_code = fct_recode(
      sp_sz,
      "cod_small" = "Cod < 35 cm",
      "had_small" = "Haddock < 30 cm",
      "whg_small" = "Whiting < 27 cm",
      "cod_large" = "Cod ≥ 35 cm",
      # "had_large" = "Haddock ≥ 30 cm",
      "whg_large" = "Whiting ≥ 27 cm"
    )
  )

# sims <- list.files(
#   path = paste0(home, "/output"),
#   pattern = "^sim.*\\.rds$",
#   full.names = TRUE
# ) |>
#   map(readRDS) |>
#   bind_rows()
# 
# sims$sp_sz <- rep(c("Cod < 35 cm", "Haddock < 30 cm", "Whiting < 27 cm", "Cod ≥ 35 cm", "Whiting ≥ 27 cm"),
#                  each = nrow(sims)/5)
# 
# sims <- sims |>
#   mutate(
#     sp_sz_code = fct_recode(
#       sp_sz,
#       "cod_small" = "Cod < 35 cm",
#       "had_small" = "Haddock < 30 cm",
#       "whg_small" = "Whiting < 27 cm",
#       "cod_large" = "Cod ≥ 35 cm",
#       # "had_large" = "Haddock ≥ 30 cm",
#       "whg_large" = "Whiting ≥ 27 cm"
#     )
#   )

for (i in unique(sims$sp_sz_code)) {
  dd <- filter(sims, sp_sz_code == i)
  
  dd_sum <- dd |> 
    filter(year_f %in% as.factor(2021:2024)) |>
    pivot_longer(V1:V30) |> 
    mutate(hotspot = ifelse(value > quantile(value, probs = 0.9), 1, 0),
           .by = c(year, quarter_merge, name)) |>
    summarise(`Hotspot stability` = mean(hotspot), 
              .by = c(X, Y, year, quarter_merge)) |> 
    mutate(quarter_merge = ifelse(quarter_merge == 1.5, "Quarter=1.5", "Quarter=2.5"))

  data_plot +
      geom_raster(
        data = dd_sum, alpha = 0.9,
        aes(x = X * 1000, y = Y * 1000, fill = `Hotspot stability`)
      ) +
      scale_fill_viridis(option = "mako") +
      facet_grid(quarter_merge~year) +
      geom_sf(data = world_lcc) +
      theme(legend.position = "bottom") +
      coord_sf(
        xlim = range(dd_sum$X) * 1000,
        ylim = range(dd_sum$Y) * 1000,
        crs = crs_lcc,
        expand = 0
      ) +
      labs(
        x = "Longitude", y = "Latitude",
        fill = "Hotspot stability",
        title = paste0(dd$sp_sz[1]),
        subtitle = "Hotspot calculated across simulations"
      )
    
    ggsave(paste0(home, "/figures/hotspots/uncertainty/sims_hotspot_recent_yrs_", i, ".png"), width = 19, height = 13, units = "cm")
    
}

# Facet by sim!
for (i in unique(sims$sp_sz_code)) {
  dd <- filter(sims, sp_sz_code == i)
  
  dd_sum <- dd |> 
    filter(year_f %in% as.factor(c(2024)),
           quarter_merge == "1.5") |>
    pivot_longer(V1:V12, names_to = "sim", values_to = "est") |> 
    mutate(`Above 90% quantile` = ifelse(est > quantile(est, probs = 0.9), "1", "0"),
           .by = c(quarter_merge, sim)) |>
    mutate(quarter_merge = ifelse(quarter_merge == 1.5, "Quarter=1.5", "Quarter=2.5"),
           sim = as.numeric(stringr::str_remove(sim, "V")))
  
  data_plot +
    geom_raster(
      data = dd_sum,
      aes(x = X * 1000, y = Y * 1000, fill = `Above 90% quantile`)
    ) +
    scale_fill_manual(values = c("grey95", "forestgreen")) +
    facet_wrap(~sim, nrow = 3) +
    geom_sf(data = world_lcc) +
    theme(legend.position = "bottom") +
    coord_sf(
      xlim = range(dd_sum$X) * 1000,
      ylim = range(dd_sum$Y) * 1000,
      crs = crs_lcc,
      expand = 0
    ) +
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      position = "inside"
    )) +
    theme(
      axis.text.x = element_text(angle = 0),
      legend.direction = "horizontal",
      legend.key.width = unit(0.3, "cm"),
      legend.position.inside = c(0.95, 0.05),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5)
    ) +
    guides(fill = guide_legend(
      title.position = "top", title.hjust = 0.5,
      position = "inside"
    )) +
    labs(
      x = "Longitude", y = "Latitude",
      fill = "Hotspot",
      title = (paste0(dd$sp_sz[1])),
      subtitle = "Hotspot calculated by sim, 2024, Q1.5"
    )
  
  ggsave(paste0(home, "/figures/hotspots/uncertainty/by_sims_hotspot_recent_yrs_", i, ".png"), width = 19, height = 17, units = "cm")
  
}




# HERE:
# CV OF PREDICTION NOT ONLY SVC
# DO THE REST OF THE PLOTS!





















# Average predictions (small)
d_avg <- preds |>
  filter(size_class == "small") |>
  filter(year_f %in% as.factor(2019:2024)) |>
  mutate(est = exp(est)) |>
  mutate(est = ifelse(est > quantile(est, 0.99), quantile(est, 0.99), est)) |>
  summarise(est = mean(est), .by = c(X, Y, sp_sz)) |>
  mutate(est_rel = est / max(est), .by = sp_sz)

data_plot +
  geom_raster(
    data = d_avg,
    aes(x = X * 1000, y = Y * 1000, fill = est_rel)
  ) +
  facet_wrap(
    ~sp_sz,
    labeller = labeller(sp_sz = c(
      "cod_small" = "Cod < 35 cm",
      "had_small" = "Haddock < 30 cm",
      "whg_small" = "Whiting < 27 cm"
    ))
  ) +
  geom_sf(data = world_lcc) +
  theme(legend.position = "bottom") +
  coord_sf(
    xlim = range(preds$X) * 1000,
    ylim = range(preds$Y) * 1000,
    crs = crs_lcc,
    expand = 0
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.direction = "horizontal",
    legend.position.inside = c(0.925, 0.085),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width = unit(0.44, "cm")
  ) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Relative density"
  )

ggsave(paste0("figures/pred_avg_recent_yrs.png"), width = 21, height = 8, units = "cm")




# Average predictions (small)
d_avg_s <- preds |>
  filter(size_class == "small") |>
  filter(year_f %in% as.factor(2019:2024)) |>
  summarise(est = mean(est, na.rm = TRUE), .by = c(X, Y, sp_sz)) |>
  mutate(
    `Above 90% quantile` = ifelse(est > quantile(est, probs = 0.9), "1", "0"),
    .by = c(sp_sz)
  )

d_avg_l <- preds |>
  filter(size_class == "large") |>
  filter(year_f %in% as.factor(2019:2024)) |>
  summarise(est = mean(est, na.rm = TRUE), .by = c(X, Y, sp_sz)) |>
  mutate(
    `Above 90% quantile` = ifelse(est > quantile(est, probs = 0.9), "1", "0"),
    .by = c(sp_sz)
  )

data_plot +
  geom_raster(
    data = d_avg_s,
    aes(x = X * 1000, y = Y * 1000, fill = `Above 90% quantile`)
  ) +
  scale_fill_manual(values = c("grey95", "forestgreen")) +
  facet_wrap(~sp_sz) +
  geom_sf(data = world_lcc) +
  theme(legend.position = "bottom") +
  coord_sf(
    xlim = range(preds$X) * 1000,
    ylim = range(preds$Y) * 1000,
    crs = crs_lcc,
    expand = 0
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.direction = "horizontal",
    legend.key.width = unit(0.3, "cm"),
    legend.position.inside = c(0.96, 0.07),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5)
  ) +
  guides(fill = guide_legend(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Hotspot"
  )

ggsave(paste0("figures/pred_avg_hotspot_small_recent_yrs.png"), width = 21, height = 8.5, units = "cm")


data_plot +
  geom_raster(
    data = d_avg_l,
    aes(x = X * 1000, y = Y * 1000, fill = `Above 90% quantile`)
  ) +
  scale_fill_manual(values = c("grey95", "forestgreen")) +
  facet_wrap(~sp_sz) +
  geom_sf(data = world_lcc) +
  theme(legend.position = "bottom") +
  coord_sf(
    xlim = range(preds$X) * 1000,
    ylim = range(preds$Y) * 1000,
    crs = crs_lcc,
    expand = 0
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.direction = "horizontal",
    legend.key.width = unit(0.3, "cm"),
    legend.position.inside = c(0.96, 0.07),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5)
  ) +
  guides(fill = guide_legend(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Hotspot"
  )

ggsave(paste0("figures/pred_avg_hotspot_large_recent_yrs.png"), width = 21, height = 8.5, units = "cm")


d_avg <- bind_rows(d_avg_s, d_avg_l)

unique(d_avg$sp_sz)

data_plot +
  geom_raster(
    data = d_avg,
    aes(x = X * 1000, y = Y * 1000, fill = `Above 90% quantile`)
  ) +
  scale_fill_manual(values = c("grey95", "forestgreen")) +
  facet_wrap(~ factor(sp_sz, levels = c(
    "Cod < 35 cm", "Haddock < 30 cm", "Whiting < 27 cm",
    "Cod ≥ 35 cm", "Haddock ≥ 30 cm", "Whiting ≥ 27 cm"
  ))) +
  geom_sf(data = world_lcc) +
  theme(legend.position = "bottom") +
  coord_sf(
    xlim = range(preds$X) * 1000,
    ylim = range(preds$Y) * 1000,
    crs = crs_lcc,
    expand = 0
  ) +
  guides(fill = guide_colorbar(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  theme(
    axis.text.x = element_text(angle = 0),
    legend.direction = "horizontal",
    legend.key.width = unit(0.3, "cm"),
    legend.position.inside = c(0.96, 0.07),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5)
  ) +
  guides(fill = guide_legend(
    title.position = "top", title.hjust = 0.5,
    position = "inside"
  )) +
  labs(
    x = "Longitude", y = "Latitude",
    fill = "Hotspot"
  )

ggsave(paste0("figures/pred_avg_hotspot_recent_yrs.png"), width = 21, height = 15, units = "cm")
