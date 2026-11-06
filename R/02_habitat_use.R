# ============================================================
# 02_habitat_use.R
# Individual vs population habitat use (depth, temperature,
# light, dissolved oxygen) relative to lake-wide availability,
# trajectories, diel shifts, OVL analysis, directional shifts, PCA.
#
# Requires derived/*.rds from 01_data_prep.R.
# Start from a fresh R session - this script does not need
# anything from 03_consequences.R.
# ============================================================

source("R/00_functions.R")
library(lme4)
library(lmerTest)
library(emmeans)
library(gt)

archival_essp_full <- read_derived("archival_essp_full")
fish_traits        <- read_derived("fish_traits")
valid_fish         <- read_derived("valid_fish")
lake_depth_dist    <- read_derived("lake_depth_dist")
hypso_layers       <- read_derived("hypso_layers")
temp_vertical      <- read_derived("temp_vertical")
temp_interp        <- read_derived("temp_interp")
light_hourly       <- read_derived("light_hourly")
light_minute_day   <- read_derived("light_minute_day")
oxygen_interp      <- read_derived("oxygen_interp")
fish_do            <- read_derived("fish_do")

####Hourly fish summaries####

#Go to hourly data for temp

archival_hourly_temp <- archival_essp_full %>%
  group_by(lake, fish_id, hour_time) %>%
  summarise(
    mean_temp = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

fish_temp_means <- archival_hourly_temp %>%
  group_by(lake, fish_id) %>%
  summarise(
    fish_mean_temp = mean(mean_temp, na.rm = TRUE),
    .groups = "drop"
  )

archival_hourly_temp <- archival_hourly_temp %>%
  left_join(fish_temp_means, by = c("lake", "fish_id")) %>%
  arrange(lake, fish_mean_temp) %>%
  mutate(fish_id = factor(fish_id, levels = unique(fish_id)))

#Hourly for depth 

#Depth

archival_hourly_depth <- archival_essp_full %>%   # use the filtered full-season dataset
  group_by(lake, fish_id, hour_time) %>%
  summarise(
    mean_depth = mean(depth_corr, na.rm = TRUE),
    .groups = "drop"
  )

fish_depth_means <- archival_hourly_depth %>%
  group_by(lake, fish_id) %>%
  summarise(
    fish_mean_depth = mean(mean_depth, na.rm = TRUE),
    .groups = "drop"
  )

archival_hourly_depth <- archival_hourly_depth %>%
  left_join(fish_depth_means, by = c("lake", "fish_id")) %>%
  arrange(lake, fish_mean_depth) %>%
  mutate(fish_id = factor(fish_id, levels = unique(fish_id)))

####Density plots of habitat use####

fish_depth_dist <- archival_hourly_depth %>%
  rename(depth_m = mean_depth)

fish_depth_pop_ES <- fish_depth_dist %>%
  dplyr::filter(lake == "ES", depth_m <= 9) %>%
  dplyr::select(depth_m)

fish_depth_pop_SP <- fish_depth_dist %>%
  dplyr::filter(lake == "SP") %>%
  dplyr::select(depth_m)

####Temperature use vs availability####

temp_with_volume <- temp_vertical %>%
  left_join(
    hypso_layers %>%
      select(lake, depth_top, depth_bottom, depth_mid, layer_volume_m3),
    by = c("lake", "depth_m" = "depth_top")
  ) %>%
  filter(!is.na(layer_volume_m3))

temp_available <- temp_with_volume %>%
  mutate(temp_bin = floor(temp_c)) %>%
  group_by(lake, temp_bin) %>%
  summarise(
    total_volume_m3 = sum(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(lake) %>%
  mutate(
    prop_available = total_volume_m3 / sum(total_volume_m3)
  ) %>%
  ungroup()

temp_min <- floor(min(c(temp_with_volume$temp_c,
                        archival_essp_full$temp_c), na.rm = TRUE))

temp_max <- ceiling(max(c(temp_with_volume$temp_c,
                          archival_essp_full$temp_c), na.rm = TRUE))

temp_breaks <- seq(temp_min, temp_max + 0.5, by = 0.5)

bin_lookup <- tibble(
  temp_bin = cut(
    temp_breaks[-length(temp_breaks)],
    breaks = temp_breaks,
    include.lowest = TRUE,
    right = FALSE
  ),
  temp_mid = temp_breaks[-length(temp_breaks)] + 0.25
)

avail_temp_ES <- temp_with_volume %>%
  filter(lake == "ES", !is.na(temp_c), !is.na(layer_volume_m3)) %>%
  mutate(
    temp_bin = cut(temp_c, breaks = temp_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(temp_bin) %>%
  summarise(
    avail_volume = sum(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(bin_lookup, by = "temp_bin") %>%
  mutate(
    avail_volume = replace_na(avail_volume, 0),
    prop = avail_volume / sum(avail_volume),
    source = "availability"
  ) %>%
  arrange(temp_mid)

fish_pop_ES <- archival_essp_full %>%
  filter(lake == "ES",
         year == 2022,
         month %in% 5:9,
         !is.na(temp_c)) %>%
  mutate(
    temp_bin = cut(temp_c, breaks = temp_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(temp_bin) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  right_join(bin_lookup, by = "temp_bin") %>%
  mutate(
    n = replace_na(n, 0),
    prop = n / sum(n),
    source = "population"
  ) %>%
  arrange(temp_mid)

fish_ind_ES <- archival_essp_full %>%
  filter(lake == "ES",
         year == 2022,
         month %in% 5:9,
         !is.na(temp_c)) %>%
  mutate(
    temp_bin = cut(temp_c, breaks = temp_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, temp_bin) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  complete(fish_id, temp_bin = bin_lookup$temp_bin, fill = list(n = 0)) %>%
  left_join(bin_lookup, by = "temp_bin") %>%
  group_by(fish_id) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup() %>%
  arrange(fish_id, temp_mid)

avail_temp_SP <- temp_with_volume %>%
  filter(lake == "SP", !is.na(temp_c), !is.na(layer_volume_m3)) %>%
  mutate(
    temp_bin = cut(temp_c, breaks = temp_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(temp_bin) %>%
  summarise(
    avail_volume = sum(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(bin_lookup, by = "temp_bin") %>%
  mutate(
    avail_volume = replace_na(avail_volume, 0),
    prop = avail_volume / sum(avail_volume),
    source = "availability"
  ) %>%
  arrange(temp_mid)

fish_pop_SP <- archival_essp_full %>%
  filter(lake == "SP",
         year == 2022,
         month %in% 5:9,
         !is.na(temp_c)) %>%
  mutate(
    temp_bin = cut(temp_c, breaks = temp_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(temp_bin) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  right_join(bin_lookup, by = "temp_bin") %>%
  mutate(
    n = replace_na(n, 0),
    prop = n / sum(n),
    source = "population"
  ) %>%
  arrange(temp_mid)

fish_ind_SP <- archival_essp_full %>%
  filter(lake == "SP",
         year == 2022,
         month %in% 5:9,
         !is.na(temp_c)) %>%
  mutate(
    temp_bin = cut(temp_c, breaks = temp_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, temp_bin) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  complete(fish_id, temp_bin = bin_lookup$temp_bin, fill = list(n = 0)) %>%
  left_join(bin_lookup, by = "temp_bin") %>%
  group_by(fish_id) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup() %>%
  arrange(fish_id, temp_mid)

####Light use vs availability (daytime 08:00-20:00)####

light_with_volume <- light_hourly %>%
  left_join(
    hypso_layers %>%
      select(lake, depth_top, layer_volume_m3),
    by = c("lake", "depth_m" = "depth_top")
  ) %>%
  filter(!is.na(layer_volume_m3))

light_breaks <- 10^seq(0, 5, by = 0.2)  # 1 lux → 100,000 lux
light_breaks <- 10^seq(
  floor(log10(min(light_with_volume$light_lux, na.rm = TRUE))),
  ceiling(log10(max(light_with_volume$light_lux, na.rm = TRUE))),
  by = 0.2
)

bin_lookup_light <- tibble(
  light_bin = cut(
    light_breaks[-length(light_breaks)],
    breaks = light_breaks,
    include.lowest = TRUE,
    right = FALSE
  ),
  light_mid = 10^(log10(light_breaks[-length(light_breaks)]) + 0.1)
)

light_with_volume_day <- light_with_volume %>%
  mutate(hour = hour(datetime_hour)) %>%
  filter(hour >= 8, hour < 20)

avail_light_ES_day <- light_with_volume %>%
  mutate(hour = hour(datetime_hour)) %>%
  filter(lake == "ES", hour >= 8, hour < 20) %>%
  mutate(
    light_bin = cut(light_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(light_bin) %>%
  summarise(
    volume = sum(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(bin_lookup_light, by = "light_bin") %>%
  mutate(
    volume = tidyr::replace_na(volume, 0),
    prop = volume / sum(volume)   # ← THIS is the key step
  )

fish_light_pop_ES_day <- archival_essp_full %>%
  mutate(hour = hour(date_time)) %>%
  filter(lake == "ES",
         month %in% 5:9,
         hour >= 8, hour < 20,
         !is.na(I_pred_lux),
         I_pred_lux > 0) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  count(light_bin, name = "n") %>%
  right_join(bin_lookup_light, by = "light_bin") %>%
  mutate(
    n = tidyr::replace_na(n, 0),
    prop = n / sum(n)   # ← recalculated AFTER filtering
  )

fish_light_ind_ES_day <- archival_essp_full %>%
  mutate(hour = hour(date_time)) %>%
  filter(lake == "ES",
         month %in% 5:9,
         hour >= 8, hour < 20,
         !is.na(I_pred_lux),
         I_pred_lux > 0,
         !is.na(fish_id)) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, light_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(fish_id, light_bin = bin_lookup_light$light_bin, fill = list(n = 0)) %>%
  left_join(bin_lookup_light, by = "light_bin") %>%
  group_by(fish_id) %>%
  mutate(
    prop = n / sum(n)   # ← recalculated per fish AFTER filtering
  ) %>%
  ungroup()

avail_light_SP_day <- light_with_volume %>%
  mutate(hour = hour(datetime_hour)) %>%
  filter(lake == "SP", hour >= 8, hour < 20) %>%
  mutate(
    light_bin = cut(light_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(light_bin) %>%
  summarise(
    volume = sum(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(bin_lookup_light, by = "light_bin") %>%
  mutate(
    volume = tidyr::replace_na(volume, 0),
    prop = volume / sum(volume)   # ← THIS is the key step
  )

fish_light_pop_SP_day <- archival_essp_full %>%
  mutate(hour = hour(date_time)) %>%
  filter(lake == "SP",
         month %in% 5:9,
         hour >= 8, hour < 20,
         !is.na(I_pred_lux),
         I_pred_lux > 0) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  count(light_bin, name = "n") %>%
  right_join(bin_lookup_light, by = "light_bin") %>%
  mutate(
    n = tidyr::replace_na(n, 0),
    prop = n / sum(n)   # ← recalculated AFTER filtering
  )

fish_light_ind_SP_day <- archival_essp_full %>%
  mutate(hour = hour(date_time)) %>%
  filter(lake == "SP",
         month %in% 5:9,
         hour >= 8, hour < 20,
         !is.na(I_pred_lux),
         I_pred_lux > 0,
         !is.na(fish_id)) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, light_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(fish_id, light_bin = bin_lookup_light$light_bin, fill = list(n = 0)) %>%
  left_join(bin_lookup_light, by = "light_bin") %>%
  group_by(fish_id) %>%
  mutate(
    prop = n / sum(n)   # ← recalculated per fish AFTER filtering
  ) %>%
  ungroup()

####Dissolved oxygen use vs availability####

oxygen_with_volume <- oxygen_interp %>%
  left_join(
    hypso_layers %>%
      select(lake, depth_top, layer_volume_m3),
    by = c("lake", "depth_m" = "depth_top")
  ) %>%
  filter(!is.na(layer_volume_m3))

do_breaks <- seq(0, 12, by = 0.5)

do_lookup <- tibble(
  do_bin = cut(
    do_breaks[-length(do_breaks)],
    breaks = do_breaks,
    include.lowest = TRUE,
    right = FALSE
  ),
  do_mid = do_breaks[-length(do_breaks)] + 0.25
)

avail_do_ES <- oxygen_with_volume %>%
  filter(lake == "ES") %>%
  mutate(
    do_bin = cut(do_mgl, breaks = do_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(do_bin) %>%
  summarise(
    volume = sum(layer_volume_m3),
    .groups = "drop"
  ) %>%
  right_join(do_lookup, by = "do_bin") %>%
  mutate(
    volume = tidyr::replace_na(volume, 0),
    prop = volume / sum(volume)
  ) %>%
  arrange(do_mid)

fish_do_pop_ES <- fish_do %>%
  filter(lake == "ES") %>%
  mutate(
    do_bin = cut(do_mgl, breaks = do_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  count(do_bin, name = "n") %>%
  right_join(do_lookup, by = "do_bin") %>%
  mutate(
    n = tidyr::replace_na(n, 0),
    prop = n / sum(n)
  ) %>%
  arrange(do_mid)

fish_do_ind_ES <- fish_do %>%
  filter(lake == "ES") %>%
  mutate(
    do_bin = cut(do_mgl, breaks = do_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, do_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(fish_id, do_bin = do_lookup$do_bin, fill = list(n = 0)) %>%
  left_join(do_lookup, by = "do_bin") %>%
  group_by(fish_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(fish_id, do_mid)

avail_do_SP <- oxygen_with_volume %>%
  filter(lake == "SP") %>%
  mutate(
    do_bin = cut(do_mgl, breaks = do_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(do_bin) %>%
  summarise(
    volume = sum(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(do_lookup, by = "do_bin") %>%
  mutate(
    volume = tidyr::replace_na(volume, 0),
    prop = volume / sum(volume)
  ) %>%
  arrange(do_mid)

fish_do_pop_SP <- fish_do %>%
  filter(lake == "SP") %>%
  mutate(
    do_bin = cut(do_mgl, breaks = do_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  count(do_bin, name = "n") %>%
  right_join(do_lookup, by = "do_bin") %>%
  mutate(
    n = tidyr::replace_na(n, 0),
    prop = n / sum(n)
  ) %>%
  arrange(do_mid)

fish_do_ind_SP <- fish_do %>%
  filter(lake == "SP") %>%
  mutate(
    do_bin = cut(do_mgl, breaks = do_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, do_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(fish_id, do_bin = do_lookup$do_bin, fill = list(n = 0)) %>%
  left_join(do_lookup, by = "do_bin") %>%
  group_by(fish_id) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup() %>%
  arrange(fish_id, do_mid)

####Selected fish (Average vs Extreme) - all four habitat dimensions####

# define groups
fish_groups <- tibble::tribble(
  ~lake, ~type,      ~fish_id,
  "ES",  "Average",  "ES_03",
  "ES",  "Extreme",  "ES_70",
  "SP",  "Average",  "SP_54",
  "SP",  "Extreme",  "SP_40"
)

fish_groups <- fish_groups %>%
  mutate(
    lake = factor(lake, levels = c("ES", "SP")),
    type = factor(type, levels = c("Average", "Extreme"))
  )

depth_plot <- ggplot() +
  
  # availability
  geom_density(
    data = lake_depth_dist %>% filter(lake %in% c("ES", "SP")),
    aes(y = depth_m),
    fill = "grey70", alpha = 0.6
  ) +
  
  # population
  geom_density(
    data = bind_rows(
      fish_depth_pop_ES %>% mutate(lake = "ES"),
      fish_depth_pop_SP %>% mutate(lake = "SP")
    ),
    aes(y = depth_m),
    color = "steelblue",
    linetype = "dashed"
  ) +
  
  # individuals (ONLY selected fish)
  geom_density(
    data = fish_depth_dist %>%
      inner_join(fish_groups, by = c("lake", "fish_id")),
    aes(y = depth_m),
    color = "black",
    linewidth = 1
  ) +
  
  facet_grid(~ lake + type) +
  scale_y_reverse() +
  theme_classic() +
  labs(x = "Density", y = "Depth (m)")

temp_plot <- ggplot() +
  
  geom_area(
    data = bind_rows(
      avail_temp_ES %>% mutate(lake = "ES"),
      avail_temp_SP %>% mutate(lake = "SP")
    ),
    aes(x = temp_mid, y = prop),
    fill = "grey70", alpha = 0.6
  ) +
  
  geom_line(
    data = bind_rows(
      fish_pop_ES %>% mutate(lake = "ES"),
      fish_pop_SP %>% mutate(lake = "SP")
    ),
    aes(x = temp_mid, y = prop),
    color = "steelblue",
    linetype = "dashed"
  ) +
  
  geom_line(
    data = bind_rows(
      fish_ind_ES %>% mutate(lake = "ES"),
      fish_ind_SP %>% mutate(lake = "SP")
    ) %>%
      inner_join(fish_groups, by = c("lake", "fish_id")),
    aes(x = temp_mid, y = prop),
    color = "black",
    linewidth = 1
  ) +
  
  facet_grid(~ lake + type) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Proportion")

light_plot <- ggplot() +
  
  geom_area(
    data = bind_rows(
      avail_light_ES_day %>% mutate(lake = "ES"),
      avail_light_SP_day %>% mutate(lake = "SP")
    ),
    aes(x = light_mid, y = prop),
    fill = "grey70", alpha = 0.6
  ) +
  
  geom_line(
    data = bind_rows(
      fish_light_pop_ES_day %>% mutate(lake = "ES"),
      fish_light_pop_SP_day %>% mutate(lake = "SP")
    ),
    aes(x = light_mid, y = prop),
    color = "steelblue",
    linetype = "dashed"
  ) +
  scale_x_log10(
    breaks = 10^(-1:5),
    labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(1e-1, 1e5)) +
  geom_line(
    data = bind_rows(
      fish_light_ind_ES_day %>% mutate(lake = "ES"),
      fish_light_ind_SP_day %>% mutate(lake = "SP")
    ) %>%
      inner_join(fish_groups, by = c("lake", "fish_id")),
    aes(x = light_mid, y = prop),
    color = "black",
    linewidth = 1
  ) +
  
  scale_x_log10() +
  facet_grid(~ lake + type) +
  theme_classic() +
  labs(x = "Light (lux)", y = "Proportion")

do_plot <- ggplot() +
  
  geom_area(
    data = bind_rows(
      avail_do_ES %>% mutate(lake = "ES"),
      avail_do_SP %>% mutate(lake = "SP")
    ),
    aes(x = do_mid, y = prop),
    fill = "grey70", alpha = 0.6
  ) +
  
  geom_line(
    data = bind_rows(
      fish_do_pop_ES %>% mutate(lake = "ES"),
      fish_do_pop_SP %>% mutate(lake = "SP")
    ),
    aes(x = do_mid, y = prop),
    color = "steelblue",
    linetype = "dashed"
  ) +
  
  geom_line(
    data = bind_rows(
      fish_do_ind_ES %>% mutate(lake = "ES"),
      fish_do_ind_SP %>% mutate(lake = "SP")
    ) %>%
      inner_join(fish_groups, by = c("lake", "fish_id")),
    aes(x = do_mid, y = prop),
    color = "black",
    linewidth = 1
  ) +
  
  facet_grid(~ lake + type) +
  theme_classic() +
  labs(x = "Dissolved Oxygen (mg/L)", y = "Proportion")



final_plot <- depth_plot /
              temp_plot /
              light_plot /
              do_plot

final_plot

depth_plot
temp_plot
light_plot
do_plot

####Per-lake plots for all individuals####

# ============================================================
# DEPTH USE VS. AVAILABILITY
# ============================================================

# Escanaba Lake (ES)
ggplot() +
  geom_density(
    data = dplyr::filter(lake_depth_dist, lake == "ES"),
    aes(y = depth_m),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_density(
    data = fish_depth_pop_ES,
    aes(y = depth_m),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_density(
    data = dplyr::filter(
      fish_depth_dist, lake == "ES", depth_m <= 9
    ),
    aes(y = depth_m, group = fish_id),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  scale_y_reverse() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Density",
    y = "Depth (m)",
    title = "Escanaba Lake – Individual Depth Use vs Population Use and Habitat Availability"
  )


# Sparkling Lake (SP)
ggplot() +
  geom_density(
    data = dplyr::filter(lake_depth_dist, lake == "SP"),
    aes(y = depth_m),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_density(
    data = fish_depth_pop_SP,
    aes(y = depth_m),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_density(
    data = dplyr::filter(fish_depth_dist, lake == "SP"),
    aes(y = depth_m, group = fish_id),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  scale_y_reverse() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Density",
    y = "Depth (m)",
    title = "Sparkling Lake – Individual Depth Use vs Population Use and Habitat Availability"
  )


# ============================================================
# TEMPERATURE USE VS. AVAILABILITY
# ============================================================

# Escanaba Lake (ES)
ggplot() +
  geom_area(
    data = avail_temp_ES,
    aes(x = temp_mid, y = prop),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_line(
    data = fish_pop_ES,
    aes(x = temp_mid, y = prop),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_line(
    data = fish_ind_ES,
    aes(x = temp_mid, y = prop, group = fish_id),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Temperature (°C)",
    y = "Proportional frequency",
    title = "Escanaba Lake – Individual Temperature Use vs Population Use and Habitat Availability"
  )


# Sparkling Lake (SP)
ggplot() +
  geom_area(
    data = avail_temp_SP,
    aes(x = temp_mid, y = prop),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_line(
    data = fish_pop_SP,
    aes(x = temp_mid, y = prop),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_line(
    data = fish_ind_SP,
    aes(x = temp_mid, y = prop, group = fish_id),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Temperature (°C)",
    y = "Proportional frequency",
    title = "Sparkling Lake – Individual Temperature Use vs Population Use and Habitat Availability"
  )


# ============================================================
# LIGHT USE VS. AVAILABILITY
# ============================================================


# Escanaba Lake daytime light use vs. availability (08:00–20:00)
ggplot() +
  # Habitat availability
  geom_area(
    data = avail_light_ES_day,
    aes(x = light_mid, y = prop),
    fill = "grey70",
    color = "grey50",
    alpha = 0.6
  ) +
  
  # Population-level use
  geom_line(
    data = fish_light_pop_ES_day,
    aes(x = light_mid, y = prop),
    color = "steelblue",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # Individual fish use
  geom_line(
    data = fish_light_ind_ES_day,
    aes(x = light_mid, y = prop, group = fish_id),
    color = "black",
    linewidth = 1
  ) +
  
  scale_x_log10(
    breaks = 10^(-1:5),
    labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(1e-1, 1e5)) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Light (lux, log scale)",
    y = "Proportional frequency",
    title = "Escanaba Lake – Individual Daytime Light Use vs Population Use and Habitat Availability"
  )
# Sparkling Lake (SP)
ggplot() +
  geom_area(
    data = avail_light_SP_day,
    aes(x = light_mid, y = prop),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_line(
    data = fish_light_pop_SP_day,
    aes(x = light_mid, y = prop),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_line(
    data = fish_light_ind_SP_day,
    aes(x = light_mid, y = prop, group = fish_id),
    color = "black", linewidth = 1
  ) +
  scale_x_log10(
    breaks = 10^(-1:5),
    labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(1e-1, 1e5)) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Light (lux, log scale)",
    y = "Proportional frequency",
    title = "Sparkling Lake – Individual Daytime Light Use vs Population Use and Habitat Availability"
  )


# ============================================================
# DISSOLVED OXYGEN USE VS. AVAILABILITY
# ============================================================

# Escanaba Lake (ES)
ggplot() +
  geom_area(
    data = avail_do_ES,
    aes(x = do_mid, y = prop),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_line(
    data = fish_do_pop_ES,
    aes(x = do_mid, y = prop),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_line(
    data = fish_do_ind_ES,
    aes(x = do_mid, y = prop, group = fish_id),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Dissolved oxygen (mg/L)",
    y = "Proportional frequency",
    title = "Escanaba Lake – Individual Dissolved Oxygen Use vs Population Use and Habitat Availability"
  )


# Sparkling Lake (SP)
ggplot() +
  geom_area(
    data = avail_do_SP,
    aes(x = do_mid, y = prop),
    fill = "grey70", color = "grey50", alpha = 0.6
  ) +
  geom_line(
    data = fish_do_pop_SP,
    aes(x = do_mid, y = prop),
    color = "steelblue", linewidth = 0.8, linetype = "dashed"
  ) +
  geom_line(
    data = fish_do_ind_SP,
    aes(x = do_mid, y = prop, group = fish_id),
    color = "black", linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Dissolved oxygen (mg/L)",
    y = "Proportional frequency",
    title = "Sparkling Lake – Individual Dissolved Oxygen Use vs Population Use and Habitat Availability"
  )

####Trajectory Plots####

target_day <- as.Date("2022-06-29")

# Light (ES only)
light_day <- light_minute_day %>%
  filter(lake == "ES") %>%
  filter(date_time >= target_day,
         date_time < target_day + 1)

light_day <- light_day %>%
  mutate(
    log10_light = log10(I_pred_lux),
    log10_light = ifelse(log10_light < 0, 0, log10_light)
  )

# Oxygen

oxygen_day <- oxygen_interp %>%
  filter(lake == "ES") %>%
  mutate(datetime = with_tz(datetime, "UTC")) %>%
  filter(datetime >= as.POSIXct("2022-06-29 00:00:00", tz = "UTC"),
         datetime <  as.POSIXct("2022-06-30 00:00:00", tz = "UTC"))

oxygen_full <- oxygen_day %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  complete(
    datetime = seq(min(datetime), max(datetime), by = "1 min"),
    depth_m
  ) %>%
  arrange(depth_m, datetime)

oxygen_full <- oxygen_full %>%
  mutate(datetime = with_tz(datetime, "America/Chicago"))


# Temperature (hourly)
temp_day <- temp_interp %>%
  filter(lake == "ES",
         datetime_hour >= target_day,
         datetime_hour < target_day + 1)



oxygen_full <- oxygen_full %>%
  group_by(depth_m) %>%
  mutate(
    do_interp = na.approx(do_mgl, x = datetime, na.rm = FALSE)
  ) %>%
  ungroup()

# Fish (ES_01)
fish_day <- archival_essp_full %>%
  filter(fish_id == "ES_01",
         date_time >= target_day,
         date_time < target_day + 1)
#Rename to not break my shit 
light_day <- light_day %>%
  rename(datetime = date_time)

temp_day <- temp_day %>%
  rename(datetime = datetime_hour)

fish_day <- fish_day %>%
  mutate(depth = depth_corr)

tz_use <- "UTC"

light_day <- light_day %>%
  mutate(datetime = with_tz(datetime, tz_use))

temp_day <- temp_day %>%
  mutate(datetime = with_tz(datetime, tz_use))

fish_day <- fish_day %>%
  mutate(date_time = with_tz(date_time, tz_use))



ggplot() +
  geom_tile(data = light_day,
            aes(x = datetime,
                y = -depth_m,
                fill = log10_light)) +
  
  geom_path(data = fish_day,
            aes(x = date_time,
                y = -depth),
            color = "black",
            size = 1) +
  
  scale_fill_viridis_c(
    option = "viridis",
  ) +
  
  theme_classic()

ggplot() +
  geom_tile(data = oxygen_full,
            aes(x = datetime,
                y = -depth_m,
                fill = do_interp)) +
  
  geom_path(data = fish_day,
            aes(x = date_time,
                y = -depth),
            color = "black",
            size = 1) +
  
  scale_fill_viridis_c(option = "inferno") +
  theme_classic()

ggplot() +
  geom_tile(data = temp_day,
            aes(x = datetime,
                y = -depth_m,
                fill = temp_c)) +
  
  geom_path(data = fish_day,
            aes(x = date_time,
                y = -depth),
            color = "black",
            size = 1) +
  
  scale_fill_viridis_c(option = "plasma") +
  theme_classic()


#Trajectory for sparkles 

target_day <- as.Date("2022-06-29")


# Light (SP only)
light_day <- light_minute_day %>%
  filter(lake == "SP") %>%
  filter(date_time >= target_day,
         date_time < target_day + 1)

light_day <- light_day %>%
  mutate(
    log10_light = log10(I_pred_lux),
    log10_light = ifelse(log10_light < 0, 0, log10_light)
  )

# Oxygen

oxygen_day <- oxygen_interp %>%
  filter(lake == "SP") %>%
  mutate(datetime = with_tz(datetime, "UTC")) %>%
  filter(datetime >= as.POSIXct("2022-06-29 00:00:00", tz = "UTC"),
         datetime <  as.POSIXct("2022-06-30 00:00:00", tz = "UTC"))

oxygen_full <- oxygen_day %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  complete(
    datetime = seq(min(datetime), max(datetime), by = "1 min"),
    depth_m
  ) %>%
  arrange(depth_m, datetime)

oxygen_full <- oxygen_full %>%
  mutate(datetime = with_tz(datetime, "UTC"))


# Temperature (hourly)
temp_day <- temp_interp %>%
  filter(lake == "SP",
         datetime_hour >= target_day,
         datetime_hour < target_day + 1)



oxygen_full <- oxygen_full %>%
  group_by(depth_m) %>%
  mutate(
    do_interp = na.approx(do_mgl, x = datetime, na.rm = FALSE)
  ) %>%
  ungroup()
fish_day <- archival_essp_full %>%
  filter(fish_id == "SP_40",
         date_time >= target_day,
         date_time < target_day + 1)

fish_day_2 <- archival_essp_full %>%
  filter(fish_id == "SP_54",
         date_time >= target_day,
         date_time < target_day + 1)
#Rename to not break my shit 
light_day <- light_day %>%
  rename(datetime = date_time)

temp_day <- temp_day %>%
  rename(datetime = datetime_hour)

fish_day <- fish_day %>%
  mutate(depth = depth_corr)

fish_day_2 <- fish_day_2 %>%
  mutate(depth = depth_corr)

tz_use <- "UTC"

light_day <- light_day %>%
  mutate(datetime = with_tz(datetime, tz_use))

temp_day <- temp_day %>%
  mutate(datetime = with_tz(datetime, tz_use))

fish_day <- fish_day %>%
  mutate(date_time = with_tz(date_time, tz_use))

fish_day_2 <- fish_day_2 %>%
  mutate(date_time = with_tz(date_time, tz_use))



ggplot() +
  geom_tile(data = light_day,
            aes(x = datetime,
                y = -depth_m,
                fill = log10_light)) +
  
  geom_path(data = fish_day,
            aes(x = date_time,
                y = -depth),
            color = "white",
            size = 1) +
  
  geom_path(data = fish_day_2,
            aes(x = date_time,
                y = -depth),
            color = "black",
            size = 1) +
  
  scale_fill_viridis_c(
    option = "viridis",
    name = "Log10 Light (Lux)"
  ) +
  labs(x = "Date & Time",
       y = "Depth (m)") +
  theme_classic()+
  
  theme(axis.title = element_text(size = 20), 
                axis.text = element_text(size = 16), 
                legend.title = element_text(size = 20),
                legend.text = element_text(size = 16))

ggplot() +
  geom_tile(data = oxygen_full,
            aes(x = datetime,
                y = -depth_m,
                fill = do_interp)) +
  
  geom_path(data = fish_day,
            aes(x = date_time,
                y = -depth),
            color = "white",
            size = 1) +
  geom_path(data = fish_day_2,
            aes(x = date_time,
                y = -depth),
            color = "black",
            size = 1) +
  
  scale_fill_viridis_c(option = "inferno",
                       name="Dissolved Oxygen (mg/L)") +
  labs(x = "Date & Time",
       y = "Depth (m)") +
  theme_classic()+
  
  theme(axis.title = element_text(size = 20), 
        axis.text = element_text(size = 16), 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))


ggplot() +
  geom_tile(data = temp_day,
            aes(x = datetime,
                y = -depth_m,
                fill = temp_c)) +
  
  geom_path(data = fish_day,
            aes(x = date_time,
                y = -depth),
            color = "white",
            size = 1) +
  geom_path(data = fish_day_2,
            aes(x = date_time,
                y = -depth),
            color = "black",
            size = 1) +
  
  scale_fill_viridis_c(option = "plasma",
                       name="Temperature (°C)") +
  labs(x = "Date & Time",
       y = "Depth (m)")  +
  theme_classic()+
  theme(axis.title = element_text(size = 20), 
      axis.text = element_text(size = 16), 
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 16))

####Diel Shifts####

diel_df <- archival_hourly_depth %>%
  inner_join(archival_hourly_temp,
             by = c("lake", "fish_id", "hour_time")) %>%
  select(lake, fish_id, hour_time, mean_depth, mean_temp)

diel_df <- diel_df %>%
  mutate(
    hour = hour(hour_time),
    
    # New diel definition
    diel = case_when(
      hour >= 8 & hour < 20 ~ "day",
      TRUE ~ "night"
    ),
    
    # Biological day: assign early morning (before 8am) to previous day
    diel_date = as.Date(hour_time) - ifelse(hour < 8, 1, 0)
  )

daily_diel <- diel_df %>%
  group_by(lake, fish_id, diel_date, diel) %>%
  summarize(
    mean_depth = mean(mean_depth, na.rm = TRUE),
    mean_temp  = mean(mean_temp, na.rm = TRUE),
    n_obs = n(),  # useful for filtering later
    .groups = "drop"
  )

daily_wide <- daily_diel %>%
  pivot_wider(
    names_from = diel,
    values_from = c(mean_depth, mean_temp, n_obs)
  )

daily_wide <- daily_wide %>%
  filter(!is.na(mean_depth_day),
         !is.na(mean_depth_night),
         !is.na(mean_temp_day),
         !is.na(mean_temp_night))

daily_shifts <- daily_wide %>%
  mutate(
    delta_depth = mean_depth_night - mean_depth_day,
    delta_temp  = mean_temp_night  - mean_temp_day
  )

daily_shifts <- daily_shifts %>%
  mutate(
    shift_type = case_when(
      delta_temp > 0 & delta_depth < 0 ~ "warmer_shallower",
      delta_temp > 0 & delta_depth > 0 ~ "warmer_deeper",
      delta_temp < 0 & delta_depth < 0 ~ "colder_shallower",
      delta_temp < 0 & delta_depth > 0 ~ "colder_deeper",
      TRUE ~ "no_change"
    )
  )

fish_summary <- daily_shifts %>%
  group_by(lake, fish_id) %>%
  summarize(
    mean_delta_depth = mean(delta_depth, na.rm = TRUE),
    mean_delta_temp  = mean(delta_temp, na.rm = TRUE),
    .groups = "drop"
  )

fish_summary_props <- daily_shifts %>%
  group_by(lake, fish_id, shift_type) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(lake, fish_id) %>%
  mutate(prop = n / sum(n))

ggplot(daily_shifts,
       aes(x = delta_depth, y = delta_temp, color = shift_type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  facet_wrap(~lake) +
  labs(x = "Δ Depth (night - day)",
       y = "Δ Temperature (night - day)",
       color = "Shift type") + 
  stat_density_2d(color = "black", alpha = 0.7, bins=6)

daily_shifts <- daily_shifts %>%
  arrange(fish_id, diel_date) %>%
  group_by(fish_id) %>%
  mutate(
    prev_shift = lag(shift_type),
    switched = shift_type != prev_shift
  )

#transitions 

transitions <- daily_shifts %>%
  arrange(lake, fish_id, diel_date) %>%
  group_by(lake, fish_id) %>%
  mutate(
    from = lag(shift_type),
    to   = shift_type
  ) %>%
  ungroup() %>%
  filter(!is.na(from))

transition_counts <- transitions %>%
  count(lake, from, to)

transition_probs <- transition_counts %>%
  group_by(lake, from) %>%
  mutate(prob = n / sum(n)) %>%
  ungroup()

ggplot(transition_probs,
       aes(x = from, y = to, fill = prob)) +
  geom_tile() +
  geom_text(aes(label = round(prob, 2))) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  facet_wrap(~lake) +
  labs(x = "Previous day strategy",
       y = "Next day strategy",
       fill = "Transition\nprobability") +
  theme_minimal()

#Probably not going to use this 

#Fish means 
dominant_strategy <- daily_shifts %>%
  group_by(lake, fish_id, shift_type) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(lake, fish_id) %>%
  slice_max(n, n = 1, with_ties = FALSE)

fish_summary <- fish_summary %>%
  left_join(dominant_strategy, by = c("lake", "fish_id"))

ggplot(fish_summary,
       aes(x = mean_delta_depth, y = mean_delta_temp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  facet_wrap(~lake) +
  labs(x = "Mean Δ Depth (night - day)",
       y = "Mean Δ Temperature (night - day)",
       color = "Dominant strategy") +
  theme_classic()


fish_means_diel <- daily_shifts %>%
  group_by(lake, fish_id) %>%
  summarize(
    mean_delta_depth = mean(delta_depth, na.rm = TRUE),
    mean_delta_temp  = mean(delta_temp, na.rm = TRUE),
    .groups = "drop"
  )

ggplot() +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  # Population cloud (all days)
  geom_point(data = daily_shifts,
             aes(x = delta_depth, y = delta_temp),
             color = "gray70", alpha = 0.4, size = 1.2) +
  
  # Individual ellipses
  stat_ellipse(data = daily_shifts,
               aes(x = delta_depth, y = delta_temp,
                   group = fish_id),
               color = "black",
               linewidth = 0.8) +
  
  geom_point(data = fish_means_diel,
             aes(x = mean_delta_depth,
                 y = mean_delta_temp,
                 color = fish_id),
             size = 3) +
  
  facet_wrap(~fish_id) +
  theme_classic() +
  labs(x = "Δ Depth (night - day)",
       y = "Δ Temperature (night - day)",
       color = "Fish ID") +
  
  theme(legend.position = "none")

#Do the monthly dance 

daily_shifts <- daily_shifts %>%
  mutate(month = lubridate::month(diel_date))

ggplot(daily_shifts,
       aes(x = delta_depth, y = delta_temp, color = shift_type)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_point(alpha = 0.5) +
  
  facet_grid(month~lake) +
  
  labs(x = "Δ Depth (night - day)",
       y = "Δ Temperature (night - day)",
       color = "Shift type") + theme(
         strip.background = element_blank(),
         panel.spacing = unit(1, "lines")
       )+ 
  stat_density_2d(color = "black", alpha = 0.7, bins=3)

#Step length? 

daily_steps <- daily_shifts %>%
  arrange(lake, fish_id, diel_date) %>%
  group_by(lake, fish_id) %>%
  mutate(
    prev_depth = lag(delta_depth),
    prev_temp  = lag(delta_temp),
    
    step_dist = sqrt(
      (delta_depth - prev_depth)^2 +
        (delta_temp  - prev_temp)^2
    )
  ) %>%
  ungroup()


ggplot(daily_steps,
       aes(x = diel_date, y = step_dist)) +
  geom_line() +
  facet_wrap(~fish_id) +
  theme_classic() +
  labs(x = "Date",
       y = "Behavioral step distance")

#Add background shading 

daily_steps_bg <- daily_steps %>%
  mutate(
    xmin = diel_date - 0.5,
    xmax = diel_date + 0.5
  )

ggplot() +
  # Background shading by shift type
  geom_rect(data = daily_steps_bg,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf,
                fill = shift_type),
            alpha = 0.2) +
  
  # Line on top
  geom_line(data = daily_steps,
            aes(x = diel_date, y = step_dist),
            color = "black") +
  
  facet_wrap(~fish_id) +
  theme_classic() +
  labs(x = "Date",
       y = "Behavioral step distance",
       fill = "Shift type")


#Now try adding counts: 

shift_counts <- daily_steps %>%
  group_by(fish_id, shift_type) %>%
  summarise(days = n(), .groups = "drop")

label_pos <- daily_steps %>%
  group_by(fish_id) %>%
  summarise(
    x = min(diel_date, na.rm = TRUE),
    y = max(step_dist, na.rm = TRUE),
    .groups = "drop"
  )

shift_labels <- daily_steps %>%
  group_by(fish_id, shift_type) %>%
  summarise(days = n(), .groups = "drop") %>%
  arrange(fish_id, shift_type) %>%
  group_by(fish_id) %>%
  summarise(
    label = paste0(shift_type, ": ", days, collapse = "\n"),
    .groups = "drop"
  )


shift_labels <- left_join(shift_labels, label_pos, by = "fish_id")

ggplot() +
  geom_rect(data = daily_steps_bg,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = 15,
                fill = shift_type),
            alpha = 0.2) +
  
  geom_line(data = daily_steps,
            aes(x = diel_date, y = step_dist),
            color = "black") +
  
  geom_text(data = shift_labels,
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 0.5,
            size = 3) +
  
  facet_wrap(~fish_id) +
  theme_classic()



##

ggplot(daily_shifts,
       aes(x = delta_depth, y = delta_temp)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_path(aes(group = fish_id),
            color = "gray60",
            alpha = 0.6) +
  
  geom_point(aes(color = fish_id), alpha = 0.6) +
  
  facet_wrap(~fish_id) +
  theme_classic()

threshold <- quantile(daily_steps$step_dist, 0.75, na.rm = TRUE)

daily_steps <- daily_steps %>%
  mutate(
    abrupt = step_dist > threshold
  )

fish_step_summary <- daily_steps %>%
  group_by(lake, fish_id) %>%
  summarize(
    mean_step = mean(step_dist, na.rm = TRUE),
    max_step  = max(step_dist, na.rm = TRUE),
    prop_abrupt = mean(step_dist > quantile(step_dist, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

#Link to quad shifting? 

# Grouped by fish so the first day of one fish is not compared
# with the last day of the previous fish
daily_steps <- daily_steps %>%
  arrange(lake, fish_id, diel_date) %>%
  group_by(lake, fish_id) %>%
  mutate(
    prev_shift = lag(shift_type),
    switched = shift_type != prev_shift
  ) %>%
  ungroup()

ggplot(daily_steps,
       aes(x = switched, y = step_dist)) +
  geom_boxplot() +
  theme_classic()

####Density OVL (leave-one-fish-out)####

# ------------------------------------------------------------
# 1. BUILD STANDARDIZED INDIVIDUAL HABITAT-USE DATA
# ------------------------------------------------------------

# Depth
ovl_depth <- fish_depth_dist %>%
  filter(
    lake %in% c("ES", "SP"),
    !is.na(depth_m),
    lake != "ES" | depth_m <= 9
  ) %>%
  select(
    fish_id,
    lake,
    value = depth_m
  ) %>%
  mutate(habitat = "Depth")


# Temperature
ovl_temp <- archival_essp_full %>%
  filter(
    lake %in% c("ES", "SP"),
    year == 2022,
    month %in% 5:9,
    !is.na(temp_c)
  ) %>%
  select(
    fish_id,
    lake,
    value = temp_c
  ) %>%
  mutate(habitat = "Temperature")


# Light
# Following your existing final light plots:
# daytime only, 08:00-20:00, positive light values
ovl_light <- archival_essp_full %>%
  mutate(hour = hour(date_time)) %>%
  filter(
    lake %in% c("ES", "SP"),
    month %in% 5:9,
    hour >= 8,
    hour < 20,
    !is.na(I_pred_lux),
    I_pred_lux > 0,
    !is.na(fish_id)
  ) %>%
  select(
    fish_id,
    lake,
    value = I_pred_lux
  ) %>%
  mutate(habitat = "Light")


# Dissolved oxygen
# fish_do is your fish-level DO dataset created by
# interpolating the measured oxygen profiles through time
ovl_do <- fish_do %>%
  filter(
    lake %in% c("ES", "SP"),
    !is.na(do_mgl)
  ) %>%
  select(
    fish_id,
    lake,
    value = do_mgl
  ) %>%
  mutate(habitat = "Oxygen")


# Combine all four
habitat_use <- bind_rows(
  ovl_temp,
  ovl_depth,
  ovl_light,
  ovl_do
) %>%
  filter(
    fish_id %in% valid_fish$fish_id
  ) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Temperature",
        "Depth",
        "Light",
        "Oxygen"
      )
    )
  )


# ------------------------------------------------------------
# 2. FUNCTION TO CALCULATE LEAVE-ONE-FISH-OUT OVL
# ------------------------------------------------------------

calculate_ovl <- function(df,
                          n_grid = 500,
                          n_min = 20) {
  
  # Remove missing/non-finite observations
  df <- df %>%
    filter(is.finite(value))
  
  # Common range for all fish in this lake/habitat
  x_min <- min(df$value, na.rm = TRUE)
  x_max <- max(df$value, na.rm = TRUE)
  
  # Common x-axis
  x_grid <- seq(
    x_min,
    x_max,
    length.out = n_grid
  )
  
  # ----------------------------------------------------------
  # Individual kernel densities
  # ----------------------------------------------------------
  
  ind_dens <- df %>%
    group_by(fish_id, lake) %>%
    group_modify(~ {
      
      values <- .x$value
      
      # Require enough observations and variation
      if (
        length(values) < n_min ||
        length(unique(values)) < 2
      ) {
        return(tibble())
      }
      
      d <- density(
        values,
        from = x_min,
        to = x_max,
        n = n_grid
      )
      
      # Explicitly normalize area to 1
      area <- sum(
        diff(d$x) *
          (
            head(d$y, -1) +
              tail(d$y, -1)
          ) / 2
      )
      
      tibble(
        x = d$x,
        density = d$y / area
      )
      
    }) %>%
    ungroup()
  
  
  # ----------------------------------------------------------
  # Population density
  #
  # Every fish receives equal weight.
  # ----------------------------------------------------------
  
  pop_dens <- ind_dens %>%
    group_by(lake, x) %>%
    summarise(
      pop_density = mean(density),
      n_fish = n(),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # Leave-one-fish-out population density
  #
  # f[-i](x) =
  # [N * f_pop(x) - f_i(x)] / (N - 1)
  #
  # This prevents each fish from contributing to its own
  # reference population distribution.
  # ----------------------------------------------------------
  
  loo_dens <- ind_dens %>%
    left_join(
      pop_dens,
      by = c("lake", "x")
    ) %>%
    mutate(
      loo_population_density =
        if_else(
          n_fish > 1,
          (n_fish * pop_density - density) /
            (n_fish - 1),
          NA_real_
        ),
      
      # Guard against tiny negative values caused by
      # floating-point precision
      loo_population_density =
        pmax(loo_population_density, 0)
    )
  
  
  # ----------------------------------------------------------
  # Calculate OVL
  #
  # OVL = integral of the minimum of the two densities
  # ----------------------------------------------------------
  
  ovl <- loo_dens %>%
    group_by(lake, fish_id) %>%
    arrange(x) %>%
    summarise(
      
      OVL = sum(
        diff(x) *
          (
            pmin(
              head(density, -1),
              head(loo_population_density, -1)
            ) +
              pmin(
                tail(density, -1),
                tail(loo_population_density, -1)
              )
          ) / 2,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  ovl
}


# ------------------------------------------------------------
# 3. CALCULATE OVL FOR ALL FOUR HABITAT DIMENSIONS
# ------------------------------------------------------------

ovl_results <- habitat_use %>%
  group_by(habitat) %>%
  group_modify(
    ~ calculate_ovl(.x)
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 4. ADD FISH INFORMATION
# ------------------------------------------------------------

ovl_results <- ovl_results %>%
  left_join(
    valid_fish %>%
      select(
        fish_id,
        lake,
        tl_pre
      ),
    by = c("fish_id", "lake")
  ) %>%
  arrange(
    lake,
    fish_id,
    habitat
  )


# ------------------------------------------------------------
# 5. CHECK THE RESULTS
# ------------------------------------------------------------

ovl_results

summary(ovl_results$OVL)


# Check for impossible OVL values
ovl_results %>%
  filter(
    OVL < 0 |
      OVL > 1 |
      is.na(OVL)
  )


# ------------------------------------------------------------
# 6. SUMMARY BY LAKE AND HABITAT
# ------------------------------------------------------------

ovl_summary <- ovl_results %>%
  group_by(
    lake,
    habitat
  ) %>%
  summarise(
    n = n(),
    mean_OVL = mean(OVL, na.rm = TRUE),
    sd_OVL = sd(OVL, na.rm = TRUE),
    median_OVL = median(OVL, na.rm = TRUE),
    min_OVL = min(OVL, na.rm = TRUE),
    max_OVL = max(OVL, na.rm = TRUE),
    .groups = "drop"
  )

ovl_summary


# ------------------------------------------------------------
# 7. PLOT THE FINAL OVL VALUES
# ------------------------------------------------------------

ggplot(
  ovl_results,
  aes(
    x = habitat,
    y = OVL
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.2
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    alpha = 0.55,
    size = 2
  ) +
  facet_wrap(~lake) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Overlap coefficient (OVL)",
    title = "Individual habitat-use overlap with population use"
  )

ovl_results
ovl_summary

# ------------------------------------------------------------
# Prepare OVL for modeling
# ------------------------------------------------------------

ovl_model <- ovl_results %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c("Temperature", "Depth", "Light", "Oxygen")
    ),
    
    # Keep values away from exactly 0 and 1
    OVL_adj = pmin(
      pmax(OVL, 0.001),
      0.999
    ),
    
    logit_OVL = qlogis(OVL_adj)
  )


# ------------------------------------------------------------
# Mixed model
# ------------------------------------------------------------

ovl_mod <- lmer(
  logit_OVL ~ habitat * lake + (1 | fish_id),
  data = ovl_model
)

summary(ovl_mod)
anova(ovl_mod)

# Variance components
VarCorr(ovl_mod)

# Extract fish and residual variance
vc <- as.data.frame(VarCorr(ovl_mod))

fish_var <- vc$vcov[vc$grp == "fish_id"]
resid_var <- vc$vcov[vc$grp == "Residual"]

# Intraclass correlation coefficient
ICC <- fish_var / (fish_var + resid_var)

ICC

ggplot(
  ovl_results,
  aes(
    x = habitat,
    y = OVL,
    group = fish_id
  )
) +
  geom_line(
    alpha = 0.25
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(~lake) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Overlap coefficient (OVL)",
    title = "Consistency in individual habitat-use overlap"
  )

# ============================================================
# ESTIMATED OVL BY HABITAT × LAKE
# ============================================================


ovl_emm <- emmeans(
  ovl_mod,
  ~ habitat | lake
)

ovl_emm


# Back-transform from logit scale to OVL scale
ovl_emm_df <- as.data.frame(ovl_emm) %>%
  mutate(
    OVL_est = plogis(emmean),
    OVL_lower = plogis(lower.CL),
    OVL_upper = plogis(upper.CL)
  )

ovl_emm_df


# Pairwise comparisons among habitats within each lake
ovl_pairs <- pairs(
  ovl_emm,
  adjust = "tukey"
)

ovl_pairs

# ============================================================
# INDIVIDUAL OVL ACROSS HABITAT DIMENSIONS
# ============================================================

ovl_wide <- ovl_results %>%
  select(
    fish_id,
    lake,
    habitat,
    OVL
  ) %>%
  pivot_wider(
    names_from = habitat,
    values_from = OVL
  )

ovl_wide

ggplot(
  ovl_results,
  aes(
    x = habitat,
    y = OVL,
    group = fish_id
  )
) +
  geom_line(
    alpha = 0.3
  ) +
  geom_point(
    size = 2
  ) +
  facet_wrap(~lake) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Overlap coefficient (OVL)",
    title = "Individual consistency in habitat-use overlap"
  )


# ============================================================
# CORRELATIONS AMONG HABITAT-SPECIFIC OVL
# ============================================================

# ES
ovl_wide %>%
  filter(lake == "ES") %>%
  select(Temperature, Depth, Light, Oxygen) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

# SP
ovl_wide %>%
  filter(lake == "SP") %>%
  select(Temperature, Depth, Light, Oxygen) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

# Function to calculate ICC from the model
get_icc <- function(model) {
  
  vc <- as.data.frame(VarCorr(model))
  
  fish_var <- vc$vcov[vc$grp == "fish_id"]
  resid_var <- vc$vcov[vc$grp == "Residual"]
  
  fish_var / (fish_var + resid_var)
}

# Observed ICC
icc_obs <- get_icc(ovl_mod)

icc_obs

set.seed(369)

icc_boot <- bootMer(
  ovl_mod,
  FUN = get_icc,
  nsim = 1000,
  type = "parametric",
  use.u = TRUE
)

icc_ci <- quantile(
  icc_boot$t,
  probs = c(0.025, 0.975),
  na.rm = TRUE
)

icc_results <- tibble(
  ICC = icc_obs,
  CI_low = icc_ci[1],
  CI_high = icc_ci[2]
)

icc_results

cor_results <- ovl_wide %>%
  group_by(lake) %>%
  group_modify(~ {
    
    vars <- c("Temperature", "Depth", "Light", "Oxygen")
    
    pairs <- combn(vars, 2, simplify = FALSE)
    
    map_dfr(pairs, function(x) {
      
      test <- cor.test(
        .x[[x[1]]],
        .x[[x[2]]],
        method = "spearman",
        exact = FALSE
      )
      
      tibble(
        habitat_1 = x[1],
        habitat_2 = x[2],
        rho = unname(test$estimate),
        p = test$p.value
      )
    })
  }) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH")
  )

cor_results

cor_plot_data <- ovl_wide %>%
  select(
    fish_id, lake,
    Temperature, Depth, Light, Oxygen
  )

ggplot() +
  
  geom_point(
    data = cor_plot_data %>% filter(lake == "ES"),
    aes(x = Depth, y = Light),
    size = 2.5,
    alpha = 0.7
  ) +
  
  geom_smooth(
    data = cor_plot_data %>% filter(lake == "ES"),
    aes(x = Depth, y = Light),
    method = "lm",
    se = TRUE
  ) +
  
  labs(
    x = "Depth OVL",
    y = "Light OVL"
  ) +
  
  theme_classic()


#Making big panel fig 


# Long-format data for pairwise comparisons
ovl_long <- ovl_wide %>%
  pivot_longer(
    cols = c(Temperature, Depth, Light, Oxygen),
    names_to = "habitat",
    values_to = "OVL"
  )

# All six unique habitat pairs
habitat_pairs <- combn(
  c("Temperature", "Depth", "Light", "Oxygen"),
  2,
  simplify = FALSE
)

# Function to calculate Spearman correlations
cor_labels <- map_dfr(
  c("ES", "SP"),
  function(lake_name) {
    
    map_dfr(
      habitat_pairs,
      function(pair) {
        
        dat <- ovl_wide %>%
          filter(lake == lake_name) %>%
          select(
            fish_id,
            all_of(pair)
          ) %>%
          drop_na()
        
        test <- cor.test(
          dat[[pair[1]]],
          dat[[pair[2]]],
          method = "spearman",
          exact = FALSE
        )
        
        tibble(
          lake = lake_name,
          habitat_x = pair[1],
          habitat_y = pair[2],
          rho = unname(test$estimate),
          p = test$p.value,
          p_adj = p.adjust(
            c(
              # correction is applied below across all comparisons
            ),
            method = "BH"
          )
        )
      }
    )
  }
)



# ---------------------------------------------------------
# 1. Correlation statistics for plotting
# ---------------------------------------------------------

cor_labels <- cor_results %>%
  mutate(
    habitat_x = habitat_1,
    habitat_y = habitat_2,
    relationship = paste(habitat_x, "vs.", habitat_y),
    label = paste0(
      "\u03c1 = ", sprintf("%.2f", rho),
      "\nP[adj] = ",
      format.pval(p_adj, digits = 2, eps = 0.001)
    )
  ) %>%
  select(
    lake,
    habitat_x,
    habitat_y,
    relationship,
    rho,
    p_adj,
    label
  )


# ---------------------------------------------------------
# 2. Create all six habitat-pair datasets
# ---------------------------------------------------------

habitat_pairs <- combn(
  c("Temperature", "Depth", "Light", "Oxygen"),
  2,
  simplify = FALSE
)

pair_plot_data <- map_dfr(
  habitat_pairs,
  function(pair) {
    
    ovl_wide %>%
      select(
        fish_id,
        lake,
        all_of(pair)
      ) %>%
      rename(
        x = !!pair[1],
        y = !!pair[2]
      ) %>%
      mutate(
        habitat_x = pair[1],
        habitat_y = pair[2],
        relationship = paste(pair[1], "vs.", pair[2])
      )
  }
)


# ---------------------------------------------------------
# 3. Add correlation statistics
# ---------------------------------------------------------

pair_plot_data <- pair_plot_data %>%
  left_join(
    cor_labels,
    by = c(
      "lake",
      "habitat_x",
      "habitat_y",
      "relationship"
    )
  )


# Check that the statistics are present
pair_plot_data %>%
  select(
    lake,
    relationship,
    rho,
    p_adj,
    label
  ) %>%
  distinct()

ggplot(
  pair_plot_data,
  aes(x = x, y = y)
) +
  
  geom_point(
    aes(color = lake),
    size = 2.5,
    alpha = 0.75
  ) +
  
  geom_text(
    data = pair_plot_data %>%
      distinct(
        lake,
        relationship,
        rho,
        p_adj,
        label
      ) %>%
      mutate(
        x = 0.42,
        y = 0.98
      ),
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.2
  ) +
  
  facet_grid(
    lake ~ relationship
  ) +
  
  scale_color_manual(
    values = c(
      "ES" = "chocolate",
      "SP" = "dodgerblue4"
    )
  ) +
  
  scale_x_continuous(
    limits = c(0.4, 1),
    breaks = seq(0.4, 1, 0.2)
  ) +
  
  scale_y_continuous(
    limits = c(0.4, 1),
    breaks = seq(0.4, 1, 0.2)
  ) +
  
  labs(
    x = "OVL",
    y = "OVL"
  ) +
  
  theme_classic() +
  
  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "black"
    ),
    strip.text = element_text(size = 9),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

# Four relationships we want to highlight
key_pairs <- tribble(
  ~lake, ~habitat_x,     ~habitat_y,
  "ES",  "Depth",        "Light",
  "ES",  "Depth",        "Oxygen",
  "SP",  "Temperature",  "Light",
  "SP",  "Depth",        "Oxygen"
)


# Create plotting data
key_plot_data <- map_dfr(
  seq_len(nrow(key_pairs)),
  function(i) {
    
    pair <- key_pairs[i, ]
    
    ovl_wide %>%
      filter(lake == pair$lake) %>%
      select(
        fish_id,
        lake,
        all_of(c(pair$habitat_x, pair$habitat_y))
      ) %>%
      rename(
        x = !!pair$habitat_x,
        y = !!pair$habitat_y
      ) %>%
      mutate(
        habitat_x = pair$habitat_x,
        habitat_y = pair$habitat_y,
        relationship = paste(
          pair$habitat_x,
          "vs.",
          pair$habitat_y
        )
      )
  }
)


# Create final panel labels
key_plot_data <- key_plot_data %>%
  left_join(
    cor_labels,
    by = c("lake", "habitat_x", "habitat_y", "relationship")
  ) %>%
  mutate(
    panel = paste(lake, "\u2014", relationship)
  )

# Plot
ggplot(
  key_plot_data,
  aes(x = x, y = y)
) +
  
  geom_point(
    aes(color = lake),
    size = 3,
    alpha = 0.8
  ) +
  
  geom_text(
    data = key_plot_data %>%
      distinct(
        lake,
        relationship,
        panel,
        rho,
        p_adj,
        label
      ) %>%
      mutate(
        x = 0.42,
        y = 0.98
      ),
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  
  facet_wrap(
    ~ panel,
    ncol = 2
  ) +
  
  scale_color_manual(
    values = c(
      "ES" = "chocolate",
      "SP" = "dodgerblue4"
    )
  ) +
  
  scale_x_continuous(
    limits = c(0.4, 1),
    breaks = seq(0.4, 1, 0.2)
  ) +
  
  scale_y_continuous(
    limits = c(0.4, 1),
    breaks = seq(0.4, 1, 0.2)
  ) +
  
  labs(
    x = "OVL",
    y = "OVL"
  ) +
  
  theme_classic() +
  
  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "black"
    ),
    strip.text = element_text(
      size = 11
    ),
    legend.position = "none"
  )


#Some result tables 

# Fixed effects
model_fixed <- summary(ovl_mod)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("Effect") %>%
  rename(
    Estimate = Estimate,
    SE = `Std. Error`,
    df = df,
    t = `t value`,
    P = `Pr(>|t|)`
  )

# Random-effect variance components
model_var <- as.data.frame(VarCorr(ovl_mod)) %>%
  transmute(
    Effect = case_when(
      grp == "fish_id" ~ "Fish variance",
      grp == "Residual" ~ "Residual variance",
      TRUE ~ grp
    ),
    Estimate = vcov
  )

# ICC
icc_row <- tibble(
  Effect = "ICC",
  Estimate = icc_results$ICC,
  SE = NA_real_,
  df = NA_real_,
  t = NA_real_,
  P = NA_real_
)

# Combine
table1 <- model_fixed %>%
  select(Effect, Estimate, SE, df, t, P) %>%
  bind_rows(icc_row)

table1

table1_gt <- table1 %>%
  gt() %>%
  fmt_number(
    columns = c(Estimate, SE, df, t),
    decimals = 3
  ) %>%
  fmt_number(
    columns = P,
    decimals = 3
  ) %>%
  cols_label(
    Effect = "Effect",
    Estimate = "Estimate",
    SE = "SE",
    df = "df",
    t = "t",
    P = "P"
  ) %>%
  tab_header(
    title = "Mixed-effects model of individual–population habitat-use overlap"
  )

table1_gt

ovl_emm_df <- summary(
  ovl_emm,
  type = "response"
) %>%
  as.data.frame()

table2 <- ovl_emm_df %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c("Temperature", "Depth", "Light", "Oxygen")
    ),
    
    OVL = plogis(emmean),
    CI_low = plogis(lower.CL),
    CI_high = plogis(upper.CL)
  ) %>%
  arrange(lake, habitat) %>%
  transmute(
    Lake = lake,
    Habitat = habitat,
    `Estimated OVL (95% CI)` = sprintf(
      "%.3f (%.3f–%.3f)",
      OVL,
      CI_low,
      CI_high
    )
  )

table2

table2_gt <- table2 %>%
  gt() %>%
  tab_header(
    title = "Estimated individual–population overlap by habitat dimension"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  )

table2_gt

ovl_pairs <- pairs(
  ovl_emm,
  adjust = "tukey"
)


table_s1 <- summary(ovl_pairs) %>%
  as.data.frame() %>%
  select(
    lake,
    contrast,
    estimate,
    SE,
    df,
    p.value
  ) %>%
  rename(
    Lake = lake,
    Comparison = contrast,
    Estimate = estimate,
    SE = SE,
    df = df,
    `P (Tukey)` = p.value
  )

table_s1_gt <- table_s1 %>%
  gt() %>%
  fmt_number(
    columns = c(Estimate, SE, df),
    decimals = 3
  ) %>%
  fmt_number(
    columns = `P (Tukey)`,
    decimals = 3
  ) %>%
  tab_header(
    title = "Pairwise comparisons of estimated OVL among habitat dimensions"
  )

table_s1_gt


table_s2 <- cor_results %>%
  mutate(
    rho = sprintf("%.3f", rho),
    p = format.pval(
      p,
      digits = 3,
      eps = 0.001
    ),
    p_adj = format.pval(
      p_adj,
      digits = 3,
      eps = 0.001
    )
  ) %>%
  select(
    Lake = lake,
    `Habitat 1` = habitat_1,
    `Habitat 2` = habitat_2,
    Spearman_rho = rho,
    P = p,
    `P (FDR)` = p_adj
  )

table_s2

table_s2_gt <- table_s2 %>%
  gt() %>%
  cols_label(
    Lake = "Lake",
    `Habitat 1` = "Habitat 1",
    `Habitat 2` = "Habitat 2",
    Spearman_rho = "Spearman \u03c1",
    P = "P",
    `P (FDR)` = "P (FDR)"
  ) %>%
  tab_header(
    title = "Pairwise correlations among individual–population overlap coefficients"
  ) %>%
  cols_align(
    align = "center",
    columns = everything()
  )

table_s2_gt


table_s3 <- ovl_results %>%
  select(
    fish_id,
    lake,
    tl_pre,
    habitat,
    OVL
  ) %>%
  arrange(
    lake,
    fish_id,
    habitat
  )

table_s3

table_s3_gt <- table_s3 %>%
  gt() %>%
  fmt_number(
    columns = OVL,
    decimals = 3
  ) %>%
  cols_label(
    fish_id = "Fish ID",
    lake = "Lake",
    tl_pre = "Length (mm)",
    habitat = "Habitat dimension",
    OVL = "OVL"
  ) %>%
  tab_header(
    title = "Individual-level habitat-use overlap coefficients"
  )

table_s3_gt


####Directional OVL####

# ------------------------------------------------------------
# 2. Function to calculate OVL + directional differences
# ------------------------------------------------------------

calculate_ovl_directional <- function(df, grid_n = 512) {
  
  fish_ids <- unique(df$fish_id)
  
  # ----------------------------------------------------------
  # Create common density grid
  # ----------------------------------------------------------
  
  x_min <- min(df$value, na.rm = TRUE)
  x_max <- max(df$value, na.rm = TRUE)
  
  x_grid <- seq(x_min, x_max, length.out = grid_n)
  
  # ----------------------------------------------------------
  # Calculate individual densities
  # ----------------------------------------------------------
  
  individual_densities <- map_dfr(fish_ids, function(id) {
    
    x <- df %>%
      filter(fish_id == id) %>%
      pull(value)
    
    # Need enough observations for density()
    if (length(x) < 2 || length(unique(x)) < 2) {
      return(NULL)
    }
    
    dens <- density(
      x,
      from = x_min,
      to = x_max,
      n = grid_n
    )
    
    tibble(
      fish_id = id,
      x = dens$x,
      density = dens$y
    )
  })
  
  # ----------------------------------------------------------
  # Equal-weight population density
  # ----------------------------------------------------------
  
  population_density <- individual_densities %>%
    group_by(x) %>%
    summarise(
      pop_density = mean(density, na.rm = TRUE),
      .groups = "drop"
    )
  
  n_fish <- n_distinct(individual_densities$fish_id)
  
  # ----------------------------------------------------------
  # Calculate results for each fish
  # ----------------------------------------------------------
  
  results <- map_dfr(
    fish_ids,
    function(id) {
      
      ind <- individual_densities %>%
        filter(fish_id == id) %>%
        arrange(x)
      
      if (nrow(ind) == 0 || n_fish < 2) {
        return(NULL)
      }
      
      # Leave-one-out population density
      pop_loo <- individual_densities %>%
        filter(fish_id != id) %>%
        group_by(x) %>%
        summarise(
          pop_density = mean(density, na.rm = TRUE),
          .groups = "drop"
        )
      
      curves <- ind %>%
        select(x, ind_density = density) %>%
        left_join(pop_loo, by = "x") %>%
        mutate(
          # Difference between individual and population
          density_diff = ind_density - pop_density,
          
          # Positive = individual overuses this part
          # Negative = individual underuses this part
          overuse = pmax(density_diff, 0),
          underuse = pmax(-density_diff, 0),
          
          # Shared density for OVL
          shared = pmin(ind_density, pop_density)
        ) %>%
        arrange(x)
      
      # --------------------------------------------------------
      # OVL
      # --------------------------------------------------------
      
      ovl <- sum(
        diff(curves$x) *
          (head(curves$shared, -1) +
             tail(curves$shared, -1)) / 2,
        na.rm = TRUE
      )
      
      # --------------------------------------------------------
      # Individual and population means
      # --------------------------------------------------------
      
      mean_individual <- sum(
        diff(curves$x) *
          (
            head(curves$x * curves$ind_density, -1) +
              tail(curves$x * curves$ind_density, -1)
          ) / 2,
        na.rm = TRUE
      )
      
      mean_population <- sum(
        diff(curves$x) *
          (
            head(curves$x * curves$pop_density, -1) +
              tail(curves$x * curves$pop_density, -1)
          ) / 2,
        na.rm = TRUE
      )
      
      delta_mean <- mean_individual - mean_population
      
      # --------------------------------------------------------
      # Amount of underuse / overuse
      # --------------------------------------------------------
      
      underuse_amount <- sum(
        diff(curves$x) *
          (
            head(curves$underuse, -1) +
              tail(curves$underuse, -1)
          ) / 2,
        na.rm = TRUE
      )
      
      overuse_amount <- sum(
        diff(curves$x) *
          (
            head(curves$overuse, -1) +
              tail(curves$overuse, -1)
          ) / 2,
        na.rm = TRUE
      )
      
      # --------------------------------------------------------
      # Location of underuse
      # --------------------------------------------------------
      
      underuse_weighted <- curves$x * curves$underuse
      
      underuse_center <- ifelse(
        underuse_amount > 0,
        sum(
          diff(curves$x) *
            (
              head(underuse_weighted, -1) +
                tail(underuse_weighted, -1)
            ) / 2,
          na.rm = TRUE
        ) / underuse_amount,
        NA_real_
      )
      
      # --------------------------------------------------------
      # Location of overuse
      # --------------------------------------------------------
      
      overuse_weighted <- curves$x * curves$overuse
      
      overuse_center <- ifelse(
        overuse_amount > 0,
        sum(
          diff(curves$x) *
            (
              head(overuse_weighted, -1) +
                tail(overuse_weighted, -1)
            ) / 2,
          na.rm = TRUE
        ) / overuse_amount,
        NA_real_
      )
      
      tibble(
        fish_id = id,
        OVL = ovl,
        mean_individual = mean_individual,
        mean_population = mean_population,
        delta_mean = delta_mean,
        underuse_amount = underuse_amount,
        overuse_amount = overuse_amount,
        underuse_center = underuse_center,
        overuse_center = overuse_center
      )
    }
  )
  
  results
}


#Change light, include two weirdos 

# ============================================================
# 1. Rebuild habitat-use data
#    Light is now log10-transformed
# ============================================================

habitat_use <- bind_rows(
  
  # Temperature
  archival_essp_full %>%
    filter(
      year == 2022,
      month %in% 5:9,
      !is.na(temp_c)
    ) %>%
    select(fish_id, lake, value = temp_c) %>%
    mutate(habitat = "Temperature"),
  
  # Depth
  fish_depth_dist %>%
    filter(
      !is.na(depth_m),
      (lake == "ES" & depth_m <= 9) | lake == "SP"
    ) %>%
    select(fish_id, lake, value = depth_m) %>%
    mutate(habitat = "Depth"),
  
  # Light
  archival_essp_full %>%
    filter(
      year == 2022,
      month %in% 5:9,
      hour(date_time) >= 8,
      hour(date_time) <= 20,
      !is.na(I_pred_lux),
      I_pred_lux > 0
    ) %>%
    transmute(
      fish_id,
      lake,
      value = log10(I_pred_lux)
    ) %>%
    mutate(habitat = "Light"),
  
  # Oxygen
  fish_do %>%
    filter(!is.na(do_mgl)) %>%
    select(fish_id, lake, value = do_mgl) %>%
    mutate(habitat = "Oxygen")
)


# ============================================================
# 2. Recalculate directional OVL results
# ============================================================

directional_results <- habitat_use %>%
  group_by(lake, habitat) %>%
  nest() %>%
  mutate(
    results = map(
      data,
      ~ calculate_ovl_directional(.x)
    )
  ) %>%
  select(lake, habitat, results) %>%
  unnest(results) %>%
  left_join(
    fish_traits %>%
      select(fish_id, tl_pre),
    by = "fish_id"
  ) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Temperature",
        "Depth",
        "Light",
        "Oxygen"
      )
    ),
    lake = factor(
      lake,
      levels = c("ES", "SP")
    )
  )


# ============================================================
# 3. Create labels for the two temperature specialists
# ============================================================

label_fish <- directional_results %>%
  filter(
    habitat == "Temperature",
    fish_id %in% c("SP_40", "SP_02")
  )


# ============================================================
# 4. OVL vs. directional habitat shift
# ============================================================

p_ovl_direction <- ggplot(
  directional_results,
  aes(
    x = OVL,
    y = delta_mean,
    color = lake
  )
) +
  
  # No directional shift
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # All individuals
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  
  # Highlight the two cold-temperature specialists
  geom_point(
    data = label_fish,
    aes(
      x = OVL,
      y = delta_mean
    ),
    size = 3.5,
    stroke = 1.2
  ) +
  
  geom_text(
    data = label_fish,
    aes(
      label = fish_id
    ),
    nudge_x = 0.015,
    nudge_y = 0.15,
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ habitat,
    scales = "free_y",
    ncol = 2
  ) +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  
  labs(
    x = "Overlap with leave-one-out population (OVL)",
    y = "Individual − population mean",
    color = "Lake"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )


p_ovl_direction

#Specifically look at SP02 and SP40


# ============================================================
# Function to extract individual + leave-one-out population
# density curves
# ============================================================

get_density_curves <- function(df, selected_fish, grid_n = 512) {
  
  fish_ids <- unique(df$fish_id)
  
  x_min <- min(df$value, na.rm = TRUE)
  x_max <- max(df$value, na.rm = TRUE)
  
  # Individual densities
  individual_densities <- map_dfr(fish_ids, function(id) {
    
    x <- df %>%
      filter(.data$fish_id == id) %>%
      pull(value)
    
    if (length(x) < 2 || length(unique(x)) < 2) {
      return(NULL)
    }
    
    dens <- density(
      x,
      from = x_min,
      to = x_max,
      n = grid_n
    )
    
    tibble(
      fish_id = id,
      x = dens$x,
      density = dens$y
    )
  })
  
  # Selected individual
  individual <- individual_densities %>%
    filter(.data$fish_id == selected_fish) %>%
    select(x, individual_density = density)
  
  # Leave-one-out population
  population <- individual_densities %>%
    filter(.data$fish_id != selected_fish) %>%
    group_by(x) %>%
    summarise(
      population_density = mean(density, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Combine
  individual %>%
    left_join(population, by = "x") %>%
    mutate(
      difference = individual_density - population_density,
      fish_id = selected_fish
    )
}


# Recalculate the two fish
sp40_curve <- get_density_curves(
  habitat_use %>% filter(habitat == "Temperature"),
  selected_fish = "SP_40"
)

sp02_curve <- get_density_curves(
  habitat_use %>% filter(habitat == "Temperature"),
  selected_fish = "SP_02"
)

temp_curves <- bind_rows(
  sp40_curve,
  sp02_curve
)

# ============================================================
# Convert to long format for plotting
# ============================================================

density_long <- temp_curves %>%
  select(
    fish_id,
    x,
    individual_density,
    population_density
  ) %>%
  pivot_longer(
    cols = c(
      individual_density,
      population_density
    ),
    names_to = "distribution",
    values_to = "density"
  ) %>%
  mutate(
    distribution = recode(
      distribution,
      individual_density = "Individual",
      population_density = "Leave-one-out population"
    )
  )


# ============================================================
# Plot distributions
# ============================================================

p_density <- ggplot(
  density_long,
  aes(
    x = x,
    y = density,
    linetype = distribution
  )
) +
  geom_line(
    linewidth = 1
  ) +
  facet_wrap(
    ~ fish_id,
    ncol = 1
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Density",
    linetype = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

p_density

# ============================================================
# Signed density difference
# ============================================================

p_difference <- ggplot(
  temp_curves,
  aes(
    x = x,
    y = difference
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_ribbon(
    aes(
      ymin = pmin(difference, 0),
      ymax = pmax(difference, 0)
    ),
    alpha = 0.4
  ) +
  geom_line(
    linewidth = 1
  ) +
  facet_wrap(
    ~ fish_id,
    ncol = 1
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Individual − population density"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_difference



#Slice and dice - heatmap 

#============================================================
  # Prepare heatmap data
  # ============================================================

heatmap_df <- directional_results %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Temperature",
        "Depth",
        "Light",
        "Oxygen"
      )
    )
  )


# ============================================================
# Standardize directional shift within habitat
# ============================================================

heatmap_direction <- heatmap_df %>%
  group_by(habitat) %>%
  mutate(
    delta_z = as.numeric(
      scale(delta_mean)
    )
  ) %>%
  ungroup()


# ============================================================
# Order fish
# ============================================================

fish_order <- heatmap_direction %>%
  group_by(fish_id, lake) %>%
  summarise(
    mean_specialization = mean(abs(delta_z), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(lake, desc(mean_specialization)) %>%
  pull(fish_id)


# ============================================================
# Directional-shift heatmap
# ============================================================

p_direction_heatmap <- ggplot(
  heatmap_direction,
  aes(
    x = habitat,
    y = factor(fish_id, levels = rev(fish_order)),
    fill = delta_z
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  facet_grid(
    ~ lake,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient2(
    midpoint = 0,
    name = "Directional shift\n(z-score)"
  ) +
  labs(
    x = NULL,
    y = "Fish"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_direction_heatmap


# ============================================================
# OVL heatmap
# ============================================================

p_ovl_heatmap <- ggplot(
  heatmap_df,
  aes(
    x = habitat,
    y = factor(fish_id, levels = rev(fish_order)),
    fill = OVL
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  facet_grid(
    ~ lake,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    name = "OVL"
  ) +
  labs(
    x = NULL,
    y = "Fish"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_ovl_heatmap


#Try a PCA? 

# ============================================================
# PCA dataset
# ============================================================

pca_df <- directional_results %>%
  ungroup() %>%
  select(
    fish_id,
    lake,
    habitat,
    delta_mean
  ) %>%
  pivot_wider(
    names_from = habitat,
    values_from = delta_mean
  ) %>%
  ungroup()


# ============================================================
# Check that everything is finite
# ============================================================

pca_df %>%
  summarise(
    across(
      c(Temperature, Depth, Light, Oxygen),
      ~ sum(!is.finite(.x))
    )
  )


# ============================================================
# Standardize each habitat dimension across ALL fish
# ============================================================

pca_scaled <- pca_df %>%
  mutate(
    across(
      c(Temperature, Depth, Light, Oxygen),
      ~ as.numeric(scale(.x))
    )
  )


# ============================================================
# PCA
# ============================================================

pca_mod <- prcomp(
  pca_scaled %>%
    select(
      Temperature,
      Depth,
      Light,
      Oxygen
    ),
  center = FALSE,
  scale. = FALSE
)


# ============================================================
# Variance explained
# ============================================================

summary(pca_mod)


# ============================================================
# PCA loadings
# ============================================================

pca_loadings <- as.data.frame(pca_mod$rotation) %>%
  mutate(
    variable = rownames(.)
  )

pca_loadings


#Making the biplot 

# ============================================================
# PCA scores
# ============================================================

pca_scores <- as.data.frame(pca_mod$x) %>%
  select(PC1, PC2) %>%
  bind_cols(
    pca_scaled %>%
      select(fish_id, lake)
  )


# ============================================================
# PCA loadings
# ============================================================

pca_loadings <- as.data.frame(pca_mod$rotation) %>%
  mutate(
    variable = rownames(.)
  )


# ============================================================
# Scale loading arrows to the score space
# ============================================================

arrow_scale <- min(
  diff(range(pca_scores$PC1)) /
    diff(range(pca_loadings$PC1)),
  
  diff(range(pca_scores$PC2)) /
    diff(range(pca_loadings$PC2))
) * 0.7


pca_loadings_plot <- pca_loadings %>%
  mutate(
    PC1_arrow = PC1 * arrow_scale,
    PC2_arrow = PC2 * arrow_scale
  )


# ============================================================
# PCA biplot
# ============================================================

p_pca_biplot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = lake
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray70"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray70"
  ) +
  
  # Fish
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  # Highlight the two cold-temperature fish
  geom_point(
    data = pca_scores %>%
      filter(fish_id %in% c("SP_02", "SP_40")),
    size = 4,
    stroke = 1.2
  ) +
  
  geom_text(
    data = pca_scores %>%
      filter(fish_id %in% c("SP_02", "SP_40")),
    aes(label = fish_id),
    nudge_y = 0.15,
    show.legend = FALSE
  ) +
  
  # Environmental vectors
  geom_segment(
    data = pca_loadings_plot,
    aes(
      x = 0,
      y = 0,
      xend = PC1_arrow,
      yend = PC2_arrow
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = unit(0.2, "cm")
    ),
    linewidth = 0.8
  ) +
  
  # Vector labels
  geom_text(
    data = pca_loadings_plot,
    aes(
      x = PC1_arrow,
      y = PC2_arrow,
      label = variable
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    nudge_x = 0.1,
    nudge_y = 0.1
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(summary(pca_mod)$importance[2, 1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_mod)$importance[2, 2] * 100, 1),
      "%)"
    ),
    color = "Lake"
  )

p_pca_biplot


#Incorporate size 

# ============================================================
# Add fish size to PCA scores
# ============================================================

pca_scores_size <- pca_scores %>%
  left_join(
    fish_traits %>%
      select(fish_id, tl_pre),
    by = "fish_id"
  )


# Check that all fish have size
pca_scores_size %>%
  summarise(
    n_fish = n(),
    n_missing_tl = sum(is.na(tl_pre)),
    min_tl = min(tl_pre, na.rm = TRUE),
    max_tl = max(tl_pre, na.rm = TRUE)
  )


# ============================================================
# PCA biplot with fish size
# ============================================================

p_pca_size <- ggplot(
  pca_scores_size,
  aes(
    x = PC1,
    y = PC2,
    color = lake,
    size = tl_pre
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray70"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray70"
  ) +
  
  # Fish
  geom_point(
    alpha = 0.8
  ) +
  
  # Highlight the two cold-temperature specialists
  geom_point(
    data = pca_scores_size %>%
      filter(fish_id %in% c("SP_02", "SP_40")),
    stroke = 1.2
  ) +
  
  geom_text(
    data = pca_scores_size %>%
      filter(fish_id %in% c("SP_02", "SP_40")),
    aes(label = fish_id),
    nudge_y = 0.15,
    show.legend = FALSE
  ) +
  
  # Environmental vectors
  geom_segment(
    data = pca_loadings_plot,
    aes(
      x = 0,
      y = 0,
      xend = PC1_arrow,
      yend = PC2_arrow
    ),
    inherit.aes = FALSE,
    arrow = arrow(
      length = unit(0.2, "cm")
    ),
    linewidth = 0.8
  ) +
  
  # Vector labels
  geom_text(
    data = pca_loadings_plot,
    aes(
      x = PC1_arrow,
      y = PC2_arrow,
      label = variable
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4,
    nudge_x = 0.1,
    nudge_y = 0.1
  ) +
  
  scale_size_continuous(
    name = "Total length (mm)",
    range = c(2.5, 7)
  ) +
  
  labs(
    x = paste0(
      "PC1 (",
      round(summary(pca_mod)$importance[2, 1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca_mod)$importance[2, 2] * 100, 1),
      "%"
    ),
    color = "Lake"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

p_pca_size

# Size relationships with the two PCA axes

summary(lm(PC1 ~ tl_pre * lake, data = pca_scores_size))
summary(lm(PC2 ~ tl_pre * lake, data = pca_scores_size))

# ============================================================
# PC2 vs fish size
# ============================================================

ggplot(
  pca_scores_size,
  aes(
    x = tl_pre,
    y = PC2,
    color = lake
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray70"
  ) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  # Highlight the two cold-temperature specialists
  geom_point(
    data = pca_scores_size %>%
      filter(fish_id %in% c("SP_02", "SP_40")),
    size = 4,
    stroke = 1.2
  ) +
  
  geom_text(
    data = pca_scores_size %>%
      filter(fish_id %in% c("SP_02", "SP_40")),
    aes(label = fish_id),
    nudge_y = 0.15,
    show.legend = FALSE
  ) +
  
  # Lake-specific exploratory trends
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  
  labs(
    x = "Total length (mm)",
    y = paste0(
      "PC2 (",
      round(summary(pca_mod)$importance[2, 2] * 100, 1),
      "%)"
    ),
    color = "Lake"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top"
  )

# ============================================================
# Directional habitat shift vs fish size
# ============================================================

size_shift_df <- directional_results %>%
  ungroup() %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Temperature",
        "Depth",
        "Light",
        "Oxygen"
      )
    )
  )


# Check
size_shift_df %>%
  select(fish_id, lake, tl_pre, habitat, delta_mean) %>%
  head()


# ============================================================
# Plot
# ============================================================

p_size_shift <- ggplot(
  size_shift_df,
  aes(
    x = tl_pre,
    y = delta_mean,
    color = lake
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray60"
  ) +
  
  geom_point(
    size = 2.8,
    alpha = 0.8
  ) +
  
  facet_wrap(
    ~ habitat,
    scales = "free_y",
    ncol = 2
  ) +
  
  labs(
    x = "Total length (mm)",
    y = "Individual − population mean",
    color = "Lake"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

p_size_shift

#Correlation bogus 

# ============================================================
# Correlation among individual directional habitat shifts
# ============================================================

cor_df <- directional_results %>%
  ungroup() %>%
  select(
    fish_id,
    lake,
    habitat,
    delta_mean
  ) %>%
  pivot_wider(
    names_from = habitat,
    values_from = delta_mean
  )


# ES correlation matrix
cor_ES <- cor_df %>%
  filter(lake == "ES") %>%
  select(
    Temperature,
    Depth,
    Light,
    Oxygen
  ) %>%
  cor(
    use = "pairwise.complete.obs",
    method = "pearson"
  )


# SP correlation matrix
cor_SP <- cor_df %>%
  filter(lake == "SP") %>%
  select(
    Temperature,
    Depth,
    Light,
    Oxygen
  ) %>%
  cor(
    use = "pairwise.complete.obs",
    method = "pearson"
  )


cor_ES
cor_SP

# Combine lake-specific correlations
cor_long <- bind_rows(
  
  as.data.frame(as.table(cor_ES)) %>%
    rename(
      variable_1 = Var1,
      variable_2 = Var2,
      correlation = Freq
    ) %>%
    mutate(lake = "ES"),
  
  as.data.frame(as.table(cor_SP)) %>%
    rename(
      variable_1 = Var1,
      variable_2 = Var2,
      correlation = Freq
    ) %>%
    mutate(lake = "SP")
) %>%
  mutate(
    variable_1 = factor(
      variable_1,
      levels = c(
        "Temperature",
        "Depth",
        "Light",
        "Oxygen"
      )
    ),
    variable_2 = factor(
      variable_2,
      levels = rev(c(
        "Temperature",
        "Depth",
        "Light",
        "Oxygen"
      ))
    )
  )


# Plot
ggplot(
  cor_long,
  aes(
    x = variable_1,
    y = variable_2,
    fill = correlation
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", correlation)
    ),
    size = 4
  ) +
  facet_wrap(~ lake) +
  scale_fill_gradient2(
    limits = c(-1, 1),
    midpoint = 0,
    name = "Pearson r"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


