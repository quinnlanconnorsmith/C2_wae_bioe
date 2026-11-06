
# ============================================================
# AUDIT: Daily volume searched and energetic requirements
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Inspect the source and derived datasets
cat("\n--- archival_essp_full ---\n")
str(archival_essp_full)

cat("\n--- df_energy_swim ---\n")
str(df_energy_swim)

cat("\n--- df_volume_energy_hourly ---\n")
str(df_volume_energy_hourly)

cat("\n--- df_volume_energy_daily ---\n")
str(df_volume_energy_daily)

# 2. Check units and distributions of the source variables
source_checks <- list(
  volume_searched_m3 = archival_essp_full$volume_searched_m3,
  swim_J_h           = df_energy_swim$swim_J_h,
  E_J_h              = df_energy_swim$E_J_h,
  E_J_g_h            = df_energy_swim$E_J_g_h
)

lapply(source_checks, function(x) {
  c(
    n = length(x),
    n_missing = sum(is.na(x)),
    n_nonfinite = sum(!is.finite(x)),
    minimum = min(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    maximum = max(x, na.rm = TRUE)
  )
})

# 3. Check hourly and daily coverage
hourly_audit <- df_volume_energy_hourly %>%
  group_by(lake, fish_id) %>%
  summarise(
    n_hour_rows = n(),
    n_unique_hours = n_distinct(hour_time),
    first_hour = min(hour_time, na.rm = TRUE),
    last_hour = max(hour_time, na.rm = TRUE),
    duplicate_hours = n() - n_distinct(hour_time),
    missing_volume = sum(is.na(mean_search)),
    missing_energy = sum(is.na(E_J_h)),
    .groups = "drop"
  )

daily_audit <- df_volume_energy_daily %>%
  group_by(lake, fish_id) %>%
  summarise(
    n_days = n(),
    duplicate_dates = n() - n_distinct(date),
    mean_daily_volume = mean(daily_volume, na.rm = TRUE),
    mean_daily_energy = mean(daily_E_J, na.rm = TRUE),
    missing_volume = sum(is.na(daily_volume)),
    missing_energy = sum(is.na(daily_E_J)),
    .groups = "drop"
  )

print(hourly_audit, n = Inf)
print(daily_audit, n = Inf)

# 4. Check whether the hourly join multiplied or lost records
join_check <- archival_hourly_search %>%
  summarise(
    search_rows = n(),
    unique_fish_hours = n_distinct(paste(fish_id, hour_time))
  )

energy_key_check <- df_energy_swim %>%
  summarise(
    energy_rows = n(),
    unique_fish_hours = n_distinct(paste(fish_id, hour_time)),
    duplicate_fish_hours = n() -
      n_distinct(paste(fish_id, hour_time))
  )

print(join_check)
print(energy_key_check)

# 5. Recalculate daily totals independently from the hourly data
daily_recalculated <- df_volume_energy_hourly %>%
  group_by(lake, fish_id, date) %>%
  summarise(
    daily_volume_check = sum(mean_search, na.rm = TRUE),
    daily_energy_check = sum(E_J_h, na.rm = TRUE),
    n_hours = n_distinct(hour_time),
    .groups = "drop"
  )

daily_comparison <- df_volume_energy_daily %>%
  select(lake, fish_id, date, daily_volume, daily_E_J) %>%
  full_join(
    daily_recalculated,
    by = c("lake", "fish_id", "date")
  ) %>%
  mutate(
    volume_difference = daily_volume - daily_volume_check,
    energy_difference = daily_E_J - daily_energy_check
  )

summary(daily_comparison$volume_difference)
summary(daily_comparison$energy_difference)

# 6. Inspect the selected fish against all fish
fish_groups_audit <- tibble::tribble(
  ~lake, ~type, ~fish_id,
  "ES", "Average", "ES_03",
  "ES", "Extreme", "ES_70",
  "SP", "Average", "SP_54",
  "SP", "Extreme", "SP_40"
)

df_volume_energy_daily %>%
  semi_join(fish_groups_audit, by = c("lake", "fish_id")) %>%
  group_by(lake, fish_id) %>%
  summarise(
    n_days = n(),
    mean_volume = mean(daily_volume, na.rm = TRUE),
    mean_energy_J = mean(daily_E_J, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = Inf)


#Alright barf let's get into this 

# Check available environmental landscape
glimpse(env_perf)
glimpse(env_pop)

# Check observed fish use
glimpse(obs_fish)
glimpse(archival_hourly_search)
glimpse(df_volume_energy_daily)

# Check the final plotting datasets, if they exist
if (exists("obs_fish_plot")) glimpse(obs_fish_plot)
if (exists("obs_hourly_plot")) glimpse(obs_hourly_plot)

# Check the available and observed ranges
env_perf %>%
  summarise(
    n = n(),
    min_temp = min(temp_c, na.rm = TRUE),
    max_temp = max(temp_c, na.rm = TRUE),
    min_light = min(light_lux, na.rm = TRUE),
    max_light = max(light_lux, na.rm = TRUE),
    n_missing_volume = sum(is.na(volume_m3_h))
  )

obs_fish %>%
  summarise(
    n = n(),
    n_fish = n_distinct(fish_id),
    n_missing_volume = sum(is.na(volume_m3_h)),
    n_missing_energy = sum(is.na(E_kJ_h))
  )

#check2 

# Check structures
str(env_hourly)
str(df_volume_energy_daily)
str(df_volume_energy_hourly)
str(archival_essp_hourly)

# Check available environmental space
env_hourly %>%
  summarise(
    n = n(),
    n_lakes = n_distinct(lake),
    n_times = n_distinct(datetime_hour),
    n_depths = n_distinct(depth_m),
    missing_temp = sum(is.na(temp_c)),
    missing_light = sum(is.na(light_lux))
  )

# Check observed fish outcomes
df_volume_energy_daily %>%
  summarise(
    n = n(),
    n_fish = n_distinct(fish_id),
    n_lakes = n_distinct(lake),
    n_dates = n_distinct(date),
    missing_volume = sum(is.na(daily_volume)),
    missing_energy = sum(is.na(daily_E_J))
  )

# Summarize observed daily outcomes by lake
df_volume_energy_daily %>%
  group_by(lake) %>%
  summarise(
    n_fish = n_distinct(fish_id),
    n_fish_days = n(),
    volume_min = min(daily_volume, na.rm = TRUE),
    volume_median = median(daily_volume, na.rm = TRUE),
    volume_max = max(daily_volume, na.rm = TRUE),
    energy_min = min(daily_E_J, na.rm = TRUE),
    energy_median = median(daily_E_J, na.rm = TRUE),
    energy_max = max(daily_E_J, na.rm = TRUE),
    .groups = "drop"
  )



#Yay finishing the auditing 


# ============================================================
# Volume searched and energetic requirements:
# Population use vs individual use
#
# Existing objects required:
#   df_volume_energy_daily
#   fish_groups
#
# This section reconstructs observed fish-use distributions.
# Lake-wide availability is added after the environmental
# performance calculation is confirmed.
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Clean the daily fish-use data
# ------------------------------------------------------------

daily_use <- df_volume_energy_daily %>%
  filter(
    lake %in% c("ES", "SP"),
    !is.na(fish_id),
    is.finite(daily_volume),
    is.finite(daily_E_J)
  ) %>%
  mutate(
    fish_id = as.character(fish_id),
    lake = as.character(lake)
  )

fish_groups_plot <- fish_groups %>%
  mutate(
    fish_id = as.character(fish_id),
    lake = as.character(lake)
  )

# Verify selected fish exist in the daily dataset
selected_fish_check <- fish_groups_plot %>%
  left_join(
    daily_use %>%
      distinct(lake, fish_id) %>%
      mutate(found = TRUE),
    by = c("lake", "fish_id")
  )

print(selected_fish_check)

# ------------------------------------------------------------
# 2. Compare observed daily distributions
# ------------------------------------------------------------

# Population: all fish in each lake.
# This pooled distribution weights each fish-day equally.
# We will also assess equal weighting by fish separately.

population_daily <- daily_use

individual_daily <- daily_use %>%
  inner_join(
    fish_groups_plot,
    by = c("lake", "fish_id")
  )

# ------------------------------------------------------------
# 3. Density normalization check
# ------------------------------------------------------------

density_area <- function(x, n = 512) {
  x <- x[is.finite(x)]
  
  if (length(x) < 2 || length(unique(x)) < 2) {
    return(NA_real_)
  }
  
  d <- density(x, n = n)
  sum(diff(d$x) * (head(d$y, -1) + tail(d$y, -1)) / 2)
}

density_check <- bind_rows(
  population_daily %>%
    group_by(lake) %>%
    summarise(
      metric = "Daily volume searched",
      n = n(),
      density_area = density_area(daily_volume),
      .groups = "drop"
    ),
  population_daily %>%
    group_by(lake) %>%
    summarise(
      metric = "Daily energy required",
      n = n(),
      density_area = density_area(daily_E_J),
      .groups = "drop"
    )
)

print(density_check)

# ------------------------------------------------------------
# 4. Plot observed population and individual distributions
# ------------------------------------------------------------

volume_plot <- ggplot() +
  geom_density(
    data = population_daily,
    aes(x = daily_volume),
    fill = "grey75",
    color = "steelblue",
    linetype = "dashed",
    alpha = 0.45,
    linewidth = 0.8
  ) +
  geom_density(
    data = individual_daily,
    aes(x = daily_volume),
    color = "black",
    linewidth = 1
  ) +
  facet_grid(~ lake + type, scales = "free_x") +
  theme_classic() +
  labs(
    x = expression("Daily volume searched (m"^3*")"),
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "bold"))

energy_plot <- ggplot() +
  geom_density(
    data = population_daily,
    aes(x = daily_E_J / 1000),
    fill = "grey75",
    color = "steelblue",
    linetype = "dashed",
    alpha = 0.45,
    linewidth = 0.8
  ) +
  geom_density(
    data = individual_daily,
    aes(x = daily_E_J / 1000),
    color = "black",
    linewidth = 1
  ) +
  facet_grid(~ lake + type, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Daily energy required (kJ)",
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "bold"))

volume_plot / energy_plot


#Oh no back to auditing 

# 1. Check the lake-wide environmental data
env_hourly %>%
  group_by(lake) %>%
  summarise(
    n_rows = n(),
    n_times = n_distinct(datetime_hour),
    n_depths = n_distinct(depth_m),
    temp_range = paste(
      round(min(temp_c, na.rm = TRUE), 2),
      round(max(temp_c, na.rm = TRUE), 2),
      sep = " to "
    ),
    light_range = paste(
      signif(min(light_lux, na.rm = TRUE), 3),
      signif(max(light_lux, na.rm = TRUE), 3),
      sep = " to "
    ),
    .groups = "drop"
  )

# 2. Check the observed population-level outcomes
df_volume_energy_daily %>%
  filter(lake %in% c("ES", "SP")) %>%
  group_by(lake) %>%
  summarise(
    n_fish = n_distinct(fish_id),
    n_fish_days = n(),
    volume_median = median(daily_volume, na.rm = TRUE),
    energy_median_kJ = median(daily_E_J, na.rm = TRUE) / 1000,
    .groups = "drop"
  )

# 3. Check whether the environmental performance dataset exists
if (exists("env_perf")) {
  print(names(env_perf))
  print(dplyr::glimpse(env_perf))
} else {
  message("env_perf does not exist in the current R session")
}



#New plot strucutre? /

library(dplyr)
library(ggplot2)
library(patchwork)

# define groups
fish_groups <- tibble::tribble(
  ~lake, ~type,      ~fish_id,
  "ES",  "Average",  "ES_04",
  "ES",  "Extreme",  "ES_03",
  "SP",  "Average",  "SP_29",
  "SP",  "Extreme",  "SP_02"
)

# Selected fish
fish_groups_plot <- fish_groups %>%
  mutate(
    lake = as.character(lake),
    fish_id = as.character(fish_id)
  )

# Population-level observed fish records
obs_fish_plot <- obs_fish %>%
  filter(
    lake %in% c("ES", "SP"),
    is.finite(volume_m3_h),
    is.finite(E_kJ_h)
  ) %>%
  mutate(
    lake = as.character(lake),
    fish_id = as.character(fish_id)
  )

# Four selected individuals
individual_fish_plot <- obs_fish_plot %>%
  inner_join(
    fish_groups_plot,
    by = c("lake", "fish_id")
  )

# Environmental availability
env_pop_plot <- env_pop %>%
  filter(
    lake %in% c("ES", "SP"),
    is.finite(volume_m3_h),
    is.finite(E_kJ_h),
    is.finite(vol_weight),
    vol_weight > 0
  ) %>%
  mutate(lake = as.character(lake))

# Volume searched
volume_plot <- ggplot() +
  geom_density(
    data = env_pop_plot,
    aes(x = volume_m3_h, weight = vol_weight),
    fill = "grey75",
    color = "grey50",
    alpha = 0.6
  ) +
  geom_density(
    data = obs_fish_plot,
    aes(x = volume_m3_h),
    color = "steelblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_density(
    data = individual_fish_plot,
    aes(x = volume_m3_h),
    color = "black",
    linewidth = 1.1
  ) +
  facet_grid(~lake + type, scales = "free_x") +
  theme_classic() +
  labs(
    x = expression("Volume searched (m"^3*"/h)"),
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "bold"))

# Energetic requirements
energy_plot <- ggplot() +
  geom_density(
    data = env_pop_plot,
    aes(x = E_kJ_h, weight = vol_weight),
    fill = "grey75",
    color = "grey50",
    alpha = 0.6
  ) +
  geom_density(
    data = obs_fish_plot,
    aes(x = E_kJ_h),
    color = "steelblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_density(
    data = individual_fish_plot,
    aes(x = E_kJ_h),
    color = "black",
    linewidth = 1.1
  ) +
  facet_grid(~lake + type, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Energetic requirement (kJ/h)",
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "bold"))

volume_plot / energy_plot



#Diagnostics 

# Compare the actual ranges and distributions behind each curve
availability_check <- env_pop %>%
  group_by(lake) %>%
  summarise(
    n = n(),
    volume_median = median(volume_m3_h, na.rm = TRUE),
    volume_95 = quantile(volume_m3_h, 0.95, na.rm = TRUE),
    energy_median = median(E_kJ_h, na.rm = TRUE),
    energy_95 = quantile(E_kJ_h, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

population_check <- obs_fish %>%
  group_by(lake) %>%
  summarise(
    n = n(),
    volume_median = median(volume_m3_h, na.rm = TRUE),
    volume_95 = quantile(volume_m3_h, 0.95, na.rm = TRUE),
    energy_median = median(E_kJ_h, na.rm = TRUE),
    energy_95 = quantile(E_kJ_h, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

availability_check
population_check

env_pop %>%
  group_by(lake) %>%
  summarise(
    max_volume = max(volume_m3_h, na.rm = TRUE),
    max_energy = max(E_kJ_h, na.rm = TRUE),
    volume_over_100 = mean(volume_m3_h > 100, na.rm = TRUE) * 100,
    volume_over_200 = mean(volume_m3_h > 200, na.rm = TRUE) * 100,
    weight_over_100 = sum(
      vol_weight[volume_m3_h > 100], na.rm = TRUE
    ),
    .groups = "drop"
  )

env_pop %>%
  filter(volume_m3_h > 100) %>%
  select(
    lake, datetime_hour, depth_m, temp_c, light_lux,
    reaction_dist_m, swim_speed_m_s, volume_m3_h,
    vol_weight, tl_pre
  ) %>%
  arrange(desc(volume_m3_h)) %>%
  print(n = 30)

env_pop %>%
  group_by(lake) %>%
  summarise(
    total_weight = sum(vol_weight, na.rm = TRUE),
    n_fish = n_distinct(fish_id),
    n_hours = n_distinct(datetime_hour),
    .groups = "drop"
  )

# Check weights before fish-size expansion
env_perf %>%
  group_by(lake, datetime_hour) %>%
  summarise(
    sum_weight = sum(vol_weight, na.rm = TRUE),
    n_depths = n_distinct(depth_m),
    .groups = "drop"
  ) %>%
  group_by(lake) %>%
  summarise(
    min_sum = min(sum_weight),
    median_sum = median(sum_weight),
    max_sum = max(sum_weight),
    .groups = "drop"
  )

# Check depth-to-hypsography assignment
env_perf %>%
  distinct(lake, depth_m, depth_bin, layer_volume_m3) %>%
  group_by(lake, depth_bin) %>%
  summarise(
    n_depths = n_distinct(depth_m),
    n_layer_volumes = n_distinct(layer_volume_m3),
    min_layer_volume = min(layer_volume_m3, na.rm = TRUE),
    max_layer_volume = max(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = 30)

#FUCK 

# Check the current env_perf structure
names(env_perf)
dplyr::glimpse(env_perf)

# Check whether the hypsography lookup still exists
exists("hypso_layers")
if (exists("hypso_layers")) {
  names(hypso_layers)
}

# Check whether env_pop still has the weighted data
names(env_pop)

# Recover one environmental record per depth and hour
env_weights_check <- env_pop %>%
  distinct(
    lake, datetime_hour, depth_m,
    depth_bin, layer_volume_m3,
    n_sublayers, sublayer_volume, vol_weight
  )

# Verify weights sum to one per lake-hour
env_weights_check %>%
  group_by(lake, datetime_hour) %>%
  summarise(
    sum_weight = sum(vol_weight, na.rm = TRUE),
    n_depths = n_distinct(depth_m),
    .groups = "drop"
  ) %>%
  group_by(lake) %>%
  summarise(
    min_sum = min(sum_weight),
    median_sum = median(sum_weight),
    max_sum = max(sum_weight),
    .groups = "drop"
  )

env_weights_check %>%
  group_by(lake, depth_bin) %>%
  summarise(
    n_depths = n_distinct(depth_m),
    n_layer_volumes = n_distinct(layer_volume_m3),
    min_layer_volume = min(layer_volume_m3, na.rm = TRUE),
    max_layer_volume = max(layer_volume_m3, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = 40)

env_weights_check %>%
  left_join(
    hypso_layers %>%
      select(lake, depth_top, depth_bottom, layer_volume_m3),
    by = c("lake", "layer_volume_m3")
  ) %>%
  group_by(lake, depth_bin) %>%
  summarise(
    min_depth = min(depth_m),
    max_depth = max(depth_m),
    min_layer_top = min(depth_top, na.rm = TRUE),
    max_layer_bottom = max(depth_bottom, na.rm = TRUE),
    n_unmatched = sum(is.na(depth_top)),
    .groups = "drop"
  ) %>%
  print(n = 40)

#Alreight lets test test test 


# Reference fish length
tl_ref <- 440

# Convert reference length to mass
mass_ref_g <- length_mm_to_mass_g(tl_ref)

# Retain one row per environmental depth and hour
env_reference <- env_pop %>%
  distinct(
    lake, datetime_hour, depth_m,
    temp_c, light_lux, vol_weight,
    thermal_perf, swim_speed_m_s,
    reaction_dist_m, distance_m_h, volume_m3_h
  ) %>%
  mutate(
    mass_g = mass_ref_g,
    R_gO2_g_h = respiration_FB4(
      mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT
    ),
    E_kJ_h = (R_gO2_g_h * 13.6) * mass_g
  )

# Check the resulting dataset
env_reference %>%
  group_by(lake) %>%
  summarise(
    n_records = n(),
    median_volume = median(volume_m3_h, na.rm = TRUE),
    median_energy = median(E_kJ_h, na.rm = TRUE),
    max_volume = max(volume_m3_h, na.rm = TRUE),
    max_energy = max(E_kJ_h, na.rm = TRUE),
    .groups = "drop"
  )

#Test figure: 


# Volume searched
volume_reference_plot <- ggplot() +
  geom_density(
    data = env_reference %>%
      filter(lake %in% c("ES", "SP"),
             is.finite(volume_m3_h),
             is.finite(vol_weight),
             vol_weight > 0),
    aes(x = volume_m3_h, weight = vol_weight),
    fill = "grey75", color = "grey50", alpha = 0.6
  ) +
  geom_density(
    data = obs_fish_plot,
    aes(x = volume_m3_h),
    color = "steelblue", linetype = "dashed", linewidth = 1
  ) +
  geom_density(
    data = individual_fish_plot,
    aes(x = volume_m3_h),
    color = "black", linewidth = 1.1
  ) +
  facet_grid(~lake + type, scales = "free_x") +
  theme_classic() +
  labs(
    x = expression("Volume searched (m"^3*"/h)"),
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "bold"))

# Energetic requirements
energy_reference_plot <- ggplot() +
  geom_density(
    data = env_reference %>%
      filter(lake %in% c("ES", "SP"),
             is.finite(E_kJ_h),
             is.finite(vol_weight),
             vol_weight > 0),
    aes(x = E_kJ_h, weight = vol_weight),
    fill = "grey75", color = "grey50", alpha = 0.6
  ) +
  geom_density(
    data = obs_fish_plot,
    aes(x = E_kJ_h),
    color = "steelblue", linetype = "dashed", linewidth = 1
  ) +
  geom_density(
    data = individual_fish_plot,
    aes(x = E_kJ_h),
    color = "black", linewidth = 1.1
  ) +
  facet_grid(~lake + type, scales = "free_x") +
  theme_classic() +
  labs(
    x = expression("Energetic requirements (kJ/h)"),
    y = "Density"
  ) +
  theme(strip.text = element_text(face = "bold"))

print(volume_reference_plot)
print(energy_reference_plot)


#Final test I hope 

# Environmental availability
env_quantiles <- env_reference %>%
  filter(lake %in% c("ES", "SP")) %>%
  pivot_longer(
    cols = c(volume_m3_h, E_kJ_h),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(is.finite(value), is.finite(vol_weight), vol_weight > 0) %>%
  group_by(lake, metric) %>%
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = quantile(value, 0.50, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    q95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distribution = "Environmental availability")

# Population-level observed use
obs_quantiles <- obs_fish_plot %>%
  pivot_longer(
    cols = c(volume_m3_h, E_kJ_h),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  group_by(lake, metric) %>%
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = quantile(value, 0.50, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    q95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distribution = "Observed population")

# Selected individuals
individual_quantiles <- individual_fish_plot %>%
  pivot_longer(
    cols = c(volume_m3_h, E_kJ_h),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  group_by(lake, type, fish_id, metric) %>%
  summarise(
    q10 = quantile(value, 0.10, na.rm = TRUE),
    median = quantile(value, 0.50, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    q95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(distribution = "Selected individual")

# Combine all summaries
quantile_comparison <- bind_rows(
  env_quantiles,
  obs_quantiles,
  individual_quantiles
) %>%
  select(
    lake, any_of(c("type", "fish_id")),
    metric, distribution, q10, median, q90, q95
  ) %>%
  arrange(metric, lake, distribution, type)

print(quantile_comparison, n = Inf)



#New fig goofinL: 
  





# --------------------------------------------------
# 1. Define the four fish to plot
# --------------------------------------------------

fish_groups <- tibble::tribble(
  ~lake, ~type,      ~fish_id,
  "ES",  "Average",  "ES_04",
  "ES",  "Extreme",  "ES_03",
  "SP",  "Average",  "SP_29",
  "SP",  "Extreme",  "SP_02"
)

# --------------------------------------------------
# 2. Prepare environmental availability
# --------------------------------------------------

env_density <- env_reference %>%
  filter(
    lake %in% c("ES", "SP"),
    is.finite(volume_m3_h),
    is.finite(E_kJ_h),
    is.finite(vol_weight),
    vol_weight > 0
  ) %>%
  mutate(lake = as.character(lake))

# --------------------------------------------------
# 3. Prepare observed fish data
# --------------------------------------------------

obs_density <- obs_fish %>%
  filter(
    lake %in% c("ES", "SP"),
    is.finite(volume_m3_h),
    is.finite(E_kJ_h),
    !is.na(fish_id)
  ) %>%
  mutate(
    lake = as.character(lake),
    fish_id = as.character(fish_id)
  )

# Selected individual observations
selected_fish <- obs_density %>%
  inner_join(
    fish_groups %>%
      mutate(
        lake = as.character(lake),
        fish_id = as.character(fish_id)
      ),
    by = c("lake", "fish_id")
  )

# --------------------------------------------------
# 4. Function to create the Average/Extreme figures
# --------------------------------------------------


plot_selected_fish <- function(metric, x_label) {
  
  # Panel order and labels
  panel_info <- fish_groups %>%
    mutate(
      lake = as.character(lake),
      fish_id = as.character(fish_id),
      type = as.character(type),
      panel = paste(lake, type, sep = "_"),
      panel_label = paste(lake, fish_id, sep = "\n")
    )
  
  panel_levels <- panel_info$panel
  
  # Population-level density, calculated once per lake
  population_density <- obs_density %>%
    group_by(lake) %>%
    group_modify(~ {
      vals <- .x[[metric]]
      vals <- vals[is.finite(vals)]
      
      d <- density(vals, n = 512)
      
      tibble(x = d$x, y = d$y)
    }) %>%
    ungroup()
  
  # Individual density for each selected fish
  individual_density <- selected_fish %>%
    group_by(lake, type, fish_id) %>%
    group_modify(~ {
      vals <- .x[[metric]]
      vals <- vals[is.finite(vals)]
      
      d <- density(vals, n = 512)
      
      tibble(x = d$x, y = d$y)
    }) %>%
    ungroup() %>%
    mutate(panel = paste(lake, type, sep = "_"))
  
  
  # Replicate population density curves into the two panels per lake
  population_panels <- dplyr::bind_rows(
    lapply(seq_len(nrow(panel_info)), function(i) {
      population_density %>%
        filter(lake == panel_info$lake[i]) %>%
        mutate(panel = panel_info$panel[i])
    })
  ) %>%
    mutate(panel = factor(panel, levels = panel_levels))
  
  
  # Replicate environmental observations into the two panels per lake
  env_panels <- dplyr::bind_rows(
    lapply(seq_len(nrow(panel_info)), function(i) {
      env_density %>%
        filter(lake == panel_info$lake[i]) %>%
        select(lake, vol_weight, all_of(metric)) %>%
        mutate(panel = panel_info$panel[i])
    })
  ) %>%
    rename(value = all_of(metric)) %>%
    mutate(panel = factor(panel, levels = panel_levels))
  
  # Draw the three distributions
  ggplot() +
    geom_density(
      data = env_panels,
      aes(x = value, weight = vol_weight),
      fill = "grey75",
      color = "grey50",
      alpha = 0.6,
      linewidth = 0.6
    ) +
    geom_line(
      data = population_panels,
      aes(x = x, y = y),
      color = "steelblue",
      linetype = "dashed",
      linewidth = 0.9
    ) +
    geom_line(
      data = individual_density,
      aes(x = x, y = y),
      color = "black",
      linewidth = 1
    ) +
    facet_grid(
      cols = vars(panel),
      scales = "free_x",
      labeller = labeller(
        panel = setNames(panel_info$panel_label, panel_info$panel)
      )
    ) +
    labs(
      x = x_label,
      y = "Density"
    ) +
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.7
      ),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11)
    )
}


# --------------------------------------------------
# 5. Create the two figures
# --------------------------------------------------
# Generate the four-panel plot for each metric
volume_average_extreme <- plot_selected_fish(
  "volume_m3_h",
  expression("Volume searched (m"^3*"/h)")
)

energy_average_extreme <- plot_selected_fish(
  "E_kJ_h",
  "Energetic requirements (kJ/h)"
)

# Remove individual plot titles so the combined figure is cleaner
volume_average_extreme <- volume_average_extreme +
  labs(title = NULL)

energy_average_extreme <- energy_average_extreme +
  labs(title = NULL)

# Combine into one figure
combined_figure <- volume_average_extreme /
  energy_average_extreme +
  plot_layout(heights = c(1, 1))

print(combined_figure)





#Some Bogus 


# Consistent panel order and labels
fish_groups_plot <- fish_groups %>%
  mutate(
    lake = as.character(lake),
    fish_id = as.character(fish_id),
    panel = factor(
      paste(lake, fish_id, sep = "_"),
      levels = paste(lake, fish_id, sep = "_")
    ),
    panel_label = paste(lake, type, sep = "\n")
  )

# Individual fish observations
individual_plot <- obs_fish %>%
  filter(
    is.finite(volume_m3_h),
    is.finite(E_kJ_h)
  ) %>%
  mutate(
    lake = as.character(lake),
    fish_id = as.character(fish_id)
  ) %>%
  inner_join(
    fish_groups_plot %>%
      select(lake, fish_id, panel, panel_label),
    by = c("lake", "fish_id")
  )

# Population-level use, repeated for each selected fish in its lake
population_plot <- obs_fish %>%
  filter(
    lake %in% c("ES", "SP"),
    is.finite(volume_m3_h),
    is.finite(E_kJ_h)
  ) %>%
  mutate(lake = as.character(lake)) %>%
  inner_join(
    fish_groups_plot %>%
      distinct(lake, panel, panel_label),
    by = "lake",
    relationship = "many-to-many"
  )

# Environmental availability, weighted by available habitat volume
available_plot <- env_pop %>%
  filter(
    lake %in% c("ES", "SP"),
    is.finite(volume_m3_h),
    is.finite(E_kJ_h),
    is.finite(vol_weight),
    vol_weight > 0
  ) %>%
  mutate(lake = as.character(lake)) %>%
  inner_join(
    fish_groups_plot %>%
      distinct(lake, panel, panel_label),
    by = "lake",
    relationship = "many-to-many"
  )

# Shared theme matching the reference layout
plot_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.6
    ),
    strip.text = element_text(
      face = "bold",
      size = 10,
      lineheight = 0.95,
      margin = margin(3, 2, 3, 2)
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(color = "black", size = 9),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5)),
    panel.spacing.x = grid::unit(0.15, "cm"),
    panel.spacing.y = grid::unit(0.35, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

# Function to plot each metric
plot_metric <- function(metric, x_label) {
  
  ggplot() +
    # Available habitat
    geom_density(
      data = available_plot,
      aes(
        x = .data[[metric]],
        weight = vol_weight
      ),
      fill = "grey75",
      color = "grey50",
      alpha = 0.65,
      linewidth = 0.6,
      adjust = 1
    ) +
    
    # Population use
    geom_density(
      data = population_plot,
      aes(x = .data[[metric]]),
      color = "steelblue",
      linetype = "dashed",
      linewidth = 0.9,
      adjust = 1
    ) +
    
    # Individual use
    geom_density(
      data = individual_plot,
      aes(x = .data[[metric]]),
      color = "black",
      linewidth = 1.0,
      adjust = 1
    ) +
    
    facet_wrap(
      ~ panel,
      nrow = 1,
      scales = "free_x",
      labeller = labeller(
        panel = stats::setNames(
          as.character(fish_groups_plot$panel_label),
          as.character(fish_groups_plot$panel)
        )
      )
    ) +
    
    labs(
      x = x_label,
      y = "Density"
    ) +
    
    plot_theme
}

# Top row: volume searched
volume_average_extreme <- plot_metric(
  "volume_m3_h",
  expression("Volume searched (m"^3*"/h)")
)

# Bottom row: energetic requirements
energy_average_extreme <- plot_metric(
  "E_kJ_h",
  "Energetic requirement (kJ/h)"
)

# Combine into one figure, without a legend
combined_average_extreme <-
  volume_average_extreme /
  energy_average_extreme +
  plot_layout(
    ncol = 1,
    heights = c(1, 1)
  )

combined_average_extreme

