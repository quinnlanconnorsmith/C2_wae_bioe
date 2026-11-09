# ============================================================
# 03_consequences.R
# Consequences of habitat use: energetic requirements (FB4) and
# effective volume searched, for individuals vs the population vs
# lake-wide availability; Pareto frontiers; trajectories.
#
# Requires derived/*.rds from 01_data_prep.R.
# Start from a fresh R session - this script does not need
# anything from 02_habitat_use.R.
# ============================================================

source("R/00_functions.R")

respiration            <- read_derived("respiration")
fish_traits            <- read_derived("fish_traits")
valid_fish             <- read_derived("valid_fish")
archival_essp_hourly   <- read_derived("archival_essp_hourly")
archival_hourly_search <- read_derived("archival_hourly_search")
hypso_layers           <- read_derived("hypso_layers")
env_hourly             <- read_derived("env_hourly")

# Observed fish positions for the Pareto frontier (saved earlier)
obs_fish_plot <- read_csv("obs_fish_plot.csv")

####Bioenergetics####

# 7. Output summary check

respiration %>%
  group_by(fish_id) %>%
  summarise(
    mean_temp = mean(mean_temp),
    mean_R = mean(R_gO2_g_h),
    mass_g = unique(round(mass_g, 1))
  )

#Running some checks: 

# 1. Check for missing values
colSums(is.na(respiration))

# 2. Quick summaries of main numeric fields
summary(respiration %>% select(mean_temp, tl_mm, mass_g, R_gO2_g_h))

# 3. Make sure hourly times are sequential for each fish
respiration %>%
  group_by(fish_id) %>%
  summarise(
    n_hours = n(),
    first_time = min(hour_time),
    last_time = max(hour_time)
  ) %>%
  arrange(desc(n_hours))

# 4. quick check plot of respiration vs temperature
ggplot(respiration, aes(x = mean_temp, y = R_gO2_g_h, color = factor(fish_id))) +
  geom_point(alpha = 0.6) +
  #geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Temperature (°C)", y = "Respiration (g O₂/g/h)", color = "Fish ID") +
  theme_minimal()

#Get some more energy flowing 

str(respiration)

df_energy <- respiration %>%
  mutate(
    # Convert O₂ to energy per g fish per hour
    E_kJ_g_h = R_gO2_g_h * 13.6,
    
    # Energy per fish per hour (kJ/h)
    E_kJ_h = E_kJ_g_h * mass_g,
    
    # Energy per fish per day (kJ/day)
    E_kJ_day = E_kJ_h * 24,
    
    # Add date and month columns
    date = as.Date(hour_time),
    month = lubridate::month(hour_time, label = TRUE)
  )

#Standarize to weight 
df_energy <- df_energy %>%
  mutate(
    E_J_day = E_kJ_day * 1000,                # convert kJ to J
    E_J_g_day = E_J_day / mass_g              # standardize by mass (J/g/day)
  )

# Daily
daily_energy <- df_energy %>%
  group_by(fish_id, date = as.Date(hour_time)) %>%
  summarise(
    mean_temp = mean(mean_temp, na.rm = TRUE),
    mean_E_kJ_day = mean(E_kJ_day, na.rm = TRUE),
    mean_E_J_g_day = mean(E_J_g_day, na.rm = TRUE),
    mass_g = mean(mass_g, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(lake = case_when(
    grepl("ES", fish_id) ~ "Escanaba",
    grepl("SP", fish_id) ~ "Sparkling",
    grepl("MC", fish_id) ~ "McDermott"
  ))

# Monthly
monthly_energy <- daily_energy %>%
  mutate(month = lubridate::month(date, label = TRUE)) %>%
  group_by(fish_id, lake, month) %>%
  summarise(
    mean_temp = mean(mean_temp, na.rm = TRUE),
    mean_E_kJ_day = mean(mean_E_kJ_day, na.rm = TRUE),
    mean_E_J_g_day = mean(mean_E_J_g_day, na.rm = TRUE),
    mass_g = mean(mass_g, na.rm = TRUE)
  ) %>%
  ungroup()

# Summer (one average per fish)
summer_energy <- daily_energy %>%
  group_by(fish_id, lake) %>%
  summarise(
    mean_temp = mean(mean_temp, na.rm = TRUE),
    mean_E_kJ_day = mean(mean_E_kJ_day, na.rm = TRUE),
    mean_E_J_g_day = mean(mean_E_J_g_day, na.rm = TRUE),
    mass_g = mean(mass_g, na.rm = TRUE)
  ) %>%
  ungroup()

daily_energy <- daily_energy %>%
  mutate(
    month = month(date)  # e.g., "Jun", "Jul", "Aug"
  )

#Add lake code
df_energy <- df_energy %>%
  mutate(
    lake = stringr::str_extract(fish_id, "^[A-Za-z]+")
  )

individual_energy <- daily_energy %>%
  group_by(fish_id) %>%
  summarise(
    mean_temp = mean(mean_temp, na.rm = TRUE),
    mean_E_J_g_day = mean(mean_E_J_g_day, na.rm = TRUE),
    mean_E_kf_day = mean(mean_E_kJ_day, na.rm = TRUE)
  ) %>%
  # Step 3: Extract lake code from wae_id
  mutate(lake = stringr::str_extract(fish_id, "^[A-Za-z]+"))

####Swim speed from field metabolism####

str(respiration$hour_time)
alpha <- 2  # field metabolic multiplier (can test 1.5–2.5 later)

df_energy_swim <- respiration %>%
  mutate(
    # ensure datetime is proper
    hour_time = as.POSIXct(hour_time, tz = "UTC"),
    
    # Energetics
    E_J_g_h = R_gO2_g_h * 13.6 * 1000,
    E_J_h = E_J_g_h * mass_g,
    E_J_day = E_J_h * 24,
    
    # Date fields
    date = as.Date(hour_time),
    month = lubridate::month(hour_time, label = TRUE),
    
    # Field metabolism
    E_J_day_field = E_J_day * alpha,
    swim_J_day = E_J_day_field - E_J_day,
    swim_J_h = swim_J_day / 24
  )


# Find mean swim energy around 20°C
ref_energy <- df_energy_swim %>%
  filter(mean_temp > 19 & mean_temp < 21) %>%
  summarise(mean_swim = mean(swim_J_h, na.rm = TRUE)) %>%
  pull(mean_swim)

U_target <- 0.4  # realistic adult walleye cruising speed (m/s)

k_drag <- ref_energy / (U_target^3)

df_energy_swim <- df_energy_swim %>%
  mutate(
    swim_speed_m_s = (swim_J_h / k_drag)^(1/3),
    
    # Cap at biological maximum
    swim_speed_m_s = pmin(swim_speed_m_s, 1)
  )

summary(df_energy_swim$swim_speed_m_s)

####Daily volume searched vs energy required####

fish_search_means <- archival_hourly_search %>%
  group_by(lake, fish_id) %>%
  summarise(
    fish_mean_search = mean(mean_search, na.rm = TRUE),
    .groups = "drop"
  )

archival_hourly_search <- archival_hourly_search %>%
  left_join(fish_search_means, by = c("lake", "fish_id")) %>%
  arrange(lake, fish_mean_search) %>%
  mutate(fish_id = factor(fish_id, levels = unique(fish_id)))

#Plotting the two together 
#Just merge hourly search - lots of little datasets  

df_volume_energy_hourly <- df_energy_swim %>%
  left_join(archival_hourly_search,
            by = c("fish_id", "hour_time"))

df_volume_energy_daily <- df_volume_energy_hourly %>%
  group_by(lake, fish_id, date) %>%
  summarise(
    daily_volume = sum(mean_search, na.rm = TRUE),
    daily_swim_J = sum(swim_J_h, na.rm = TRUE),
    daily_E_J = sum(E_J_h, na.rm = TRUE),
    daily_E_J_g = sum(E_J_g_h, na.rm = TRUE),  # convert hourly to daily - sum or mean*24??
    .groups = "drop"
  )

#Plot faily 
ggplot(  na.omit(df_volume_energy_daily[, c("daily_volume",
                                            "daily_E_J",
                                            "fish_id",
                                            "lake")]),
         aes(x = daily_E_J,
             y = daily_volume, 
             color = fish_id)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~ lake) +
  theme_classic() +
  labs(x = "Daily Energy Required (J)",
       y = "Daily Volume Searched (m³)",
       title = "Daily Search Effort vs Energy Required")


#Compute centroids and plot: 

df_daily_sve <- df_volume_energy_daily %>%
  select(daily_volume, daily_E_J, fish_id, lake) %>%
  na.omit() %>%
  filter(daily_volume > 0)

centroids <- df_daily_sve %>%
  group_by(lake) %>%
  summarise(
    centroid_E = mean(daily_E_J, na.rm = TRUE),
    centroid_V = mean(daily_volume, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_daily_sve,
       aes(x = daily_E_J / 1000,
           y = daily_volume,
           color = fish_id)) +
  
  geom_point(alpha = 0.4) +
  
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  
  geom_point(data = centroids,
             aes(x = centroid_E / 1000, y = centroid_V),
             inherit.aes = FALSE,
             color = "black",
             size = 4,
             shape =19,
             stroke = 1.5) +
  
  facet_wrap(~ lake) +
  theme_classic() +
  theme(legend.position = "none") +
  
  labs(x = "Daily Energy Required (kJ)",
       y = "Daily Volume Searched (m³)",
       title = "Daily Search Effort vs Energy Required")





df_volume_energy_daily_clean <- df_volume_energy_daily %>%
  dplyr::filter(
    !is.na(daily_volume),
    !is.na(daily_E_J)
  )

df_volume_energy_cv <- df_volume_energy_daily_clean %>%
  group_by(lake, fish_id) %>%
  summarise(
    
    mean_volume = mean(daily_volume),
    sd_volume = sd(daily_volume),
    cv_volume = sd_volume / mean_volume,
    
    mean_energy = mean(daily_E_J),
    sd_energy = sd(daily_E_J),
    cv_energy = sd_energy / mean_energy,
    
    mean_efficiency = mean(daily_volume / daily_E_J),
    sd_efficiency = sd(daily_volume / daily_E_J),
    cv_efficiency = sd_efficiency / mean_efficiency,
    
    .groups = "drop"
  )

#4-way bogus? 

ggplot(df_volume_energy_cv,
       aes(x = mean_energy,
           y = mean_volume)) +
  
  geom_point(aes(size = cv_volume,
                 color = cv_energy),
             alpha = 0.9) +
  
  scale_color_viridis_c() +
  
  theme_classic() +
  facet_wrap(~lake)+
  
  labs(x = "Mean Daily Energy Required (J)",
       y = "Mean Daily Volume Searched (m³)",
       size = "CV Search Volume",
       color = "CV Energy",
       title = "Energetic Costs, Search Effort, and Behavioral Variability")

#Behavioral strategy map 

df_volume_energy_cv <- df_volume_energy_daily %>%
  dplyr::filter(!is.na(daily_volume),
                !is.na(daily_E_J)) %>%
  group_by(lake, fish_id) %>%
  summarise(
    mean_volume = mean(daily_volume, na.rm = TRUE),
    sd_volume = sd(daily_volume, na.rm = TRUE),
    cv_volume = sd_volume / mean_volume,
    
    mean_energy = mean(daily_E_J, na.rm = TRUE),
    sd_energy = sd(daily_E_J, na.rm = TRUE),
    cv_energy = sd_energy / mean_energy,
    
    mean_ratio = mean(daily_volume / daily_E_J, na.rm = TRUE),
    sd_ratio = sd(daily_volume / daily_E_J, na.rm = TRUE),
    cv_ratio = sd_ratio / mean_ratio,
    .groups = "drop"
  )

df_strategy <- df_volume_energy_cv %>%
  group_by(lake) %>%
  mutate(
    rel_energy = mean_energy / mean(mean_energy, na.rm = TRUE),
    rel_volume = mean_volume / mean(mean_volume, na.rm = TRUE),
    rel_ratio  = mean_ratio / mean(mean_ratio, na.rm = TRUE)
  ) %>%
  ungroup()

ggplot(df_strategy,
       aes(x = rel_energy,
           y = rel_volume)) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "grey40") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey40") +
  geom_point(aes(size = cv_volume,
                 color = cv_energy),
             alpha = 0.9) +
  facet_wrap(~ lake) +
  scale_color_viridis_c() +
  theme_classic() +
  labs(x = "Relative Mean Energy Required",
       y = "Relative Mean Volume Searched",
       size = "CV Volume",
       color = "CV Energy",
       title = "Behavioral Strategy Map")

df_strategy_clean <- df_strategy %>%
  dplyr::filter(!is.na(lake))

# define quadrant rectangles
quad_df <- expand.grid(
  lake = unique(df_strategy_clean$lake),
  quadrant = c("UL", "LR", "LL", "UR")
) %>%
  mutate(
    xmin = case_when(
      quadrant %in% c("UL", "LL") ~ -Inf,
      TRUE ~ 1
    ),
    xmax = case_when(
      quadrant %in% c("UL", "LL") ~ 1,
      TRUE ~ Inf
    ),
    ymin = case_when(
      quadrant %in% c("LL", "LR") ~ -Inf,
      TRUE ~ 1
    ),
    ymax = case_when(
      quadrant %in% c("LL", "LR") ~ 1,
      TRUE ~ Inf
    ),
    fill_col = case_when(
      quadrant == "UL" ~ "green",
      quadrant == "LR" ~ "red",
      TRUE ~ "yellow"
    )
  )

ggplot(df_strategy_clean,
       aes(x = rel_energy,
           y = rel_volume)) +
  
  # Quadrant shading (NOW WORKS)
  geom_rect(data = quad_df,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = fill_col),
            alpha = 0.2,
            inherit.aes = FALSE) +
  
  scale_fill_identity() +
  
  # Reference lines
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "grey40") +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey40") +
  
  # Points
  geom_point(aes(size = cv_energy,
                 color = cv_volume),
             alpha = 0.9) +
  
  facet_wrap(~ lake) +
  scale_color_viridis_c(option="mako") +
  
  labs(x = "Relative Mean Energy Required",
       y = "Relative Mean Volume Searched",
       size = "CV Energy",
       color = "CV Volume") +
  
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

####Lake-wide environmental performance (availability)####

# Every hour x 0.1 m cell in each lake, with swim speed, reaction
# distance and volume searched. Volume searched does not depend on
# fish size, so this one table serves every fish.

env_perf <- env_hourly %>%
  filter(!is.na(lake)) %>%
  add_search_performance() %>%
  arrange(lake, depth_m)

# Attach hypsography: each depth gets the volume of its 1 m layer
hypso_lookup <- hypso_layers %>%
  select(lake, depth_top, depth_bottom, layer_volume_m3)

env_perf <- env_perf %>%
  group_split(lake) %>%
  lapply(function(df) {
    hyp <- hypso_lookup %>% filter(lake == unique(df$lake))
    if (nrow(hyp) == 0) stop(paste("No hypsography for lake:", unique(df$lake)))
    idx <- findInterval(df$depth_m, hyp$depth_top)
    idx[idx == 0] <- 1
    idx[idx > nrow(hyp)] <- nrow(hyp)
    df$layer_volume_m3 <- hyp$layer_volume_m3[idx]
    df
  }) %>%
  bind_rows()

# Split each 1 m layer volume among its 0.1 m cells, then
# normalise to weights that sum to 1 per lake-hour
env_perf <- env_perf %>%
  mutate(depth_bin = floor(depth_m)) %>%
  group_by(datetime_hour, lake, depth_bin) %>%
  mutate(n_sublayers = n(),
         sublayer_volume = layer_volume_m3 / n_sublayers) %>%
  group_by(datetime_hour, lake) %>%
  mutate(vol_weight = sublayer_volume / sum(sublayer_volume, na.rm = TRUE)) %>%
  ungroup()

# Checks: weights should sum to 1 within each lake-hour
env_perf %>%
  group_by(lake, datetime_hour) %>%
  summarise(sum_weight = sum(vol_weight, na.rm = TRUE), .groups = "drop") %>%
  group_by(lake) %>%
  summarise(min_sum = min(sum_weight),
            median_sum = median(sum_weight),
            max_sum = max(sum_weight),
            .groups = "drop")

env_perf %>% count(lake)

# ------------------------------------------------------------
# Availability of energetic requirements
#
# Energetic requirement depends on fish mass, so availability is
# computed for each fish. The old version crossed every hour x
# depth cell with every fish (env_pop), which multiplied the
# largest table by the number of fish - the main reason the
# workspace got so big. Energy depends only on temperature and
# mass, so we first sum the volume weights by lake x temperature
# (to 0.01 C), then cross with the fish. The weighted
# distributions are the same, with a far smaller table.
# ------------------------------------------------------------

env_pop <- env_perf %>%
  mutate(temp_c = round(temp_c, 2)) %>%
  group_by(lake, temp_c) %>%
  summarise(vol_weight = sum(vol_weight, na.rm = TRUE), .groups = "drop") %>%
  inner_join(valid_fish, by = "lake", relationship = "many-to-many") %>%
  mutate(mass_g = length_mm_to_mass_g(tl_pre),
         E_kJ_h = energy_kJ_h(mass_g, temp_c))

# Availability for a single reference-size fish (440 mm)
tl_ref <- 440
mass_ref_g <- length_mm_to_mass_g(tl_ref)

env_reference <- env_perf %>%
  select(lake, datetime_hour, depth_m, temp_c, light_lux, vol_weight,
         thermal_perf, swim_speed_m_s, reaction_dist_m,
         distance_m_h, volume_m3_h) %>%
  mutate(mass_g = mass_ref_g,
         R_gO2_g_h = respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT),
         E_kJ_h = (R_gO2_g_h * OXY_KJ) * mass_g)

env_reference %>%
  group_by(lake) %>%
  summarise(n_records = n(),
            median_volume = median(volume_m3_h, na.rm = TRUE),
            median_energy = median(E_kJ_h, na.rm = TRUE),
            max_volume = max(volume_m3_h, na.rm = TRUE),
            max_energy = max(E_kJ_h, na.rm = TRUE),
            .groups = "drop")

# ------------------------------------------------------------
# Observed (used) energetic requirement and volume searched
# ------------------------------------------------------------

obs_fish <- archival_essp_hourly %>%
  filter(fish_id %in% valid_fish$fish_id) %>%
  left_join(valid_fish, by = c("fish_id", "lake")) %>%
  mutate(mass_g = length_mm_to_mass_g(tl_pre),
         R_gO2_g_h = respiration_FB4(mass_g, mean_temp, RA, RB, RQ, RTO, RTM, ACT),
         E_kJ_h = (R_gO2_g_h * OXY_KJ) * mass_g,
         volume_m3_h = mean_volume)

####Pareto Frontier####

# 1. Unique fish sizes

fish_sizes <- fish_traits %>%
  select(fish_id, lake, tl_pre) %>%
  distinct()

# 2. Build env space per fish
# Volume searched comes straight from env_perf (it does not depend
# on fish size); only the energetic requirement is fish-specific.
build_env_space_fish <- function(fish_id_val, tl_val, lake_val) {

  env_perf %>%
    filter(lake == lake_val) %>%
    mutate(
      fish_id = fish_id_val,
      mass_g = length_mm_to_mass_g(tl_val),
      E_kJ_h = energy_kJ_h(mass_g, pmin(temp_c, RTM - 0.5)),
      volume_searched_m3_h = volume_m3_h
    ) %>%
    filter(
      !is.na(E_kJ_h),
      !is.na(volume_searched_m3_h),
      E_kJ_h > 0,
      E_kJ_h < 50,
      light_lux > 0
    ) %>%
    # reduce size early
    transmute(
      fish_id,
      lake,
      E_kJ_h = round(E_kJ_h, 2),
      volume_searched_m3_h = round(volume_searched_m3_h, 2)
    ) %>%
    distinct()
}

# 3. Run for all fish
env_space_fish <- purrr::pmap_dfr(
  fish_sizes,
  ~ build_env_space_fish(..1, ..3, ..2)
)


# 4. Build envelope per fish

envelope_fish <- env_space_fish %>%
  group_by(fish_id, lake, E_kJ_h) %>%
  summarise(
    vol_max = quantile(volume_searched_m3_h, 0.95),
    vol_min = quantile(volume_searched_m3_h, 0.05),
    .groups = "drop"
  )

# 5. Plot (facet by fish)

p_fish_individual <- ggplot() +
  geom_line(
    data = envelope_fish,
    aes(x = E_kJ_h, y = vol_max, color = fish_id),
    linewidth = 1
  ) +
  geom_point(
    data = obs_fish_plot,
    aes(x = E_kJ_h, y = volume_searched_m3_h, color = fish_id),
    size = 2
  ) +
  # `envelope_df` was never created in the script; smoothing the
  # per-fish envelope instead
  geom_smooth(
    data = envelope_fish,
    aes(x = E_kJ_h, y = vol_max),
    method = "loess",
    se = FALSE
  ) + 
  coord_cartesian(xlim = c(0, 30)) +
  facet_wrap(~ fish_id) +
  theme_classic() +
  labs(
    x = "Energy requirement (kJ/h)",
    y = "Volume searched (m3/h)",
    title = "Individual-specific feasible space and observed position"
  )

print(p_fish_individual)



#Calculating diatnace to frontier :


# 1. Match each fish to its frontier

distance_df <- obs_fish_plot %>%
  left_join(
    envelope_fish,
    by = c("fish_id", "lake")
  ) %>%
  
  # find closest energy match
  group_by(fish_id) %>%
  mutate(
    energy_diff = abs(E_kJ_h.x - E_kJ_h.y)
  ) %>%
  
  slice_min(energy_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  rename(
    E_obs = E_kJ_h.x,
    V_obs = volume_searched_m3_h,
    E_frontier = E_kJ_h.y,
    V_max = vol_max
  ) %>%
  
# 2. Calculate distance

mutate(
  dist_to_frontier = V_max - V_obs,
  
  # optional normalized version
  rel_efficiency = V_obs / V_max
)

#Vis
ggplot(distance_df, aes(x = fish_id, y = rel_efficiency, fill = lake)) +
  geom_col() +
  theme_classic() +
  labs(
    y = "Relative efficiency (V_obs / V_max)",
    x = "Fish ID"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#More vis
ggplot(distance_df, aes(x = fish_id, y = dist_to_frontier, fill = lake)) +
  geom_col() +
  theme_classic() +
  labs(
    y = "Distance to frontier (m³/h)",
    x = "Fish ID"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Try to build out the same kind of metric for horizontal and vertical 

vertical_df <- obs_fish_plot %>%
  group_by(fish_id, lake) %>%
  summarise(
    E_obs = first(E_kJ_h),
    V_obs = first(volume_searched_m3_h),
    
    V_max = {
      sub <- envelope_fish[
        envelope_fish$fish_id == first(fish_id) &
          envelope_fish$lake == first(lake),
      ]
      
      if (nrow(sub) == 0) {
        NA_real_
      } else {
        sub$vol_max[which.min(abs(sub$E_kJ_h - first(E_kJ_h)))]
      }
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    dist_vertical = V_max - V_obs,
    rel_efficiency = V_obs / V_max
  )

horizontal_frontier <- env_space_fish %>%
  group_by(fish_id, lake, volume_searched_m3_h) %>%
  summarise(
    E_min = quantile(E_kJ_h, 0.05, na.rm = TRUE),
    .groups = "drop"
  )

horizontal_df <- obs_fish_plot %>%
  group_by(fish_id, lake) %>%
  summarise(
    E_obs = first(E_kJ_h),
    V_obs = first(volume_searched_m3_h),
    
    E_min = {
      sub <- horizontal_frontier[
        horizontal_frontier$fish_id == first(fish_id) &
          horizontal_frontier$lake == first(lake),
      ]
      
      sub <- sub[sub$E_min <= first(E_kJ_h), ]
      
      if (nrow(sub) == 0) {
        NA_real_
      } else {
        sub$E_min[which.min(abs(sub$volume_searched_m3_h - first(volume_searched_m3_h)))]
      }
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    dist_horizontal = E_obs - E_min,
    rel_energy_eff = E_min / E_obs
  )

pareto_df <- vertical_df %>%
  left_join(
    horizontal_df %>%
      select(fish_id, dist_horizontal, rel_energy_eff),
    by = "fish_id"
  ) %>%
  mutate(
    joint_efficiency = rel_efficiency * rel_energy_eff
  )

pareto_df <- pareto_df %>%
  group_by(lake) %>%
  mutate(fish_order = fct_reorder(fish_id, rel_efficiency)) %>%
  ungroup()

ggplot(pareto_df, aes(x = fish_order, y = dist_horizontal, fill = lake)) +
  geom_col() +
  facet_wrap(~ lake, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Plot the b's 


pareto_df <- pareto_df %>%
  group_by(lake) %>%
  arrange(rel_efficiency, .by_group = TRUE) %>%
  mutate(fish_order = factor(fish_id, levels = unique(fish_id))) %>%
  ungroup()


ggplot(pareto_df, aes(x = fish_order, y = rel_efficiency, fill = lake)) +
  geom_col() +
  facet_wrap(~ lake, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Fish ID",
    y = "Relative search efficiency",
    title = "Fish  search efficiency relative to Pareto frontier"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pareto_df <- pareto_df %>%
  group_by(lake) %>%
  arrange(rel_energy_eff, .by_group = TRUE) %>%
  mutate(fish_order = factor(fish_id, levels = unique(fish_id))) %>%
  ungroup()

ggplot(pareto_df, aes(x = fish_order, y = rel_energy_eff, fill = lake)) +
  geom_col() +
  facet_wrap(~ lake, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Fish ID",
    y = "Relative energy efficiency",
    title = "Energetic inefficiency across fish"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pareto_df <- pareto_df %>%
  group_by(lake) %>%
  arrange(joint_efficiency, .by_group = TRUE) %>%
  mutate(fish_order_joint = factor(fish_id, levels = unique(fish_id))) %>%
  ungroup()

ggplot(pareto_df, aes(x = fish_order_joint, y = joint_efficiency, fill = lake)) +
  geom_col() +
  facet_wrap(~ lake, scales = "free_x") +
  theme_classic() +
  labs(
    x = "Fish ID",
    y = "Joint efficiency",
    title = "Combined energetic + search performance efficiency"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p_2d <- ggplot(pareto_df, aes(x = dist_horizontal, y = dist_vertical, color = lake)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic() +
  labs(
    x = "Excess energy (kJ/h)",
    y = "Lost search volume (m³/h)",
    title = "Distance from Pareto frontier (search & energy)"
  )

print(p_2d)

p_2d_labels <- p_2d +
  geom_text(aes(label = fish_id), size = 3, vjust = -0.5) +
  geom_vline(xintercept = median(pareto_df$dist_horizontal, na.rm = TRUE), linetype = "dashed") +
  geom_hline(yintercept = median(pareto_df$dist_vertical, na.rm = TRUE), linetype = "dashed")

print(p_2d_labels)


#Size class pareto frontiers: 

pareto_bin_df <- pareto_df %>%
  left_join(fish_traits, by = c("fish_id", "lake")) %>%
  mutate(
    size_class = case_when(
      tl_pre >= 310 & tl_pre <= 380 ~ "small",
      tl_pre >= 381 & tl_pre <= 456 ~ "medium",
      tl_pre > 456 ~ "large"
    ),
    size_class = factor(size_class, levels = c("small", "medium", "large"))
  )

size_bounds <- tibble::tibble(
  size_class = c("small", "medium", "large"),
  min_len = c(310, 381, 457),
  max_len = c(380, 456, 631)
)

#Pareto building function: 

build_frontier <- function(tl_val) {
  env_perf %>%
    mutate(
      mass_g = length_mm_to_mass_g(tl_val),
      E_kJ_h = energy_kJ_h(mass_g, pmin(temp_c, RTM - 0.5)),
      volume_searched_m3_h = volume_m3_h
    ) %>%
    filter(
      !is.na(E_kJ_h),
      !is.na(volume_searched_m3_h),
      E_kJ_h > 0,
      E_kJ_h < 50,
      light_lux > 0
    ) %>%
    group_by(lake, E_bin = round(E_kJ_h, 1)) %>%
    summarise(
      E_kJ_h = mean(E_kJ_h),
      V_max = quantile(volume_searched_m3_h, 0.95),
      .groups = "drop"
    )
}


frontier_min_df <- size_bounds %>%
  mutate(frontier = map(min_len, build_frontier)) %>%
  tidyr::unnest(frontier) %>%
  mutate(bound = "min")

frontier_max_df <- size_bounds %>%
  mutate(frontier = map(max_len, build_frontier)) %>%
  tidyr::unnest(frontier) %>%
  mutate(bound = "max")

frontier_bin_df <- bind_rows(frontier_min_df, frontier_max_df)

# Ribbon between the smallest and largest fish in each size class.
# Matched on the 0.1 kJ/h energy bin (E_kJ_h itself differs between
# the two sizes, so it cannot be used as the key).
ribbon_df <- frontier_bin_df %>%
  select(size_class, lake, E_bin, bound, V_max) %>%
  pivot_wider(names_from = bound, values_from = V_max)

ggplot() +
  
  geom_line(
    data = frontier_bin_df,
    aes(x = E_kJ_h, y = V_max, linetype = bound),
    linewidth = 1.2
  ) +
  
  geom_point(
    data = pareto_bin_df,
    aes(x = E_obs, y = V_obs, color = lake),
    size = 3,
    alpha = 0.8
  ) +
  
  facet_grid(size_class ~ lake, scales = "free") +
  
  xlim(0,30)+
  geom_ribbon(
    data = ribbon_df,
    aes(x = E_bin, ymin = min, ymax = max),
    alpha = 0.15
  ) +
  
  theme_classic() +
  labs(
    x = "Energy (kJ/h)",
    y = "Search volume (m³/h)",
    title = "Size-class Pareto frontiers"
  )

####Trajectories for consequences####

#Try again 


# 1. SELECT FISH

fish_id_sel <- "SP_40"
day_sel     <- as.Date("2022-06-29")

fish_info <- fish_traits %>%
  filter(fish_id == fish_id_sel) %>%
  slice(1)

lake_sel <- fish_info$lake

# 2. GET FISH MASS
mass_sel_g <- length_mm_to_mass_g(fish_info$tl_pre)


# 3. SUBSET ENVIRONMENT

env_day <- env_hourly %>%
  mutate(date = as.Date(datetime_hour)) %>%
  filter(lake == lake_sel, date == day_sel)


# 4. BUILD FULL DEPTH GRID

depth_seq <- round(seq(
  floor(min(env_day$depth_m, na.rm = TRUE)*10)/10,
  ceiling(max(env_day$depth_m, na.rm = TRUE)*10)/10,
  by = 0.1
), 1)

env_full <- env_day %>%
  complete(datetime_hour, depth_m = depth_seq)


# 5. INTERPOLATE TEMP + LIGHT

env_interp <- env_full %>%
  # interpolate across depth
  arrange(datetime_hour, depth_m) %>%
  group_by(datetime_hour) %>%
  mutate(
    temp_c    = zoo::na.approx(temp_c, na.rm = FALSE, rule = 2),
    light_lux = zoo::na.approx(light_lux, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup() %>%
  # interpolate across time
  arrange(depth_m, datetime_hour) %>%
  group_by(depth_m) %>%
  mutate(
    temp_c    = zoo::na.approx(temp_c, na.rm = FALSE, rule = 2),
    light_lux = zoo::na.approx(light_lux, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup()


# 6. COMPUTE PERFORMANCE
# (named env_perf_day so it does not overwrite the full-season env_perf)

env_perf_day <- env_interp %>%
  add_search_performance() %>%
  mutate(
    mass_g = mass_sel_g,
    E_kJ_h = energy_kJ_h(mass_g, temp_c),
    efficiency = volume_m3_h / E_kJ_h
  )


# 7. GET FISH TRAJECTORY

fish_day <- archival_essp_hourly %>%
  mutate(date = as.Date(hour_time)) %>%
  filter(
    fish_id == fish_id_sel,
    lake == lake_sel,
    date == day_sel
  )


# 8. PLOTS


# Volume plot 
ggplot() +
  geom_raster(
    data = env_perf_day,
    aes(datetime_hour, depth_m, fill = volume_m3_h),
    interpolate = TRUE
  ) +
  geom_line(
    data = fish_day,
    aes(hour_time, mean_depth),
    color = "white",
    linewidth = 1.2
  ) +
  scale_y_reverse() +
  scale_fill_viridis_c(name = "Volume (m³/h)") +
  theme_classic() +
  labs(title = paste("Volume searched |", fish_id_sel))

# nergy plot 
ggplot() +
  geom_raster(
    data = env_perf_day,
    aes(datetime_hour, depth_m, fill = E_kJ_h),
    interpolate = TRUE
  ) +
  geom_line(
    data = fish_day,
    aes(hour_time, mean_depth),
    color = "white",
    linewidth = 1.2
  ) +
  scale_y_reverse() +
  scale_fill_viridis_c(name = "Energy (kJ/h)") +
  theme_classic() +
  labs(title = paste("Energy requirement |", fish_id_sel))

# Efficiency plot (bonus)
ggplot() +
  geom_raster(
    data = env_perf_day,
    aes(datetime_hour, depth_m, fill = efficiency),
    interpolate = TRUE
  ) +
  geom_line(
    data = fish_day,
    aes(hour_time, mean_depth),
    color = "white",
    linewidth = 1.2
  ) +
  scale_y_reverse() +
  scale_fill_viridis_c(name = "Efficiency (m³/kJ)") +
  theme_classic() +
  labs(title = paste("Search efficiency |", fish_id_sel))

####Distributions: availability vs population vs individuals####

#Nicer plots methinks: 


# --------------------------------------------------
# 1. Prepare environmental availability
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
# 2. Prepare ALL observed fish
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


# --------------------------------------------------
# 3. Plotting function
# --------------------------------------------------

plot_all_fish <- function(lake_name, metric, x_label) {
  
  # Environmental availability
  env_lake <- env_reference %>%
    filter(
      lake == lake_name,
      is.finite(.data[[metric]]),
      is.finite(vol_weight),
      vol_weight > 0
    )
  
  # All observed fish in this lake
  obs_lake <- obs_fish %>%
    filter(
      lake == lake_name,
      is.finite(.data[[metric]]),
      !is.na(fish_id)
    ) %>%
    mutate(fish_id = as.character(fish_id))
  
  # Calculate the population-level density independently
  pop_density <- density(
    obs_lake[[metric]],
    na.rm = TRUE,
    n = 512
  )
  
  pop_df <- data.frame(
    x = pop_density$x,
    y = pop_density$y
  )
  
  # Precompute each individual fish's density
  individual_df <- obs_lake %>%
    group_by(fish_id) %>%
    group_modify(~ {
      d <- density(.x[[metric]], na.rm = TRUE, n = 512)
      data.frame(x = d$x, y = d$y)
    }) %>%
    ungroup()
  
  ggplot() +
    # Environmental availability
    geom_density(
      data = env_lake,
      aes(
        x = .data[[metric]],
        weight = vol_weight
      ),
      fill = "grey75",
      color = "grey50",
      alpha = 0.6,
      linewidth = 0.5
    ) +
    
    # Population-level observed use
    geom_line(
      data = pop_df,
      aes(x = x, y = y),
      color = "steelblue",
      linetype = "dashed",
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    
    # Individual fish
    geom_line(
      data = individual_df,
      aes(x = x, y = y, group = fish_id),
      color = "black",
      linewidth = 0.65
    ) +
    
    facet_wrap(~fish_id, ncol = 6) +
    labs(
      title = paste(lake_name, "-", x_label),
      x = x_label,
      y = "Density"
    ) +
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      strip.background = element_rect(
        fill = "white", color = "black", linewidth = 0.5
      ),
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      panel.spacing = grid::unit(0.5, "lines")
    )
}

# Four plots
ES_volume_all <- plot_all_fish(
  "ES", "volume_m3_h",
  expression("Volume searched (m"^3*"/h)")
)

SP_volume_all <- plot_all_fish(
  "SP", "volume_m3_h",
  expression("Volume searched (m"^3*"/h)")
)

ES_energy_all <- plot_all_fish(
  "ES", "E_kJ_h",
  "Energetic requirements (kJ/h)"
)

SP_energy_all <- plot_all_fish(
  "SP", "E_kJ_h",
  "Energetic requirements (kJ/h)"
)

print(ES_volume_all)
print(SP_volume_all)
print(ES_energy_all)
print(SP_energy_all)

# Selected fish: Average vs Extreme in each lake

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

# Environmental availability, weighted by available habitat volume.
# Volume searched comes from env_perf (size-independent); energetic
# requirement comes from env_pop (all valid fish in the lake).
add_panels <- function(df) {
  df %>%
    filter(lake %in% c("ES", "SP"),
           is.finite(vol_weight),
           vol_weight > 0) %>%
    mutate(lake = as.character(lake)) %>%
    inner_join(
      fish_groups_plot %>%
        distinct(lake, panel, panel_label),
      by = "lake",
      relationship = "many-to-many"
    )
}

available_plot <- list(
  volume_m3_h = env_perf %>%
    filter(is.finite(volume_m3_h)) %>%
    select(lake, volume_m3_h, vol_weight) %>%
    add_panels(),
  E_kJ_h = env_pop %>%
    filter(is.finite(E_kJ_h)) %>%
    select(lake, E_kJ_h, vol_weight) %>%
    add_panels()
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
      data = available_plot[[metric]],
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

####Quantile summary of the three distributions####

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
obs_quantiles <- obs_density %>%
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
individual_quantiles <- selected_fish %>%
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

