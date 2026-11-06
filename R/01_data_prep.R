# ============================================================
# 01_data_prep.R
# C2 - population/individual walleye habitat use & energetics
# ES & SP, no McD
#
# Run this ONCE (or whenever raw data change). It reads the raw
# CSVs, builds every shared intermediate dataset, saves each one
# to derived/*.rds, and keeps nothing else. 02_habitat_use.R and
# 03_consequences.R only read the .rds files they need, so neither
# has to hold the raw archival + minutely light data in memory.
#
# Run from the project root (the folder with C2_wae_bioe.Rproj).
# ============================================================

source("R/00_functions.R")
dir.create(derived_dir, showWarnings = FALSE)

# ============================================================
# 1. Fish archival data and fish traits
# ============================================================

wae_summer22_final <- read_csv("wae_summer22_final.csv") %>%
  mutate(
    # Midnight timestamps are stored as a bare date ("2022-05-24"),
    # which plain ymd_hms() turns into NA; truncated = 3 keeps them
    datetime = ymd_hms(datetime, truncated = 3),
    month = lubridate::month(datetime, label = TRUE, abbr = TRUE),
    length_bin = case_when(
      tl_pre >= 310 & tl_pre <= 380 ~ "small",
      tl_pre >= 381 & tl_pre <= 456 ~ "medium",
      tl_pre >= 457                 ~ "large",
      TRUE                          ~ NA_character_
    ),
    length_bin = factor(length_bin,
                        levels = c("small", "medium", "large"),
                        ordered = TRUE)
  ) %>%
  # Only fish present in each of the 5 summer months
  group_by(fish_id) %>%
  filter(all(c("May", "Jun", "Jul", "Aug", "Sep") %in% month)) %>%
  ungroup()

# One row per fish (first observation)
archival_length <- wae_summer22_final %>%
  group_by(fish_id) %>%
  slice(1) %>%
  ungroup()

fish_traits <- archival_length %>%
  select(fish_id, lake, sex, tl_pre) %>%
  distinct() %>%
  mutate(sex = factor(sex))

# ------------------------------------------------------------
# Hourly respiration for each fish (FB4)
# ------------------------------------------------------------

respiration <- wae_summer22_final %>%
  mutate(
    hour_time = floor_date(datetime, "hour")
  ) %>%
  group_by(fish_id, hour_time) %>%
  summarise(
    mean_temp = mean(temp_c, na.rm = TRUE),
    tl_mm = first(tl_pre),
    .groups = "drop"
  ) %>%
  mutate(
    mass_g = length_mm_to_mass_g(tl_mm),
    R_gO2_g_h = respiration_FB4(mass_g, mean_temp, RA, RB, RQ, RTO, RTM, ACT)
  )

save_derived(fish_traits, "fish_traits")
save_derived(respiration, "respiration")
rm(wae_summer22_final, archival_length); gc()

# ------------------------------------------------------------
# Archival data already joined to light, swim speed and
# volume searched (see the commented "Volume searched" code in
# c2_code_clean.R for how these columns were made)
# ------------------------------------------------------------

archival_essp <- bind_rows(
  read_csv("es_archival_joined.csv"),
  read_csv("sp_archival_joined.csv")
) %>%
  mutate(month = month(date))

# Fish present in all of May-Sept (drops the 4 McD fish)
fish_full_season <- archival_essp %>%
  filter(month %in% 5:9) %>%
  group_by(fish_id) %>%
  summarise(n_months = n_distinct(month), .groups = "drop") %>%
  filter(n_months == 5)

archival_essp_full <- archival_essp %>%
  filter(fish_id %in% fish_full_season$fish_id,
         month %in% 5:9)
rm(archival_essp); gc()

# Fish used in the OVL and environmental-availability analyses:
# full-season fish with a length measurement.
# NOTE: `valid_fish` was used in c2_code_clean.R but never created
# there; this is a reconstruction - check that it matches what you
# had in the old workspace.
valid_fish <- fish_traits %>%
  filter(fish_id %in% fish_full_season$fish_id,
         !is.na(tl_pre)) %>%
  select(fish_id, lake, tl_pre)

# Hourly fish summaries used by the consequences analysis
archival_essp_hourly <- archival_essp_full %>%
  group_by(fish_id, lake, year, date, hour_time) %>%
  summarise(
    mean_temp = mean(temp_c, na.rm = TRUE),
    mean_depth = mean(depth_corr, na.rm = TRUE),
    mean_light = mean(I_pred_lux, na.rm = TRUE),
    mean_speed = mean(swim_speed_m_s, na.rm = TRUE),
    mean_volume = mean(volume_searched_m3, na.rm = TRUE),
    total_distance = sum(distance_m, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

archival_hourly_search <- archival_essp_full %>%
  group_by(lake, fish_id, hour_time) %>%
  summarise(mean_search = sum(volume_searched_m3, na.rm = TRUE),
            .groups = "drop")

save_derived(archival_essp_full, "archival_essp_full")
save_derived(valid_fish, "valid_fish")
save_derived(archival_essp_hourly, "archival_essp_hourly")
save_derived(archival_hourly_search, "archival_hourly_search")

# ============================================================
# 2. Hypsography
# ============================================================

essp_hypso <- bind_rows(
  read_csv("sparkling_hypso.csv") %>%
    mutate(lake = "SP",
           depth_m = depth_ft * 0.3048,
           area_m2 = area_ft2 * 0.092903),
  read_csv("escanaba_hypso.csv") %>%
    mutate(lake = "ES")
) %>%
  select(lake, depth_m, proportion_area, area_m2)

# Depth availability as a distribution (for density plots)
lake_depth_dist <- essp_hypso %>%
  uncount(weights = round(proportion_area * 1000))

# Area interpolated to 1 m, then 1 m layer volumes
hypso_layers <- essp_hypso %>%
  group_by(lake) %>%
  group_modify(~ {
    df <- arrange(.x, depth_m)
    grid_depths <- seq(0, floor(max(df$depth_m, na.rm = TRUE)), by = 1)
    out <- approx(df$depth_m, df$area_m2, xout = grid_depths,
                  method = "linear", rule = 2)
    tibble(depth_m = out$x, area_m2 = out$y)
  }) %>%
  ungroup() %>%
  arrange(lake, depth_m) %>%
  group_by(lake) %>%
  mutate(
    depth_top = depth_m,
    depth_bottom = lead(depth_m),
    area_top = area_m2,
    area_bottom = lead(area_m2),
    layer_thickness_m = depth_bottom - depth_top,
    layer_volume_m3 = ((area_top + area_bottom) / 2) * layer_thickness_m,
    depth_mid = depth_top + layer_thickness_m / 2
  ) %>%
  ungroup() %>%
  filter(!is.na(layer_volume_m3))

# 1 m depth grid per lake (for profile interpolation)
depth_grid <- essp_hypso %>%
  group_by(lake) %>%
  summarise(max_depth_m = max(depth_m, na.rm = TRUE), .groups = "drop") %>%
  rowwise() %>%
  mutate(depth_seq = list(seq(0, floor(max_depth_m), by = 1))) %>%
  unnest(depth_seq) %>%
  rename(depth_m = depth_seq) %>%
  ungroup()

save_derived(essp_hypso, "essp_hypso")
save_derived(lake_depth_dist, "lake_depth_dist")
save_derived(hypso_layers, "hypso_layers")

# ============================================================
# 3. Temperature availability (HOBO + LTER), 2022 May-Sep
# ============================================================

hobo_22 <- read_csv("hobo_all.csv") %>%
  mutate(lake = case_when(lake == "Escanaba" ~ "ES",
                          lake == "Sparkling" ~ "SP",
                          TRUE ~ lake)) %>%
  filter(lake %in% c("ES", "SP"),
         year(datetime) == 2022,
         month(datetime) %in% 5:9) %>%
  transmute(lake, string_id = hobo_df_name, datetime, depth_m, temp_c)

lter_22 <- read_csv("ntl5_2022_v1.csv") %>%
  filter(month(sampledate) %in% 5:9) %>%
  mutate(
    lake = "SP",
    string_id = "lter_sp_string",
    datetime = as.POSIXct(sampledate) + as.numeric(sampletime),
    depth_m = depth,
    temp_c = wtemp
  ) %>%
  filter(is.na(flag_wtemp) | flag_wtemp == FALSE) %>%
  transmute(lake, string_id, datetime, depth_m, temp_c)

# Hourly mean per string/depth, then across strings
temp_obs <- bind_rows(hobo_22, lter_22) %>%
  filter(!is.na(temp_c), !is.na(depth_m), !is.na(datetime),
         temp_c <= 26) %>%
  mutate(datetime_hour = floor_date(datetime, unit = "hour")) %>%
  group_by(lake, string_id, datetime_hour, depth_m) %>%
  summarise(temp_c = mean(temp_c, na.rm = TRUE), .groups = "drop") %>%
  group_by(lake, datetime_hour, depth_m) %>%
  summarise(temp_c = mean(temp_c, na.rm = TRUE),
            n_strings = n_distinct(string_id),
            .groups = "drop") %>%
  arrange(lake, datetime_hour, depth_m)
rm(hobo_22, lter_22); gc()

# Hourly profiles interpolated to 1 m (~3 min to run)
temp_vertical <- temp_obs %>%
  group_by(lake, datetime_hour) %>%
  group_modify(~ {
    grid_i <- depth_grid %>% filter(lake == .y$lake[[1]]) %>% pull(depth_m)
    interp_vertical_profile(.x, grid_i)
  }) %>%
  ungroup()

# ... and to 0.1 m (used for the environmental performance grid).
# depth_m is rounded so it joins exactly to the 0.1 m light data
# (seq(..., by = 0.1) gives values like 0.30000000000000004).
temp_interp <- temp_vertical %>%
  arrange(lake, datetime_hour, depth_m) %>%
  group_by(lake, datetime_hour) %>%
  group_modify(~ {
    prof <- .x %>% distinct(depth_m, .keep_all = TRUE) %>% arrange(depth_m)
    new_depths <- seq(min(prof$depth_m), max(prof$depth_m), by = 0.1)
    out <- approx(prof$depth_m, prof$temp_c, xout = new_depths,
                  method = "linear", rule = 1)
    tibble(depth_m = round(out$x, 1), temp_c = out$y)
  }) %>%
  ungroup()

save_derived(temp_vertical, "temp_vertical")
save_derived(temp_interp, "temp_interp")
rm(temp_obs); gc()

# ============================================================
# 4. Light availability
# Minutely interpolation of the 0.1 m light grids (from
# archive/light_grids_one_time.R), then collapsed to hourly.
# The minutely data are the biggest objects in the project, so
# only one day is kept (for the trajectory plots).
# ============================================================

interp_light_minutely <- function(light) {
  light <- light %>%
    mutate(date_time = ymd_hms(date_time, truncated = 3)) %>%
    filter(!is.na(date_time))

  num_cols <- light %>% select(where(is.numeric), -depth_m) %>% names()

  minute_seq <- seq(min(light$date_time), max(light$date_time), by = "1 min")

  light %>%
    group_by(depth_m) %>%
    arrange(date_time) %>%
    complete(date_time = minute_seq) %>%
    mutate(across(all_of(num_cols),
                  ~ approx(date_time, ., date_time, rule = 2)$y)) %>%
    ungroup()
}

traj_day <- as.Date("2022-06-29")

light_hourly_list <- list()
light_check_list <- list()
light_day_list <- list()

for (lk in c("ES", "SP")) {
  file_lk <- if (lk == "ES") "escanaba_light_sum22.csv" else "sparkling_light_sum22.csv"

  light_minute <- interp_light_minutely(read_csv(file_lk)) %>%
    mutate(lake = lk)

  # One day of minutely light for trajectory plots
  light_day_list[[lk]] <- light_minute %>%
    filter(date_time >= traj_day, date_time < traj_day + 1)

  # Hourly, 1 m bins, positive light only (availability plots)
  light_hourly_list[[lk]] <- light_minute %>%
    filter(!is.na(I_pred_lux), I_pred_lux > 0) %>%
    mutate(datetime_hour = floor_date(date_time, "hour"),
           depth_m = floor(depth_m)) %>%
    group_by(lake, datetime_hour, depth_m) %>%
    summarise(light_lux = mean(I_pred_lux, na.rm = TRUE), .groups = "drop")

  # Hourly, 0.1 m (environmental performance grid)
  light_check_list[[lk]] <- light_minute %>%
    mutate(datetime_hour = floor_date(date_time, unit = "hour"),
           depth_m = round(depth_m, 1)) %>%
    group_by(lake, datetime_hour, depth_m) %>%
    summarise(
      n_obs = n(),
      n_nonmiss = sum(!is.na(I_pred_lux)),
      light_lux = mean(I_pred_lux, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(light_lux = ifelse(n_nonmiss == 0, NA_real_, light_lux))

  rm(light_minute); gc()
}

light_hourly <- bind_rows(light_hourly_list)
light_hourly_check <- bind_rows(light_check_list)
light_minute_day <- bind_rows(light_day_list)
rm(light_hourly_list, light_check_list, light_day_list); gc()

save_derived(light_hourly, "light_hourly")
save_derived(light_minute_day, "light_minute_day")

# ============================================================
# 5. Hourly temperature x light environment (0.1 m), restricted
#    to depths covered by both temperature and light each hour
# ============================================================

depth_limits <- temp_interp %>%
  group_by(lake, datetime_hour) %>%
  summarise(min_temp_depth = min(depth_m, na.rm = TRUE),
            max_temp_depth = max(depth_m, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(
    light_hourly_check %>%
      group_by(lake, datetime_hour) %>%
      summarise(min_light_depth = min(depth_m, na.rm = TRUE),
                max_light_depth = max(depth_m, na.rm = TRUE),
                .groups = "drop"),
    by = c("lake", "datetime_hour")
  ) %>%
  mutate(min_shared_depth = pmax(min_temp_depth, min_light_depth),
         max_shared_depth = pmin(max_temp_depth, max_light_depth))

env_hourly <- temp_interp %>%
  left_join(light_hourly_check %>% select(lake, datetime_hour, depth_m, light_lux),
            by = c("lake", "datetime_hour", "depth_m")) %>%
  left_join(depth_limits %>% select(lake, datetime_hour,
                                    min_shared_depth, max_shared_depth),
            by = c("lake", "datetime_hour")) %>%
  filter(depth_m >= min_shared_depth,
         depth_m <= max_shared_depth,
         !is.na(light_lux)) %>%
  select(lake, datetime_hour, depth_m, temp_c, light_lux)

save_derived(env_hourly, "env_hourly")
rm(light_hourly_check, depth_limits); gc()

# ============================================================
# 6. Dissolved oxygen
#    Profiles -> 1 m grid -> interpolated to every fish timestamp
# ============================================================

oxygen_clean <- read_csv("oxygen_df.csv") %>%
  mutate(datetime = mdy_hm(date_time),
         lake = case_when(lake == "Escanaba" ~ "ES",
                          lake == "Sparkling" ~ "SP",
                          TRUE ~ lake)) %>%
  filter(!is.na(datetime), !is.na(depth_m), !is.na(do_mgl))

oxygen_vertical <- oxygen_clean %>%
  group_by(lake, datetime) %>%
  group_modify(~ {
    grid_i <- depth_grid %>% filter(lake == .y$lake[[1]]) %>% pull(depth_m)
    interp_do_profile(.x, grid_i)
  }) %>%
  ungroup()

fish_times <- archival_essp_full %>%
  filter(!is.na(date_time)) %>%
  distinct(lake, datetime = date_time)

oxygen_interp <- oxygen_vertical %>%
  group_by(lake, depth_m) %>%
  group_modify(~ {
    target_times <- fish_times %>% filter(lake == .y$lake[[1]]) %>% pull(datetime)
    interp_time_series(.x, target_times)
  }) %>%
  ungroup()

# DO experienced by each fish (only the columns the analyses use)
fish_do <- archival_essp_full %>%
  filter(!is.na(depth_m), !is.na(date_time)) %>%
  select(fish_id, lake, date_time, depth_m) %>%
  left_join(oxygen_interp, by = c("lake", "depth_m", "date_time" = "datetime")) %>%
  filter(!is.na(do_mgl))

save_derived(oxygen_interp, "oxygen_interp")
save_derived(fish_do, "fish_do")

# ============================================================
# Done - everything needed downstream is in derived/.
# Restart R (Session > Restart R) before running 02 or 03.
# ============================================================
list.files(derived_dir)
