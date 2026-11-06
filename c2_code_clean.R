#C2 - population/individual wae energetics (very theroetical)
#ES & SP, no McD

####Load libraries####
library(tidyverse)
library(readr)
library(zoo)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(lubridate)
library(tibble)
library(gt)

####Load & Manipulate Data####
#Load in archival data for summer 2022 
wae_summer22_final <- read_csv("wae_summer22_final.csv")

sp_archival <- wae_summer22_final %>%
  filter(lake == "SP")

sp_archival <- sp_archival %>%
  rename(date_time = datetime)

es_archival <- wae_summer22_final %>%
  filter(lake == "ES")

es_archival <- es_archival %>%
  rename(date_time = datetime)

oxygen_df <- read_csv("oxygen_df.csv")

#Hobo bogus 
hobo_all <- read_csv("hobo_all.csv") #This is what should be used methinks 

es_surface <- read_csv("Escanaba_surface.csv")
sp_surface <- read_csv("Sparkling_surface.csv")

#Light stuff that def won't break the comp- 

light_es <- read_csv("escanaba_light_sum22.csv")
light_sp <- read_csv("sparkling_light_sum22.csv")

#Load hypsos

sparkling_hypso <- read_csv("sparkling_hypso.csv") %>%
  mutate(
    lake = "SP",
    depth_m = depth_ft * 0.3048,
    area_m2 = area_ft2 * 0.092903
  ) %>%
  select(lake, depth_m, proportion_area, area_m2)

escanaba_hypso <- read_csv("escanaba_hypso.csv") %>%
  mutate(
    lake = "ES",
    depth_m = depth_m
  ) %>%
  select(lake, depth_m, proportion_area, area_m2)

essp_hypso <- bind_rows(
  sparkling_hypso,
  escanaba_hypso
)

lake_depth_dist <- essp_hypso %>%
  uncount(weights = round(proportion_area * 1000))

#Big data joined to light 

es_archival_joined <- read_csv("es_archival_joined.csv")
sp_archival_joined <- read_csv("sp_archival_joined.csv")

#Pareto frontier from before 

obs_fish_plot <- read_csv("obs_fish_plot.csv") #This is what should be used methinks 
obs_hourly_plot <- read_csv("obs_hourly_plot.csv") #This is what should be used methinks 

#Manipulate some things 
# Escanaba daily surface averages
es_surface_daily <- es_surface %>%
  mutate(doy = yday(date_time)) %>%
  group_by(doy) %>%
  summarise(surface_temp = mean(temp_c, na.rm = TRUE), .groups = "drop")

#Sparkling daily surface averages
sp_surface_daily <- sp_surface %>%
  mutate(doy = yday(date_time)) %>%
  group_by(doy) %>%
  summarise(surface_temp = mean(temp_c, na.rm = TRUE), .groups = "drop")

#Add some columns to wae_archival data 
wae_summer22_final <- wae_summer22_final %>%
  mutate(
    month = lubridate::month(datetime, label = TRUE, abbr = TRUE),
    
    length_bin = case_when(
      tl_pre >= 310 & tl_pre <= 380 ~ "small",
      tl_pre >= 381 & tl_pre <= 456 ~ "medium",
      tl_pre >= 457                ~ "large",
      TRUE                         ~ NA_character_
    )
  )

#Order those columns 

wae_summer22_final <- wae_summer22_final %>%
  mutate(
    length_bin = factor(
      length_bin,
      levels = c("small", "medium", "large"),
      ordered = TRUE
    )
  )

#Further restrict this dataset to fish that appear in all summer months

wae_summer22_final <- wae_summer22_final %>%
  group_by(fish_id) %>%
  filter(all(c("May", "Jun", "Jul", "Aug", "Sep") %in% month)) %>%
  ungroup()

#Make it so only fish present in each of the 5 months are included 

# Keep one tl_mm value per unique fish ID (e.g., first observation per fish)
archival_length <- wae_summer22_final %>%
  group_by(fish_id) %>%
  slice(1) %>%     # or summarise(mean_tl = mean(tl_mm), .groups = "drop") if you want an average
  ungroup()

archival_daily <- wae_summer22_final %>%
  group_by(lake, fish_id, doy) %>%
  summarise(daily_temp = mean(temp_c, na.rm = TRUE), .groups = "drop")

####Bioenergetics####

# FB4 respiration model for Walleye (Sander vitreus)
# Based on Kitchell et al. 1977 (Equation 2)


#prep dataframe 

# Ensure date_time is POSIXct and summarize hourly
bioe_archival <- wae_summer22_final %>%
  mutate(
    datetime = ymd_hms(datetime),          # ensure datetime format
    hour_time = floor_date(datetime, "hour")  # round down to the nearest hour
  ) %>%
  group_by(fish_id, hour_time) %>%
  summarise(
    mean_temp = mean(temp_c, na.rm = TRUE),
    tl_mm = first(tl_pre)  # assuming total length doesn't change
  ) %>%
  ungroup()

# 1. Input data 

# Required columns:
#   id        = unique fish identifier
#   datetime  = timestamp (POSIXct)
#   temp_C    = temperature in °C
#   length_mm = fish total length (mm)
## I'm an idiot and changed column names below so be warned 




# 2. FB4 parameters for Walleye (adult) 

#Int of allometric function for smr
RA  <- 0.0108
#Slope of allometric function for smr
RB  <- -0.2
#Q10 rate rate at which fun increases 
RQ  <- 2.1
#Optimum temp for resp. 
RTO <- 27
#Lethal water temp
RTM <- 32
#Activity constant - can be edited later and edited below  
ACT <- 1


# 3. Length–weight conversion

# Using FishBase Bayesian parameters for wae:
# W = a * L^b  (L in cm, W in g)
#Fix this later with data from GS, MW, JTM?
lw_a <- 0.00631
lw_b <- 3.15

length_mm_to_mass_g <- function(tl_mm, a = lw_a, b = lw_b) {
  L_cm <- tl_mm / 10
  W_g <- a * (L_cm ^ b)
  return(W_g)
}


# 4. Temperature function (TF) the fuck

TF_calc <- function(T_C, RTO, RTM) {
  TF <- ((T_C - RTO) / (RTM - RTO)) / (1 - ((T_C - RTO) / (RTM - RTO)))
  return(TF)
}


# 5. Respiration equation 

respiration_FB4 <- function(W_g, T_C, RA, RB, RQ, RTO, RTM, ACT) {
  TF <- TF_calc(T_C, RTO, RTM)
  R <- RA * (W_g ^ RB) * exp(RQ * TF) * ACT
  return(R)
}


# 6. Apply model to dataset

respiration <- bioe_archival%>%
  mutate(
    mass_g = length_mm_to_mass_g(tl_mm),
    R_gO2_g_h = respiration_FB4(mass_g, mean_temp, RA, RB, RQ, RTO, RTM, ACT)
  )

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






####Let there be light####

#The big equation 
process_time_step_midpoints <- function(df_time,
                                        delta_z = 0.1,
                                        small_light = 0.1,
                                        min_sensors = 3) {
  
  df_time <- df_time %>%
    arrange(depth_m) %>%
    mutate(light_lux = ifelse(is.na(light_lux) | light_lux <= 0,
                              small_light, light_lux))
  
  # Require enough sensors
  if (sum(df_time$light_lux > 0) < min_sensors) {
    return(list(midpoint_k = tibble(), light_grid = tibble()))
  }
  
  d <- df_time$depth_m
  I <- df_time$light_lux
  
  # ---------- 1. Bulk k ----------
  lm_fit <- lm(log(I) ~ d)
  k_bulk <- -coef(lm_fit)[2]
  I0_reg <- exp(coef(lm_fit)[1])
  
  # ---------- 2. Local midpoint k-values (NO interpolation) ----------
  k_local <- -log(I[-1] / I[-length(I)]) / diff(d)
  mid_depth <- (d[-1] + d[-length(d)]) / 2
  
  midpoint_df <- tibble(depth_mid = mid_depth,
                        k_mid = k_local)
  
  # Fix zeros: if k == 0 replace by previous k
  for (i in seq_len(nrow(midpoint_df))) {
    if (midpoint_df$k_mid[i] <= 0 &&
        i > 1) {
      midpoint_df$k_mid[i] <- midpoint_df$k_mid[i - 1]
    }
  }
  
  # ---------- 3. Interpolate k(z) only for light prediction grid ----------
  max_z <- max(d)
  grid_z <- seq(0, max_z, by = delta_z)
  
  if (nrow(midpoint_df) == 1) {
    k_interp <- rep(midpoint_df$k_mid[1], length(grid_z))
  } else {
    k_interp <- approx(x = midpoint_df$depth_mid,
                       y = midpoint_df$k_mid,
                       xout = grid_z,
                       rule = 2)$y
  }
  
  # Prevent negative / zero k
  k_interp[k_interp <= 0] <- 0.0001
  
  # ---------- 4. Integrate k(z) and compute I(z) ----------
  dz <- delta_z
  tau <- numeric(length(k_interp))
  for (i in 2:length(k_interp)) {
    tau[i] <- tau[i - 1] + 0.5*(k_interp[i - 1] + k_interp[i])*dz
  }
  
  I0 <- I0_reg
  if (is.na(I0) || I0 <= 0) {
    I0 <- I[1] * exp(midpoint_df$k_mid[1] * d[1])
  }
  
  I_pred <- I0 * exp(-tau)
  
  light_grid <- tibble(
    depth_m = grid_z,
    k_interp = k_interp,
    tau = tau,
    I_pred_lux = I_pred
  )
  
  return(list(
    midpoint_k = midpoint_df %>% mutate(k_bulk = k_bulk,
                                        I0_reg = I0_reg),
    light_grid = light_grid %>% mutate(k_bulk = k_bulk,
                                       I0_reg = I0_reg)
  ))
}
####Light things to not run ####
#Sparkling - done (don't run this again)
#May

hobo_sparkling_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Sparkling",
    datetime >= as.Date("2022-05-01"),
    datetime <  as.Date("2022-06-01")              # strictly before Oct 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_sparkling_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_sp_may<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))


#June 
hobo_sparkling_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Sparkling",
    datetime >= as.Date("2022-06-01"),
    datetime <  as.Date("2022-07-01")              # strictly before july 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_sparkling_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_sp_june<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))
#July

hobo_sparkling_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Sparkling",
    datetime >= as.Date("2022-07-01"),
    datetime <  as.Date("2022-08-01")              # strictly before aug 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_sparkling_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_sp_july<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))

#August
hobo_sparkling_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Sparkling",
    datetime >= as.Date("2022-08-01"),
    datetime <  as.Date("2022-09-01")              # strictly before sept1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_sparkling_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_sp_aug<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))

#september

hobo_sparkling_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Sparkling",
    datetime >= as.Date("2022-09-01"),
    datetime <  as.Date("2022-10-01")              # strictly before Oct 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_sparkling_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_sp_sep<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))



#Combine all the bogus 

# Combine all monthly datasets
results_sp_all <- bind_rows(
  results_sp_may,
  results_sp_june,
  results_sp_july,
  results_sp_aug,
  results_sp_sep
) %>%
  # Make sure date_time is POSIXct
  mutate(date_time = ymd_hms(date_time))

#Extract k values 
midpoints_sp <- results_sp_all %>%
  mutate(midpoint_k = map(proc, "midpoint_k")) %>%
  select(date_time, midpoint_k) %>%
  unnest(midpoint_k)

write.csv(midpoints_sp, file = "sparkling_k_mids.csv", row.names=FALSE)

#Extract light grids 
light_sp <- results_sp_all %>%
  mutate(light_grid = map(proc, "light_grid")) %>%
  select(date_time, light_grid) %>%
  unnest(light_grid)

write.csv(light_sp, file = "sparkling_light_sum22.csv", row.names=FALSE)

#For safekeeping and for later 
#light_sp is my light file 

#ES light grids 

#May - mcd

hobo_es_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Escanaba",
    datetime >= as.Date("2022-05-01"),
    datetime <  as.Date("2022-06-01")              # strictly before Oct 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_es_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_es_may<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))


#June 
hobo_es_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Escanaba",
    datetime >= as.Date("2022-06-01"),
    datetime <  as.Date("2022-07-01")              # strictly before july 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_es_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_es_june<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))
#July

hobo_es_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Escanaba",
    datetime >= as.Date("2022-07-01"),
    datetime <  as.Date("2022-08-01")              # strictly before aug 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_es_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_es_july<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))

#August
hobo_es_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Escanaba",
    datetime >= as.Date("2022-08-01"),
    datetime <  as.Date("2022-09-01")              # strictly before sept1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_es_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_es_aug<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))

#september

hobo_es_subset <- all_hobo %>%
  mutate(datetime = ymd_hms(datetime)) %>%         # ensure datetime is POSIXct
  filter(
    lake == "Escanaba",
    datetime >= as.Date("2022-09-01"),
    datetime <  as.Date("2022-10-01")              # strictly before Oct 1
  ) %>%
  select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m)

future::plan(multisession, workers = parallel::detectCores() - 1)

nested_hobo <- hobo_es_subset %>%
  mutate(date_time = datetime) %>%
  arrange(date_time, depth_m) %>%
  group_by(date_time) %>%
  nest()

results_es_sep<- nested_hobo %>%
  mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                           .progress = TRUE))



#Combine all the bogus 

# Combine all monthly datasets
results_es_all <- bind_rows(
  results_es_may,
  results_es_june,
  results_es_july,
  results_es_aug,
  results_es_sep
) %>%
  # Make sure date_time is POSIXct
  mutate(date_time = ymd_hms(date_time))



#Extract k values 
midpoints_es <- results_es_all %>%
  mutate(midpoint_k = map(proc, "midpoint_k")) %>%
  select(date_time, midpoint_k) %>%
  unnest(midpoint_k)

write.csv(midpoints_es, file = "es_k_mids.csv", row.names=FALSE)

#Extract light grids 
light_es <- results_es_all %>%
  mutate(light_grid = map(proc, "light_grid")) %>%
  select(date_time, light_grid) %>%
  unnest(light_grid)

write.csv(light_es, file = "escanaba_light_sum22.csv", row.names=FALSE)

#Done with light grids - joining to archival data by lake below. 

####Joining light to archival data ####

#Data loaded above - skip the calculations besides minutely light interpolations! 

# Ensure datetime is POSIXct
light_sp <- light_sp %>%
  mutate(date_time = ymd_hms(date_time))

summary(light_sp$date_time)
any(is.na(light_sp$date_time))

light_sp <- light_sp %>% 
  filter(!is.na(date_time))

# Identify numeric columns to interpolate (excluding depth)
num_cols <- light_sp %>%
  select(where(is.numeric), -depth_m) %>%
  names()

# Create full minutely time sequence
minute_seq <- tibble(date_time = seq(min(light_sp$date_time),
                                     max(light_sp$date_time),
                                     by = "1 min"))

# Interpolate by depth (very important!)
light_sp_minute <- light_sp %>%
  group_by(depth_m) %>%              # interpolate each depth independently
  arrange(date_time) %>%
  # Complete to minutely grid
  complete(date_time = minute_seq$date_time) %>%
  # Linear interpolation for numeric columns
  mutate(across(all_of(num_cols), ~ approx(date_time, ., date_time, rule = 2)$y)) %>%
  ungroup()



#################ES

# Ensure datetime is POSIXct
Light_es <- light_es %>%
  mutate(date_time = ymd_hms(date_time))

summary(light_es$date_time)
any(is.na(light_es$date_time))

light_es <- light_es %>% 
  filter(!is.na(date_time))

# Identify numeric columns to interpolate (excluding depth)
num_cols <- light_es %>%
  select(where(is.numeric), -depth_m) %>%
  names()

# Create full minutely time sequence
minute_seq <- tibble(date_time = seq(min(light_es$date_time),
                                     max(light_es$date_time),
                                     by = "1 min"))

# Interpolate by depth (very important!)
light_es_minute <- light_es %>%
  group_by(depth_m) %>%              # interpolate each depth independently
  arrange(date_time) %>%
  # Complete to minutely grid
  complete(date_time = minute_seq$date_time) %>%
  # Linear interpolation for numeric columns
  mutate(across(all_of(num_cols), ~ approx(date_time, ., date_time, rule = 2)$y)) %>%
  ungroup()

#All this is important for light density plots! 
#There's additional stuff for filtering and day/night - might need for new C2

####Volume searched - also to not run####
#Use resp instead? 
#Yea idiot start here 



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

#SP volume searched below - STOP HERE

#sp_archival_joined <- sp_archival_joined %>%
#  mutate(
#    hour_time = floor_date(date_time, "hour")
#  ) %>%
#  left_join(
#    df_energy_swim %>%
#      select(fish_id, hour_time, swim_speed_m_s),
#    by = c("fish_id", "hour_time")
#  )

#sp_archival_joined <- sp_archival_joined %>%
#  arrange(fish_id, date_time) %>%
#  group_by(fish_id) %>%
#  mutate(
#    dt_sec = as.numeric(difftime(date_time,
#                                 lag(date_time),
#                                 units = "secs")),
#    
#    distance_m = swim_speed_m_s * dt_sec
#  ) %>%
#  ungroup()

#sp_archival_joined <- sp_archival_joined %>%
#  mutate(
#    reaction_dist_m = case_when(
#      I_pred_lux <= 0.1 ~ 0.10,
#      I_pred_lux <= 1   ~ 0.20,
#      I_pred_lux > 1    ~ 0.30,
#      TRUE ~ NA_real_
#    )
#  )

#sp_archival_joined <- sp_archival_joined %>%
#  mutate(
#    volume_searched_m3 = pi * (reaction_dist_m^2) * distance_m
#  )

summary(sp_archival_joined$distance_m)
summary(sp_archival_joined$volume_searched_m3)

#Run this 
daily_volume_sp <- sp_archival_joined %>%
  mutate(date = as.Date(date_time)) %>%
  group_by(fish_id, date) %>%
  summarise(
    daily_volume_m3 = sum(volume_searched_m3, na.rm = TRUE),
    mean_temp = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

summary(daily_volume_sp$daily_volume_m3)

#ES volume searched 

#str(es_archival_joined)

#setdiff(unique(es_archival_joined$fish_id),
#        unique(df_energy_swim$fish_id))

#setdiff(unique(es_archival_joined$fish_id),
#        unique(es_archival$fish_id))

#setdiff(unique(es_archival_joined$fish_id),
#        unique(es_archival$fish_id))

#Some ES fish are surviving the pruning process leading to NA's - get rid of these 

#es_archival_joined <- es_archival_joined %>%
#  mutate(
#    hour_time = floor_date(date_time, "hour")
#  ) %>%
#  left_join(
#    df_energy_swim %>%
#      select(fish_id, hour_time, swim_speed_m_s),
#    by = c("fish_id", "hour_time")
#  )

#es_archival_joined <- es_archival_joined %>%
#  arrange(fish_id, date_time) %>%
#  group_by(fish_id) %>%
#  mutate(
#    dt_sec = as.numeric(difftime(date_time,
#                                 lag(date_time),
#                                 units = "secs")),
#    
#    distance_m = swim_speed_m_s * dt_sec
#  ) %>%
#  ungroup()

#es_archival_joined <- es_archival_joined %>%
#  mutate(
#    reaction_dist_m = case_when(
#      I_pred_lux <= 0.1 ~ 0.10,
#      I_pred_lux <= 1   ~ 0.20,
#      I_pred_lux > 1    ~ 0.30,
#      TRUE ~ NA_real_
#    )
#  )

#es_archival_joined <- es_archival_joined %>%
#  mutate(
#    volume_searched_m3 = pi * (reaction_dist_m^2) * distance_m
#  )

summary(es_archival_joined$distance_m)
summary(es_archival_joined$volume_searched_m3)

#And run this! 
daily_volume_es <- es_archival_joined %>%
  mutate(date = as.Date(date_time)) %>%
  group_by(fish_id, date) %>%
  summarise(
    daily_volume_m3 = sum(volume_searched_m3, na.rm = TRUE),
    mean_temp = mean(temp_c, na.rm = TRUE),
    .groups = "drop"
  )

summary(daily_volume_es$daily_volume_m3)

####Combined dataset for foraging efficiency####
#Make the move to hourly data to compare volume searched and energetic costs: 

archival_essp <- bind_rows(es_archival_joined,
                           sp_archival_joined)

str(archival_essp)

# Add month column
archival_essp <- archival_essp %>%
  mutate(month = month(date))

# Identify fish present in May–Sept (5–9)
fish_full_season <- archival_essp %>%
  filter(month %in% 5:9) %>%
  group_by(fish_id) %>%
  summarise(n_months = n_distinct(month),
            .groups = "drop") %>%
  filter(n_months == 5)

#4 fish missing - McD fish
archival_essp_full <- archival_essp %>%
  filter(fish_id %in% fish_full_season$fish_id,
         month %in% 5:9)

archival_essp_full %>%
  group_by(fish_id) %>%
  summarise(months = n_distinct(month))

#Vol searched

archival_hourly_search <- archival_essp_full %>%   # use your filtered full-season dataset
  group_by(lake, fish_id, hour_time) %>%
  summarise(
    mean_search = sum(volume_searched_m3, na.rm = TRUE),
    .groups = "drop"
  )

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


####Density plots of habitat use####

fish_depth_dist <- archival_hourly_depth %>%
  rename(depth_m = mean_depth)

fish_depth_pop_ES <- fish_depth_dist %>%
  dplyr::filter(lake == "ES", depth_m <= 9) %>%
  dplyr::select(depth_m)

fish_depth_pop_SP <- fish_depth_dist %>%
  dplyr::filter(lake == "SP") %>%
  dplyr::select(depth_m)

#ES depth use vs availability
ggplot() +
  
  # habitat availability
  geom_density(
    data = dplyr::filter(lake_depth_dist, lake == "ES"),
    aes(y = depth_m),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  # population depth use
  geom_density(
    data = fish_depth_pop_ES,
    aes(y = depth_m),
    fill = "steelblue",
    alpha = 0.6,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  
  # individual fish use
  geom_density(
    data = dplyr::filter(fish_depth_dist,
                         lake == "ES",
                         depth_m <= 9),
    aes(y = depth_m, color = fish_id),
    linewidth = 1,
    alpha = 1
  ) +
  
  facet_wrap(~fish_id) +
  scale_y_reverse() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Density",
       y = "Depth (m)",
       title = "Escanaba Lake – Individual Depth Use vs Population and Availability")

#Sparkling depth use vs availability
ggplot() +
  
  # habitat availability
  geom_density(
    data = dplyr::filter(lake_depth_dist, lake == "SP"),
    aes(y = depth_m),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  # population depth use
  geom_density(
    data = fish_depth_pop_SP,
    aes(y = depth_m),
    fill = "steelblue",
    alpha = 0.6,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  
  # individual fish use
  geom_density(
    data = dplyr::filter(fish_depth_dist,
                         lake == "SP"),
    aes(y = depth_m, color = fish_id),
    linewidth = 1,
    alpha = 1
  ) +
  
  facet_wrap(~fish_id) +
  scale_y_reverse() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Density",
       y = "Depth (m)",
       title = "SP Lake – Individual Depth Use vs Population and Availability")

#Temp available vs. use 
lter_sp_22 <- read_csv("ntl5_2022_v1.csv")

str(lter_sp_22)

# 1. HOBO data: ES + SP, 2023, May-Sep

hobo_22 <- hobo_all %>%
  mutate(
    lake = case_when(
      lake == "Escanaba" ~ "ES",
      lake == "Sparkling" ~ "SP",
      TRUE ~ lake
    )
  ) %>%
  filter(
    lake %in% c("ES", "SP"),
    year(datetime) == 2022,
    month(datetime) %in% 5:9
  ) %>%
  transmute(
    lake,
    string_id = hobo_df_name,
    datetime,
    depth_m = depth_m,
    temp_c = temp_c
  )

# 2. LTER SP data: create datetime and standardize
lter_22 <- lter_sp_22 %>%
  filter(month(sampledate) %in% 5:9) %>%
  mutate(
    lake = "SP",
    string_id = "lter_sp_string",
    datetime = as.POSIXct(sampledate) + as.numeric(sampletime),
    depth_m = depth,
    temp_c = wtemp
  ) %>%
  transmute(
    lake,
    string_id,
    datetime,
    depth_m,
    temp_c,
    flag_wtemp
  ) %>%
  filter(is.na(flag_wtemp) | flag_wtemp == FALSE) %>%
  select(-flag_wtemp)

# 3. Combine

temp_raw_22 <- bind_rows(hobo_22, lter_22) %>%
  filter(!is.na(temp_c), !is.na(depth_m), !is.na(datetime))

#Building temp grid 
temp_hourly <- temp_raw_22 %>%
  filter(!is.na(temp_c), temp_c <= 26) %>%
  mutate(datetime_hour = floor_date(datetime, unit = "hour")) %>%
  group_by(lake, string_id, datetime_hour, depth_m) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

temp_obs <- temp_hourly %>%
  group_by(lake, datetime_hour, depth_m) %>%
  summarise(
    temp_c = mean(temp_c, na.rm = TRUE),
    n_strings = n_distinct(string_id),
    .groups = "drop"
  ) %>%
  arrange(lake, datetime_hour, depth_m)

max_depths <- essp_hypso %>%
  group_by(lake) %>%
  summarise(max_depth_m = max(depth_m, na.rm = TRUE), .groups = "drop")

depth_grid <- max_depths %>%
  rowwise() %>%
  mutate(depth_seq = list(seq(0, floor(max_depth_m), by = 1))) %>%
  unnest(depth_seq) %>%
  rename(depth_m = depth_seq) %>%
  ungroup()

interp_vertical_profile <- function(df_profile, grid_depths) {
  df_profile <- df_profile %>%
    filter(!is.na(depth_m), !is.na(temp_c)) %>%
    arrange(depth_m) %>%
    distinct(depth_m, .keep_all = TRUE)
  
  if (nrow(df_profile) < 2) {
    return(tibble(
      depth_m = numeric(0),
      temp_c = numeric(0)
    ))
  }
  
  out <- approx(
    x = df_profile$depth_m,
    y = df_profile$temp_c,
    xout = grid_depths,
    method = "linear",
    rule = 1
  )
  
  tibble(
    depth_m = out$x,
    temp_c = out$y
  ) %>%
    filter(!is.na(temp_c))
}

temp_obs %>%
  group_by(lake, datetime_hour) %>%
  summarise(n_depths = n_distinct(depth_m), .groups = "drop") %>%
  count(n_depths)

temp_vertical <- temp_obs %>%
  group_by(lake, datetime_hour) %>%
  group_modify(~ {
    lake_i <- .y$lake[[1]]
    
    grid_i <- depth_grid %>%
      filter(lake == lake_i) %>%
      pull(depth_m)
    
    interp_vertical_profile(.x, grid_i)
  }) %>%
  ungroup()

#Check check check - 3 min to run 

glimpse(temp_vertical)

temp_vertical %>%
  count(lake)

temp_vertical %>%
  group_by(lake, datetime_hour) %>%
  summarise(n_depths = n(), .groups = "drop") %>%
  count(n_depths)

temp_vertical %>%
  filter(lake == "ES") %>%
  filter(datetime_hour == max(datetime_hour)) %>%
  ggplot(aes(x = temp_c, y = depth_m)) +
  geom_path() +
  scale_y_reverse()

#Hypso to 1m? 

hypso_1m <- essp_hypso %>%
  group_by(lake) %>%
  group_modify(~{
    df <- .x %>%
      arrange(depth_m)
    
    grid_depths <- seq(0, floor(max(df$depth_m, na.rm = TRUE)), by = 1)
    
    out <- approx(
      x = df$depth_m,
      y = df$area_m2,
      xout = grid_depths,
      method = "linear",
      rule = 2
    )
    
    tibble(
      depth_m = out$x,
      area_m2 = out$y
    )
  }) %>%
  ungroup()

hypso_layers <- hypso_1m %>%
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

#Check 

hypso_layers %>%
  group_by(lake) %>%
  summarise(
    total_volume_m3 = sum(layer_volume_m3, na.rm = TRUE),
    max_depth_layer = max(depth_bottom, na.rm = TRUE),
    .groups = "drop"
  )

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

#First pass at plotting 

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

#Big ol plots 
ggplot() +
  geom_area(
    data = avail_temp_ES,
    aes(x = temp_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_pop_ES,
    aes(x = temp_mid, y = prop),
    color = "steelblue",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_ind_ES,
    aes(x = temp_mid, y = prop, color = fish_id, group = fish_id),
    linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Temperature (°C)",
    y = "Proportional frequency",
    title = "Escanaba Lake – Individual Temperature Use vs Population and Availability"
  )

#Now sparkles 

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


ggplot() +
  geom_area(
    data = avail_temp_SP,
    aes(x = temp_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_pop_SP,
    aes(x = temp_mid, y = prop),
    color = "steelblue",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_ind_SP,
    aes(x = temp_mid, y = prop, color = fish_id, group = fish_id),
    linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Temperature (°C)",
    y = "Proportional frequency",
    title = "Sparkling Lake – Individual Temperature Use vs Population and Availability"
  )

# Light density stuff - 4/23

light_es_minute$lake <- "ES"
light_sp_minute$lake <- "SP"

light_essp <- bind_rows(light_es_minute, light_sp_minute)

str(light_essp)


light_hourly <- light_essp %>%
  filter(!is.na(I_pred_lux), I_pred_lux > 0) %>%   # remove zeros for log scale later
  mutate(
    datetime_hour = floor_date(date_time, "hour"),
    depth_m_round = floor(depth_m)   # collapse to 1 m bins
  ) %>%
  group_by(lake, datetime_hour, depth_m_round) %>%
  summarise(
    light_lux = mean(I_pred_lux, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(depth_m = depth_m_round)

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

avail_light_ES <- light_with_volume %>%
  filter(lake == "ES") %>%
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
    prop = volume / sum(volume)
  ) %>%
  arrange(light_mid)

fish_light_pop_ES <- archival_essp_full %>%
  filter(lake == "ES",
         month %in% 5:9,
         !is.na(I_pred_lux),
         I_pred_lux > 0) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  count(light_bin, name = "n") %>%
  right_join(bin_lookup_light, by = "light_bin") %>%
  mutate(
    n = tidyr::replace_na(n, 0),
    prop = n / sum(n)
  )

fish_light_ind_ES <- archival_essp_full %>%
  filter(lake == "ES",
         month %in% 5:9,
         !is.na(I_pred_lux),
         I_pred_lux > 0) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  group_by(fish_id, light_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::complete(fish_id, light_bin = bin_lookup_light$light_bin, fill = list(n = 0)) %>%
  left_join(bin_lookup_light, by = "light_bin") %>%
  group_by(fish_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot() +
  geom_area(
    data = avail_light_ES,
    aes(x = light_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_light_pop_ES,
    aes(x = light_mid, y = prop),
    color = "steelblue",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_light_ind_ES,
    aes(x = light_mid, y = prop, color = fish_id, group = fish_id),
    linewidth = 1
  ) +
  scale_x_log10(
    breaks = 10^(-1:5),
    labels = scales::label_number()
  ) +
  coord_cartesian(xlim = c(1e-2, 1e5)) +
  facet_wrap(~fish_id) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Light (lux, log scale)",
    y = "Proportional frequency",
    title = "Escanaba Lake – Light Use vs Availability"
  )


#Sparkles light

avail_light_SP <- light_with_volume %>%
  filter(lake == "SP") %>%
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
    prop = volume / sum(volume)
  ) %>%
  arrange(light_mid)

fish_light_pop_SP <- archival_essp_full %>%
  filter(lake == "SP",
         month %in% 5:9,
         !is.na(I_pred_lux),
         I_pred_lux > 0) %>%
  mutate(
    light_bin = cut(I_pred_lux, breaks = light_breaks, include.lowest = TRUE, right = FALSE)
  ) %>%
  count(light_bin, name = "n") %>%
  right_join(bin_lookup_light, by = "light_bin") %>%
  mutate(
    n = tidyr::replace_na(n, 0),
    prop = n / sum(n)
  ) %>%
  arrange(light_mid)

fish_light_ind_SP <- archival_essp_full %>%
  filter(lake == "SP",
         month %in% 5:9,
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
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(fish_id, light_mid)

ggplot() +
  geom_area(
    data = avail_light_SP,
    aes(x = light_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_light_pop_SP,
    aes(x = light_mid, y = prop),
    color = "steelblue",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_light_ind_SP,
    aes(x = light_mid, y = prop, color = fish_id, group = fish_id),
    linewidth = 1
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
    title = "Sparkling Lake – Light Use vs Availability"
  )


#Just daytime light - 

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

ggplot() +
  geom_area(
    data = avail_light_ES_day,
    aes(x = light_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_light_pop_ES_day,
    aes(x = light_mid, y = prop),
    color = "steelblue",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_light_ind_ES_day,
    aes(x = light_mid, y = prop, color = fish_id, group = fish_id),
    linewidth = 1
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
    title = "Escanaba Lake (08:00–20:00) – Light Use vs Availability"
  )

#Onto sparkles 

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

ggplot() +
  geom_area(
    data = avail_light_SP_day,
    aes(x = light_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_light_pop_SP_day,
    aes(x = light_mid, y = prop),
    color = "steelblue",
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_light_ind_SP_day,
    aes(x = light_mid, y = prop, color = fish_id, group = fish_id),
    linewidth = 1
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
    title = "Sparkling Lake (08:00–20:00) – Light Use vs Availability"
  )

#Oxygen density plots 
str(oxygen_df)


oxygen_clean <- oxygen_df %>%
  mutate(
    datetime = mdy_hm(date_time),
    lake = case_when(
      lake == "Escanaba" ~ "ES",
      lake == "Sparkling" ~ "SP",
      TRUE ~ lake
    )
  ) %>%
  filter(!is.na(datetime), !is.na(depth_m), !is.na(do_mgl))

max_depths <- essp_hypso %>%
  group_by(lake) %>%
  summarise(max_depth = max(depth_m), .groups = "drop")

depth_grid <- max_depths %>%
  rowwise() %>%
  mutate(depth_seq = list(seq(0, floor(max_depth), by = 1))) %>%
  tidyr::unnest(depth_seq) %>%
  rename(depth_m = depth_seq)

interp_do_profile <- function(df, grid_depths) {
  df <- df %>%
    arrange(depth_m) %>%
    distinct(depth_m, .keep_all = TRUE)
  
  if (nrow(df) < 2) {
    return(tibble(depth_m = numeric(0), do_mgl = numeric(0)))
  }
  
  out <- approx(
    x = df$depth_m,
    y = df$do_mgl,
    xout = grid_depths,
    rule = 1
  )
  
  tibble(
    depth_m = out$x,
    do_mgl = out$y
  ) %>%
    filter(!is.na(do_mgl))
}

oxygen_vertical <- oxygen_clean %>%
  group_by(lake, datetime) %>%
  group_modify(~{
    lake_i <- .y$lake[[1]]
    
    grid_i <- depth_grid %>%
      filter(lake == lake_i) %>%
      pull(depth_m)
    
    interp_do_profile(.x, grid_i)
  }) %>%
  ungroup()

interp_time_series <- function(df, target_times) {
  df <- df %>%
    arrange(datetime) %>%
    distinct(datetime, .keep_all = TRUE)
  
  if (nrow(df) < 2) {
    return(tibble(datetime = as.POSIXct(character()), do_mgl = numeric()))
  }
  
  out <- approx(
    x = as.numeric(df$datetime),
    y = df$do_mgl,
    xout = as.numeric(target_times),
    rule = 1
  )
  
  tibble(
    datetime = as.POSIXct(out$x, origin = "1970-01-01"),
    do_mgl = out$y
  ) %>%
    filter(!is.na(do_mgl))
}

fish_times <- archival_essp_full %>%
  filter(!is.na(date_time)) %>%
  distinct(lake, datetime = date_time)

oxygen_interp <- oxygen_vertical %>%
  group_by(lake, depth_m) %>%
  group_modify(~{
    lake_i <- .y$lake[[1]]
    
    target_times <- fish_times %>%
      filter(lake == lake_i) %>%
      pull(datetime)
    
    interp_time_series(.x, target_times)
  }) %>%
  ungroup()

fish_do <- archival_essp_full %>%
  filter(!is.na(depth_m), !is.na(date_time)) %>%
  left_join(
    oxygen_interp,
    by = c("lake", "depth_m", "date_time" = "datetime")
  ) %>%
  filter(!is.na(do_mgl))

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

ggplot() +
  geom_area(
    data = avail_do_ES,
    aes(x = do_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_do_pop_ES,
    aes(x = do_mid, y = prop),
    color = "steelblue",
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_do_ind_ES,
    aes(x = do_mid, y = prop, color = fish_id),
    linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  labs(
    x = "Dissolved Oxygen (mg/L)",
    y = "Density"
  )


#SP stuff for oxy 

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

ggplot() +
  geom_area(
    data = avail_do_SP,
    aes(x = do_mid, y = prop),
    fill = "grey70",
    alpha = 0.6
  ) +
  geom_line(
    data = fish_do_pop_SP,
    aes(x = do_mid, y = prop),
    color = "steelblue",
    linetype = "dashed"
  ) +
  geom_line(
    data = fish_do_ind_SP,
    aes(x = do_mid, y = prop, color = fish_id),
    linewidth = 1
  ) +
  facet_wrap(~fish_id) +
  theme_classic() +
  labs(
    x = "Dissolved Oxygen (mg/L)",
    y = "Density"
  )

#Ze big one


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


#Figs for all feesh 



#ze big one for consequences - hopefully don't need conseq_code.R

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



#Nicer plots for all habitat use of all individuals 

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
light_day <- light_es_minute %>%
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
  ungroup()# Fish (ES_01)
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
light_day <- light_sp_minute %>%
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

daily_steps <- daily_steps %>%
  mutate(
    prev_shift = lag(shift_type),
    switched = shift_type != prev_shift
  )

ggplot(daily_steps,
       aes(x = switched, y = step_dist)) +
  geom_boxplot() +
  theme_classic()


####Pareto Frontier####

fish_traits <- archival_length %>%
  select(fish_id, lake, sex, tl_pre) %>%
  distinct()
fish_traits <- fish_traits %>%
  mutate(sex = factor(sex))

#Build light - temp grid 

#Interpolate temperature to 0.1m 

#temp_interp <- temp_vertical %>%
#  arrange(lake, datetime_hour, depth_m) %>%
#  group_by(lake, datetime_hour) %>%
##  group_modify(~{
#    prof <- .x %>%
#      distinct(depth_m, .keep_all = TRUE) %>%
#      arrange(depth_m)
#    
#    new_depths <- seq(min(prof$depth_m), max(prof$depth_m), by = 0.1)
#    
#    interp <- approx(
#      x = prof$depth_m,
#      y = prof$temp_c,
#      xout = new_depths,
#      method = "linear",
#      rule = 1
#    )
    
#    tibble(
#      depth_m = interp$x,

#      temp_c = interp$y
 #   )
#  }) %>%
#  ungroup()

#light_es_minute$lake <- "ES"
#light_sp_minute$lake <- "SP"

#light_essp <- bind_rows(light_es_minute, light_sp_minute)

#Check and find na's 
sum(is.na(light_essp$I_pred_lux))
mean(is.na(light_essp$I_pred_lux))

light_essp %>%
  group_by(lake) %>%
  summarise(
    n = n(),
    n_na = sum(is.na(I_pred_lux)),
    prop_na = mean(is.na(I_pred_lux))
  )

light_hourly_check <- light_essp %>%
  mutate(
    datetime_hour = floor_date(date_time, unit = "hour"),
    depth_m = round(depth_m, 1)
  ) %>%
  group_by(lake, datetime_hour, depth_m) %>%
  summarise(
    n_obs = n(),
    n_nonmiss = sum(!is.na(I_pred_lux)),
    n_miss = sum(is.na(I_pred_lux)),
    light_lux = mean(I_pred_lux, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(light_lux = ifelse(n_nonmiss == 0, NA_real_, light_lux))

#More checks 
light_hourly_check %>%
  summarise(groups_with_no_data = sum(n_nonmiss == 0))

#light_hourly <- light_essp %>%
#  mutate(datetime_hour = floor_date(date_time, unit = "hour")) %>%
#  group_by(lake, datetime_hour, depth_m) %>%
#  summarise(
#    light_lux = mean(I_pred_lux, na.rm = TRUE),
#    .groups = "drop"
#  )

#The joining 
env_hourly <- temp_interp %>%
  left_join(light_hourly, by = c("lake", "datetime_hour", "depth_m"))

#Check - build around NA's? 

env_hourly %>%
  filter(is.na(light_lux)) %>%
  count(lake, depth_m) %>%
  arrange(lake, depth_m)

env_hourly %>%
  filter(is.na(light_lux)) %>%
  count(lake, datetime_hour) %>%
  arrange(lake, datetime_hour)

light_hourly_check %>%
  group_by(lake, datetime_hour) %>%
  summarise(
    min_light_depth = min(depth_m, na.rm = TRUE),
    max_light_depth = max(depth_m, na.rm = TRUE),
    n_light_depths = n(),
    .groups = "drop"
  )


#Set some limits 
light_limits <- light_hourly_check %>%
  group_by(lake, datetime_hour) %>%
  summarise(
    min_light_depth = min(depth_m, na.rm = TRUE),
    max_light_depth = max(depth_m, na.rm = TRUE),
    .groups = "drop"
  )

temp_limits <- temp_interp %>%
  group_by(lake, datetime_hour) %>%
  summarise(
    min_temp_depth = min(depth_m, na.rm = TRUE),
    max_temp_depth = max(depth_m, na.rm = TRUE),
    .groups = "drop"
  )

depth_limits <- temp_limits %>%
  left_join(light_limits, by = c("lake", "datetime_hour")) %>%
  mutate(
    min_shared_depth = pmax(min_temp_depth, min_light_depth),
    max_shared_depth = pmin(max_temp_depth, max_light_depth)
  )

env_hourly <- temp_interp %>%
  left_join(light_hourly_check, by = c("lake", "datetime_hour", "depth_m")) %>%
  left_join(depth_limits, by = c("lake", "datetime_hour")) %>%
  filter(
    depth_m >= min_shared_depth,
    depth_m <= max_shared_depth
  ) %>%
  filter(!is.na(light_lux))

env_coverage <- env_hourly %>%
  group_by(lake, datetime_hour) %>%
  summarise(
    min_depth = min(depth_m, na.rm = TRUE),
    max_depth = max(depth_m, na.rm = TRUE),
    n_depths = n(),
    .groups = "drop"
  )


# 1. Unique fish sizes

fish_sizes <- fish_traits %>%
  select(fish_id, lake, tl_pre) %>%
  distinct()

# 2. Build env space per fish
build_env_space_fish <- function(fish_id_val, tl_val, lake_val) {
  
  env_hourly %>%
    filter(lake == lake_val) %>%
    mutate(
      fish_id = fish_id_val,
      tl_mm = tl_val,
      
      temp_c = pmin(temp_c, RTM - 0.5),
      mass_g = length_mm_to_mass_g(tl_mm),
      
      # energy
      R_gO2_g_h = respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT),
      E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
      
      # swimming
      opt_temp = 22,
      U_opt = 0.4,
      thermal_perf = exp(-((temp_c - opt_temp)^2) / (2 * 6^2)),
      swim_speed_m_s = U_opt * thermal_perf,
      
      # light
      reaction_dist_m = 0.05 + 0.08 * log10(light_lux + 1),
      reaction_dist_m = pmin(reaction_dist_m, 0.30),
      reaction_dist_m = pmax(reaction_dist_m, 0.05),
      
      # volume
      distance_m_h = swim_speed_m_s * 3600,
      volume_searched_m3_h = pi * (reaction_dist_m^2) * distance_m_h
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

#Gotta rebuild this hourly bogus 
str(archival_essp)

archival_essp_hourly <- archival_essp %>%
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


build_env_space_size <- function(tl_val) {
  env_hourly %>%
    mutate(
      tl_mm = tl_val,
      
      # ---- Fix temperature blow-ups ----
      temp_c = pmin(temp_c, RTM - 0.5),
      
      # ---- Mass ----
      mass_g = length_mm_to_mass_g(tl_mm),
      
      # ---- Energy ----
      R_gO2_g_h = respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT),
      E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
      
      # ---- Swimming performance ----
      opt_temp = 22,
      U_opt = 0.4,
      thermal_perf = exp(-((temp_c - opt_temp)^2) / (2 * 6^2)),
      swim_speed_m_s = U_opt * thermal_perf,
      
      # ---- Reaction distance ----
      reaction_dist_m = 0.05 + 0.08 * log10(light_lux + 1),
      reaction_dist_m = pmin(reaction_dist_m, 0.30),
      reaction_dist_m = pmax(reaction_dist_m, 0.05),
      
      # ---- Volume ----
      distance_m_h = swim_speed_m_s * 3600,
      volume_searched_m3_h = pi * (reaction_dist_m^2) * distance_m_h
    ) %>%
    filter(
      !is.na(E_kJ_h),
      !is.na(volume_searched_m3_h),
      E_kJ_h > 0,
      E_kJ_h < 50,
      light_lux > 0
    )
}


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
  geom_smooth(
    data = envelope_df,
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
ggplot(distance_df, aes(x = fish_id, y = dist_to_frontier /or/ rel_efficiency, fill = lake)) +
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
  env_hourly %>%
    mutate(
      tl_mm = tl_val,
      temp_c = pmin(temp_c, RTM - 0.5),
      mass_g = length_mm_to_mass_g(tl_mm),
      
      R_gO2_g_h = respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT),
      E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
      
      opt_temp = 22,
      U_opt = 0.4,
      thermal_perf = exp(-((temp_c - opt_temp)^2) / (2 * 6^2)),
      swim_speed_m_s = U_opt * thermal_perf,
      
      reaction_dist_m = 0.05 + 0.08 * log10(light_lux + 1),
      reaction_dist_m = pmin(reaction_dist_m, 0.30),
      reaction_dist_m = pmax(reaction_dist_m, 0.05),
      
      distance_m_h = swim_speed_m_s * 3600,
      volume_searched_m3_h = pi * (reaction_dist_m^2) * distance_m_h
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

ribbon_df <- frontier_bin_df %>%
  pivot_wider(names_from = bound, values_from = V_max)

frontier_bin_df <- bind_rows(frontier_min_df, frontier_max_df)

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
    aes(x = E_kJ_h, ymin = min, ymax = max),
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

fish_info <- archival_length %>%
  filter(fish_id == fish_id_sel) %>%
  slice(1)

lake_sel <- fish_info$lake
tl_mm    <- fish_info$tl_pre

# 2. GET FISH MASS
mass_g <- length_mm_to_mass_g(tl_mm)


# 3. SUBSET ENVIRONMENT

env_day <- env_hourly %>%
  mutate(date = as.Date(datetime_hour)) %>%
  filter(lake == lake_sel, date == day_sel)


# 4. BUILD FULL DEPTH GRID

depth_seq <- seq(
  floor(min(env_day$depth_m, na.rm = TRUE)*10)/10,
  ceiling(max(env_day$depth_m, na.rm = TRUE)*10)/10,
  by = 0.1
)

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

env_perf <- env_interp %>%
  mutate(
    mass_g = mass_g,
    
    # ENERGY
    R_gO2_g_h = respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT),
    E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
    
    # SWIMMING
    thermal_perf = exp(-((temp_c - 22)^2) / (2 * 6^2)),
    swim_speed_m_s = 0.4 * thermal_perf,
    
    # REACTION DISTANCE
    reaction_dist_m = 0.05 + 0.08 * log10(light_lux + 1),
    reaction_dist_m = pmin(reaction_dist_m, 0.30),
    reaction_dist_m = pmax(reaction_dist_m, 0.05),
    
    # VOLUME
    distance_m_h = swim_speed_m_s * 3600,
    volume_m3_h = pi * (reaction_dist_m^2) * distance_m_h,
    
    # EFFICIENCY
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
    data = env_perf,
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
    data = env_perf,
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
    data = env_perf,
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


#Density plots for vol searched and energetic requirements? 

hypso_lookup <- hypso_layers %>%
  select(lake, depth_top, depth_bottom, layer_volume_m3)

env_perf <- env_perf %>%
  filter(!is.na(lake))

env_perf <- env_perf %>%
  arrange(lake, depth_m)

env_perf <- env_perf %>%
  group_split(lake) %>%
  lapply(function(df) {
    
    lk <- unique(df$lake)
    hyp <- hypso_lookup %>% filter(lake == lk)
    
    idx <- findInterval(df$depth_m, hyp$depth_top)
    
    idx[idx == 0] <- 1
    idx[idx > nrow(hyp)] <- nrow(hyp)
    
    df$layer_volume_m3 <- hyp$layer_volume_m3[idx]
    
    df
  }) %>%
  bind_rows()

env_perf <- env_perf %>%
  mutate(depth_bin = floor(depth_m)) %>%
  group_by(datetime_hour, lake, depth_bin) %>%
  mutate(
    n_sublayers = n(),
    sublayer_volume = layer_volume_m3 / n_sublayers
  ) %>%
  ungroup()

env_perf <- env_perf %>%
  group_by(datetime_hour, lake) %>%
  mutate(
    vol_weight = sublayer_volume / sum(sublayer_volume, na.rm = TRUE)
  ) %>%
  ungroup()

fish_sub <- fish_traits %>%
  filter(lake == unique(env_perf$lake)) %>%
  select(fish_id, tl_pre)

env_pop <- env_perf %>%
  crossing(fish_sub) %>%
  mutate(
    mass_g = length_mm_to_mass_g(tl_pre),
    
    R_gO2_g_h = respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT),
    E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
    
    efficiency = volume_m3_h / E_kJ_h
  )

obs_fish <- archival_essp_hourly %>%
  left_join(fish_traits, by = c("fish_id", "lake")) %>%
  mutate(
    mass_g = length_mm_to_mass_g(tl_pre),
    
    R_gO2_g_h = respiration_FB4(mass_g, mean_temp, RA, RB, RQ, RTO, RTM, ACT),
    E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
    
    volume_m3_h = mean_volume
  )

ggplot() +
  geom_density(
    data = env_pop,
    aes(x = volume_m3_h, weight = vol_weight),
    fill = "gray70",
    alpha = 0.6
  ) +
  geom_density(
    data = obs_fish,
    aes(x = volume_m3_h),
    color = "blue",
    linewidth = 1
  ) +
  facet_wrap(~ fish_id, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Volume searched (m³/h)",
    y = "Density",
    title = "Available vs used volume searched"
  )

ggplot() +
  geom_density(
    data = env_pop,
    aes(x = E_kJ_h, weight = vol_weight),
    fill = "gray70",
    alpha = 0.6
  ) +
  geom_density(
    data = obs_fish,
    aes(x = E_kJ_h),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ fish_id, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Energy (kJ/h)",
    y = "Density",
    title = "Available vs used energy"
  )


#The above works, but let's see if we can spice this up 


# Attach hypsography to env_perf

hypso_lookup <- hypso_layers %>%
  select(lake, depth_top, depth_bottom, layer_volume_m3)

env_perf <- env_perf %>%
  group_split(lake) %>%
  lapply(function(df) {
    
    lk <- unique(df$lake)
    hyp <- hypso_lookup %>% filter(lake == lk)
    
    if (nrow(hyp) == 0) {
      stop(paste("No hypsography for lake:", lk))
    }
    
    # match each depth to its layer
    idx <- findInterval(df$depth_m, hyp$depth_top)
    
    # handle edge cases
    idx[idx == 0] <- 1
    idx[idx > nrow(hyp)] <- nrow(hyp)
    
    df$layer_volume_m3 <- hyp$layer_volume_m3[idx]
    
    df
  }) %>%
  bind_rows()

env_perf %>%
  summarise(missing = sum(is.na(layer_volume_m3)))

env_perf <- env_perf %>%
  mutate(depth_bin = floor(depth_m))

env_perf <- env_perf %>%
  group_by(datetime_hour, lake, depth_bin) %>%
  mutate(
    n_sublayers = n(),
    sublayer_volume = layer_volume_m3 / n_sublayers
  ) %>%
  ungroup()

env_perf <- env_perf %>%
  group_by(datetime_hour, lake) %>%
  mutate(
    vol_weight = sublayer_volume / sum(sublayer_volume, na.rm = TRUE)
  ) %>%
  ungroup()

#Checks 
env_perf %>%
  group_by(datetime_hour, lake) %>%
  summarise(sum_w = sum(vol_weight)) %>%
  summary()

env_perf %>% count(lake)

env_pop <- env_perf %>%
  inner_join(valid_fish, by = "lake", multiple = "all") %>%
  mutate(
    mass_g = length_mm_to_mass_g(tl_pre),
    
    R_gO2_g_h = respiration_FB4(
      mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT
    ),
    
    E_kJ_h = (R_gO2_g_h * 13.6) * mass_g
  )

obs_fish <- archival_essp_hourly %>%
  filter(fish_id %in% valid_fish$fish_id) %>%
  left_join(valid_fish, by = c("fish_id", "lake")) %>%
  mutate(
    mass_g = length_mm_to_mass_g(tl_pre),
    
    R_gO2_g_h = respiration_FB4(
      mass_g, mean_temp, RA, RB, RQ, RTO, RTM, ACT
    ),
    
    E_kJ_h = (R_gO2_g_h * 13.6) * mass_g,
    volume_m3_h = mean_volume
  )


ggplot() +
  
  geom_density(
    data = env_pop,
    aes(x = volume_m3_h, weight = vol_weight),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  geom_density(
    data = obs_fish,
    aes(x = volume_m3_h),
    color = "blue",
    linewidth = 1
  ) +
  
  facet_wrap(~ lake, scales = "free") +
  
  theme_classic()

ggplot() +
  
  # AVAILABLE
  geom_density(
    data = env_pop,
    aes(y = volume_m3_h, weight = vol_weight),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  # POPULATION
  geom_density(
    data = obs_fish,
    aes(y = volume_m3_h),
    fill = "steelblue",
    alpha = 0.6,
    linetype = "dashed"
  ) +
  
  # INDIVIDUAL
  geom_density(
    data = obs_fish,
    aes(y = volume_m3_h, color = fish_id),
    linewidth = 1
  ) +
  
  facet_grid(lake ~ fish_id, scales = "free_y") +
  
  theme_classic() +
  theme(legend.position = "none")

env_ES <- env_pop %>% filter(lake == "ES")
env_SP <- env_pop %>% filter(lake == "SP")

obs_ES <- obs_fish %>% filter(lake == "ES")
obs_SP <- obs_fish %>% filter(lake == "SP")

ggplot() +
  
  geom_density(
    data = env_ES,
    aes(x = volume_m3_h, weight = vol_weight),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  geom_density(
    data = obs_ES,
    aes(x = volume_m3_h),
    alpha = 0.6,
    linetype = "dashed"
  ) +
  
  geom_density(
    data = obs_ES,
    aes(x = volume_m3_h, color = fish_id),
    linewidth = 1
  ) +
  
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~fish_id)+
  
  labs(
    x = "Volume searched (m³/h)",
    y = "Density",
    title = "ES – Volume searched (Availability vs Population vs Individuals)"
  )

ggplot() +
  
  geom_density(
    data = env_SP,
    aes(x = volume_m3_h, weight = vol_weight),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  geom_density(
    data = obs_SP,
    aes(x = volume_m3_h),
    alpha = 0.6,
    linetype = "dashed"
  ) +
  
  geom_density(
    data = obs_SP,
    aes(x = volume_m3_h, color = fish_id),
    linewidth = 1
  ) +
  
  theme_classic() +
  theme(legend.position = "none") +
   facet_wrap(~fish_id)+
  
  labs(
    x = "Volume searched (m³/h)",
    y = "Density",
    title = "SP – Volume searched (Availability vs Population vs Individuals)"
  )

ggplot() +
  
  geom_density(
    data = env_ES,
    aes(x = E_kJ_h, weight = vol_weight),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  geom_density(
    data = obs_ES,
    aes(x = E_kJ_h),
    alpha = 0.6,
    linetype = "dashed"
  ) +
  
  geom_density(
    data = obs_ES,
    aes(x = E_kJ_h, color = fish_id),
    linewidth = 1
  ) +
  
  theme_classic() +
  theme(legend.position = "none") +
   facet_wrap(~fish_id)+
  
  labs(
    x = "Energy (kJ/h)",
    y = "Density",
    title = "ES – Energy use (Availability vs Population vs Individuals)"
  )

ggplot() +
  
  geom_density(
    data = env_SP,
    aes(x = E_kJ_h, weight = vol_weight),
    fill = "grey70",
    alpha = 0.6
  ) +
  
  geom_density(
    data = obs_SP,
    aes(x = E_kJ_h),
    alpha = 0.6,
    linetype = "dashed"
  ) +
  
  geom_density(
    data = obs_SP,
    aes(x = E_kJ_h, color = fish_id),
    linewidth = 1
  ) +
  
  theme_classic() +
  theme(legend.position = "none") +
   facet_wrap(~fish_id)+
  
  labs(
    x = "Energy (kJ/h)",
    y = "Density",
    title = "SP – Energy use (Availability vs Population vs Individuals)"
  )


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



####Density OVL fanciness####

# ============================================================
# INDIVIDUAL HABITAT-USE DATA FOR OVL ANALYSIS
# ============================================================

# Depth
ovl_depth <- fish_depth_dist %>%
  filter(
    lake %in% c("ES", "SP"),
    !is.na(depth_m),
    lake != "ES" | depth_m <= 9
  ) %>%
  select(fish_id, lake, value = depth_m)

# Temperature
ovl_temp <- archival_essp_full %>%
  filter(
    lake %in% c("ES", "SP"),
    year == 2022,
    month %in% 5:9,
    !is.na(temp_c)
  ) %>%
  select(fish_id, lake, value = temp_c)

# Light
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
  select(fish_id, lake, value = I_pred_lux)

# Dissolved oxygen
ovl_do <- fish_do %>%
  filter(
    lake %in% c("ES", "SP"),
    !is.na(do_mgl)
  ) %>%
  select(fish_id, lake, value = do_mgl)


# Combine
habitat_use <- bind_rows(
  ovl_temp %>% mutate(habitat = "Temperature"),
  ovl_depth %>% mutate(habitat = "Depth"),
  ovl_light %>% mutate(habitat = "Light"),
  ovl_do %>% mutate(habitat = "Oxygen")
) %>%
  filter(fish_id %in% valid_fish$fish_id) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c("Temperature", "Depth", "Light", "Oxygen")
    )
  )

# ============================================================
# FUNCTION TO CALCULATE OVL
# ============================================================

calculate_ovl <- function(df, n_grid = 500, n_min = 20) {

  # Remove missing values
  df <- df %>%
    filter(is.finite(value))

  # Common range
  x_min <- min(df$value, na.rm = TRUE)
  x_max <- max(df$value, na.rm = TRUE)

  # Common grid
  x_grid <- seq(x_min, x_max, length.out = n_grid)

  # ----------------------------------------------------------
  # Individual densities
  # ----------------------------------------------------------

  ind_dens <- df %>%
    group_by(fish_id, lake) %>%
    group_modify(~{

      x <- .x$value

      if (
        length(x) < n_min ||
        length(unique(x)) < 2
      ) {
        return(tibble())
      }

      d <- density(
        x,
        from = x_min,
        to = x_max,
        n = n_grid
      )

      # Normalize density to area = 1
      area <- sum(
        diff(d$x) *
          (head(d$y, -1) + tail(d$y, -1)) / 2
      )

      tibble(
        x = d$x,
        density = d$y / area
      )

    }) %>%
    ungroup()

  # ----------------------------------------------------------
  # Equal-weight population density
  # ----------------------------------------------------------

  pop_dens <- ind_dens %>%
    group_by(lake, x) %>%
    summarise(
      pop_density = mean(density),
      n_fish = n(),
      .groups = "drop"
    )

  # ----------------------------------------------------------
  # Individual vs population overlap
  # ----------------------------------------------------------

  ovl <- ind_dens %>%
    left_join(
      pop_dens,
      by = c("lake", "x")
    ) %>%
    group_by(lake, fish_id) %>%
    arrange(x) %>%
    summarise(
      OVL = sum(
        diff(x) *
          (
            pmin(
              head(density, -1),
              head(pop_density, -1)
            ) +
            pmin(
              tail(density, -1),
              tail(pop_density, -1)
            )
          ) / 2,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  ovl
}


# ============================================================
# RUN FOR ALL FOUR HABITAT VARIABLES
# ============================================================

ovl_results <- habitat_use %>%
  group_by(habitat) %>%
  group_modify(~calculate_ovl(.x)) %>%
  ungroup()

ovl_results <- ovl_results %>%
  left_join(
    valid_fish %>%
      select(fish_id, lake, tl_pre),
    by = c("fish_id", "lake")
  )


ggplot(
  ovl_results,
  aes(x = habitat, y = OVL)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.2
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.5
  ) +
  facet_wrap(~lake) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Overlap coefficient",
    title = "Individual habitat use overlap with population use"
  )



# ============================================================
# INDIVIDUAL HABITAT-USE OVERLAP
# LEAVE-ONE-FISH-OUT OVL
# ============================================================


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


#Directionality with OVL analysis 

# ------------------------------------------------------------
# 1. Standardize habitat-use data
# ------------------------------------------------------------

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
    select(fish_id, lake, value = I_pred_lux) %>%
    mutate(habitat = "Light"),
  
  # Oxygen
  fish_do %>%
    filter(!is.na(do_mgl)) %>%
    select(fish_id, lake, value = do_mgl) %>%
    mutate(habitat = "Oxygen")
)


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


# ------------------------------------------------------------
# 3. Run for each lake × habitat
# ------------------------------------------------------------

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
  unnest(results)


# ------------------------------------------------------------
# 4. Add fish traits
# ------------------------------------------------------------

directional_results <- directional_results %>%
  left_join(
    fish_traits %>%
      select(fish_id, tl_pre),
    by = "fish_id"
  ) %>%
  relocate(
    habitat,
    lake,
    fish_id,
    tl_pre
  )


# ------------------------------------------------------------
# 5. Look at the results
# ------------------------------------------------------------

directional_results


# ------------------------------------------------------------
# 6. Quick check: OVL should match your existing results
# ------------------------------------------------------------

directional_results %>%
  select(
    habitat,
    lake,
    fish_id,
    OVL,
    mean_individual,
    mean_population,
    delta_mean,
    underuse_amount,
    overuse_amount,
    underuse_center,
    overuse_center
  ) %>%
  arrange(habitat, lake, fish_id)


#Making figs 

# ------------------------------------------------------------
# Prepare directional OVL results
# ------------------------------------------------------------

plot_df <- directional_results %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c("Temperature", "Depth", "Light", "Oxygen")
    ),
    
    lake = factor(
      lake,
      levels = c("ES", "SP")
    )
  )


# ------------------------------------------------------------
# OVL vs. directional habitat shift
# ------------------------------------------------------------

p_ovl_direction <- ggplot(
  plot_df,
  aes(
    x = OVL,
    y = delta_mean,
    color = lake
  )
) +
  
  # Zero = no directional shift from population
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # Individual fish
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  
  # Labels for particularly divergent individuals
  geom_text(
    aes(label = fish_id),
    vjust = -0.8,
    size = 3,
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

#or try: 

p_ovl_direction_clean <- ggplot(
  plot_df,
  aes(
    x = OVL,
    y = delta_mean,
    color = lake
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_point(
    size = 2.5,
    alpha = 0.8
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

p_ovl_direction_clean

#Checking out individuals: 

# ------------------------------------------------------------
# Identify the most directionally divergent fish
# ------------------------------------------------------------

directional_results %>%
  group_by(habitat) %>%
  arrange(delta_mean) %>%
  slice_head(n = 5) %>%
  select(
    habitat, lake, fish_id, OVL,
    mean_individual,
    mean_population,
    delta_mean,
    underuse_center,
    overuse_center
  )

directional_results %>%
  group_by(habitat) %>%
  arrange(desc(delta_mean)) %>%
  slice_head(n = 5) %>%
  select(
    habitat, lake, fish_id, OVL,
    mean_individual,
    mean_population,
    delta_mean,
    underuse_center,
    overuse_center
  )

directional_results %>%
  group_by(habitat) %>%
  arrange(OVL) %>%
  slice_head(n = 5) %>%
  select(
    habitat, lake, fish_id, OVL,
    mean_individual,
    mean_population,
    delta_mean,
    underuse_center,
    overuse_center
  )

directional_results %>%
  filter(habitat == "Temperature") %>%
  arrange(delta_mean) %>%
  select(
    lake,
    fish_id,
    OVL,
    mean_individual,
    mean_population,
    delta_mean,
    underuse_center,
    overuse_center
  )


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

