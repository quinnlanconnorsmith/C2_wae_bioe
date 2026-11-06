# ============================================================
# 00_functions.R
# Shared libraries, parameters and helper functions.
# Sourced at the top of every other script - nothing here
# creates data, so it is cheap to run.
# ============================================================

library(tidyverse)
library(lubridate)
library(zoo)
library(patchwork)

# Where intermediate .rds files live (written by 01_data_prep.R)
derived_dir <- "derived"

read_derived <- function(name) {
  readRDS(file.path(derived_dir, paste0(name, ".rds")))
}

save_derived <- function(obj, name) {
  saveRDS(obj, file.path(derived_dir, paste0(name, ".rds")))
}

# ------------------------------------------------------------
# FB4 respiration model for Walleye (Sander vitreus)
# Based on Kitchell et al. 1977 (Equation 2)
# ------------------------------------------------------------

RA  <- 0.0108  # Intercept of allometric function for SMR
RB  <- -0.2    # Slope of allometric function for SMR
RQ  <- 2.1     # Q10-type rate at which the function increases
RTO <- 27      # Optimum temperature for respiration
RTM <- 32      # Lethal water temperature
ACT <- 1       # Activity multiplier

# Oxycalorific coefficient (kJ per g O2)
OXY_KJ <- 13.6

# Length-weight conversion: W = a * L^b (L in cm, W in g)
# FishBase Bayesian parameters for walleye
# TODO: replace with data from GS, MW, JTM?
lw_a <- 0.00631
lw_b <- 3.15

length_mm_to_mass_g <- function(tl_mm, a = lw_a, b = lw_b) {
  a * ((tl_mm / 10) ^ b)
}

TF_calc <- function(T_C, RTO, RTM) {
  ((T_C - RTO) / (RTM - RTO)) / (1 - ((T_C - RTO) / (RTM - RTO)))
}

respiration_FB4 <- function(W_g, T_C, RA, RB, RQ, RTO, RTM, ACT) {
  RA * (W_g ^ RB) * exp(RQ * TF_calc(T_C, RTO, RTM)) * ACT
}

# Hourly energetic requirement (kJ/h) for a fish of a given mass
energy_kJ_h <- function(mass_g, temp_c) {
  respiration_FB4(mass_g, temp_c, RA, RB, RQ, RTO, RTM, ACT) *
    OXY_KJ * mass_g
}

# ------------------------------------------------------------
# Environmental search performance
#
# Same formulas that were repeated in build_env_space_fish(),
# build_env_space_size(), build_frontier() and env_perf:
#   swim speed  = 0.4 m/s scaled by a Gaussian thermal
#                 performance curve (optimum 22 C, sd 6)
#   reaction    = 0.05 + 0.08 * log10(lux + 1), bounded 0.05-0.30 m
#   volume      = pi * reaction^2 * distance swum per hour
#
# Adds: thermal_perf, swim_speed_m_s, reaction_dist_m,
#       distance_m_h, volume_m3_h
# ------------------------------------------------------------

add_search_performance <- function(df,
                                   opt_temp = 22,
                                   U_opt = 0.4,
                                   perf_sd = 6) {
  df %>%
    mutate(
      thermal_perf    = exp(-((temp_c - opt_temp)^2) / (2 * perf_sd^2)),
      swim_speed_m_s  = U_opt * thermal_perf,
      reaction_dist_m = 0.05 + 0.08 * log10(light_lux + 1),
      reaction_dist_m = pmin(pmax(reaction_dist_m, 0.05), 0.30),
      distance_m_h    = swim_speed_m_s * 3600,
      volume_m3_h     = pi * (reaction_dist_m^2) * distance_m_h
    )
}

# ------------------------------------------------------------
# Interpolation helpers (temperature, oxygen)
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Light attenuation: k(z) and predicted light grid for one
# HOBO time step. Only used by archive/light_grids_one_time.R
# ------------------------------------------------------------

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
