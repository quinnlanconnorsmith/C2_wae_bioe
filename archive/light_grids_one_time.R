# ============================================================
# ONE-TIME SCRIPT - DO NOT RE-RUN UNLESS THE HOBO DATA CHANGE
#
# Builds the k(z) midpoints and 0.1 m predicted light grids from
# the HOBO strings, month by month, for Sparkling and Escanaba.
# Outputs (already in the project folder):
#   sparkling_k_mids.csv, sparkling_light_sum22.csv
#   es_k_mids.csv,        escanaba_light_sum22.csv
#
# This is a condensed version of the "Light things to not run"
# section of c2_code_clean.R (same calculations; the ten copied
# monthly blocks are replaced by one loop). The original
# referred to `all_hobo`; the loaded object is `hobo_all`.
# Run from the project root.
# ============================================================

source("R/00_functions.R")
library(furrr)

hobo_all <- read_csv("hobo_all.csv")

future::plan(multisession, workers = parallel::detectCores() - 1)

month_starts <- as.Date(c("2022-05-01", "2022-06-01", "2022-07-01",
                          "2022-08-01", "2022-09-01", "2022-10-01"))

run_light_month <- function(lake_name, start, end) {
  hobo_all %>%
    mutate(datetime = ymd_hms(datetime)) %>%
    filter(lake == lake_name, datetime >= start, datetime < end) %>%
    select(datetime, temp_c, light_lux, lake, hobo_df_name, depth_m) %>%
    mutate(date_time = datetime) %>%
    arrange(date_time, depth_m) %>%
    group_by(date_time) %>%
    nest() %>%
    mutate(proc = future_map(data, ~ process_time_step_midpoints(.x),
                             .progress = TRUE))
}

run_light_lake <- function(lake_name) {
  map2(head(month_starts, -1), tail(month_starts, -1),
       ~ run_light_month(lake_name, .x, .y)) %>%
    bind_rows() %>%
    mutate(date_time = ymd_hms(date_time))
}

extract_part <- function(results, part) {
  results %>%
    mutate(out = map(proc, part)) %>%
    select(date_time, out) %>%
    unnest(out)
}

# Sparkling
results_sp_all <- run_light_lake("Sparkling")
write.csv(extract_part(results_sp_all, "midpoint_k"),
          "sparkling_k_mids.csv", row.names = FALSE)
write.csv(extract_part(results_sp_all, "light_grid"),
          "sparkling_light_sum22.csv", row.names = FALSE)
rm(results_sp_all); gc()

# Escanaba
results_es_all <- run_light_lake("Escanaba")
write.csv(extract_part(results_es_all, "midpoint_k"),
          "es_k_mids.csv", row.names = FALSE)
write.csv(extract_part(results_es_all, "light_grid"),
          "escanaba_light_sum22.csv", row.names = FALSE)
rm(results_es_all); gc()
