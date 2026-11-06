# C2_wae_bioe
Hit memory issues with the last one - this is definitely going to stay clean, right? 

## How to run

Open `C2_wae_bioe.Rproj` (workspace saving is turned off for this project) and run, each in a fresh R session:

1. `R/01_data_prep.R` - once, or whenever raw data change. Writes `derived/*.rds`.
2. `R/02_habitat_use.R` - depth, temperature, light and DO use vs availability; trajectories; diel shifts; OVL, directional shifts, PCA.
3. `R/03_consequences.R` - energetic requirements and volume searched; availability vs population vs individuals; Pareto frontiers.

`R/00_functions.R` holds the shared parameters and functions and is sourced by the others.
`archive/light_grids_one_time.R` rebuilds the light grids from the HOBO strings (already done; no need to re-run).
`c2_code_clean.R` and `conseq_code.R` are the original combined scripts, kept for reference.
