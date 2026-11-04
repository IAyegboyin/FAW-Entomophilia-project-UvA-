# ============================================================================
# EDA: Effect of Serratia (low/high dose) on Spodoptera frugiperda (FAW)
# Project: FAW-Entomophilia-project-UvA
# Author: Ismail A. Ayegboyin (ayegboyinismail@gmail.com)
# Date: Spetmeber 23rd, 2026
# ============================================================================

# libraries and work directory settings  --------------------------------------

library(tidyverse)
library(readxl)
library(writexl)

#pathways
project_dir <- "~/Desktop/Data Science Library/Data for play/FAW-Entomophilia-project-UvA"
data_path   <- file.path(project_dir, "Data/Serratia effect on SF_31-8-26.xlsx")

fig_dir   <- file.path(project_dir, "Figures and Notes")
table_dir <- file.path(project_dir, "Figures and Notes/Tables")
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# colour palletes
treat_colors <- c(
  "PBS"      = "#4D4D4D",
  "lowdose"  = "#1B9E77",
  "highdose" = "#D95F02"
)


# Load raw data --------------------------------------------------------

raw <- readxl::read_excel(data_path, sheet = "combined treatments")
# view(raw) commented out just to see what the data looks like first finally in R
glimpse(raw)

# Rename columns (positional) -----------------------------------------
# The original Excel headers contain stray spaces/parentheses (e.g.
# "larvalwt_day3_ 26-7-26)") which makes them fragile to reference by name.
# We rename by POSITION instead, using the column order confirmed by
# glimpse(). If your sheet's column order ever changes, this stopifnot()
# will fail loudly rather than silently mis-labelling a column.

new_names <- c(
  "insect", "insect_family", "treatment",
  "wt_day0", "wt_day3", "wt_day6",
  "pup_day8", "pup_day9", "pup_day11", 
  "pup_day12", "pup_day13", "days_to_pupation", 
  "mortality_day", "larvae_died", "adult_emergence", 
  "day_emergence", "sex", "mating_pair", 
  "mating_status", "egg_24h", "egg_48h",
  "egg_72h", "adult_mortality_dai", "total_egg_count"
)
stopifnot(length(new_names) == ncol(raw))
names(raw) <- new_names

# glimpse(raw) good! names have changed here already
raw <- raw %>%
  filter(!is.na(treatment)) %>%
  mutate(treatment = factor(treatment, levels = c("PBS", "lowdose", "highdose")))