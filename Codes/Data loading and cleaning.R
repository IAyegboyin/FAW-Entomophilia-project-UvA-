# ============================================================================
# EDA: Effect of Serratia (low/high dose) on Spodoptera frugiperda (FAW)
# Project: FAW-Entomophilia-project-UvA
# Author: Ismail A. Ayegboyin (iayegboyin@aucegypt.edu)
# Date: Spetmeber 23rd, 2026
# ============================================================================

# librarie and work directory settings  --------------------------------------

library(tidyverse)
library(readxl)
library(writexl)

#pathways
project_dir <- "~/Desktop/Data Science Library/Data for play/FAW-Entomophilia-project-UvA"
data_path   <- file.path(project_dir, "Data/Serratia effect on SF_31-8-26.xlsx")

fig_dir   <- file.path(project_dir, "Figures")
table_dir <- file.path(project_dir, "Tables")
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 12))

# one consistent color per treatment, used in every plot
treat_colors <- c(
  "PBS"      = "#4D4D4D",
  "lowdose"  = "#1B9E77",
  "highdose" = "#D95F02"
)


# Load raw data --------------------------------------------------------

raw <- readxl::read_excel(data_path, sheet = "combined treatments")
view(raw)
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
  "pup_day8", "pup_day9", "pup_day11", "pup_day12", "pup_day13",
  "days_to_pupation", "mortality_day", "larvae_died",
  "adult_emergence", "day_emergence", "sex", "mating_pair", "mating_status",
  "egg_24h", "egg_48h", "egg_72h", "adult_mortality_dai", "total_egg_count"
)
stopifnot(length(new_names) == ncol(raw))
names(raw) <- new_names

raw <- raw %>%
  filter(!is.na(treatment)) %>%
  mutate(treatment = factor(treatment, levels = c("PBS", "lowdose", "highdose")))


# Helper functions -----------------------------------------------------

# Converts a character column to numeric, but FLAGS any value that:
#   - is not a recognised missing token ("NA", "N/A"), and
#   - still fails to parse as a number
# This gives explicit flags instead of silently guessing what a weird
# entry meant (e.g. a stray "no" typed into a numeric column by mistake).
safe_numeric <- function(x) {
  x_chr <- as.character(x)
  is_missing_tok <- is.na(x_chr) | x_chr %in% c("NA", "N/A")
  val  <- suppressWarnings(as.numeric(x_chr))
  flag <- is.na(val) & !is_missing_tok
  list(value = val, flag = flag)
}

# Converts "yes"/"no" text to TRUE/FALSE, anything else (N/A, blank) -> NA
yn_to_logical <- function(x) {
  x <- str_trim(str_to_lower(as.character(x)))
  case_when(
    x == "yes" ~ TRUE,
    x == "no"  ~ FALSE,
    TRUE       ~ NA
  )
}

# Parses egg-batch cells like "8x", "3xx", "1xxx", "0", "N/A" into a count
# and a size category, per the project notes:
#   x = small egg mass, xx = medium egg mass, xxx = big egg mass
parse_egg_batch <- function(x) {
  x <- str_trim(as.character(x))
  count <- rep(NA_real_, length(x))
  size  <- rep(NA_character_, length(x))
  
  is_zero <- x == "0"
  count[is_zero] <- 0
  
  m <- str_match(x, "^([0-9]*)(x{1,3})$")
  has_match <- !is.na(m[, 1]) & !is_zero
  count[has_match] <- ifelse(m[has_match, 2] == "", 1, as.numeric(m[has_match, 2]))
  size[has_match]  <- recode(m[has_match, 3], x = "small", xx = "medium", xxx = "big")
  
  is_missing_tok <- is.na(x) | x %in% c("NA", "N/A")
  flag <- is.na(count) & !is_missing_tok   # value present but didn't match pattern
  
  list(count = count,
       size  = factor(size, levels = c("small", "medium", "big")),
       flag  = flag)
}
