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

#  Cleaning variable groups ------------------------------------------------

# For larval weight data -- "x" means the larva was already dead at that
# weigh-in (per project notes). as.numeric() turns "x" into NA
wt0 <- safe_numeric(raw$wt_day0)
wt3 <- safe_numeric(raw$wt_day3)
wt6 <- safe_numeric(raw$wt_day6)

sf_data <- raw %>%
  mutate(
    died_by_day3 = wt_day3 == "x",
    died_by_day6 = wt_day6 == "x",
    wt_day0 = wt0$value,  wt_day0_flag = wt0$flag,
    wt_day3 = wt3$value,  wt_day3_flag = wt3$flag,
    wt_day6 = wt6$value,  wt_day6_flag = wt6$flag
  )

# Pupation-check columns, larval death, adult emergence, mating status
sf_data <- sf_data %>%
  mutate(
    across(c(pup_day8, pup_day9, pup_day11, pup_day12, pup_day13,
             adult_emergence, larvae_died),
           yn_to_logical),
    mating_status = yn_to_logical(mating_status),
    sex = na_if(str_trim(str_to_lower(sex)), "n/a"),
    sex = factor(sex, levels = c("male", "female"))
  )

# view(sf_data)
# Numeric "days"/"mortality" columns
dtp   <- safe_numeric(df$days_to_pupation)
mday  <- safe_numeric(df$mortality_day)
demg  <- safe_numeric(df$day_emergence)
amort <- safe_numeric(df$adult_mortality_dai)

sf_data <- sf_data %>%
  mutate(
    days_to_pupation     = dtp$value,   days_to_pupation_flag     = dtp$flag,
    mortality_day        = mday$value,  mortality_day_flag        = mday$flag,
    day_emergence        = demg$value,  day_emergence_flag        = demg$flag,
    adult_mortality_dai  = amort$value, adult_mortality_dai_flag  = amort$flag
  ) %>%
  mutate(days_to_pupation_flag = days_to_pupation_flag | (days_to_pupation == 0))
# Egg-batch columns cleaning (count + size category at 24h/48h/72h)
egg24 <- parse_egg_batch(df$egg_24h)
egg48 <- parse_egg_batch(df$egg_48h)
egg72 <- parse_egg_batch(df$egg_72h)

sf_data <- sf_data %>%
  mutate(
    egg_24h_count = egg24$count, egg_24h_size = egg24$size, egg_24h_flag = egg24$flag,
    egg_48h_count = egg48$count, egg_48h_size = egg48$size, egg_48h_flag = egg48$flag,
    egg_72h_count = egg72$count, egg_72h_size = egg72$size, egg_72h_flag = egg72$flag
  )










# ---- 5. Derived life-stage variables ----------------------------------------

sf_data <-sf_data %>%
  mutate(pupated_ever = pup_day13)   # day13 is the final pupation check

survival_funnel <- df %>%
  group_by(treatment) %>%
  summarise(
    n_start             = n(),
    n_larvae_died       = sum(larvae_died, na.rm = TRUE),
    n_pupated           = sum(pupated_ever, na.rm = TRUE),
    n_emerged           = sum(adult_emergence, na.rm = TRUE),
    pct_larvae_died     = 100 * n_larvae_died / n_start,
    pct_survived_larva  = 100 - pct_larvae_died,
    pct_pupated         = 100 * n_pupated / n_start,
    pct_emerged         = 100 * n_emerged / n_start,
    .groups = "drop"
  )


# ---- 6. Data-quality report ---------------------------------------------------
# "Explicit flagging rather than silent correction": every row where a
# conversion above produced an unexpected result is captured here so you
# can go back to the raw Excel file and check it.

qc_summary <- df %>%
  summarise(
    n_total = n(),
    n_wt_day0_flag             = sum(wt_day0_flag, na.rm = TRUE),
    n_wt_day3_flag             = sum(wt_day3_flag, na.rm = TRUE),
    n_wt_day6_flag             = sum(wt_day6_flag, na.rm = TRUE),
    n_days_to_pupation_flag    = sum(days_to_pupation_flag, na.rm = TRUE),
    n_mortality_day_flag       = sum(mortality_day_flag, na.rm = TRUE),
    n_day_emergence_flag       = sum(day_emergence_flag, na.rm = TRUE),
    n_adult_mortality_dai_flag = sum(adult_mortality_dai_flag, na.rm = TRUE),
    n_egg_24h_flag             = sum(egg_24h_flag, na.rm = TRUE),
    n_egg_48h_flag             = sum(egg_48h_flag, na.rm = TRUE),
    n_egg_72h_flag             = sum(egg_72h_flag, na.rm = TRUE)
  )

flagged_rows <- df %>%
  filter(if_any(ends_with("_flag"), isTRUE)) %>%
  select(insect, treatment, everything())

write_xlsx(
  list(qc_summary = qc_summary, flagged_rows = flagged_rows),
  file.path(table_dir, "00_data_quality_flags.xlsx")
)


# ---- 7. Summary tables --------------------------------------------------------

weight_summary <- df %>%
  select(treatment, wt_day0, wt_day3, wt_day6) %>%
  pivot_longer(-treatment, names_to = "day", values_to = "weight") %>%
  group_by(treatment, day) %>%
  summarise(n = sum(!is.na(weight)),
            mean = mean(weight, na.rm = TRUE),
            sd = sd(weight, na.rm = TRUE),
            se = sd / sqrt(n),
            .groups = "drop")

mortality_summary <- df %>%
  group_by(treatment) %>%
  summarise(pct_died = 100 * mean(larvae_died, na.rm = TRUE))

pupation_summary <- df %>%
  group_by(treatment) %>%
  summarise(pct_pupated = 100 * mean(pupated_ever, na.rm = TRUE))

days_to_pupation_summary <- df %>%
  filter(!is.na(days_to_pupation), days_to_pupation > 0) %>%
  group_by(treatment) %>%
  summarise(n = n(), mean = mean(days_to_pupation), sd = sd(days_to_pupation))

emergence_summary <- df %>%
  group_by(treatment) %>%
  summarise(pct_emerged = 100 * mean(adult_emergence, na.rm = TRUE))

day_emergence_summary <- df %>%
  filter(!is.na(day_emergence)) %>%
  group_by(treatment) %>%
  summarise(n = n(), mean = mean(day_emergence), sd = sd(day_emergence))

sex_ratio_table <- df %>%
  filter(!is.na(sex)) %>%
  count(treatment, sex)

# cumulative larval mortality by day, per treatment (used by the trend plot below)
n_start_tbl <- df %>% count(treatment, name = "n_start")

mortality_trend <- df %>%
  filter(larvae_died, !is.na(mortality_day)) %>%
  count(treatment, mortality_day, name = "n_died") %>%
  complete(treatment, mortality_day = full_seq(mortality_day, 1), fill = list(n_died = 0)) %>%
  arrange(treatment, mortality_day) %>%
  group_by(treatment) %>%
  mutate(cum_died = cumsum(n_died)) %>%
  ungroup() %>%
  left_join(n_start_tbl, by = "treatment") %>%
  mutate(pct_cum_mortality = 100 * cum_died / n_start)

write_xlsx(
  list(
    weight_summary            = weight_summary,
    mortality_summary         = mortality_summary,
    mortality_trend           = mortality_trend,
    pupation_summary          = pupation_summary,
    days_to_pupation_summary  = days_to_pupation_summary,
    emergence_summary         = emergence_summary,
    day_emergence_summary     = day_emergence_summary,
    sex_ratio                 = sex_ratio_table,
    survival_funnel           = survival_funnel
  ),
  file.path(table_dir, "01_summary_tables.xlsx")
)


# ---- 8. Plots -----------------------------------------------------------------

## -- Larval growth --

p01 <- ggplot(df, aes(treatment, wt_day0, fill = treatment)) +
  geom_boxplot() + geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = treat_colors) +
  labs(title = "Larval weight at day 0 (pre-injection)", x = NULL, y = "Weight (g)")

p02 <- ggplot(df, aes(treatment, wt_day3, fill = treatment)) +
  geom_boxplot() + geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = treat_colors) +
  labs(title = "Larval weight at day 3", x = NULL, y = "Weight (g)")

p03 <- ggplot(df, aes(treatment, wt_day6, fill = treatment)) +
  geom_boxplot() + geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_manual(values = treat_colors) +
  labs(title = "Larval weight at day 6", x = NULL, y = "Weight (g)")

weight_long <- df %>%
  select(insect, treatment, wt_day0, wt_day3, wt_day6) %>%
  pivot_longer(starts_with("wt_day"), names_to = "day", values_to = "weight") %>%
  mutate(day = parse_number(day))

weight_trend <- weight_long %>%
  group_by(treatment, day) %>%
  summarise(mean_wt = mean(weight, na.rm = TRUE),
            se_wt = sd(weight, na.rm = TRUE) / sqrt(sum(!is.na(weight))),
            .groups = "drop")

p04 <- ggplot(weight_trend, aes(day, mean_wt, color = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_wt - se_wt, ymax = mean_wt + se_wt), width = 0.3) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Mean larval weight over time (\u00B1 SE)", x = "Day post-injection", y = "Weight (g)")

p05 <- ggplot(weight_long, aes(day, weight, group = insect)) +
  geom_line(alpha = 0.25, color = "grey40") +
  facet_wrap(~treatment) +
  labs(title = "Individual larval growth trajectories", x = "Day post-injection", y = "Weight (g)")

## -- Larval mortality --

p06 <- ggplot(mortality_summary, aes(treatment, pct_died, fill = treatment)) +
  geom_col(position = "dodge") + scale_fill_manual(values = treat_colors) +
  labs(title = "% larval mortality by treatment", x = NULL, y = "% died")

p07 <- df %>%
  filter(larvae_died) %>%
  ggplot(aes(treatment, mortality_day, fill = treatment)) +
  geom_boxplot() + scale_fill_manual(values = treat_colors) +
  labs(title = "Day of larval death post-injection", x = NULL, y = "Day of death")

p08 <- df %>%
  filter(larvae_died) %>%
  ggplot(aes(mortality_day, fill = treatment)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  scale_fill_manual(values = treat_colors) +
  labs(title = "Distribution of day of larval death", x = "Day post-injection", y = "Count")

# NEW: cumulative larval mortality trend over time (line per treatment)
p08b <- ggplot(mortality_trend, aes(mortality_day, pct_cum_mortality, color = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Cumulative larval mortality over time",
       x = "Day post-injection", y = "% cumulative mortality (of starting n)")

## -- Pupation --

p09 <- ggplot(pupation_summary, aes(treatment, pct_pupated, fill = treatment)) +
  geom_col(position = "dodge") + scale_fill_manual(values = treat_colors) +
  labs(title = "% larvae reaching pupation", x = NULL, y = "% pupated")

pupation_days_df <- df %>%
  select(insect, treatment, pup_day8, pup_day9, pup_day11, pup_day12, pup_day13) %>%
  pivot_longer(starts_with("pup_day"), names_to = "day", values_to = "pupated") %>%
  mutate(day = parse_number(day)) %>%
  group_by(treatment, day) %>%
  summarise(pct_pupated = 100 * mean(pupated, na.rm = TRUE), .groups = "drop")

p10 <- ggplot(pupation_days_df, aes(day, pct_pupated, color = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Cumulative pupation over time", x = "Day post-injection", y = "% pupated")

p11 <- df %>%
  filter(!is.na(days_to_pupation), days_to_pupation > 0) %>%
  ggplot(aes(treatment, days_to_pupation, fill = treatment)) +
  geom_boxplot() + scale_fill_manual(values = treat_colors) +
  labs(title = "Days to pupation", x = NULL, y = "Days")

p12 <- df %>%
  filter(!is.na(days_to_pupation), days_to_pupation > 0) %>%
  ggplot(aes(days_to_pupation, fill = treatment)) +
  geom_density(alpha = 0.4) + scale_fill_manual(values = treat_colors) +
  labs(title = "Distribution of days to pupation", x = "Days", y = "Density")

## -- Adult emergence --

p13 <- ggplot(emergence_summary, aes(treatment, pct_emerged, fill = treatment)) +
  geom_col(position = "dodge") + scale_fill_manual(values = treat_colors) +
  labs(title = "% adult emergence", x = NULL, y = "% emerged")

p14 <- df %>%
  filter(!is.na(day_emergence)) %>%
  ggplot(aes(treatment, day_emergence, fill = treatment)) +
  geom_boxplot() + scale_fill_manual(values = treat_colors) +
  labs(title = "Day of adult emergence", x = NULL, y = "Day post-injection")

p15 <- df %>%
  filter(!is.na(day_emergence)) %>%
  ggplot(aes(day_emergence, fill = treatment)) +
  geom_density(alpha = 0.4) + scale_fill_manual(values = treat_colors) +
  labs(title = "Distribution of day of adult emergence", x = "Day", y = "Density")

## -- Survival funnel --

funnel_long <- survival_funnel %>%
  transmute(
    treatment,
    `Injected`               = 100,
    `Survived larval stage`  = pct_survived_larva,
    `Pupated`                = pct_pupated,
    `Emerged as adult`       = pct_emerged
  ) %>%
  pivot_longer(-treatment, names_to = "stage", values_to = "pct") %>%
  mutate(stage = factor(stage, levels = c("Injected", "Survived larval stage",
                                          "Pupated", "Emerged as adult")))

p16 <- ggplot(funnel_long, aes(stage, pct, color = treatment, group = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Survival funnel across life stages", x = NULL, y = "% of insects") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

## -- Sex ratio (dodged) --

p17 <- df %>%
  filter(!is.na(sex)) %>%
  count(treatment, sex) %>%
  group_by(treatment) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(treatment, pct, fill = sex)) +
  geom_col(position = "dodge") +
  labs(title = "Sex ratio of emerged adults", x = NULL, y = "% of emerged adults")

## -- Mating --

p18 <- df %>%
  filter(!is.na(mating_status)) %>%
  group_by(treatment) %>%
  summarise(pct_mated = 100 * mean(mating_status, na.rm = TRUE)) %>%
  ggplot(aes(treatment, pct_mated, fill = treatment)) +
  geom_col(position = "dodge") + scale_fill_manual(values = treat_colors) +
  labs(title = "% mating success", x = NULL, y = "% mated")

# dodged (was position = "fill"): % within each sex x treatment that mated
p19 <- df %>%
  filter(!is.na(mating_status), !is.na(sex)) %>%
  count(treatment, sex, mating_status) %>%
  group_by(treatment, sex) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(sex, pct, fill = mating_status)) +
  geom_col(position = "dodge") +
  facet_wrap(~treatment) +
  labs(title = "Mating status by sex and treatment", y = "% within sex", x = NULL, fill = "Mated")

## -- Egg laying --

egg_long <- df %>%
  select(insect, treatment, egg_24h_count, egg_48h_count, egg_72h_count) %>%
  pivot_longer(-c(insect, treatment), names_to = "timepoint", values_to = "count") %>%
  mutate(timepoint = parse_number(timepoint)) %>%
  group_by(treatment, timepoint) %>%
  summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

p20 <- ggplot(egg_long, aes(timepoint, mean_count, color = treatment)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Mean number of egg masses over time", x = "Hours post-mating", y = "Mean egg mass count")

# dodged (was position = "stack")
p21 <- df %>%
  filter(!is.na(egg_72h_size)) %>%
  count(treatment, egg_72h_size) %>%
  group_by(treatment) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(treatment, pct, fill = egg_72h_size)) +
  geom_col(position = "dodge") +
  labs(title = "Egg mass size distribution at 72h", x = NULL, y = "% of egg masses", fill = "Size")

p22 <- df %>%
  filter(sex == "female", !is.na(mating_status)) %>%
  mutate(laid_eggs = egg_72h_count > 0) %>%
  group_by(treatment) %>%
  summarise(pct_laid = 100 * mean(laid_eggs, na.rm = TRUE)) %>%
  ggplot(aes(treatment, pct_laid, fill = treatment)) +
  geom_col(position = "dodge") + scale_fill_manual(values = treat_colors) +
  labs(title = "% of females that laid eggs by 72h", x = NULL, y = "% laid eggs")

## -- Adult mortality --

p23 <- df %>%
  filter(!is.na(adult_mortality_dai)) %>%
  ggplot(aes(treatment, adult_mortality_dai, fill = treatment)) +
  geom_boxplot() + scale_fill_manual(values = treat_colors) +
  labs(title = "Adult mortality (days after injection)", x = NULL, y = "Days after injection")

## -- Block / family check --

p24 <- ggplot(df, aes(factor(insect_family), wt_day0, fill = treatment)) +
  geom_boxplot() + scale_fill_manual(values = treat_colors) +
  labs(title = "Day-0 larval weight by insect family (block check)",
       x = "Insect family", y = "Weight (g)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## -- Relationships between traits --

p25 <- df %>%
  filter(!is.na(days_to_pupation), days_to_pupation > 0) %>%
  ggplot(aes(wt_day0, days_to_pupation, color = treatment)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Day-0 weight vs. days to pupation", x = "Day-0 weight (g)", y = "Days to pupation")

p26 <- df %>%
  filter(!is.na(day_emergence)) %>%
  ggplot(aes(wt_day0, day_emergence, color = treatment)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = treat_colors) +
  labs(title = "Day-0 weight vs. day of adult emergence", x = "Day-0 weight (g)", y = "Day of emergence")


# ---- 9. Save all plots to Figures/ --------------------------------------------

plots <- list(
  p01_weight_day0 = p01, p02_weight_day3 = p02, p03_weight_day6 = p03,
  p04_weight_trajectory_mean = p04, p05_weight_trajectory_individual = p05,
  p06_pct_larval_mortality = p06, p07_day_of_death_box = p07, p08_day_of_death_hist = p08,
  p08b_larval_mortality_trend = p08b,
  p09_pct_pupated = p09, p10_cumulative_pupation = p10,
  p11_days_to_pupation_box = p11, p12_days_to_pupation_density = p12,
  p13_pct_emergence = p13, p14_day_emergence_box = p14, p15_day_emergence_density = p15,
  p16_survival_funnel = p16, p17_sex_ratio = p17,
  p18_pct_mating = p18, p19_mating_by_sex = p19,
  p20_egg_count_trend = p20, p21_egg_size_distribution = p21, p22_pct_females_laid_eggs = p22,
  p23_adult_mortality_dai = p23, p24_weight_by_insect_family = p24,
  p25_weight_vs_pupation = p25, p26_weight_vs_emergence = p26
)

walk2(plots, names(plots), ~ ggsave(
  filename = file.path(fig_dir, paste0(.y, ".png")),
  plot = .x, width = 7, height = 5, dpi = 300
))

cat("Done.\n",
    "-", nrow(flagged_rows), "rows flagged for data-quality review -> Tables/00_data_quality_flags.xlsx\n",
    "- Summary tables written to Tables/01_summary_tables.xlsx\n",
    "-", length(plots), "plots written to", fig_dir, "\n")