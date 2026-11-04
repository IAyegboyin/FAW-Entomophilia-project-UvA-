# Helper functions -----------------------------------------------------
# Converts a character column to numeric
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

view(sf_data)