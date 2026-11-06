# Statistics: Effect of Serratia on Spodoptera frugiperda (FAW)

#packages - newly added to the ones in the load data era
library(survival)
library(lme4)
library(car)
library(multcomp)
library(multcompView)
library(broom)

# data mutation to get more variables and also check for consistency here ----

censor_day <- 13   # last day larvae were checked (day 13 pupation check)

sf_data <- sf_data %>%
  mutate(
    pupated_ever = pup_day13,
    # growth rate in mg/day (weights are in grams -> x1000)
    gr_0_3 = (wt_day3 - wt_day0) / 3 * 1000,
    gr_3_6 = (wt_day6 - wt_day3) / 3 * 1000,
    gr_0_6 = (wt_day6 - wt_day0) / 6 * 1000,
    # pupation days of 0 are flagged as errors in cleaning -> treat as NA here
    days_to_pupation = if_else(days_to_pupation > 0, days_to_pupation, NA_real_)
  )

# Survival dataset, time = day of death; survivors censored at censor_day
surv_data <- sf_data %>%
  filter(!is.na(larvae_died)) %>%
  mutate(time  = if_else(larvae_died, mortality_day, censor_day),
         event = as.integer(larvae_died)) %>%
  filter(!is.na(time))

# Consistency checks -- these print rows to look at in Excel, nothing is deleted
cat("\n--- Consistency checks (rows to double-check in the Excel file) ---\n")
chk <- list(
  "weight 'x' (dead) but larvae_died = no" =
    sf_data %>% filter((died_by_day3 %in% TRUE | died_by_day6 %in% TRUE),
                       larvae_died %in% FALSE) %>% pull(insect),
  "larvae_died = yes but pupated" =
    sf_data %>% filter(larvae_died %in% TRUE, pupated_ever %in% TRUE) %>% pull(insect),
  "adult emerged but never pupated" =
    sf_data %>% filter(adult_emergence %in% TRUE, pupated_ever %in% FALSE) %>% pull(insect),
  "adult emerged but day_emergence missing" =
    sf_data %>% filter(adult_emergence %in% TRUE, is.na(day_emergence)) %>% pull(insect),
  "negative growth day 0-3 (weight dropped)" =
    sf_data %>% filter(gr_0_3 < 0) %>% pull(insect)
)
iwalk(chk, ~ cat(sprintf("%-45s : %s\n", .y,
                         if (length(.x) == 0) "none" else paste(.x, collapse = ", "))))

# helper functions ----
# Compare a CONTINUOUS variable across treatments.
# Returns: descriptives, normality, Levene, ANOVA, Kruskal-Wallis,
#          Tukey contrasts, pairwise Wilcoxon (BH), compact letters (Tukey)
compare_continuous <- function(data, var, label = var) {
  d <- data %>% filter(!is.na(.data[[var]])) %>% mutate(y = .data[[var]])
  
  desc <- d %>%
    group_by(treatment) %>%
    summarise(n = n(), mean = mean(y), sd = sd(y), se = sd / sqrt(n),
              median = median(y), min = min(y), max = max(y), .groups = "drop") %>%
    mutate(variable = label, .before = 1)
  
  normality <- d %>%
    group_by(treatment) %>%
    summarise(shapiro_p = if (n() >= 3 && n() <= 5000 && sd(y) > 0)
      shapiro.test(y)$p.value else NA_real_, .groups = "drop") %>%
    mutate(variable = label, .before = 1)
  
  lev <- car::leveneTest(y ~ treatment, data = d)
  fit <- aov(y ~ treatment, data = d)
  an  <- summary(fit)[[1]]
  kw  <- kruskal.test(y ~ treatment, data = d)
  
  overall <- tibble(
    variable = label,
    levene_p = lev[["Pr(>F)"]][1],
    anova_F  = an[["F value"]][1],
    anova_p  = an[["Pr(>F)"]][1],
    kruskal_chisq = unname(kw$statistic),
    kruskal_p = kw$p.value
  )
  emm    <- emmeans(fit, ~ treatment)
  tukey  <- as.data.frame(pairs(emm, adjust = "tukey")) %>% mutate(variable = label, .before = 1)
  cld_df <- as.data.frame(multcomp::cld(emm, Letters = letters, adjust = "tukey")) %>%
    transmute(treatment, letters = str_replace_all(.group, " ", ""))
  
  wil <- suppressWarnings(pairwise.wilcox.test(d$y, d$treatment, p.adjust.method = "BH"))
  wil <- broom::tidy(wil) %>% mutate(variable = label, .before = 1)
  
  list(desc = desc, normality = normality, overall = overall,
       tukey = tukey, wilcox = wil, letters = cld_df)
}

# Compare a TRUE/FALSE variable across treatments (Fisher exact + pairwise)
compare_proportion <- function(data, var, label = var) {
  d <- data %>%
    filter(!is.na(.data[[var]])) %>%
    mutate(y = factor(.data[[var]], levels = c(FALSE, TRUE)))
  
  desc <- d %>%
    group_by(treatment) %>%
    summarise(n = n(), n_yes = sum(y == "TRUE"), pct = 100 * n_yes / n, .groups = "drop") %>%
    mutate(variable = label, .before = 1)
  
  overall <- tibble(variable = label, test = "Fisher exact (all groups)",
                    p = tryCatch(fisher.test(table(d$treatment, d$y))$p.value,
                                 error = function(e) NA_real_))
  pairs_df <- combn(levels(d$treatment), 2, simplify = FALSE) %>%
    map_dfr(function(p) {
      sub <- d %>% filter(treatment %in% p)
      tb  <- table(droplevels(sub$treatment), sub$y)
      ft  <- tryCatch(fisher.test(tb), error = function(e) NULL)
      tibble(variable = label, comparison = paste(p, collapse = " vs "),
             p = if (is.null(ft)) NA_real_ else ft$p.value,
             odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate))
    }) %>%
    mutate(p_adj_BH = p.adjust(p, method = "BH"))
  
  pv <- pairs_df$p_adj_BH
  pv[is.na(pv)] <- 1
  lets <- multcompView::multcompLetters(
    setNames(pv, str_replace(pairs_df$comparison, " vs ", "-")))$Letters
  letters_df <- tibble(treatment = factor(names(lets), levels = levels(d$treatment)),
                       letters = unname(lets))
  
  list(desc = desc, overall = overall, pairwise = pairs_df, letters = letters_df)
}

# Larval survival (Kaplan-Meier, log-rank, Cox) -----
km_fit <- survfit(Surv(time, event) ~ treatment, data = surv_data)
print(km_fit)

logrank <- survdiff(Surv(time, event) ~ treatment, data = surv_data)
logrank_p <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
cat("\nLog-rank test across treatments: p =", signif(logrank_p, 3), "\n")

# pairwise log-rank (BH-adjusted)
pairwise_logrank <- combn(levels(surv_data$treatment), 2, simplify = FALSE) %>%
  map_dfr(function(p) {
    sub <- surv_data %>% filter(treatment %in% p) %>% droplevels()
    sd  <- survdiff(Surv(time, event) ~ treatment, data = sub)
    tibble(comparison = paste(p, collapse = " vs "),
           chisq = sd$chisq, p = 1 - pchisq(sd$chisq, df = 1))
  }) %>%
  mutate(p_adj_BH = p.adjust(p, method = "BH"))

# Cox proportional hazards: hazard ratio vs PBS
cox_fit <- coxph(Surv(time, event) ~ treatment, data = surv_data)
cox_table <- broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(hazard_ratio = estimate)
ph_check <- cox.zph(cox_fit)       # proportional-hazards assumption check
print(ph_check)

# Does starting weight predict death?  (covariate check)
cox_wt <- coxph(Surv(time, event) ~ treatment + wt_day0, data = surv_data)
cox_wt_table <- broom::tidy(cox_wt, exponentiate = TRUE, conf.int = TRUE) %>%
  rename(hazard_ratio = estimate)

# median survival + survival at final check
km_summary <- summary(km_fit, times = c(3, 6, 9, 13))$table
km_surv_at <- tibble(
  treatment = str_remove(summary(km_fit, times = c(3, 6, 9, 13))$strata, "treatment="),
  day       = summary(km_fit, times = c(3, 6, 9, 13))$time,
  survival  = summary(km_fit, times = c(3, 6, 9, 13))$surv,
  lower     = summary(km_fit, times = c(3, 6, 9, 13))$lower,
  upper     = summary(km_fit, times = c(3, 6, 9, 13))$upper
)
# Proportions: larval death, pupation, adult emergence ----
prop_res <- list(
  larvae_died     = compare_proportion(sf_data, "larvae_died",     "Larval mortality"),
  pupated_ever    = compare_proportion(sf_data, "pupated_ever",    "Reached pupation"),
  adult_emergence = compare_proportion(sf_data, "adult_emergence", "Adult emergence (survival to adulthood)")
)
# Emergence among insects that pupated (pupal/adult-stage failure)
emerg_given_pup <- sf_data %>% filter(pupated_ever %in% TRUE)
prop_res$emerg_given_pup <- compare_proportion(emerg_given_pup, "adult_emergence",
                                               "Adult emergence | pupated")
# Growth -----
cont_res <- list(
  wt_day0 = compare_continuous(sf_data, "wt_day0", "Weight day 0 (g)"),
  wt_day3 = compare_continuous(sf_data, "wt_day3", "Weight day 3 (g)"),
  wt_day6 = compare_continuous(sf_data, "wt_day6", "Weight day 6 (g)"),
  gr_0_3  = compare_continuous(sf_data, "gr_0_3",  "Growth rate day 0-3 (mg/day)"),
  gr_3_6  = compare_continuous(sf_data, "gr_3_6",  "Growth rate day 3-6 (mg/day)"),
  gr_0_6  = compare_continuous(sf_data, "gr_0_6",  "Growth rate day 0-6 (mg/day)")
)
# Repeated-measures mixed model on log weight (treatment by day)
# random intercepts: individual insect nested in family
long_w <- sf_data %>%
  select(insect, insect_family, treatment, wt_day0, wt_day3, wt_day6) %>%
  pivot_longer(starts_with("wt_day"), names_to = "day", values_to = "wt") %>%
  mutate(day = factor(parse_number(day)), log_wt = log(wt)) %>%
  filter(!is.na(wt), wt > 0)

m_wt <- tryCatch(
  lme4::lmer(log_wt ~ treatment * day + (1 | insect_family) + (1 | insect), data = long_w),
  error = function(e) NULL)

if (!is.null(m_wt)) {
  m_wt_anova <- as.data.frame(anova(m_wt)) %>% tibble::rownames_to_column("term")
  m_wt_contrasts <- as.data.frame(pairs(emmeans(m_wt, ~ treatment | day), adjust = "tukey"))
  print(m_wt_anova)
}

# Insect development time ----
cont_res$days_to_pupation <- compare_continuous(sf_data, "days_to_pupation", "Days to pupation")
cont_res$day_emergence    <- compare_continuous(sf_data, "day_emergence",    "Day of adult emergence")

# Sex ratio ----
sex_tab <- sf_data %>% filter(!is.na(sex)) %>% 
  count(treatment, sex) %>%
  group_by(treatment) %>% 
  mutate(pct = 100 * n / sum(n)) %>% 
  ungroup()

sex_test_overall <- fisher.test(table(sf_data$treatment[!is.na(sf_data$sex)],
                                      sf_data$sex[!is.na(sf_data$sex)]))

# Is each treatment different from a 1:1 sex ratio?
sex_binom <- sex_tab %>%
  group_by(treatment) %>%
  summarise(n_total = sum(n),
            n_female = sum(n[sex == "female"]),
            pct_female = 100 * n_female / n_total,
            binom_p_vs_50pct = binom.test(n_female, n_total, p = 0.5)$p.value,
            .groups = "drop")

cat("\nSex ratio ~ treatment (Fisher): p =", signif(sex_test_overall$p.value, 3), "\n")