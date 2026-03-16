suppressPackageStartupMessages({
  library(tidycensus)
  library(dplyr)
  library(sf)
})

year <- 2021
state <- "SD"
moe_mult <- qnorm(0.95)

income21.ed <-
  get_acs(
    geography = "tract",
    variables = "B19013_001",
    summary_var = "B01003_001",
    state = state,
    geometry = FALSE,
    year = year,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  rename(population = summary_est, medium_income = estimate, income_moe = moe) %>%
  select(GEOID, NAME, population, medium_income, income_moe) %>%
  filter(population > 0, !is.na(medium_income), !is.na(income_moe), medium_income > 0) %>%
  mutate(
    var_income = (income_moe / moe_mult)^2,
    lg_income = log(medium_income),
    lg_var = var_income / medium_income^2
  )

povrate21 <-
  get_acs(
    geography = "tract",
    variables = "B17001_002",
    summary_var = "B01003_001",
    state = state,
    geometry = TRUE,
    year = year,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  rename(population = summary_est, pov_count = estimate, pov_moe = moe) %>%
  filter(population > 0, !is.na(pov_count), !is.na(pov_moe)) %>%
  mutate(
    pov_rate = pov_count / population,
    var_pov_count = (pov_moe / moe_mult)^2,
    var_pov = (pov_moe / (population * moe_mult))^2,
    m = pov_rate * (1 - pov_rate) / pmax(var_pov, 1e-12),
    ystar = m * pov_rate
  ) %>%
  select(GEOID, NAME, population, pov_rate, m, var_pov, var_pov_count, pov_count, ystar, geometry)

white_df <-
  get_acs(
    geography = "tract",
    variables = "B01001A_001",
    summary_var = "B01003_001",
    state = state,
    geometry = FALSE,
    year = year,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  transmute(GEOID, white = estimate)

aligned <-
  povrate21 %>%
  inner_join(income21.ed, by = "GEOID", suffix = c("", "_inc")) %>%
  inner_join(white_df, by = "GEOID") %>%
  select(-NAME_inc, -population_inc) %>%
  mutate(X_white = white / population) %>%
  filter(
    is.finite(lg_income),
    is.finite(lg_var),
    is.finite(pov_rate),
    is.finite(var_pov),
    is.finite(m),
    m > 0,
    !is.na(X_white)
  ) %>%
  arrange(GEOID)

lg_mean <- mean(aligned$lg_income)
lg_sd <- sd(aligned$lg_income)

aligned <-
  aligned %>%
  mutate(
    income_scaled = (lg_income - lg_mean) / lg_sd,
    lg_var = lg_var / lg_sd^2
  )

MO_county_sf <- aligned %>% select(GEOID, NAME, geometry)

income21.ed <-
  aligned %>%
  st_drop_geometry() %>%
  select(GEOID, NAME, population, medium_income, income_moe, var_income, lg_income, lg_var, income_scaled)

povrate21 <-
  aligned %>%
  select(GEOID, NAME, population, pov_rate, m, var_pov, var_pov_count, pov_count, ystar, geometry)

X_white <- aligned$X_white

save(MO_county_sf, income21.ed, povrate21, X_white, file = "SD cleaned.RData")
