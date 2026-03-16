library(tidycensus)
library(dplyr)
library(sf)

state_code <- "SD"
year_use <- 2021
z_90 <- qnorm(0.95)
eps <- 1e-06

income21.ed <-
  get_acs(
    geography = "tract",
    variables = "B19013_001",
    summary_var = "B01003_001",
    state = state_code,
    geometry = FALSE,
    year = year_use
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    NAME,
    population = summary_est,
    medium_income = estimate,
    moe,
    var_income = (moe / z_90)^2,
    lg_income = log(medium_income),
    lg_var = (moe / z_90)^2 / (medium_income^2)
  ) %>%
  filter(
    population > 0,
    !is.na(medium_income),
    !is.na(moe),
    medium_income > 0,
    is.finite(lg_income),
    is.finite(lg_var)
  )

povrate21 <-
  get_acs(
    geography = "tract",
    variables = "B17001_002",
    summary_var = "B01003_001",
    state = state_code,
    geometry = TRUE,
    year = year_use
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    NAME,
    population = summary_est,
    pov_count = estimate,
    pov_moe = moe,
    geometry
  ) %>%
  filter(population > 0, !is.na(pov_count), !is.na(pov_moe)) %>%
  mutate(
    var_pov_count = (pov_moe / z_90)^2,
    var_pov = (pov_moe / (population * z_90))^2,
    pov_rate = pmin(pmax(pov_count / population, eps), 1 - eps),
    m = pov_rate * (1 - pov_rate) / pmax(var_pov, 1e-12)
  ) %>%
  filter(is.finite(pov_rate), is.finite(var_pov), is.finite(m), m > 0)

white_df <-
  get_acs(
    geography = "tract",
    variables = "B01001A_001",
    state = state_code,
    geometry = FALSE,
    year = year_use
  ) %>%
  arrange(GEOID) %>%
  transmute(GEOID, white = estimate)

aligned_sf <-
  povrate21 %>%
  inner_join(income21.ed, by = "GEOID", suffix = c("", "_inc")) %>%
  inner_join(white_df, by = "GEOID") %>%
  mutate(
    NAME = coalesce(NAME, NAME_inc),
    X_white = white / population
  ) %>%
  filter(
    !is.na(X_white),
    is.finite(X_white)
  ) %>%
  arrange(GEOID)

income21.ed <-
  aligned_sf %>%
  st_drop_geometry() %>%
  transmute(
    GEOID,
    NAME,
    population,
    medium_income,
    moe,
    var_income,
    lg_income,
    lg_var,
    income_scaled = as.numeric(scale(lg_income))
  )

povrate21 <-
  aligned_sf %>%
  select(GEOID, NAME, population, pov_rate, m, var_pov, var_pov_count, pov_count, geometry)

MO_county_sf <- povrate21 %>% select(GEOID, NAME, geometry)
X_white <- aligned_sf$X_white

save(MO_county_sf, income21.ed, povrate21, X_white, file = "SD cleaned.RData")
