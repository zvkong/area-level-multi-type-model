library(tidycensus)
library(dplyr)
library(sf)
library(purrr)
library(spdep)
library(spatialreg)
states_use <- c("IA", "KS", "MN", "MO", "NE", "ND", "SD")
year_use <- 2021
z_90 <- qnorm(0.95)
eps <- 1e-06

fetch_acs_states <- function(variables, summary_var = NULL, geometry = FALSE, survey = "acs5") {
  map_dfr(
    states_use,
    ~ get_acs(
      geography = "tract",
      variables = variables,
      summary_var = summary_var,
      state = .x,
      geometry = geometry,
      year = year_use,
      survey = survey
    )
  )
}

income21.ed <-
  fetch_acs_states(
    variables = "B19013_001",
    summary_var = "B01003_001",
    geometry = FALSE
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
  fetch_acs_states(
    variables = "B17001_002",
    summary_var = "B01003_001",
    geometry = TRUE
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
  fetch_acs_states(
    variables = "B01001A_001",
    geometry = FALSE
  ) %>%
  arrange(GEOID) %>%
  transmute(GEOID, white = estimate)

bachelor_df <-
  fetch_acs_states(
    variables = "DP02_0068PE",
    geometry = FALSE,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  transmute(GEOID, X_bachelor = estimate / 100)

snap_df <-
  fetch_acs_states(
    variables = "DP03_0074PE",
    geometry = FALSE,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  transmute(GEOID, X_snap = estimate / 100)

aligned_sf <-
  povrate21 %>%
  inner_join(income21.ed, by = "GEOID", suffix = c("", "_inc")) %>%
  inner_join(white_df, by = "GEOID") %>%
  inner_join(bachelor_df, by = "GEOID") %>%
  inner_join(snap_df, by = "GEOID") %>%
  mutate(
    NAME = coalesce(NAME, NAME_inc),
    X_white = white / population
  ) %>%
  filter(
    !is.na(X_white),
    !is.na(X_bachelor),
    !is.na(X_snap),
    is.finite(X_white),
    is.finite(X_bachelor),
    is.finite(X_snap)
  ) %>%
  arrange(GEOID)

povrate21 <-
  aligned_sf %>%
  select(GEOID, NAME, population, pov_rate, m, var_pov, var_pov_count, pov_count, geometry)

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

MO_county_sf <- povrate21 %>% select(GEOID, NAME, geometry)
X_white <- aligned_sf$X_white
X_bachelor <- aligned_sf$X_bachelor
X_snap <- aligned_sf$X_snap
X_1 <- as.matrix(cbind(1, X_bachelor, X_snap, X_white))
X_2 <- as.matrix(cbind(1, X_white, X_snap, X_bachelor))

if (file.exists("functions.r")) {
  source("functions.r")
}
if (!exists("large_abs", mode = "function")) {
  large_abs <- function(x, c) which(abs(x) > quantile(abs(x), probs = c, na.rm = TRUE))
}

MO_nb <- poly2nb(MO_county_sf, queen = TRUE)
W_list <- nb2listw(MO_nb, style = "B", zero.policy = TRUE)
B <- as.matrix(as_dgRMatrix_listw(W_list))
eig <- eigen(B, symmetric = TRUE)
idx <- large_abs(eig$values, 0.25)
S <- eig$vectors[, idx, drop = FALSE]

save(MO_county_sf, income21.ed, povrate21, X_white, X_bachelor, X_snap, X_1, X_2, B, S, file = "westnorthc.Rdata")

