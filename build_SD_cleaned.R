library(tidycensus)
library(dplyr)
library(sf)

state_code <- "SD"
year_use <- 2021
z_90 <- 1.645
eps <- 1e-6

total_pop_df <-
  get_acs(
    geography = "tract",
    variables = "B01003_001",
    state = state_code,
    geometry = FALSE,
    year = year_use
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    tract_pop = estimate,
    tract_pop_moe = moe
  ) %>%
  filter(
    !is.na(tract_pop),
    tract_pop > 0
  )

income21.ed <-
  get_acs(
    geography = "tract",
    variables = "B19013_001",
    state = state_code,
    geometry = FALSE,
    year = year_use
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    NAME,
    median_income = estimate,
    income_moe = moe
  ) %>%
  filter(
    !is.na(median_income),
    !is.na(income_moe),
    median_income > 0
  ) %>%
  mutate(
    var_income = (income_moe / z_90)^2,
    lg_income = log(median_income),
    lg_var = var_income / (median_income^2)
  ) %>%
  filter(
    is.finite(lg_income),
    is.finite(lg_var),
    lg_var > 0
  )

# povrate21 <-
#   get_acs(
#     geography = "tract",
#     variables = "B17001_002",
#     summary_var = "B17001_001",
#     state = state_code,
#     geometry = TRUE,
#     year = year_use
#   ) %>%
#   arrange(GEOID) %>%
#   transmute(
#     GEOID,
#     NAME,
#     pov_universe = summary_est,
#     pov_universe_moe = summary_moe,
#     pov_count = estimate,
#     pov_moe = moe,
#     geometry
#   ) %>%
#   filter(
#     !is.na(pov_universe),
#     !is.na(pov_universe_moe),
#     !is.na(pov_count),
#     !is.na(pov_moe),
#     pov_universe > 0
#   ) %>%
#   mutate(
#     pov_rate_raw = pmin(pmax(pov_count / pov_universe, eps), 1 - eps),
#     pov_rate = pmin(pmax((pov_count + 0.5) / (pov_universe + 1), eps), 1 - eps),
#     pov_rate_moe = moe_prop(pov_count, pov_universe, pov_moe, pov_universe_moe),
#     var_pov = (pov_rate_moe / z_90)^2,
#     m = pov_rate * (1 - pov_rate) / pmax(var_pov, 1e-12)
#   ) %>%
#   filter(
#     is.finite(pov_rate),
#     is.finite(var_pov),
#     is.finite(m),
#     var_pov > 0,
#     m > 0
#   )

povrate21 <-
  get_acs(
    geography = "tract",
    variables = c(
      pov_rate_pct = "S1701_C03_001",
      pov_count = "S1701_C02_001"
    ),
    state = state_code,
    geometry = TRUE,
    year = year_use,
    survey = "acs5"
  ) %>%
  select(GEOID, NAME, variable, estimate, moe, geometry) %>%
  tidyr::pivot_wider(
    names_from = variable,
    values_from = c(estimate, moe)
  ) %>%
  transmute(
    GEOID,
    NAME,
    pov_count = estimate_pov_count,
    pov_count_moe = moe_pov_count,
    pov_rate_raw = pmin(pmax(estimate_pov_rate_pct / 100, eps), 1 - eps),
    pov_rate = pmin(pmax(estimate_pov_rate_pct / 100, eps), 1 - eps),
    pov_rate_moe = moe_pov_rate_pct / 100,
    var_pov = (pov_rate_moe / z_90)^2,
    m = pov_rate * (1 - pov_rate) / pmax(var_pov, 1e-12),
    geometry
  ) %>%
  filter(
    !is.na(pov_rate),
    !is.na(pov_rate_moe),
    !is.na(var_pov),
    !is.na(m),
    is.finite(pov_rate),
    is.finite(var_pov),
    is.finite(m),
    var_pov > 0,
    m > 0
  ) %>%
  arrange(GEOID)

# white_df <-
#   get_acs(
#     geography = "tract",
#     variables = "B01001A_001",
#     state = state_code,
#     geometry = FALSE,
#     year = year_use
#   ) %>%
#   arrange(GEOID) %>%
#   transmute(
#     GEOID,
#     white = estimate
#   ) %>%
#   filter(!is.na(white))

white_df <-
  get_acs(
    geography = "tract",
    variables = "DP05_0037PE",
    state = state_code,
    geometry = FALSE,
    year = year_use,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    X_white = estimate / 100
  ) %>%
  filter(!is.na(X_white))

aligned_sf <-
  povrate21 %>%
  inner_join(total_pop_df, by = "GEOID") %>%
  inner_join(
    income21.ed %>%
      select(
        GEOID,
        NAME_income = NAME,
        median_income,
        income_moe,
        var_income,
        lg_income,
        lg_var
      ),
    by = "GEOID"
  ) %>%
  inner_join(white_df, by = "GEOID") %>%
  mutate(
    NAME = coalesce(NAME, NAME_income),
    X_white = X_white
  ) %>%
  filter(
    !is.na(NAME),
    !is.na(X_white),
    is.finite(X_white),
    X_white >= 0,
    tract_pop > 0
  ) %>%
  arrange(GEOID)

income21.ed <-
  aligned_sf %>%
  st_drop_geometry() %>%
  transmute(
    GEOID,
    NAME,
    population = tract_pop,
    median_income,
    income_moe,
    var_income,
    lg_income,
    lg_var
  )

# povrate21 <-
#   aligned_sf %>%
#   transmute(
#     GEOID,
#     NAME,
#     population = tract_pop,
#     pov_universe,
#     pov_universe_moe,
#     pov_count,
#     pov_moe,
#     pov_rate,
#     pov_rate_moe,
#     var_pov,
#     m,
#     geometry
#   )
povrate21 <-
  aligned_sf %>%
  transmute(
    GEOID,
    NAME,
    population = tract_pop,
    pov_count,
    pov_count_moe,
    pov_rate,
    pov_rate_moe,
    var_pov,
    m,
    geometry
  )

tract_sf <- aligned_sf %>% select(GEOID, NAME, geometry)
X_white <- aligned_sf$X_white

save(
  tract_sf,
  income21.ed,
  povrate21,
  X_white,
  file = "SD cleaned.RData"
)
