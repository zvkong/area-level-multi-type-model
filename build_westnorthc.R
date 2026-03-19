source("packages.r")
source("functions.r")

states_use <- c("IA", "KS", "MN", "MO", "NE", "ND", "SD")
year_use <- 2021
z_90 <- 1.645
eps <- 1e-6
basis_q <- 0.25

tract_population_df <-
  fetch_acs_region(
    variables = "B01003_001",
    geometry = FALSE
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    tract_population = estimate,
    tract_population_moe = moe
  ) %>%
  filter(
    !is.na(tract_population),
    !is.na(tract_population_moe),
    tract_population > 0
  )

income_df <-
  fetch_acs_region(
    variables = "B19013_001",
    geometry = FALSE
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
    income_variance = (income_moe / z_90)^2,
    log_income = log(median_income),
    log_income_variance = income_variance / (median_income^2)
  ) %>%
  filter(
    is.finite(log_income),
    is.finite(log_income_variance),
    log_income_variance > 0
  )

poverty_df <-
  fetch_acs_region(
    variables = "B17001_002",
    summary_var = "B17001_001",
    geometry = TRUE
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    NAME,
    poverty_universe = summary_est,
    poverty_universe_moe = summary_moe,
    poverty_count = estimate,
    poverty_count_moe = moe,
    geometry
  ) %>%
  filter(
    !is.na(poverty_universe),
    !is.na(poverty_universe_moe),
    !is.na(poverty_count),
    !is.na(poverty_count_moe),
    poverty_universe > 0
  ) %>%
  mutate(
    poverty_rate = pmin(pmax(poverty_count / poverty_universe, eps), 1 - eps),
    poverty_rate_moe = tidycensus::moe_prop(
      poverty_count,
      poverty_universe,
      poverty_count_moe,
      poverty_universe_moe
    ),
    poverty_count_variance = (poverty_count_moe / z_90)^2,
    poverty_rate_variance = (poverty_rate_moe / z_90)^2,
    effective_n = poverty_rate * (1 - poverty_rate) / pmax(poverty_rate_variance, 1e-12)
  ) %>%
  filter(
    is.finite(poverty_rate),
    is.finite(poverty_rate_variance),
    is.finite(effective_n),
    poverty_rate_variance > 0,
    effective_n > 0
  )

white_df <-
  fetch_acs_region(
    variables = "B01001A_001",
    geometry = FALSE
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    white_population = estimate
  ) %>%
  filter(!is.na(white_population))

bachelor_df <-
  fetch_acs_region(
    variables = "DP02_0068PE",
    geometry = FALSE,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    prop_bachelor = estimate / 100
  ) %>%
  filter(!is.na(prop_bachelor))

snap_df <-
  fetch_acs_region(
    variables = "DP03_0074PE",
    geometry = FALSE,
    survey = "acs5"
  ) %>%
  arrange(GEOID) %>%
  transmute(
    GEOID,
    prop_snap = estimate / 100
  ) %>%
  filter(!is.na(prop_snap))

region_tract_sf <-
  poverty_df %>%
  inner_join(tract_population_df, by = "GEOID") %>%
  inner_join(
    income_df %>%
      dplyr::select(
        GEOID,
        NAME_income = NAME,
        median_income,
        income_moe,
        income_variance,
        log_income,
        log_income_variance
      ),
    by = "GEOID"
  ) %>%
  inner_join(white_df, by = "GEOID") %>%
  inner_join(bachelor_df, by = "GEOID") %>%
  inner_join(snap_df, by = "GEOID") %>%
  mutate(
    NAME = dplyr::coalesce(NAME, NAME_income),
    prop_white = white_population / tract_population
  ) %>%
  filter(
    !is.na(NAME),
    !is.na(prop_white),
    !is.na(prop_bachelor),
    !is.na(prop_snap),
    is.finite(prop_white),
    is.finite(prop_bachelor),
    is.finite(prop_snap),
    prop_white >= 0,
    tract_population > 0
  ) %>%
  arrange(GEOID)

income_df <-
  region_tract_sf %>%
  sf::st_drop_geometry() %>%
  transmute(
    GEOID,
    NAME,
    tract_population,
    median_income,
    income_moe,
    income_variance,
    log_income,
    log_income_variance,
    log_income_scaled = as.numeric(scale(log_income))
  )

poverty_df <-
  region_tract_sf %>%
  transmute(
    GEOID,
    NAME,
    tract_population,
    poverty_universe,
    poverty_universe_moe,
    poverty_count,
    poverty_count_moe,
    poverty_count_variance,
    poverty_rate,
    poverty_rate_moe,
    poverty_rate_variance,
    effective_n,
    geometry
  )

x_white <- region_tract_sf$prop_white
x_bachelor <- region_tract_sf$prop_bachelor
x_snap <- region_tract_sf$prop_snap

X_income <- as.matrix(cbind(1, x_bachelor, x_snap, x_white))
X_poverty <- as.matrix(cbind(1, x_white, x_snap, x_bachelor))

region_nb <- spdep::poly2nb(region_tract_sf, queen = TRUE)
region_weights <- spdep::nb2listw(region_nb, style = "B", zero.policy = TRUE)
adjacency_matrix <- as.matrix(spatialreg::as_dgRMatrix_listw(region_weights))

eigen_decomp <- eigen(adjacency_matrix, symmetric = TRUE)
basis_index <- large_abs(eigen_decomp$values, basis_q)
basis_matrix <- eigen_decomp$vectors[, basis_index, drop = FALSE]

save(
  region_tract_sf,
  income_df,
  poverty_df,
  x_white,
  x_bachelor,
  x_snap,
  X_income,
  X_poverty,
  adjacency_matrix,
  basis_matrix,
  basis_q,
  z_90,
  file = "westnorth_region_clean.RData"
)
