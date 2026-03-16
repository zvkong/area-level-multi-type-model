suppressPackageStartupMessages({
  library(tidycensus)
  library(dplyr)
  library(sf)
  library(purrr)
  library(spdep)
  library(Matrix)
})

year <- 2021
states <- c("IA", "KS", "MN", "MO", "NE", "ND", "SD")
moe_mult <- qnorm(0.95)

large_abs <- function(x, threshold = 0.25) {
  which(abs(x) > threshold * max(abs(x), na.rm = TRUE))
}

get_income_state <- function(state) {
  get_acs(
    geography = "tract",
    variables = "B19013_001",
    summary_var = "B01003_001",
    state = state,
    geometry = FALSE,
    year = year,
    survey = "acs5"
  ) %>%
    transmute(
      GEOID,
      NAME,
      population = summary_est,
      medium_income = estimate,
      income_moe = moe
    )
}

get_poverty_state <- function(state) {
  get_acs(
    geography = "tract",
    variables = "B17001_002",
    summary_var = "B01003_001",
    state = state,
    geometry = TRUE,
    year = year,
    survey = "acs5"
  ) %>%
    transmute(
      GEOID,
      NAME,
      population = summary_est,
      pov_count = estimate,
      pov_moe = moe,
      geometry
    )
}

get_white_state <- function(state) {
  get_acs(
    geography = "tract",
    variables = "B01001A_001",
    summary_var = "B01003_001",
    state = state,
    geometry = FALSE,
    year = year,
    survey = "acs5"
  ) %>%
    transmute(GEOID, white = estimate)
}

get_bachelor_state <- function(state) {
  get_acs(
    geography = "tract",
    variables = c(
      total25 = "B15003_001",
      bachelor = "B15003_022",
      master = "B15003_023",
      professional = "B15003_024",
      doctorate = "B15003_025"
    ),
    state = state,
    geometry = FALSE,
    year = year,
    survey = "acs5"
  ) %>%
    select(GEOID, variable, estimate) %>%
    tidyr::pivot_wider(names_from = variable, values_from = estimate) %>%
    mutate(bachelor_plus = bachelor + master + professional + doctorate) %>%
    transmute(GEOID, adults25 = total25, bachelor_plus)
}

get_snap_state <- function(state) {
  get_acs(
    geography = "tract",
    variables = c(total_hh = "B22001_001", snap_hh = "B22001_002"),
    state = state,
    geometry = FALSE,
    year = year,
    survey = "acs5"
  ) %>%
    select(GEOID, variable, estimate) %>%
    tidyr::pivot_wider(names_from = variable, values_from = estimate) %>%
    transmute(GEOID, total_hh, snap_hh)
}

income21.ed <- map_dfr(states, get_income_state) %>%
  arrange(GEOID) %>%
  filter(population > 0, !is.na(medium_income), !is.na(income_moe), medium_income > 0) %>%
  mutate(
    var_income = (income_moe / moe_mult)^2,
    lg_income = log(medium_income),
    lg_var = var_income / medium_income^2
  )

povrate21 <- map_dfr(states, get_poverty_state) %>%
  st_as_sf() %>%
  arrange(GEOID) %>%
  filter(population > 0, !is.na(pov_count), !is.na(pov_moe)) %>%
  mutate(
    pov_rate = pov_count / population,
    var_pov_count = (pov_moe / moe_mult)^2,
    var_pov = (pov_moe / (population * moe_mult))^2,
    m = pov_rate * (1 - pov_rate) / pmax(var_pov, 1e-12),
    ystar = m * pov_rate
  ) %>%
  select(GEOID, NAME, population, pov_rate, m, var_pov, var_pov_count, pov_count, ystar, geometry)

white_df <- map_dfr(states, get_white_state) %>% arrange(GEOID)
bachelor_df <- map_dfr(states, get_bachelor_state) %>% arrange(GEOID)
snap_df <- map_dfr(states, get_snap_state) %>% arrange(GEOID)

aligned <-
  povrate21 %>%
  inner_join(income21.ed, by = "GEOID", suffix = c("", "_inc")) %>%
  inner_join(white_df, by = "GEOID") %>%
  inner_join(bachelor_df, by = "GEOID") %>%
  inner_join(snap_df, by = "GEOID") %>%
  select(-NAME_inc, -population_inc) %>%
  mutate(
    X_white = white / population,
    X_bachelor = ifelse(adults25 > 0, bachelor_plus / adults25, NA_real_),
    X_snap = ifelse(total_hh > 0, snap_hh / total_hh, NA_real_)
  ) %>%
  filter(
    is.finite(lg_income),
    is.finite(lg_var),
    is.finite(pov_rate),
    is.finite(var_pov),
    is.finite(m),
    m > 0,
    is.finite(X_white),
    is.finite(X_bachelor),
    is.finite(X_snap)
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

nb_obj <- poly2nb(MO_county_sf, queen = TRUE)
W_list <- nb2listw(nb_obj, style = "B", zero.policy = TRUE)
B <- as.matrix(as_dgRMatrix_listw(W_list))
eig <- eigen(B, symmetric = TRUE)
index <- large_abs(eig$values, 0.25)
S <- eig$vectors[, index, drop = FALSE]

X_white <- aligned$X_white
X_bachelor <- aligned$X_bachelor
X_snap <- aligned$X_snap
X_1 <- as.matrix(cbind(1, X_bachelor, X_snap, X_white))
X_2 <- as.matrix(cbind(1, X_white, X_snap, X_bachelor))

income21.ed <-
  aligned %>%
  st_drop_geometry() %>%
  select(GEOID, NAME, population, medium_income, income_moe, var_income, lg_income, lg_var, income_scaled)

povrate21 <-
  aligned %>%
  select(GEOID, NAME, population, pov_rate, m, var_pov, var_pov_count, pov_count, ystar, geometry)

save(
  MO_county_sf,
  income21.ed,
  povrate21,
  B,
  X_white,
  X_bachelor,
  X_snap,
  X_1,
  X_2,
  S,
  file = "westnorthc.Rdata"
)
