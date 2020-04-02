setwd("~/Documents/lshtm/github repos/CFR_calculation/global_estimates/")

zmeanHDT <- 13
zsdHDT <- 12.7
zmedianHDT <- 9.1
muHDT <- log(zmedianHDT)
sigmaHDT <- sqrt(2*(log(zmeanHDT) - muHDT))
cCFRBaseline <- 1.38
cCFREstimateRange <- c(1.23, 1.53)
#cCFRIQRRange <- c(1.3, 1.4)

# Hospitalisation to death distribution
hospitalisation_to_death_truncated <- function(x) {
  plnorm(x + 1, muHDT, sigmaHDT) - plnorm(x, muHDT, sigmaHDT)
}

# Function to work out correction CFR
scale_cfr <- function(data_1_in, delay_fun){
  case_incidence <- data_1_in$new_cases
  death_incidence <- data_1_in$new_deaths
  cumulative_known_t <- 0 # cumulative cases with known outcome at time tt
  # Sum over cases up to time tt
  for(ii in 1:nrow(data_1_in)){
    known_i <- 0 # number of cases with known outcome at time ii
    for(jj in 0:(ii - 1)){
      known_jj <- (case_incidence[ii - jj]*delay_fun(jj))
      known_i <- known_i + known_jj
    }
    cumulative_known_t <- cumulative_known_t + known_i # Tally cumulative known
  }
  
  total_deaths <- sum(death_incidence)
  total_cases <- sum(case_incidence)
  
  # naive CFR value
  b_tt <- total_deaths / total_cases
  
  # underestimation of CFR due to unknown outcomes
  u <- cumulative_known_t / total_cases
  
  # correct the CFR
  p_tt <- b_tt / u
  
  data.frame(
    nCFR = b_tt,
    cCFR = p_tt,
    underestimation = u,
    total_deaths = sum(death_incidence),
    cum_known_t = round(cumulative_known_t),
    total_cases = sum(case_incidence)
  )
}



# Get data
httr::GET("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", httr::authenticate(":", ":", type="ntlm"), httr::write_disk(tf <- tempfile(fileext = ".csv")))
allDat <- read.csv(tf)

library(dplyr)

allDatDesc <- allDat %>% 
  dplyr::arrange(countriesAndTerritories, dateRep) %>% 
  dplyr::mutate(dateRep = lubridate::dmy(dateRep))%>% 
  dplyr::rename(date = dateRep, new_cases = cases, new_deaths = deaths, country = countriesAndTerritories) %>%
  dplyr::select(date, country, new_cases, new_deaths) %>%
  dplyr::filter(country != "CANADA", 
                country != "Cases_on_an_international_conveyance_Japan")
# Do analysis
allTogetherClean2 <- allDatDesc %>%
  dplyr::group_by(country) %>%
  padr::pad() %>%
  dplyr::mutate(new_cases = tidyr::replace_na(new_cases, 0),
                new_deaths = tidyr::replace_na(new_deaths, 0)) %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(cum_deaths = sum(new_deaths)) %>%
  dplyr::filter(cum_deaths > 0) %>%
  dplyr::select(-cum_deaths) %>%
  dplyr::do(scale_cfr(., delay_fun = hospitalisation_to_death_truncated)) %>%
  dplyr::filter(cum_known_t > 0 & cum_known_t >= total_deaths)  %>%
  dplyr::filter(total_deaths > 10)

# Now estimate the reporting rate of symptomatics for all countries
# simultaneously, using a Bayesian model with uninformative priors fit by MCMC.
# Note: there is currently no information shared between the country estimates
# so this is the same as modelling each country separately. The advantage of
# modelling them all simultaneously is computational only.

# Using the following model:
#
#   total_deaths ~ Binomial(total_cases, nCFR)
#   nCFR = cCFR * underestimation
#   cCFR = baseline_CFR / psi
#
# where psi is the reporting rate. Underestimation is known, and the (prior)
# distribution of baseline_CFR is known, so we can estimate psi from this model,
# incorporating uncertainty in the baseline CFR estimate into the estimates.

library(greta)

n <- nrow(allTogetherClean2)
deaths <- allTogetherClean2$total_deaths
cases <- allTogetherClean2$total_cases
underestimation <- allTogetherClean2$underestimation

# Distribution over plausible baseline CFR values from China study. The 95% CIs
# are symmetric around the estimate, so we assume it's a an approximately
# Gaussian distribution, truncated to allowable values. 
baseline_cfr <- normal(1.38, 0.077, dim = n, truncation = c(0, 1))

# A separate reporting rate for each country, with all reporting rates a priori
# equally as likely.
psi <- uniform(0, 1, dim = n)

# Observation model. Equivalent to, but more numerically stable than:
#   nCFR <- (baseline_cfr / 100) * underestimation / psi
log_nCFR <- log(baseline_cfr) - log(100) + log(underestimation) - log(psi)
nCFR <- exp(log_nCFR)
distribution(deaths) <- binomial(cases, nCFR)

set.seed(2020-04-02)
m <- model(psi)
draws <- mcmc(m, chains = 10, n_samples = 3000)

# check convergence before continuing
coda::gelman.diag(draws)
sry <- summary(draws)

# summarise estimates
allTogetherClean2$underreporting_estimate <- sry$statistics[, "Mean"]
allTogetherClean2$lower <- sry$quantiles[, "2.5%"]
allTogetherClean2$upper <- sry$quantiles[, "97.5%"]

reportDataFinal <- allTogetherClean2 %>%
  dplyr::select(country, total_cases, total_deaths, underreporting_estimate, lower,
                upper) %>%
  #dplyr::mutate(is.numeric, signif, digits=2)  %>%
  dplyr::mutate(underreporting_estimate = ifelse(underreporting_estimate <= 1, underreporting_estimate, 1)) %>%
  dplyr::mutate(upper = ifelse(upper <= 1, upper, 1)) %>%
  #dplyr::mutate(top = ifelse(top <= 1, top, 1)) %>%
  dplyr::mutate(underreporting_estimate = signif(underreporting_estimate, 2)) %>%
  dplyr::mutate(lower = signif(lower, 2)) %>%
  dplyr::mutate(upper = signif(upper, 2)) %>%
  #dplyr::mutate(bottom = signif(bottom, 2)) %>%
  #dplyr::mutate(top = signif(top, 2)) %>%
  dplyr::ungroup(country) %>%
  dplyr::mutate(country = country %>% stringr::str_replace_all("_", " ")) %>% 
  dplyr::mutate(underreporting_estimate_clean = paste0(underreporting_estimate*100,
                                                "% (",lower*100,"% - ",upper*100,"%)"))

saveRDS(reportDataFinal, "data/reportDataFinal.rds")
