# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)

# Load and clean data
eq.data <- read.csv("wildfires.csv", header=TRUE)
eq.data$DISCOVERY_DATE <- as.Date(eq.data$DISCOVERY_DATE, format="%m/%d/%Y")

# Filter for 2018 and California
eq.data <- subset(eq.data, FIRE_YEAR == 2018 & STATE == "CA")

# Drop rows with missing or malformed time
eq.data <- eq.data[grepl("^\\d{4}$", eq.data$DISCOVERY_TIME), ]

# Parse time
eq.data$DISCOVERY_HOUR <- as.numeric(substr(eq.data$DISCOVERY_TIME, 1, 2))
eq.data$DISCOVERY_MIN <- as.numeric(substr(eq.data$DISCOVERY_TIME, 3, 4))
eq.data$datetime <- as.POSIXct(
  paste(eq.data$DISCOVERY_DATE, sprintf("%02d:%02d:00", eq.data$DISCOVERY_HOUR, eq.data$DISCOVERY_MIN)),
  format="%Y-%m-%d %H:%M:%S"
)

# Sort by datetime
eq.data <- eq.data[order(eq.data$datetime), ]

# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)
eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))

# Construct cumulative counts for fitting
cum_counts <- seq_along(eq.data$time_in_days)
time_days <- eq.data$time_in_days

# Load minpack.lm for robust nonlinear fitting
if(!require(minpack.lm)) install.packages("minpack.lm")
library(minpack.lm)

# 3-parameter logistic function
logistic_cdf <- function(t, L, k, x0) {
  L / (1 + exp(-k * (t - x0)))
}

# Starting parameter guesses
start_logistic <- list(
  L = max(cum_counts) * 1.1,
  k = 0.1,
  x0 = median(time_days)
)

# Fit the logistic curve
logistic_fit <- nlsLM(
  cum_counts ~ logistic_cdf(time_days, L, k, x0),
  start = start_logistic,
  control = nls.lm.control(maxiter = 1000, ftol = 1e-10)
)

# Extract parameters
params_logistic <- coef(logistic_fit)
print(params_logistic)

# Define cumulative intensity function Lambda(t)
Lambda.fn <- function(t) {
  params_logistic["L"] / (1 + exp(-params_logistic["k"] * (t - params_logistic["x0"])))
}

# Derivative (intensity function lambda(t))
lambda.fn <- function(t) {
  L <- params_logistic["L"]
  k <- params_logistic["k"]
  x0 <- params_logistic["x0"]
  exp_term <- exp(-k * (t - x0))
  (L * k * exp_term) / (1 + exp_term)^2
}

# Compute rescaled times using Lambda(t)
rescaled_times <- Lambda.fn(time_days)

# Compute rescaled interarrival times
rescaled_interarrivals <- diff(c(0, rescaled_times))

# Add tiny jitter to avoid ties in KS test
set.seed(123)
rescaled_interarrivals_jittered <- rescaled_interarrivals + runif(length(rescaled_interarrivals), 0, 1e-10)

# Perform KS test against Exp(1)
ks_result <- ks.test(rescaled_interarrivals_jittered, "pexp", 1)
print(ks_result)

# Empirical CDF of elapsed times
elapsed_times <- sort(time_days)
n <- length(elapsed_times)
empirical_cdf <- (1:n) / n

# Theoretical CDF normalized to [0,1]
theoretical_cdf <- Lambda.fn(elapsed_times) / params_logistic["L"]

# Plot empirical vs theoretical CDF
plot(elapsed_times, empirical_cdf, type = "s", lwd = 2, col = "black",
     xlab = "Days since January 1, 2018", ylab = "CDF",
     main = "Empirical vs Theoretical CDF")
lines(elapsed_times, theoretical_cdf, col = "blue", lwd = 2)
legend("bottomright", legend = c("Empirical CDF", "Theoretical CDF"),
       col = c("black", "blue"), lwd = 2)
       
#Important results
nrow(eq.data)
lambda.fn
Lambda.fn
