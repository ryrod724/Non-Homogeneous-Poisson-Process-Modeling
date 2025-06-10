# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)
library(readr)

# Load data
eq.data <- read_csv("tornadoes.csv") %>%
  filter(
    st %in% c("OK", "TX", "KS", "NE"),
    mo %in% 1:12,
    yr %in% c(2020)
  )

eq.data$datetime <- as.POSIXct(eq.data$datetime_utc, format="%Y-%m-%d %H:%M:%S")

# Lag and interarrival filtering
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))
eq.data <- eq.data[-1, ]
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600
eq.data <- eq.data[eq.data$elapsed.time > 1, ]

# Time in days since start
start_time <- min(eq.data$datetime, na.rm = TRUE)
eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units = "days"))

# Gumbel PDF
gumbel_pdf <- function(t, mu, beta) {
  z <- (t - mu) / beta
  (1 / beta) * exp(-(z + exp(-z)))
}

# Objective: minimize KS statistic with scaling to observed event count
objective_fn <- function(params) {
  mu <- params[1]
  beta <- params[2]

  # Raw lambda without scaling
  raw_lambda <- function(t) {
    gumbel_pdf(t, mu, beta)
  }

  # Calculate total integral of raw_lambda over observed time
  total_lambda <- tryCatch(integrate(raw_lambda, 0, max(eq.data$time_in_days))$value,
                           error = function(e) NA)
  if (!is.finite(total_lambda) || total_lambda == 0) return(Inf)

  # Scale factor to match expected events with observed number
  scale_factor <- nrow(eq.data) / total_lambda

  # Scaled lambda function
  lambda.fn <- function(t) scale_factor * raw_lambda(t)

  # Compute Lambda(t) for each observed time
  Lambda <- function(t) {
    sapply(t, function(x) {
      tryCatch(integrate(lambda.fn, 0, x)$value,
               error = function(e) NA)
    })
  }

  rescaled_times <- Lambda(eq.data$time_in_days)
  interarrivals <- diff(c(0, rescaled_times))

  # If any non-finite values, return large penalty
  if (any(!is.finite(interarrivals))) return(Inf)

  # KS test of rescaled interarrival times against Exp(1)
  ks <- suppressWarnings(ks.test(interarrivals, "pexp", 1))
  return(ks$statistic)
}

# Optimize mu and beta (adjust initial guess as needed)
opt_result <- optim(par = c(40, 10), fn = objective_fn, method = "L-BFGS-B", lower = c(0.01, 0.01))
mu_opt <- opt_result$par[1]
beta_opt <- opt_result$par[2]
cat("Estimated mu:", mu_opt, "\n")
cat("Estimated beta:", beta_opt, "\n")

# Final scaled lambda function
raw_lambda <- function(t) gumbel_pdf(t, mu_opt, beta_opt)
total_lambda <- integrate(raw_lambda, 0, max(eq.data$time_in_days))$value
scale_factor <- nrow(eq.data) / total_lambda
lambda.fn <- function(t) scale_factor * raw_lambda(t)

# Rescale times for goodness of fit testing
rescaled_times <- sapply(eq.data$time_in_days, function(t) {
  if (is.finite(t) && !is.na(t)) {
    integrate(lambda.fn, 0, t)$value
  } else {
    NA
  }
})
rescaled_interarrivals <- diff(c(0, rescaled_times))

# KS test on rescaled interarrival times
ks_result <- ks.test(rescaled_interarrivals, "pexp", 1)
print(ks_result)

# Empirical CDF of observed times
elapsed_times <- sort(eq.data$time_in_days)
n <- length(elapsed_times)
empirical_cdf <- (1:n)/n

# Theoretical CDF based on scaled Lambda(t)
Lambda.fn <- function(t) {
  sapply(t, function(x) {
    tryCatch(integrate(lambda.fn, 0, x)$value,
             error = function(e) NA)
  })
}
theoretical_cdf <- Lambda.fn(elapsed_times)
theoretical_cdf <- theoretical_cdf / max(theoretical_cdf, na.rm = TRUE)  # Normalize to 1

# Plot empirical vs theoretical CDF
plot(elapsed_times, empirical_cdf, type = "s", lwd = 2, col = "black",
     xlab = "Days since January 1, 2020", ylab = "CDF")
lines(elapsed_times, theoretical_cdf, col = "blue", lwd = 2)
legend("bottomright", legend = c("Empirical CDF", "Theoretical CDF"),
       col = c("black", "blue"), lwd = 2)
