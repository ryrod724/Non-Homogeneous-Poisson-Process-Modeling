# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)

# Load data
eq.data <- read.csv("oilpipelineaccidents.csv", header=TRUE, sep=",")
eq.data$datetime <- as.POSIXct(eq.data$Accident.Date.Time, format="%m/%d/%Y %I:%M %p")

# Computing lag
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))

# Remove first row
eq.data <- eq.data[-1, ]

# Interarrival times (hours)
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600

### MODEL lambda(t) as constant ###

# Calculate total events and total observation time (days)
total_events <- nrow(eq.data)
total_days <- as.numeric(difftime(max(eq.data$datetime), min(eq.data$datetime), units = "days"))

# Constant intensity per day
lambda_const <- total_events / total_days

# Define lambda function as constant
lambda.fn <- function(t) {
  lambda_const
}

# Define cumulative intensity function Lambda(t)
Lambda.fn <- function(t) {
  lambda_const * t
}

### Time-rescaling ###

# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)

eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))

# Rescale each event time by integrating lambda from 0 to t (simple with constant lambda)
rescaled_times <- sapply(eq.data$time_in_days, function(t) {
  if (is.finite(t) && !is.na(t)) {
    return(lambda_const * t)
  } else {
    return(NA)
  }
})

# Compute rescaled interarrival times
rescaled_interarrivals <- diff(c(0, rescaled_times))  # Add 0 to start

### KS Goodness of fit test for Exp(1) ###

ks_result <- ks.test(jitter(rescaled_interarrivals), "pexp", 1)

# Show results
print(ks_result)

### Construct empirical CDF ###
elapsed_times <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))
elapsed_times <- sort(elapsed_times)
n <- length(elapsed_times)
empirical_cdf <- (1:n)/n

# Construct theoretical CDF
theoretical_cdf <- Lambda.fn(elapsed_times)
theoretical_cdf <- theoretical_cdf / max(theoretical_cdf)  # Normalize to [0,1]

# Plot empirical vs theoretical CDF
plot(elapsed_times, empirical_cdf, type="s", lwd=2, col="black", 
     xlab="Days since January 1, 2016", ylab="CDF")
lines(elapsed_times, theoretical_cdf, col="blue", lwd=2)
legend("bottomright", legend=c("Empirical CDF", "Theoretical CDF"), 
       col=c("black", "blue"), lwd=2)
