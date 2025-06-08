# Load libraries
library(Hmisc)
library(gtools)
library(dplyr)
library(lubridate)

# Load data
eq.data <- read.csv("/Users/ryanrodrigue/Downloads/earthquakes.csv", header=TRUE, sep=",")
eq.data$datetime <- as.POSIXct(eq.data$datetime, format="%Y-%m-%d %H:%M:%S")

# Filter to October-December 2024 only (most recent 3 month period)
eq.data <- eq.data[290257:294855, ]

# Computing lag
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))

# Remove first row
eq.data <- eq.data[-1, ]

# Interarrival times (hours)
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600

# Remove immediate aftershocks (within 2 hours)
eq.data <- eq.data[eq.data$elapsed.time > 2, ]

### MODEL lambda(t) ###

# Create year-month variable
eq.data$year.month <- format(as.Date(eq.data$datetime), "%Y-%m")

# Number of earthquakes per month
freq.month <- data.frame(table(eq.data$year.month))
year.month.unique <- freq.month[,1]
neq.month <- freq.month[,2]

# Number of days per month
day1.month <- ymd(paste(year.month.unique, "01", sep="-"))
ndays.month <- monthDays(as.Date(day1.month, "%Y-%m-%d"))

# Estimate intensity per day
lambda <- neq.month / ndays.month

# Cumulative number of days to each month
median.time <- c()
ndays.total <- c()
median.time[1] <- ndays.month[1] / 2
ndays.total[1] <- ndays.month[1]
for (i in 2:length(ndays.month)) {
  median.time[i] <- ndays.total[i-1] + ndays.month[i]/2
  ndays.total[i] <- ndays.total[i-1] + ndays.month[i]
}
median.time <- as.numeric(median.time)

# Polynomial regression
median.time.re <- median.time / 1000
median.time.sq <- median.time.re^2
median.time.cu <- median.time.re^3
median.time.qd <- median.time.re^4
median.time.qu <- median.time.re^5
median.time.sx <- median.time.re^6

model <- lm(lambda ~ median.time.re + median.time.sq + median.time.cu)
coefs <- coef(model)

# Define estimated lambda function
lambda.fn <- function(t) {
  coefs[1] + coefs[2]*(t/1000) + coefs[3]*(t/1000)^2 + 
    coefs[4]*(t/1000)^3
}

### Time-rescaling ###

# Create "time in days" since start
start_time <- min(eq.data$datetime, na.rm = TRUE)

eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))

# Now rescale each event time by integrating lambda from 0 to t
# Ensure t is valid before integrating
rescaled_times <- sapply(eq.data$time_in_days, function(t) {
  if (is.finite(t) && !is.na(t)) {
    return(integrate(lambda.fn, lower = 0, upper = t)$value)
  } else {
    return(NA)  # Return NA if t is not valid
  }
})


# Compute rescaled interarrival times
rescaled_interarrivals <- diff(c(0, rescaled_times))  # Add 0 to start

### KS Goodness of fit via test for Exp(1) ###

ks_result <- ks.test(rescaled_interarrivals, "pexp", 1)

# Show results
print(ks_result)

### Construct empirical CDF ###
elapsed_times <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))
elapsed_times <- sort(elapsed_times)
n <- length(elapsed_times)
empirical_cdf <- (1:n)/n

# Construct theoretical CDF
# Define cumulative intensity function Lambda(t) by integrating lambda(t)
Lambda.fn <- function(t) {
  sapply(t, function(x) {
    tryCatch({
      integrate(lambda.fn, 0, x)$value
    }, error = function(e) {
      NA  # Return NA if integration fails
    })
  })
}

theoretical_cdf <- Lambda.fn(elapsed_times)
theoretical_cdf <- theoretical_cdf / max(theoretical_cdf)  # Normalize to [0,1]

# Plot empirical vs theoretical CDF
plot(elapsed_times, empirical_cdf, type="s", lwd=2, col="black", 
     xlab="Days since October 1, 2024", ylab="CDF")
lines(elapsed_times, theoretical_cdf, col="blue", lwd=2)
legend("bottomright", legend=c("Empirical CDF", "Theoretical CDF"), 
       col=c("black", "blue"), lwd=2)
