library(lubridate)

# Load and clean data
eq.data <- read.csv("/Users/ryanrodrigue/Downloads/wildfireDATA.csv", header=TRUE, sep=",")
eq.data$DISCOVERY_DATE <- as.Date(eq.data$DISCOVERY_DATE, format="%m/%d/%Y")
eq.data <- subset(eq.data, FIRE_YEAR == 2018 & STATE == "OR")

# Convert time
eq.data$DISCOVERY_HOUR <- as.numeric(substr(eq.data$DISCOVERY_TIME, 1, 2))
eq.data$DISCOVERY_MIN <- as.numeric(substr(eq.data$DISCOVERY_TIME, 3, 4))
eq.data$datetime <- as.POSIXct(
  paste(eq.data$DISCOVERY_DATE, sprintf("%02d:%02d:00", eq.data$DISCOVERY_HOUR, eq.data$DISCOVERY_MIN)),
  format="%Y-%m-%d %H:%M:%S"
)

# Remove first row and compute interarrival time
eq.data$datetime.lag <- c(0, head(eq.data$datetime, -1))
eq.data <- eq.data[-1, ]
eq.data$elapsed.time <- (as.numeric(eq.data$datetime) - as.numeric(eq.data$datetime.lag)) / 3600
eq.data <- eq.data[eq.data$elapsed.time > 2, ]

# Time in days since first fire
start_time <- min(eq.data$datetime, na.rm = TRUE)
eq.data$time_in_days <- as.numeric(difftime(eq.data$datetime, start_time, units="days"))

# Sort and build empirical CDF
elapsed_times <- sort(eq.data$time_in_days)
n <- length(elapsed_times)
empirical_cdf <- (1:n) / n
cdf_df <- data.frame(
  t = elapsed_times,
  cdf = empirical_cdf
)

# Fit logistic CDF via nls
logistic_model <- nls(
  cdf ~ 1 / (1 + exp(-a * (t - b))),
  data = cdf_df,
  start = list(a = 0.04, b = 180),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

# Extract fitted parameters
coefs <- coef(logistic_model)
a_fit <- coefs["a"]
b_fit <- coefs["b"]
cat(sprintf("Fitted logistic parameters: a = %.4f, b = %.2f\n", a_fit, b_fit))

# Define fitted CDF
Lambda.fn <- function(t) {
  1 / (1 + exp(-a_fit * (t - b_fit)))
}

# Evaluate theoretical CDF at same time points
theoretical_cdf <- Lambda.fn(elapsed_times)

# Kolmogorov-Smirnov test between empirical and theoretical CDFs
ks_result <- ks.test(empirical_cdf, theoretical_cdf)
print(ks_result)

# Plot
plot(elapsed_times, empirical_cdf, type="s", lwd=2, col="black",
     xlab="Days since Jan 1, 2018", ylab="CDF")
lines(elapsed_times, theoretical_cdf, col="blue", lwd=2)
legend("bottomright", legend=c("Empirical CDF", "Fitted Logistic CDF"), 
       col=c("black", "blue"), lwd=2)
