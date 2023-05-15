################################################################################
#######  This code defines a random walk using the random variables in 
#######  log-normal distribution.
####### 

library(ggplot2)
library(tseries)

set.seed(10)

# 1. Defining Function to sample random variables from the log normal distribution
lognormal_from_uniform <- function(n, mean_dist, sd_dist) {
  u1 <- runif(n)
  u2 <- runif(n)
  z <- sqrt(-2*log(u1))*cos(2*pi*u2)
  x <- mean_dist + sd_dist*z
  x <- exp(x)
  return(x)
}

mean_dist = 3
sd_dist = 1

n <- 10000
N <- 5

ZN <- rep(NA,n)

# Combining N log-normal distributions
for (i in 1:n) {
  ZN[i] <- sum(lognormal_from_uniform(N, mean_dist, sd_dist))
}

# 2. Theoretically compute the cumulative distribution function of ZN

# Define the range of values for the CDF
x <- seq(0, max(ZN), length=n)

# Calculate the theoretical CDF using dlnorm
lognorm_cdf <- function(x, mu, sigma, N) {
  p <- rep(0, length(x))
  for (i in 1:N) {
    p <- p + dlnorm(x, meanlog = mu, sdlog = sigma, log = FALSE)
  }
  p <- p / N
  print(p)
  cdf <- cumsum(p * mean(diff(x)))
  return(cdf)
}


cdf <- lognorm_cdf(x, mean_dist, sd_dist, N)

# 3. Write a Monte Carlo program to generate samples of ZN from above CDF

# Define x values for CDF plot
x <- seq(min(ZN), max(ZN), length.out = length(ZN))

# Define Monte Carlo function to compute ECDF
monte_carlo_ecdf <- function(x, ZN, n_sims = 1000) {
  ecdf_vals <- numeric(length(x))
  for (i in 1:n_sims) {
    # Sample from ZN with replacement
    ZN_sample <- sample(ZN, replace = TRUE)
    # Compute empirical CDF at each x value
    ecdf_vals <- ecdf_vals + (ZN_sample <= x)
  }
  # Normalize ECDF values by number of simulations
  ecdf_vals <- ecdf_vals / n_sims
  return(ecdf_vals)
}

# Compute empirical CDF using Monte Carlo simulation
ecdf_ZN <- monte_carlo_ecdf(x, ZN)


# 4. Plot a histogram to compare the above two steps
# Plot the ECDF and the theoretical CDF
plot(ecdf_ZN, col="blue", main="Empirical vs Theoretical CDF", xlab="ZN", ylab="Cumulative Probability", lwd=1)
lines(x, cdf, col="red", lwd=2)
legend("bottomright", legend=c("Empirical CDF", "Theoretical CDF"), col=c("blue", "red"), lwd=2)


# 5. Central Limit Theorem. Show that UN converges to standard theorem

n <- 10000
N <- 100000

ZN <- rep(NA,n)

for (i in 1:n) {
  ZN[i] <- sum(lognormal_from_uniform(N, mean_dist, sd_dist))
}

# Compute the mean and standard deviation of the log-normal distribution
E_X = exp(mean_dist + 0.5*sd_dist^2)
S_X = sqrt((exp(sd_dist^2)-1)*exp(2*mean_dist+sd_dist^2))

UN = (ZN-N*E_X)/S_X

hist(UN,100,probability='TRUE')

library(tseries)

jarque.bera.test(UN)


