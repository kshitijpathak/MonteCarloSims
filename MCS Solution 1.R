###############################################################################
####### This code generates n samples from LOG-NORMAL Distribution ############
#######       and Compares with the Theoretical distribution       ############
###############################################################################
set.seed(10)

# Importing libraries
library(ggplot2)

# Setting Number of samples
n <- 10000

# Set parameters of the distribution
mean_dist <- 0
sd_dist <- 0.5

# Function to Generate normally distributed random numbers using uniform distribution
# We are using the Box-Muller method to generate normally distributed random 
# variables using uniform distribution

normal_from_uniform <- function(n, mean_dist, sd_dist) {
  u1 <- runif(n)
  u2 <- runif(n)
  z <- sqrt(-2*log(u1))*cos(2*pi*u2)
  x <- mean_dist + sd_dist*z
  return(x)
}

# Generating Normally distributed random variables
calc_gen = normal_from_uniform(n, mean_dist, sd_dist)

# Checking mean is approximately equal to the set parameter
mean(calc_gen)

# Generating Log-Normally distributed random variables from normal random variables
x <- exp(calc_gen)

# Plotting histogram of random numbers with ggplot
ggplot(data = data.frame(x = x), aes(x = x)) +
  geom_histogram(binwidth = 0.1, aes(y = ..density..), fill = "gray", color = "blue") +
  labs(title = "Histogram of Log-Normal Distribution", x = "X", y = "PDF") +
  # Add theoretical curve
  stat_function(fun = dlnorm, args = list(meanlog = mean_dist, sdlog = sd_dist), color = "red", size = 1.5) +
  # Add labels
  geom_text(aes(x = 6, y = 0.9, label = "Generated Samples"), color = "blue", size = 3) +
  geom_text(aes(x = 6, y = 0.8, label = "Theoretical curve"), color = "red", size = 3) +
  # Setting Theme of GGplot
  theme_bw()


# Comparing the generated log-normal distribution with the theoretical log-normal 
# distribution using the chi-squared goodness-of-fit test and the Kolmogorov-Smirnov (KS) test

# Setting number of bins in the histogram
num_bins <- 50

# Calculating observed frequencies
observed_freq <- hist(x, breaks = seq(0, max(x), length = num_bins+1), plot = FALSE)$counts

# Calculating expected frequencies from theoretical distribution
theoretical_freq <- n * diff(plnorm(seq(0, max(x), length = num_bins+1), meanlog = mean_dist, sdlog = sd_dist))

# Conducting Chi-squared goodness-of-fit test
chisq_test <- chisq.test(observed_freq, p = theoretical_freq, rescale.p = TRUE)

# Printing Chi-squared statistic and p-value
cat("Chi-squared statistic:", chisq_test$statistic, "\n")
cat("p-value:", chisq_test$p.value, "\n")
print("Rejecting the  Null Hypothesis, that states that no relationship exists between observed and expected distributions, at 95% confidence level.")

# Conducting KS test
ks_test <- ks.test(x, "plnorm", meanlog = mean_dist, sdlog = sd_dist)
p_value_ks <- ks_test$p.value

cat("Kolmogorov-Smirnov test results:\n")
cat(sprintf("  Test statistic: %.3f\n", ks_test$statistic))
cat(sprintf("  P-value: %.3e\n", p_value_ks))
print("Can not reject the null hypothesis of KS test that states no relationship exists between observed and expected distributions.")


