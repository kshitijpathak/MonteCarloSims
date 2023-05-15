################################################################################
####### This Code compares the Theoretically calculated invariant distribution 
####### of a Markov Chain to a Monte Carlo simulated Markov chain

library(markovchain)

# Define the transition probability matrix
# Row 1 of Transition probability Matrix
a = 0.55
b = 0.36
c = 1-(a+b)

# Row 2 of Transition probability Matrix
d = 0.42
e = 0.29
f = 1-(d+e)

# Row 3 of Transition probability Matrix
g = 0.08
h = 0.22
i = 1-(g+h)

# Verify row sum should be equal to 1
sum(a,b,c)
sum(d,e,f)
sum(g,h,i)

P <- matrix(c(a,b,c,d,e,f,g,h,i), nrow = 3, byrow = TRUE)
P

# Calculating theoretical Invariant distribution
eig <- eigen(t(P))
pie <- eig$vectors[, which(abs(eig$values - 1) < 1e-8)]
pie <- pie / sum(pie)  # Normalize the vector so that it sums to 1

# Define the number of iterations
n <- 100000

# Set the initial state
x <- sample(c(1,2,3), 1)

# Define the inverse CDF function for the first state
inv_cdf_state1 <- function(u) {
  qnorm(u, mean = 2, sd = 1)
}

# Define the inverse CDF function for the second state
inv_cdf_state2 <- function(u) {
  qnorm(u, mean = 0, sd = 1)
}

# Define the inverse CDF function for the third state
inv_cdf_state3 <- function(u) {
  qnorm(u, mean = -2, sd = 1)
}

# Define the function to perform the Monte Carlo simulation with quantile transformation
simulate_markov_chain <- function(x, transition_probs, inv_cdf_state1, inv_cdf_state2, inv_cdf_state3) {
  # Generate a uniform random number
  u <- runif(1)
  
  # Use the inverse CDF function to transform the uniform random number to the next state
  if (x == 1) {
    x_new <- ifelse(u < transition_probs[1, 1], 1, ifelse(u < sum(transition_probs[1, 1:2]), 2, 3))
    y <- inv_cdf_state1(u)
  } else if (x == 2) {
    x_new <- ifelse(u < transition_probs[2, 1], 1, ifelse(u < sum(transition_probs[2, 1:2]), 2, 3))
    y <- inv_cdf_state2(u)
  } else {
    x_new <- ifelse(u < transition_probs[3, 1], 1, ifelse(u < sum(transition_probs[3, 1:2]), 2, 3))
    y <- inv_cdf_state3(u)
  }
  
  # Return the new state and the transformed random number
  list(x_new, y)
}

# Initialize the output vectors
chain <- numeric(n)
samples <- numeric(n)

# Generate the Markov chain using the Monte Carlo simulation with quantile transformation
for (i in 1:n) {
  # Perform the Monte Carlo simulation
  result <- simulate_markov_chain(x, P, inv_cdf_state1, inv_cdf_state2, inv_cdf_state3)
  
  # Store the new state and the transformed random number
  x <- result[[1]]
  y <- result[[2]]
  
  # Store the new state and the transformed random number in the output vectors
  chain[i] <- x
  samples[i] <- y
}

# Compute the empirical distribution
empirical_dist <- table(chain) / n


# Simulating the Markov chain with n transitions using library function
#######################################################################

mc <- new("markovchain", states = c("1","2","3"), transitionMatrix = P, name = "Three-state Markov chain")
class(mc)

# Simulate n from the Markov chain
current_state = sample(c("1","2","3"), 1)
simulated_states <- c(current_state)

for (i in 1:n) {
  current_state <- rmarkovchain(n = 1, object = mc, t0 = current_state)#$state
  simulated_states <- c(simulated_states, current_state)
}

# Compute the empirical distribution of states
empirical_dist_lib <- table(simulated_states) / n

# Compare the empirical distributions with the theoretical invariant distribution
print(empirical_dist)
print(empirical_dist_lib)
print(pie)


# Compare the empirical distribution with the theoretical invariant distribution in a plot
barplot(rbind(empirical_dist, empirical_dist_lib, pie), beside = TRUE, col = c("red", "blue", "green"), names.arg = c("State 1", "State 2", "State 3"), ylim = c(0, 0.5), ylab = "Probability", main = "Empirical vs Theoretical Invariant Distribution")
legend("top", c("Empirical", "Empirical Lib", "Theoretical"), fill = c("red", "blue", "green"), border = NA)

###############################################################################
# We can also see the convergence of the empirical distribution and theoretical
# as a function of n, number of iterations. ################################### 
###############################################################################

# Initialize vectors to store the results
n_vec <- seq(100, 50000, 100)
empirical_dist_vec <- matrix(0, nrow = length(n_vec), ncol = 3)
theoretical_dist_vec <- matrix(0, nrow = length(n_vec), ncol = 3)

# Looping through the difference number of transitions

for (n in n_vec) {
  print(n)
  chain <- numeric(n)
  
  # Randomly setting the initial state
  x <- sample(1:3, 1)
  
  # Generate the Markov chain using the Monte Carlo simulation with quantile transformation
  for (i in 1:n) {
    # Perform the Monte Carlo simulation
    result <- simulate_markov_chain(x, P, inv_cdf_state1, inv_cdf_state2, inv_cdf_state3)
    
    # Store the new state and the transformed random number
    x <- result[[1]]
    y <- result[[2]]
    
    # Store the new state and the transformed random number in the output vectors
    chain[i] <- x
    samples[i] <- y
  }
  
  # Compute the empirical distribution
  empirical_dist <- table(chain) / n
  empirical_dist_vec[which(n==n_vec),] <- empirical_dist
  
  # Calculating theoretical Invariant distribution
  eig <- eigen(t(P))
  pie <- eig$vectors[, which(abs(eig$values - 1) < 1e-8)]
  pie <- pie / sum(pie)  # Normalize the vector so that it sums to 1
  theoretical_dist_vec[which(n==n_vec),] <- pie
}

# Plot the results
plot(n_vec, empirical_dist_vec[,1], type = "l", col = "blue",
     xlab = "Number of Iterations (n)", ylab = "Proportion of Time Spent in State",
     ylim = c(0.2, 0.6))
lines(n_vec, empirical_dist_vec[,2], type = "l", col = "green")
lines(n_vec, empirical_dist_vec[,3], type = "l", col = "red")
lines(n_vec, theoretical_dist_vec[,1], type = "l", col = "blue", lty = 2)
lines(n_vec, theoretical_dist_vec[,2], type = "l", col = "green", lty = 2)
lines(n_vec, theoretical_dist_vec[,3], type = "l", col = "red", lty = 2)
legend("top", legend = c("Empirical (State 1)", "Empirical (State 2)",
                              "Empirical (State 3)", "Theoretical (State 1)",
                              "Theoretical (State 2)", "Theoretical (State 3)"),
       col = c("blue", "green", "red", "blue", "green", "red"),
       lty = c(1, 1, 1,2,2,2),ncol=2)


