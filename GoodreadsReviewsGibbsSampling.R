############################
# Gibbs Sampling Example: Comparison for Two Populations
# Goodreads review scores of NY Times Bestsellers: Fiction vs. Nonfiction
# Nov 14, 2023
############################

#load the library for the inverse gamma distribution 
library(invgamma)

# PRIOR PARAMETERS: We'll use the same prior for both populations
# How would I change the code to allow for different priors for the two populations?
# Prior parameters for mu: 
lambda <- 4 #Bestsellers will have higher reviews
tau2 <- 10 #leave room for uncertainty
# Prior parameters for sigma2: 
# possible range for reviews is ~1-5, 
# then an approximate standard deviation would be 4/6
# so an approximate variance would be 0.4444
gamma <- 2.01
phi <- 0.4444
#prior expected value of variance: 
phi/(gamma-1)

# Plot the prior distributions to make sure they seem reasonable
par(mfrow=c(1,2))
curve(dnorm(x, lambda, sqrt(tau2)), xlim=c(-4, 12), ylab="prior density", main=expression(pi(mu)), xlab=expression(mu))
curve(dinvgamma(x, gamma, phi), xlim=c(0, 2), ylab="prior density", main=expression(pi(sigma^2)), xlab=expression(sigma^2))


# COLLECT DATA
#Goodreads average scores from NY Times Bestsellers dated Nov 19, 2023
fiction <- c(4.44, 3.84, 4.65, 4.2, 4.18, 4.33, 4.31, 4.32, 3.97, 4.44, 4, 4.34, 4.17, 4.22, 3.84)
n.f <- length(fiction)
nonfiction <- c(4.22, 3.78, 4.15, 4.37, 4.43, 4.56, 4.45, 4.4, 4.42, 4.49, 3.83, 3.85, 4.43, 4.23, 4.55)
n.n <- length(nonfiction)

#Plot data
hist(fiction)
hist(nonfiction)
#qq plot to check for normality 
#(not great but not too bad for only 15 observations)
qqnorm(fiction)
qqline(fiction)
qqnorm(nonfiction)
qqline(nonfiction)

# POSTERIOR DISTRIBUTIONS: Must use Gibbs Sampling Algorithm to approximate

#Starting values 
mu.f <- 4
sigma2.f <- 0.5
mu.n <- 4
sigma2.n <- 0.5

# initializations for the Gibbs Sampling Algorithm
iters <- 10000
muf.save <- mun.save <- rep(0, iters)
muf.save[1] <- mu.f
mun.save[1] <- mu.n
sigma2f.save <- sigma2n.save <- rep(0, iters)
sigma2f.save[1] <- sigma2.f
sigma2n.save[1] <- sigma2.n

#Gibbs Sampling Algorithm
for(t in 2:iters){
  
  # Full conditional of mu.f (update the value of the parameters)
  lambdaf.p <- (tau2*sum(fiction) + sigma2.f*lambda)/(tau2*n.f + sigma2.f)
  tau2f.p <- sigma2.f*tau2/(tau2*n.f + sigma2.f)
  #sample a new value of mu.f
  mu.f <- rnorm(1, lambdaf.p, sqrt(tau2f.p))
  #save the value of mu.f
  muf.save[t] <- mu.f
  
  #Repeat for mu.n
  lambdan.p <- (tau2*sum(nonfiction) + sigma2.n*lambda)/(tau2*n.n + sigma2.n)
  tau2n.p <- sigma2.n*tau2/(tau2*n.n + sigma2.n)
  mu.n <- rnorm(1, lambdan.p, sqrt(tau2n.p))
  mun.save[t] <- mu.n
  
  # full conditional of sigma2 (update the value of the parameters)
  gammaf.p <- gamma + n.f/2
  phif.p <- phi + sum((fiction - mu.f)^2 )/2
  #sample new value of sigma2.f
  sigma2.f <- rinvgamma(1, gammaf.p, phif.p) 
  #save the value of sigma2.f
  sigma2f.save[t] <- sigma2.f
  
  #Repeat for sigma2.n
  gamman.p <- gamma + n.n/2
  phin.p <- phi + sum((nonfiction - mu.n)^2 )/2
  sigma2.n <- rinvgamma(1, gamman.p, phin.p) 
  sigma2n.save[t] <- sigma2.n
  
}


# Trace plots (decide if we need to throw out the first few values)
par(mfrow=c(1,2))
plot(muf.save, type='l')
plot(sigma2f.save, type='l')

plot(mun.save, type='l')
plot(sigma2n.save, type='l')

#throw out the first few values
burn <- 100
muf.use <- muf.save[-(1:burn)]
sigma2f.use <- sigma2f.save[-(1:burn)]
mun.use <- mun.save[-(1:burn)]
sigma2n.use <- sigma2n.save[-(1:burn)]
plot(muf.use, type='l')
plot(sigma2f.use, type='l')
plot(mun.use, type='l')
plot(sigma2n.use, type='l')

#SUMMARIZE THE POSTERIOR DISTRIBUTION(S)

# posterior distributions of mu.f and mu.n 
#plot 
plot(density(muf.use), xlab=expression(mu), ylab="density", main=expression(pi(mu~"|"~data)), lwd=2)
lines(density(mun.use), col='gray', lwd=2)
legend("topleft", c("Fiction", "Nonfiction"), col=c("black", "gray"), lty=1, lwd=2)

# posterior distributions of sigma2.f and sigma2.n
plot(density(sigma2f.use), xlab=expression(sigma^2), ylab="density", main=expression(pi(sigma^2~"|"~data)), lwd=2)
lines(density(sigma2n.use), col='gray', lwd=2)
legend("topright", c("Fiction", "Nonfiction"), col=c("black", "gray"), lty=1, lwd=2)

#Compare mu.f to mu.n
mean(muf.use > mun.use)
#Given our data and prior knowledge, there is about a 30% chance that 
#NYTimes Bestseller fiction books have a higher goodreads rating than nonfiction books

# Credible interval
diff <- muf.use - mun.use
plot(density(diff), xlab=expression(mu[f]-mu[n]), main="Posterior of Average Fiction vs Nonfiction")
abline(v=0, lty=2)

quantile(diff, c(.025, .975))
#Given our data and prior knowledge, there is a 95% chance that the 
#NYTimes bestseller Fiction goodreads ratings are between 0.3 points 
#lower and 0.18 points higher than the nonfiction books
