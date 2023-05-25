#Assignment 1

#Let y1, ..., yn|θ ∼ Bern(θ), and assume that you have obtained a sample
#with s = 22 successes in n = 70 trials. Assume a Beta(α0, β0) prior for θ and let α0 = β0 = 8.


###Starting assignment###


###Variables###
s=22
n=70
a0=8
b0=8
###############

#Assignment 1
meanBeta <- function(a,b){
  return(a/(a+b))
}

varBeta <-function(a,b) {
  return(a*b/(((a+b)**2)*(a+b+1)))
}

BetaPlot <- function(a,b){
  xGrid <- seq(0.001, 0.999, by=0.001)
  prior = dbeta(xGrid, a, b)
  maxDensity <- max(prior) # Use to make the y-axis high enough
  plot(xGrid, prior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
       ylab = 'Density', main = 'Beta(a,b) density')
}

#a) Draw 10000 random values (nDraws = 10000) from the posterior θ|y ∼ Beta(α0+ s, β0 + f), where y = (y1, . . . , yn),
#and verify graphically that the posterior mean E [θ|y] and standard deviation SD [θ|y] converges to the true values
#as the number of random draws grows large. [Hint: use rbeta() to draw random values and make graphs of the sample means 
#and standard deviations of θ as a function of the accumulating number of drawn values].

BetaPlot(a0,b0) #Pdf of prior

BetaPlot(a0+s,b0+n-s) #Pdf of posterior
meanPost = meanBeta(a0+s,b0+n-s) #Calculate the actual mean of the posterior
meanPost #actual is 0.3488
varPost = varBeta(a0+s,b0+n-s) #Calculate the actual variance of the posterior
varPost #actual variance is 0.00261
stdPost = sqrt(varPost) #Standard deviation of actual
stdPost #actual standard deviation is 0.0511

Ndraws = 10000 #Number of samples to draw
samples = rbeta(Ndraws, shape1=a0+s, shape2=b0+n-s) #Draws 10000 random samples from the posterior

mean_vals = rep(0,10000) #List to store mean
std = rep(0,10000) #List to store standard deviation

for (i in 1:10000) {
  mean = mean(samples[1:i])
  mean_vals[i] = mean
  std[i] = sd(samples[1:i])
}

plot(mean_vals, type ="l", col="red") #Plotting accum. mean vals based on increased Ndraws
abline(meanPost,0) #plotting tangent line of actual mean

plot(std, type="l", col="red", ylim=c(0,0.06)) #Plotting accum. standard deviation based on increased Ndraws
abline(stdPost,0) #plotting tangent line of actual std

#As seen in the plotted graphs, the mean and std dev converges towards the true values for big n:s

#b) Draw 10000 random values from the posterior to compute the posterior probability
#Pr(θ > 0.3|y) and compare with the exact value from the Beta posterior. [Hint: use pbeta()].


threshold = 0.3
true_probability = 1 - pbeta(threshold,shape1=a0+s, shape2=b0+n-s) #Calculating the true probability of the posterior
true_probability #True value is 0.83

samples = rbeta(Ndraws, shape1=a0+s, shape2=b0+n-s) #Drawing 10000 samples from the beta distribution
p_greater = mean(samples > threshold) #Calculating the mean of all random samples greater than 0.3
p_greater #Calculated value is 0.8286

#We can see that the probability is ~83% for both, indicating we get good 
#approximations for large data sets of random samples

#c) Draw 10000 random values from the posterior of the odds φ = θ / 1−θ by using
#the previous random draws from the Beta posterior for θ and plot the posterior distribution of φ.
#[Hint: hist() and density() can be utilized].


samples = rbeta(Ndraws, shape1=a0+s, shape2=b0+n-s) #Sample 10000 random samples
phi = samples/(1-samples) #Calculate phi
hist(phi) #Plot histogram
plot(density(phi)) #Plot density function
