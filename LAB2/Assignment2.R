# Assignment 2 Posterior approximation for classification with logistic regression
#The dataset WomenAtWork.dat contains n = 132 observations on the following eight variables related to women:

#2.a) 
#Consider the logistic regression model: Pr(y = 1|x, β) =  exp(xTβ) / (1 + exp (xT β))
# where y equals 1 if the woman works and 0 if she does not.
#x is a 7-dimensional vector containing the seven features (including a 1 to model the intercept). 
#The goal is to approximate the posterior distribution of the parameter vector β with a multivariate normal distribution
#  β|y, x ∼ N(B_hat, J^-1(B_hat)) where B_hat is the posterior mode and J(B_hat) = −∂2 lnp(β|y)|
#  Note that   ∂ β ∂ β T β = B_hat is the negative of ∂2 ln p(β|y) the observed Hessian evaluated at the posterior mode.
#Note that ∂β∂βT is a 7 × 7 matrix with second derivatives on the diagonal and cross-derivatives ∂2 ln p(β|y) ∂βi∂βj on the off-diagonal. 
#You can compute this derivative by hand, but we will let the computer do it numerically for you. 
#Calculate both B_hat and J(B_hat) by using the optim function in R. 
#[Hint: You may use code snippets from my demo of logistic regression in Lecture 6.] Use the prior β ∼ N(0,τ2I), where τ = 2.
#  Present the numerical values of β ̃ and J−1(β ̃) for the WomenAtWork dat
#a. Compute an approximate 95% equal tail posterior probability interval for the regression coefficient to the variable NSmallChild. 
#Would you say that this feature is of importance for the probability that a woman works?
#[Hint: You can verify that your estimation results are reasonable by comparing the posterior means
# to the maximum likelihood estimates, given by: glmModel <- glm(Work ~ 0 + ., data = WomenAtWork, family = binomial).]


#Import packages
library(mvtnorm)

#Read the data
data <-  (read.table("WomenAtWork.dat", header=TRUE))
Xnames <- names(data[,2:ncol(data)]) #Names of columns
head(data)
n_rows <- dim(data)[1]  #number of observations


LogPostLogistic <- function(betas,y,X,mu,Sigma){
  linPred <- X%*%betas;
  logLik <- sum( linPred*y - log(1 + exp(linPred)) );
  #if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  logPrior <- dmvnorm(betas, mu, Sigma, log=TRUE);
  
  return(logLik + logPrior)
}

#initialize variables
set.seed(123)
y <- as.numeric(data[,1]) #target output
X <- as.matrix(data[,-1]) #training data
X <- apply(X, c(1, 2), as.numeric) #convert from characters to doubles
n_col <- dim(X)[2] #number of columns in training data
tao = 2 #prior std.dev for Beta
initVal <- matrix(0,n_col,1) #initial values for Beta
mu <- as.matrix(rep(0,n_col)) # Prior mean vector for Beta
Sigma <- tao^2*diag(n_col) # Prior covariance matrix for Beta
logPost <- LogPostLogistic; #function to be minimized

OptimRes <- optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

# Printing the results of mean and variance of normally distributed Beta 
OptimBeta <- OptimRes$par
names(OptimBeta) <- Xnames # Naming the coefficient by covariates
PostCov <- solve(-OptimRes$hessian) #get covariance matrix by inverting the negative hessian
approxPostStd <- sqrt(diag(PostCov)) # Computing approximate standard deviations.
names(approxPostStd) <- Xnames # Naming the coefficient by covariates
print('The posterior mode is:')
print(OptimBeta)
print('The approximate posterior standard deviation is:')
print(approxPostStd)

#Draw samples from distribution of NSmallChild
my <- OptimBeta[6] #mean of Beta[NSmallChild]
sigma <-approxPostStd[6] #std.dev of Beta[NSmallChild]
samples <- rnorm(1000, my, sigma)
quantiles <-  quantile(samples, c(0.025, 0.975))
plot(density(samples), main="Posterior probability interval")
abline(v=quantiles)
print("Upper and lower bounds are:")
print(quantiles)

#Check that estimation results are reasonable
glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)
print(OptimBeta)
print(glmModel$coefficients)
#values are similar, estimation seems reasonable. Since the boundary is negative, it seems like a higher number 
#of small children makes a woman less likely to work. 

#b) Use your normal approximation to the posterior from (a). 
#Write a function that simulate draws from the posterior predictive distribution of Pr(y = 0|x), 
#where the values of x corresponds to a 40-year-old woman, with two children (4 and 7 years old), 
#11 years of education, 7 years of experience, and a husband with an income of 18. 
#Plot the posterior predictive distribution of Pr(y = 0|x) for this woman.
#[Hints: The R package mvtnorm will be useful. Remember that Pr(y = 0|x) can be calculated for each posterior draw of β.]

#sigmoid function
sigmoid <- function(x, beta) {
  s <- beta %*% x
  return (exp(s) / (1+exp(s)))
}

set.seed(123)
PostStd <- c(approxPostStd) #convert to vector
beta_samples <- rmvnorm(n=1000, mean=OptimBeta, sigma=PostCov) #draw samples of Beta from the distribution
x <- c(1, 18, 11, 7, 40, 1, 1) #values for a single sample of x

y <- rep(0,1000)
for(i in 1:1000) {
  y[i] = sigmoid(x, beta_samples[i,]) 
}

#draw random binomial samples given the probabilities
predictions <- rep(0,1000)
for(i in seq_along(y)) {
  p = y[i] #probability
  predictions[i] <- rbinom(1,1,p) 
}
plot(density(predictions), main="Posterior distribution")
barplot(table(predictions), main="Posterior predictive distribution") 
#We can see that the prediction is 0 for approximately 80% of the samples, which means it is
#likely that this woman does not work.

#c) Now, consider 13 women which all have the same features as the woman in (b). 
#Rewrite your function and plot the posterior predictive distribution for the number of women, 
#out of these 13, that are not working. 
#[Hint: Simulate from the binomial distribution, which is the distribution for a sum of Bernoulli random variables.]
set.seed(123)
n_women <- 13
y <- rep(0,1000)
predictions <- c()

#samples from Beta distribution 1000 times and make predictions. Simulate from binomial distribution for the 13 women. 
beta_samples <- rmvnorm(n=1000, mean=OptimBeta, sigma=PostCov) #draw samples of Beta from the distribution
for (i in 1:1000) {
  y[i] <- sigmoid(x, beta_samples[i,]) #probabilities
  predictions <- c(predictions, n_women-rbinom(n=1, size=n_women, y[i])) #simulate from binomial distribution where each result is prediction of number of women working
}

barplot(table(predictions), main=paste("Post. pre. dist. for", n_women, "women not working"), xlab="Number of women")
#The conclusion is that according to the posterior predictive distribution, it seems that 11 women most likely does not work, with high numbers
#for 8-12 women not working. 

