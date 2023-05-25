#The dataset Linkoping2022.xlsx contains daily average temperatures (in degree Celcius) 
#in Linköping over the course of the year 2022. Use the function read_xlsx(), which is included
#in the R package readxl (install.packages("readxl")), to import the dataset in R. The response variable is 
#temp and the covariate time that you need to create yourself is de􏰁ned by

#time =(#days since beginning of year)/365

#A Bayesian analysis of the following quadratic regression model is to be performed:
#temp=β0+β1·time+β2·time +ε,ε∼N(0,σ).

#Read data and import library
library(readxl)
data = read_xlsx("Linkoping2022.xlsx")

#Calculate variable time
data$time = as.integer(difftime(data$datetime,"2022-01-01", units = "day"))/365

#a)
# Use the conjugate prior for the linear regression model. The prior hyperparameters μ0, Ω0, ν0 and σ02
#shall be set to sensible values. Start with μ0 = (0,100,−100)T, Ω0 = 0.01 · I3, ν0 = 1 and σ02 = 1.
#Check if this prior agrees with your prior opinions by simulating draws from the joint prior of all parameters 
#and for every draw compute the regression curve. This gives a collection of regression curves; one for each draw 
#from the prior. Does the collection of curves look reasonable? If not, change the prior hyperparame- ters until the 
#collection of prior regression curves agrees with your prior beliefs about the regression curve.
#[Hint: R package mvtnorm can be used and your Inv-χ2 simulator of random draws from Lab 1.]

library(mvtnorm)

### LOAD DATA ###
my_zero = c(0,100,-100)
omega_zero = 0.1*diag(3)
v_zero = 100 #Increased from 1 since we got bad regression curves when using low values
sigma2_zero = 1

#################


###Simulate draws from prior for sigma squared
Ndraw = 1000 #Integer for how many draws we want to make
set.seed(12345)
chi_vals = rchisq(Ndraw,v_zero) #Drawing 10000 random samples from chi-squared distribution
sigmaSquared = v_zero*sigma2_zero/chi_vals #Calculating prior sigma squared

###Simulate draws from prior beta
beta=matrix(0,length(sigmaSquared),3) #Create matrix to store beta draws
for (i in 1:length(sigmaSquared)) {
beta[i,] = rmvnorm(1, my_zero, sigmaSquared[i]*solve(omega_zero)) #Drawing beta values from prior
}

reg = matrix(0,365,Ndraw) #Creating matrix to store the set of regressions
for (i in 1:Ndraw) {
reg[,i]=beta[i,1]+data$time*beta[i,2]+(data$time**2)*beta[i,3] #Calculating regressions
}
###Creating plot where all regressions can be plotted
plot.new()
plot.window(xlim=c(0,1), ylim=c(-10,40))
axis(side = 1)
axis(side = 2)
for (i in 1:Ndraw) {
  lines(data$time, reg[,i], type = "l") #Plotting regression
}
#Conclusion: By increasing the prior degrees of freedom, we got better regression curves. 

#b)
#Write a function that simulate draws from the joint posterior distribution of β0, β1,β2 and σ2.
#i) Plot a histogram for each marginal posterior of the parameters.
#ii) Make a scatter plot of the temperature data and overlay a curve for the posterior median of the 
#regression function f(time) = E[temp|time] = β0 + β1 · time + β2 · time2, i.e. the median of f (time)
#is computed for every value of time. In addition, overlay curves for the 90% equal tail posterior probability 
#intervals of f(time), i.e. the 5 and 95 posterior percentiles of f (time) is computed for every value of time. 
#Does the posterior probability intervals contain most of the data points? Should they?


###DATA###
n = length(data$temp)
v_n = v_zero + n
X = cbind(1, data$time, data$time**2)
beta_hat=solve(t(X)%*%X)%*%t(X)%*%data$temp
my_n=solve(t(X)%*%X+omega_zero)%*%(t(X)%*%X%*%beta_hat+omega_zero%*%my_zero)
omega_n = t(X)%*%X+omega_zero
sigma2_n = as.numeric((sigma2_zero*v_zero + (t(data$temp)%*%data$temp+t(my_zero)%*%omega_zero%*%my_zero-t(my_n)%*%omega_n%*%my_n))/v_n)
###########

###Draw posterior values for sigma squared
Ndraw = 1000 #Integer for how many draws we want to make
set.seed(12345)
chi_vals = rchisq(Ndraw,v_n) #Drawing 10000 random samples from chi-squared distribution
sigmaSquared = v_n*sigma2_n/chi_vals #Calculating posterior sigma squared

###Simulate draws for posterior beta
beta=matrix(0,length(sigmaSquared),3) #Create matrix to store beta simulations
for (i in 1:length(sigmaSquared)) {
  set.seed(12345)
  beta[i,] = rmvnorm(1, my_n, sigmaSquared[i]*solve(omega_n)) #Drawing beta values from posterior
}

###Plotting histogram for posterior draws
hist(sigmaSquared) #Histogram for sigma_squared
hist(beta[,1]) #Histogram for beta_zero
hist(beta[,2]) #Histogram for beta_one
hist(beta[,3]) #Histogram for beta_two

###Scatter plot
plot.new()
plot.window(xlim=c(0,1), ylim=c(-10,30))
axis(side = 1)
axis(side = 2)
points(data$time, data$temp, pch=19)

reg = matrix(0,365,Ndraw) #Creating matrix to store the set of regressions
for (i in 1:Ndraw) {
  reg[,i]=X%*%beta[i,] #Calculating regressions
}

temp_medians = rep(0,dim(reg)[1]) #Creating list to store median values of temps each day
for(i in 1:dim(reg)[1]) {
 temp_medians[i] = median(reg[i,]) 
}

#Plotting median values
lines(data$time, temp_medians, pch=19, col="red")

#Creating 90% crediible interval
min = 0.05*dim(reg)[2]+1
max = 0.95*dim(reg)[2]

under_limit = rep(0,dim(reg)[1])
upper_limit = rep(0,dim(reg)[1])
for(i in 1:dim(reg)[1]) {
  sorted = sort(reg[i,])
  under_limit[i]=quantile(sorted,probs=0.05)
  upper_limit[i]=quantile(sorted,probs=0.95)
}

lines(data$time, under_limit, col="green") 
lines(data$time, upper_limit, col="green")  

###Conclusion:The equal tail intervals are very narrow as we can see in the plot. As we saw in the previously
### presented histograms, the posterior draws are narrow which yields a very slim upper and under limit
### (they are very similar, but not the same however)

#c)
#It is of interest to locate the time with the highest expected temperature (i.e. the time where f(time) is maximal). 
#Let's call this value x ̃. Use the simulated draws in (b) to simulate from the posterior distribution of x ̃.
#[Hint: the regression curve is a quadratic polynomial. Given each posterior draw of β0, β1 and β2,
#you can find a simple formula for x ̃.]


#As it is a polynomial function with a maximum point, the time for the highest temp
#should be located where the derivative is 0

###Simulating the expected time where the temperature is at its peak
peak_temp = rep(0,Ndraw)
for (i in 1:Ndraw) {
  peak_temp[i]=(-1)*(beta[i,2]/(2*beta[i,3])) #Estimating time when temperature reaches its peak
}

hist(peak_temp) #Plotting histogram of when the temperature is expected to peak

###Conclusion: Once again, the posterior distribution is very narrow.
###Furthermore, posterior distribution that temperature is likely to peak in mid July
### which is reasonable

#d)
#Say now that you want to estimate a polynomial regression of order 10, but you suspect that higher 
#order terms may not be needed, and you worry about overfitting the data. Suggest a suitable prior that mitigates 
#this potential problem. You do not need to compute the posterior. 
#Just write down your prior. [Hint: the task is to specify μ0 and Ω0 in a suitable way.]

#Answer:
#The prior my_zero can be set to previous values. Prior omega_zero should be X*I (X times identity matrix).
#Rest of coefficients should be close to 0 to avoid overfitting from higher order terms.
#Furthermore, the X value for the omega_zero should be set to a large value in order to obtain
# a smaller spread for the beta values, making them more likely to stay close to zero.
#On the other hand, you could also argue that the X value should be low as this would give the data
#less impact and therefore decrease the risk of overfitting to the data.
