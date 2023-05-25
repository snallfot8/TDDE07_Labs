#Assignment 2

#Assume that you have asked 8 randomly selected persons about their monthly income 
#(in thousands Swedish Krona) and obtained the following eight observations: 33, 24, 48, 32, 55, 74, 23, and 17. 
#A common model for non-negative continuous variables is the log-normal distribution. 
#The log-normal distribution logN(μ,σ2) has density function

#where y > 0, −∞ < μ < ∞ and σ2 > 0. The log-normal distribution is related
#to the normal distribution as follows: if y ∼ logN(μ,σ2) then logy ∼ N(μ,σ2).
#2iid 2 2 Let y1,...,yn|μ,σ ∼ logN(μ,σ ), where μ = 3.6 is assumed to be known but σ
#is unknown with non-informative prior p(σ2) ∝ 1/σ2. The posterior for σ2 is the Inv−χ2(n,τ2) distribution, where

###Starting assignment###

###Variables###
income = c(33,24,48,32,55,74,23,17)
my = 3.6
Ndraw = 10000
###############

#Calculating tao
taoFunc = function(y) {
  return((sum((log(y)-my)**2))/length(y))
}

#Function calculate sigma squared
sigmaSquaredFunc = function(chi_vals,tao,n) {
 return(((n)*tao)/chi_vals)
}


#a) Draw 10000 random values from the posterior of σ2 by assuming μ = 3.6 and plot the posterior distribution.


tao = taoFunc(income) #Calculating tao value
n = length(income) #Defining number of data points
set.seed(12345)
chi_vals = rchisq(Ndraw,n) #Drawing 10000 random samples from chi-squared distribution
sigmaSquared = sigmaSquaredFunc(chi_vals, tao, n) #Calculating approximation for sigmaSquared
#Plotting the posterior distribution of simgasquared
plot(density(sigmaSquared))


#b) The most common measure of income inequality is the Gini coefficient, G, 
#where 0 ≤ G ≤ 1. G = 0 means a completely equal income distribution, whereas G = 1 means complete 
#income inequality (see e.g. Wikipedia for more information about the Gini coefficient). 
#It can be shown that G = 2Φ 􏰃σ/√2􏰄−1 when incomes follow a log N (μ, σ2) distribution.
#Φ(z) is the cumulative distribution function (CDF) for the standard normal distribution with mean zero 
#and unit variance. Use the posterior draws in a) to compute the posterior distribution of the Gini coe􏰂cient G 
#for the current data set.


G = 2*pnorm(sigmaSquared/sqrt(2))-1 #Calculating gini values
plot(density(G), main="Density distribution of Gini coefficient") #Density distribution of Gini values

#Centered around low values (~0.03) which shows equal distribution?

#c) Use the posterior draws from b) to compute a 95% equal tail credible interval for G. 
#A 95% equal tail credible interval (a,b) cuts o􏰀 2.5% percent of the posterior probability mass to the 
#left of a, and 2.5% to the right of b.


sorted_G = sort(G) #Sorting values in G
min = length(G)*0.025 #Calculating lower bound
max = length(G)*0.975 #Calculating upper bound
filt_G = sorted_G[min:max] #Filtering out tails
plot(density(G), main="Density distribution of 95% tail credible interval")
abline(v=sorted_G[min]) #Lower end of credible interval
abline(v=sorted_G[max]) #Upper end of credible interval
plot(density(filt_G), main="Density distribution of 95% tail credible ineterval") #Plot with only values in credible interval

#d) Use the posterior draws from b) to compute a 95% Highest Posterior Density Interval (HPDI) for G. 
#Compare the two intervals in (c) and (d). [Hint: do a kernel density estimate of the posterior of G using the 
#density function in R with default settings, and use that kernel density estimate to compute the HPDI. 
#Note that you need to order/sort the estimated density values to obtain the HPDI.].


density_vals = density(G) #Using density function to load info for Gini coefficients density
df = data.frame(x=density_vals$x, y=density_vals$y) #Storing coordinates and indexes in a data frame
df_ordered=df[order(df$y),] #Sorting the data frame in ascending order
df_ordered$index=seq(1,length(density_vals$x),1) #Creating index to keep track of x & y coordinates from G:s density
ind_to_remove = round(length(df$x)*0.05) #Calculating how many indexes to remove to get the 95% HPDI


library(dplyr)
df_ordered_filt = df_ordered %>% slice(-c(1:ind_to_remove)) #Removing lowest y values
df_ordered_final = df_ordered_filt[order(df_ordered_filt$x),] #Filtering in order of x values
plot(df_ordered_final$x, df_ordered_final$y, type="l") #Plotting result
