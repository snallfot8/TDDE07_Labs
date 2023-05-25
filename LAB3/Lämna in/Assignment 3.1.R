#Gibbs sampler for a normal model
#The dataset Precipitation.rds consists of daily records of weather with rain or snow (in units of mm)
#from the beginning of 1962 to the end of 2008 in a certain area. 
#Assume the natural log of the daily precipitation {y1,...,yn} to be indepent normally distributed, 
#ln y1, ..., ln yn|μ, σ ∼ N (μ, σ ), where both μ and σ are unknown. 
#Let μ ∼ N (μ0, τ02) independently of σ2 ∼ Inv-χ2(ν0, σ02).

#(a) Implement (code!) a Gibbs sampler that simulates from the joint posterior p(μ,σ2|lny1,...,lnyn). 
#The full conditional posteriors are given on the slides from Lecture 7. 
#Evaluate the convergence of the Gibbs sampler by calculating the Ineficiency Factors (IFs) 
#and by plotting the trajectories of the sampled Markov chains.

data=readRDS("Precipitation.rds") #Reading data
data = log(data)

#a)

# Initial setting (prior settings)
my0 = 1
tao0=1
v0=1
sigma0sq=1
nDraws = 1000
n = nDraws
true_mean = mean(data)
vn = v0 + n

#Lists for storing draws
my_vals = rep(0,nDraws)
sigma_vals = rep(0,nDraws)
sigma_vals[1] = 1 #Initial sigmasq value to be used

#Drawing posteriors
for(i in 1:nDraws) {

  w = (n/sigma_vals[i])/((n/sigma_vals[i])+(1/tao0**2)) #Calculating weights
  myn = w*true_mean+(1-w)*my0 #Calculating myn
  taonsq = 1/((n/sigma_vals[i])+(1/tao0)) #Calculating taonsq
  my_vals[i] = rnorm(1,mean=myn, taonsq) #Calculating posterior my_val
 
  if (i <= (nDraws-1)) {
  chi_val = rchisq(1,vn) #Drawing random value from chi distribution
  sigmansq = (v0*sigma0sq+sum((data-my_vals[i])**2))/(n+v0) #Calculating approximation sigmasq
  sigma_vals[i+1] = vn*sigmansq/chi_val #Calculating posterior sigma_val
  }
}

#Plotting my convergence
plot(my_vals[10:nDraws], ,main = "Plot of gibbs sampled distribution for my",type = "l")
abline(h=mean(my_vals), col="red") #Expected value

#Plotting sigma convergence
plot(sqrt(sigma_vals[10:nDraws]), ,main = "Plot of gibbs sampled distribution for sigma",type = "l")
abline(h=mean(sqrt(sigma_vals)), col="red") #Expected value

#Calculating inefficienty factor (IF)

#Autocorrelation for my
a_my = acf(my_vals)
inefficiency_factor_my = 1+2*sum(a_my$acf[-1])
inefficiency_factor_my

#Autocorrelation for my
a_sigma = acf(sqrt(sigma_vals))
inefficiency_factor_sigma = 1+2*sum(a_sigma$acf[-1])
inefficiency_factor_sigma

#Conclusion: As the plots suggest, we can't distinguish any convergence for the mean and standard deviation
#as we clearly can see that the values oscillate around the means (red line)
#Moreover, the inefficiency factors shows how correlated the draws are, as we can see, the results are
#relatively close to 1 which indicates low correlation

#b)
#1) Histogram or kernel density estimate of the daily preception y1,...,yn

hist(data) #Histogram with true data

#2) The resulting posterior predictive density p(y ̃|y1,...,yn) using the simulated posterior draws from (a)
predictions = rep(0,nDraws)
for (i in 1:nDraws) {
predictions[i] = rnorm(1,my_vals[i],sigma_vals[i])
}

hist(predictions)
#Conclusion: As we can see, the draws from the predictions differ slightly from the true distribution.
#The main difference is that the normal distribution have a bit longer tails on each end, however
#the distributions seem to be centered around the same values in both plots

