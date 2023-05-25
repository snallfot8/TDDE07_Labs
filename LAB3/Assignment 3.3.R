#Time series models in Stan
#Write a function in R that simulates data from the AR(1)-process
#xt=μ+φ(xt−1−μ)+εt,εt∼N0 , for given values of μ, φ and σ2. 
#Start the process at x1 = μ and then simulate values for xt for t = 2, 3 . . . , T 
#and return the vector x1:T containing all time points. 
#Use μ = 13, σ2 = 3 and T = 300 and look at some different realizations (simulations) of x1:T 
#for values of φ between −1 and 1 (this is the interval of φ where the AR(1)-process is stationary). 
#Include a plot of at least one realization in the report. What effect does the value of φ have on x1:T ?
#a)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Initial values
my=13
sigmasq=3
T=300
x1 = my
phi = seq(-1,1,0.25)

#Autoregression function
ar_function = function(my,sigmasq,T,x1,phi){
  result = rep(0,T)
  result[1] = x1
  for(i in 2:T) {
  epsilon = rnorm(1,0,sigmasq)
  result[i] = my + phi*(result[i-1]-my) + epsilon
  }
  return(result)
}

#Calculating output
result = matrix(0,T,length(phi))
for(i in 1:length(phi)) {
  result[,i] = ar_function(my,sigmasq,T,x1,phi[i])
}

#Plotting result for phi=-1 and phi=1
plot(seq(1,300,1),result[,1], xlab="Order", ylab="Realization", main="Realization for phi = -1",type="l")
plot(seq(1,300,1),result[,3], xlab="Order", ylab="Realization", main="Realization for phi = -0.5",type="l")
plot(seq(1,300,1),result[,7], xlab="Order", ylab="Realization", main="Realization for phi = 0.5",type="l")
plot(seq(1,300,1),result[,9], xlab="Order", ylab="Realization", main="Realization for phi = 1",type="l")

#Conclusion: Increased phi values reduces oscilations. We can see clear oscilations for phi=-1
#and that the output becomes more and more correlated with increased phi values

#b)
#Use your function from a) to simulate two AR(1)-processes, x1:T with φ = 0.2 and y1:T with φ = 0.95. 
#Now, treat your simulated vectors as synthetic data, and treat the values of μ, φ and σ2 as unknown parameters. 
#Implement Stan code that samples from the posterior of the three parameters, 
#using suitable non-informative priors of your choice.
#[Hint: Look at the time-series models examples in the Stan user's guide/reference manual, 
#and note the dierent parameterization used here.]

#i. Report the posterior mean, 95% credible intervals and the number of effective posterior samples
#for the three inferred parameters for each of the simulated AR(1)-process. Are you able to estimate the true values?
#ii. For each of the two data sets, evaluate the convergence of the samplers
#and plot the joint posterior of μ and φ.Comments?

#Drawing x1:T values using phi=0.2
x=rep(0,T)
x=ar_function(my,sigmasq,T,x1,0.2)

#Drawing y1:T values using phi=0.95
y=rep(0,T)
y=ar_function(my,sigmasq,T,x1,0.95)

#Implementing Stan model

StanModel = '
data {
  int<lower=0> N;
  vector[N] z;
}
parameters {
  real mu;
  real phi;
  real<lower=0> sigma2;
}
model {
  for(i in 2:N){
    z[i] ~ normal(mu+phi*(z[i-1]-mu),sqrt(sigma2));
  }
}'

#For X
data <- list(N=T, z=x)
fit_x <- stan(model_code=StanModel,data=data)

#For Y
data <- list(N=T, z=y)
fit_y <- stan(model_code=StanModel,data=data)

#Extract posterior draws
post_x = extract(fit_x)
post_y = extract(fit_y)

#i)
#AR-process for X
#Calculating mean and credible interval for all parameters
#Can be achived with print(fit_x) too
mean_mu = mean(post_x$mu)
cred_int_mu = quantile(post_x$mu,probs=c(0.025,0.975))
print(mean_mu)
print(cred_int_mu)

mean_phi= mean(post_x$phi)
cred_int_phi = quantile(post_x$phi,probs=c(0.025,0.975))
print(mean_phi)
print(cred_int_phi)

mean_sigma2 = mean(post_x$sigma2)
cred_int_sigma2 = quantile(post_x$sigma2,probs=c(0.025,0.975))
print(mean_sigma2)
print(cred_int_sigma2)

#AR-process for Y
#Calculating mean and credible interval for all parameters
mean_mu = mean(post_y$mu)
cred_int_mu = quantile(post_y$mu,probs=c(0.025,0.975))
print(mean_mu)
print(cred_int_mu)

mean_phi= mean(post_y$phi)
cred_int_phi = quantile(post_y$phi,probs=c(0.025,0.975))
print(mean_phi)
print(cred_int_phi)

mean_sigma2 = mean(post_y$sigma2)
cred_int_sigma2 = quantile(post_y$sigma2,probs=c(0.025,0.975))
print(mean_sigma2)
print(cred_int_sigma2)

#Plotting using print
print(fit_x)
print(fit_y)

#Conclusion: As the values in the credible interval are more narrow for the one where phi=0.2
#were used, we get a better understanding of the true values for this one compared to the one where
#phi=0.95 were used. For this one, the values are really wide in the credible intervals

#ii)
#Joint posterior for X
plot(post_x$mu[1000:2000],post_x$phi[1000:2000],xlab="My vals",ylab="Phi vals", main="Joint posterior of mu and phi (X)")

#Joint posterior for X
plot(post_y$mu[1000:2000],post_y$phi[1000:2000],xlab="My vals",ylab="Phi vals", main="Joint posterior of mu and phi (Y)")


#Conclusion: We can distinguish a convergence for the first sampler. However it is more difficult to distinguish a
#convergence for the second sampler as shown in the figure. This follows the previous reasoning, that the 
#second sampler have a wider span for its values compared to the first sampler