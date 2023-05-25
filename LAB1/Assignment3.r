#Assignment 3

#This exercise is concerned with directional data. The point is to show you that the posterior distribution
#for somewhat weird models can be obtained by plotting it over a grid of values. 
#The data points are observed wind directions at a given location on ten di􏰀erent days. 
#The data are recorded in degrees: (20, 314, 285, 40, 308, 314, 299, 296, 303, 326),
#where North is located at zero degrees (see Figure 1 on the next page, where the angles are measured clockwise). 
#To 􏰁t with Wikipedia's description of probability distributions for circular data we convert 
#the data into radians −π ≤ y ≤ π . The 10 observations in 
#radians are (−2.79, 2.33, 1.83, −2.44, 2.23, 2.33, 2.07, 2.02, 2.14, 2.54).
#Assume that these data points conditional on (μ,κ) are 
#independent observations from the following von Mises distribution:

#where I0(κ) is the modi􏰁ed Bessel function of the 􏰁rst kind of order zero [see ?besselI in R]. 
#The parameter μ (−π ≤ μ ≤ π) is the mean direction and κ > 0 is called the concentration parameter.
#Large κ gives a small variance around μ, and vice versa. Assume that μ is known to be 2.4. Let κ ∼ Exponential(λ = 0.5)
#a priori, where λ is the rate parameter of the exponential distribution (so that the mean is 1/λ).


###Starting assignment###

#a) Derive the expression for what the posterior p(κ|y, μ) is proportional to. Hence, derive the function f (κ) such
#that p(κ|y, μ) ∝ f (κ). Then, plot the posterior distribution of κ for the wind direction data over a 􏰁ne grid of κ
#values. [Hint: you need to normalize the posterior distribution of κ so that it integrates to one.]


#Derive the posterior p(k|y,my) is proportional to p(y|k,my)*p(k|my) is proportional to p(y|k,my)*p(k)
#p(k|my) is proportional to p(k) since k is independent of my

#Data and know parameters
degrees = c(-2.79,2.33,1.83,-2.44,2.23,2.33,2.07,2.02,2.14,2.54) #Directions in radians
my = 2.4 #Given my value
lambda = 0.5 #Given lambda value

#Function for von Mises distribution
vonMises = function(kappa, y, my) {
  likelihood = 1
  for (data in y) {
  likelihood = likelihood*exp(kappa*cos(data-my))/(2*pi*besselI(kappa,0))
  }
  return(likelihood)
}

#Function for exponential distribution
exponential = function(tetha, data) {
  return(tetha*exp(-tetha*data))
}

kappa = seq(0.1,10,0.1) #Kappa values to be used in distribution
VMProbs = rep(0, length(kappa)) #List to store probs from von Mis.
expProbs = rep(0, length(kappa)) #List to store probs from exp. distr

#Loop to generate all probs from given distribution
for(i in 1:length(kappa)) {
VMProbs[i] = vonMises(kappa[i],degrees,my)
expProbs[i] = exponential(lambda,kappa[i])
}

posterior = expProbs*VMProbs #Calculating posterior
postIntegral = sum(posterior*0.1) #Calculating the integral for the posterior
normPosterior = posterior/postIntegral #Normalizing the posterior to integrate to 1
plot(kappa,normPosterior, type="l", xlab="Kappa values", ylab="Likelihood") #Plotting the normalized posterior

#b) #Find the (approximate) posterior mode of κ from the information in a).


#As seen in the plot, the approximate posterior mode for kappa is when kappa is between 2 & 4
#A good estimation judging by the graph seems to be around ~2.8
