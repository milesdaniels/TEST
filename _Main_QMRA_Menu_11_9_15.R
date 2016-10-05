######################################################################
#****************        QMRA MODEL                  ****************#
#****************   WITH BAYES AND MCMC              ****************#
#****************  Developed by Miles Daniels        ****************#
#****************     Copyright (C) 2015.            ****************#
#****************   DO NOT ** DISTRUBTUE WITHOUT     ****************#
#****************     AUTHOR'S CONSESNT             ****************#
######################################################################


################################################################################
# remove past objects, set your working directory
# and file export path, and load desired libraries:
# file format = "C:\\Location_A\\Location_B\\... File.csv"
################################################################################
rm(list=ls())
library(emdbook) # this library holds the beta-binomial functions
library(boa) # this library holds the convergence diagnostics, i.e. gelman and rubin test for convergence
setwd("/Users/mdaniels/Documents/Orissa_Project/QMRA/Modules")
output_path=( "/Users/mdaniels/Documents/Orissa_Project/QMRA/Outputs/" )

################################################################################
# Load code for the core model and functions:
################################################################################
{
# this function imports a .csv file to be used as the data set:
source("recov_data_fun_11_9_15.r",echo=F)

# this fuction proposes the next parameter value in the markov chain:
source("propose_fun_11_9_15.r",echo = F)

# this function sets up priors for recovery assuming a given distribtion (beata now):
source("prior_recov_fun_11_9_15.r",echo = F)

# this function calculated likelohood assuming a betabinomial distribution:
source("betabinom_likelihood_fun_11_9_15.r",echo = F)

# this function sets up priors for recovery assuming a given distribtion (beata now):
source("betabinom_likelihood_fun_11_9_15.r",echo = F)

# this function sets up priors for recovery assuming a given distribtion (gamma now):
source("gamma_likelihood_fun_11_9_15.r",echo = F)

# this function runs the code to calculate the distribution of the percent recovery from spiking trials:
source("fit_recov_dist_fun_11_9_15.r",echo = F)

# this function runs the code to calculate the distribution of observations:
source("fit_obs_dist_fun_11_9_15.r",echo = F)

# this function export results to a .csv:
source("export_fun_11_9_15.r",echo = F)

# this function export results to a .csv:
source("export_chain_fun_11_9_15.r",echo = F)

# this function tests of chains converge:
source("chain_converge_fun_11_20_15.r",echo = F)
}

################################################################################
# What type of system and pathogen you going to model?
# Select Recovery, Private, or Government.
# Select Cryptosporidium or Giardia:
# This will determine what data to use for modeling.
################################################################################
Type = "Recovery"
Pathogen = "Giardia"

# First fit spike data to fit percent recovery data ######
# Assuming a betabinomail distribution as in Ashbolt et al. 2007
# Incorporating method recovery uncertainties in stochastic estimates of raw water protozoan concentrations for QMRA
# Journal of water and health

# load the spike data set:
Recovery_Dat = f.import()

# set how jiggly the proposal distributions should be (by trial and error):
# very important to check trace plots to make sure parameter space
# is searched enough.
ishape1 = 3 # index of parameter
shape1.jiggle = 0.5 # search step 
shape1.lo = .5 # lower bound
shape1.hi = 10 # upper bound

ishape2 = 4 # index of parameter
shape2.jiggle = 0.5 # search step 
shape2.lo = .5 # lower bound
shape2.hi = 20 # uppper bound

# set the probability bounds of percent recovery to search:
pp.lo = 0.1 # lower bound
pp.hi = 0.9 # upper bound

# run the chain of a given length
# do this three times to determine convergecne:
recov_chain_1  = f.chain_recovery_fit(50000)
recov_chain_2  = f.chain_recovery_fit(50000)
recov_chain_3  = f.chain_recovery_fit(50000)

#Burn in chains by 10K
recov_chain_1  = recov_chain_1[1000:50000,]
recov_chain_2  = recov_chain_2[1000:50000,]
recov_chain_3  = recov_chain_3[1000:50000,]

BOA_Result = f.converge_chain(recov_chain_1,recov_chain_2, recov_chain_3 )

if(BOA_Result$mpsrf < 1.2)
{
  print("CHAINS HAVE CONVERGED")
  CONVERGE_TRUE = TRUE
}
if(BOA_Result$mpsrf > 1.2)
{
  CONVERGE_TRUE = FALSE
  stop()
}
# run g and r diagnostic:
# Note: mpsrf should be < 1.2 for convergence to be acheived.


# export results for records:
f.export(recov_chain_1)
# f.export(recov_chain_2)
# f.export(recov_chain_3)

# If covergence is achieved we can plot the distribution using one of the chains:
if(CONVERGE_TRUE == TRUE)
{
par(mfrow = c(2,2))
plot(seq(0,10),Recovery_Dat/200, pch = 16,col  = addTrans("purple",150), yaxt = 'n', xaxt = "n",
     ylab ="", xlab = "", axes = F, ylim = c(0, 1))
mean_re = mean(Recovery_Dat/200) 
lines(c(0,20), c(mean_re,mean_re), col  = addTrans("purple",150), lty = 2)
axis(side = 1, at = c(seq(0, 10)), labels = c(1, 2, 3, 4, 5, 6,7, 8, 9, 10,11))
axis(side = 2, at = c(0, .50, 1), labels = c(0, 50, 100))
mtext("Spike Experiment Trial", side = 1, line = 2.5, cex = 1)
mtext("% Spike Recovered", side = 2, line = 2.5, cex = 1)

alpha1 = density(recov_chain_1[,3],adjust = .25)
beta1 = density(recov_chain_1[,4],adjust = .25)
plot(alpha1$x, alpha1$y/max(alpha1$y), col  = addTrans("grey",200) , lwd =2, type = 'l')  

plot(beta1$x, beta1$y/max(beta$y),adjust = .25, col  = addTrans("grey",200) , lwd =2, type = 'l')  
plot_index_max_LL_recov = which.max( recov_chain_1[,2] )
plot(density(rbetabinom(n = 100000, size = 100, theta = 1, shape1 = mean(recov_chain_1[plot_index_max_LL_recov,3]), shape2= mean(recov_chain_1[plot_index_max_LL_recov,4]))))
}
# PART 2: Fit Observational Data: ####

# Now fit observational data to a gamma distribution. 
# Select what water source type you would like to work with:
# choose "Private_Well" or "Government_Well"
# choose "Cryptosporidium" or "Giardia"
#par(mfrow = c(3,2))

Pathogen = "Giardia"

Type = "Government_Well"

if(Type == "Government_Well" )
{
  Type_Plot = "Governments Wells"
}
if(Type == "Private_Well" )
{
  Type_Plot = "Private Wells"
}



# Select what year of data you would like to work with:
# choose 2012 or 2013.
Water_Year = 2012

# load the obervational data set:
Obs_Dat_n_20_L = f.import()

# for non-detects, assign half SLOD:

Vol_IMS_ml = 25
Limit_of_Detect = 1
Obs_Dat_n_20_L$Adjust_Zero_Low = Obs_Dat_n_20_L[,2]
Obs_Dat_n_20_L$Adjust_Zero_High = Obs_Dat_n_20_L[,2]
for(i in 1:length(Obs_Dat_n_20_L[,1]))
{
  if(Obs_Dat_n_20_L$Adjust_Zero_High[i] == 0)
  {
    Obs_Dat_n_20_L$Adjust_Zero_High[i] = (Obs_Dat_n_20_L[i,5]/Vol_IMS_ml)*(Limit_of_Detect/2)
  }
  if(Obs_Dat_n_20_L$Adjust_Zero_Low[i] == 0)
  {
    Obs_Dat_n_20_L$Adjust_Zero_Low[i] = (Limit_of_Detect/2)
  }
}

#extract just counts from dataframe
Obs_Dat_n_20_L = Obs_Dat_n_20_L$Adjust_Zero_High

Obs_Dat_n_20_L = subset(Obs_Dat_n_20_L, Obs_Dat_n_20_L < 500)

imu = 3
isd = 4

mu.lo = 1
mu.hi = 40

sd.lo = 0
sd.hi = 40

# run the chain of a given length

obs_chain_1  = f.chain_obs_fit(50000)
obs_chain_2  = f.chain_obs_fit(50000)
obs_chain_3  = f.chain_obs_fit(50000)


#Burn in chains by 10K
obs_chain_1  = obs_chain_1[1000:50000,]
obs_chain_2  = obs_chain_2[1000:50000,]
obs_chain_3  = obs_chain_3[1000:50000,]

BOA_Result = f.converge_chain(obs_chain_1,obs_chain_2, obs_chain_3 )

if(BOA_Result$mpsrf < 1.2)
{
  print("CHAINS HAVE CONVERGED")
  CONVERGE_TRUE = TRUE
}
if(BOA_Result$mpsrf > 1.2)
{
  CONVERGE_TRUE = FALSE
  stop()
}
# run g and r diagnostic:
# Note: mpsrf should be < 1.2 for convergence to be acheived.

max_index = which.max(obs_chain_1[,2])
obs_chain_1[max_index,3]
obs_chain_1[max_index,4]

plot(density(rgamma(n = 100000,  shape = obs_chain_1[max_index,4], scale= obs_chain_1[max_index,3]/obs_chain_1[max_index,4] )))
points(Obs_Dat_n_20_L, rep(0, length(Obs_Dat_n_20_L)), col = 'red')
#f.export.chain(obs_chain_1)

if(CONVERGE_TRUE == TRUE)
{
  par(mfrow = c(3,1))
  par(mai = c(0.6, .6, .1, .2))

max_index_ch1 = which.max(obs_chain_1[,2])
obs_chain_1_den_scale = density(obs_chain_1[,3]/obs_chain_1[,4], adjust= .5)
obs_chain_1_den_shape = density(obs_chain_1[,4], adjust= .5)
max_index_ch2 = which.max(obs_chain_2[,2])
obs_chain_2_den_scale = density(obs_chain_2[,3]/obs_chain_2[,4], adjust= .5)
obs_chain_2_den_shape = density(obs_chain_2[,4], adjust= .5)
max_index_ch3 = which.max(obs_chain_3[,2])
obs_chain_3_den_scale = density(obs_chain_3[,3]/obs_chain_3[,4], adjust= .5)
obs_chain_3_den_shape = density(obs_chain_3[,4], adjust= .5)


plot_index_max_LL_obs = which.max( obs_chain_1[,2] )
Dat = rgamma(n = 1000000,  shape = obs_chain_1[max_index,4], scale= obs_chain_1[max_index,3]/obs_chain_1[max_index,4] )
Dat_Den = density(Dat)
plot(Dat_Den$x, Dat_Den$y/max(Dat_Den$y), 
     main = '', xlim = c(0, max(Obs_Dat_n_20_L)), type = 'l',ylab  = "Probability", 
     xlab = print(paste("Conc. ", Pathogen, " in ",Type_Plot, " during ", Water_Year, "(# per 20L)")),yaxt = 'n', cex.lab =1.5)
axis(side = 2, at = c(0, .5, 1), labels = c(0,.5, 1))
points(Obs_Dat_n_20_L, rep(0, length(Obs_Dat_n_20_L)), col = 'black', pch =16)

legend(
  x      = max(Obs_Dat_n_20_L) - (.5*max(Obs_Dat_n_20_L)),
  y      = .8,
  legend = 'Observations',
  pch   = c( 16),lty=c(0), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1.5, lwd = c(2,2), merge = F)
legend(
  x      = max(Obs_Dat_n_20_L)- (.5*max(Obs_Dat_n_20_L)),
  y      = .6,
  legend = 'Model fit',
  ,lty=c(1), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1.5, lwd = c(2,2), merge = F)



plot(obs_chain_1_den_scale$x,obs_chain_1_den_scale$y/max(obs_chain_1_den_scale$y) , main = '', 
     xlab = "Posterior PDF of scale parameter for gamma distribution", type = 'l', lty = 1, ylab  = "Probability", yaxt = 'n', cex.lab =1.5)
lines(obs_chain_2_den_scale$x,obs_chain_2_den_scale$y/max(obs_chain_1_den_scale$y), lty = 2)
lines(obs_chain_3_den_scale$x,obs_chain_3_den_scale$y/max(obs_chain_1_den_scale$y), lty = 3)
axis(side = 2, at = c(0, .5, 1), labels = c(0,.5, 1))

legend(
  x      = mean(obs_chain_1_den_scale$x),
  y      = .9,
  legend = print(paste("Chain 1 = ", round(obs_chain_1[max_index_ch1,3]/obs_chain_1[max_index_ch1,4],3)))
  ,lty=c(1), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1, lwd = c(2,2), merge = F)

legend(
  x      = mean(obs_chain_1_den_scale$x),
  y      = .8,
  legend = print(paste("Chain 2 = ", round(obs_chain_2[max_index_ch2,3]/obs_chain_2[max_index_ch2,4],3)))
  ,lty=c(2), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1, lwd = c(2,2), merge = F)

legend(
  x      = mean(obs_chain_1_den_scale$x),
  y      = .7,
  legend = print(paste("Chain 3 = ", round(obs_chain_3[max_index_ch3,3]/obs_chain_3[max_index_ch3,4],3)))
  ,lty=c(3), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1, lwd = c(2,2), merge = F)




plot(obs_chain_1_den_shape$x,obs_chain_1_den_shape$y/max(obs_chain_1_den_shape$y) , main = '', 
     xlab = "Posterior PDF of shape parameter for gamma distribution", type = 'l', lty = 1,ylab  = "Probability", yaxt = 'n', cex.lab =1.5)
lines(obs_chain_2_den_shape$x,obs_chain_2_den_shape$y/max(obs_chain_2_den_shape$y), lty = 2)
lines(obs_chain_3_den_shape$x,obs_chain_3_den_shape$y/max(obs_chain_3_den_shape$y), lty = 3)
axis(side = 2, at = c(0, .5, 1), labels = c(0,.5, 1))

legend(
  x      = mean(obs_chain_1_den_shape$x)+(.1*mean(obs_chain_1_den_shape$x)),
  y      = .9,
  legend = print(paste("Chain 1 = ", round(obs_chain_1[max_index_ch1,4],3)))
  ,lty=c(1), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1, lwd = c(2,2), merge = F)

legend(
  x      = mean(obs_chain_1_den_shape$x)+(.1*mean(obs_chain_1_den_shape$x)),
  y      = .8,
  legend = print(paste("Chain 2 = ", round(obs_chain_2[max_index_ch2,4],3)))
  ,lty=c(2), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1, lwd = c(2,2), merge = F)

legend(
  x      = mean(obs_chain_1_den_shape$x)+(.1*mean(obs_chain_1_den_shape$x)),
  y      = .7,
  legend = print(paste("Chain 3 = ", round(obs_chain_3[max_index_ch3,4],3)))
  ,lty=c(3), col = "black" ,
  , bty = 'n', ncol = 1, cex = 1, lwd = c(2,2), merge = F)



# export results for records:
path = "/Users/mdaniels/Documents/Orissa_Project/QMRA/Outputs/Chains_Fit_Gamma_Conc/"
write.csv( obs_chain_1,  paste( path, Pathogen, "_", Type,"_", Water_Year, "_Ran_2_4_16_","Remove_Outlier_Chain_1_Fit_Gamma_Dist_from_",  length(obs_chain_1[,1]), "_Sims.csv", sep="" ))
write.csv( obs_chain_2,  paste( path, Pathogen, "_", Type,"_", Water_Year, "_Ran_2_4_16_","Remove_Outlier_Chain_2_Fit_Gamma_Dist_from_",  length(obs_chain_1[,1]), "_Sims.csv", sep="" ))
write.csv( obs_chain_3,  paste( path, Pathogen, "_", Type,"_", Water_Year, "_Ran_2_4_16_","Remove_Outlier_Chain_3_Fit_Gamma_Dist_from_",  length(obs_chain_1[,1]), "_Sims.csv", sep="" ))


}


plot(density(obs_chain_1[,3]/obs_chain_1[,4]), main = '')



f.export.chain(obs_chain_1)
f.export.chain(obs_chain_2)
f.export.chain(obs_chain_3)

}




n= 10000
data = matrix(ncol = 3, nrow = n)
#adjust_conc = rep(0, n)
for(i in 1:n)
{
  data[i,1] = rgamma(n = 1,  shape = mean(obs_chain_1[plot_index_max_LL_obs,4]), scale= mean(obs_chain_1[plot_index_max_LL_obs,3])/ mean(obs_chain_1[plot_index_max_LL_obs,4]))
  
  data[i,2] = rbetabinom(n = 1, size = 100, theta = 1, shape1 = (recov_chain_1[plot_index_max_LL_recov,3]), shape2= (recov_chain_1[plot_index_max_LL_recov,4]))
  data[i,3] = data[i,1] * (100/data[i,2])
}

par(mfrow = c(1,1))
plot(density( data[,1] ))
lines(density( data[,3] ), col= 'red')


# now make uniform distribution for volume of water
# injested per day. 
# using data from: 
# estimates for children 1-5 years old is 1.3 - 1.7 Liters per day.

Water_Intake 

n= 365
Infection = matrix(ncol = 7, nrow = n)
for(i in 1:length(Infection[,1]))
{
conc_20_L = rgamma(n = 1,  shape = mean(obs_chain_1[plot_index_max_LL_obs,4]), scale= mean(obs_chain_1[plot_index_max_LL_obs,3])/ mean(obs_chain_1[plot_index_max_LL_obs,4]))
Infection[i,1] = conc_20_L/20
Infection[i,2] = rbetabinom(n = 1, size = 100, theta = 1, shape1 = (recov_chain_1[plot_index_max_LL_recov,3]), shape2= (recov_chain_1[plot_index_max_LL_recov,4]))
Infection[i,3] = Infection[i,1] * (100/Infection[i,2])
Infection[i,4] = runif(n = 1, min = 1.3, max = 1.7)
Infection[i,5] = Infection[i,3]*Infection[i,4]
#Infection[i,4] = pexp(q = Infection[i,1]*Infection[i,2], rate = runif(n = 1, min = 0.00419, max = 0.09))
Infection[i,6] = pexp(q = Infection[i,1]*Infection[i,2], rate = 0.09)
Infection[i,7] = rbinom(n = 1, size = 1, prob = Infection[i,6])
}

round(Infection, 2)
hist(Infection[,3])
as.integer(Infection)


# first attempt at making infection model
# date 10/4/15

mat = matrix(ncol = 100, nrow = 100)
mat[1:100,1:100] = 0 
prob = .1

for(i in 1:100)
{
  for(j in 1:95)
  {
    if(j == 1)
    {
      p = as.integer(rbinom(1, 1, runif(n = 1, min = .1, max = .8)))
      if(p == 1)
      {
        mat[i,j] = 1
        mat[i,j+1] = 1
        mat[i,j+2] = 1
      }
      else  mat[i,j] = 0
    }
    if(j > 1)
    {
      p = as.integer(rbinom(1, 1, runif(n = 1, min = .1, max = .4)))
      if(mat[i,j] == 1)
      {
        
      }
      else
        if(p == 1)
        {
          mat[i,j] = 1
          mat[i,j+1] = 1
          mat[i,j+2] = 1
        }
      else mat[i,j] = p
    }
  }
}



sums  = matrix(ncol = 100, nrow = 1)
for(i in 1:100)
{
  sums[,i] = sum(mat[,i])
}

plot(seq(1, length(sums[1,]), 1),sums, type = "l")


