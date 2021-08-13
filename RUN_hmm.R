# Run the simulation and the model

# Created: May 28, 2021
# Last modified: August 13, 2021 by EPD

# Load packages
library(lattice)
library(jagsUI)
library(rje)

# Set working directory
setwd(paste(mypath, "HMM", sep = ""))

# Contents (ctrl-f):
#	I. Make some fish and let swim in your river
#	II. Make some habitat for your fish
#	III. Which sub-state were your fish in?
#	IV. Time to look for your fish
#	V. Run HMM
#	VI. Run SMM


########## I. Make some fish and let them swim in your river ##########

# Load simulation and parameter scripts
source("Code/simulation_hmm.R")
source("Code/parameters_hmm.R")

# Set seed to reproduce results from Duskey et al 2021
set.seed(8074)
#set.seed(8235)
#set.seed(2457)
#set.seed(7808)
#set.seed(3197)
#set.seed(2352)
#set.seed(3980)
#set.seed(2581)
#set.seed(4172)
#set.seed(9586)

# Get a sequence of locations for your fish
statseq = move(init, trn, trnhm, nfish, ndays)
write.table(statseq$stat, "Data/statseq_hmm.txt")

# Plot sequence of locations for your fish
#for(i in 1:nfish) {
#	xl = "Day"
#	yl = "State"
#	ml = paste("Fish", i, " - click anywhere to advance to the next fish")
#	plot(statseq$stat[i,], type = 'l', lwd = 2, ylim = c(1,nstates), xlab = xl, ylab = yl, main = ml)
#	locator(1)
#}

# Plot sequence of behavior for your fish
#for(i in 1:nfish) {
#	xl = "Day"
#	yl = "Behavior (1 = up-migration; 2 = down-migration"
#	ml = paste("Fish", i, " - click anywhere to advance to the next fish")
#	plot(statseq$hmm[i,], type = 'l', lwd = 2, ylim = c(1,2), xlab = xl, ylab = yl, main = ml)
#	locator(1)
#}


########## II. Make some habitat for your fish ##########

# Puts a random proportion of sand in each sub-state
sand = hack(pcs, nstates)

# Plots the proportion sand in each sub-state
levelplot(sand, xlab = "Sub-state", ylab = "State", main = "Proportion sand in each sub-state", col.regions = hcl.colors(100, "BrwnYl", rev = F))

# Based on how much your fish like sand, they will have a certain probability of occupying each sub-state
occ = apply(sand, 2, pocc, a0 = a0, a1 = a1)	

# Plot the probability of occupying each sub-state
levelplot(occ, xlab = "Sub-state", ylab = "State", main = "Probability of occupying each sub-state", col.regions = hcl.colors(100, "Magenta", rev = T))


########## III. Which sub-state were your fish in? ##########

# Places your fish in sub-states based on the occupancy probability calculated in step II
substat = where(statseq$stat, occ)


########## IV. Time to look for your fish ########## 

# Would your fish have been detected by your equipment?  If not, cut those observations out
statdet = see(statseq$stat, p)
hypobs = cut(statseq$stat, statdet)	

# Where did you search with your equipment?  Can't find fish where you don't look
obs = search(rg, hypobs, substat, statdet)

# Store the final results in separate matrices
actobs = obs$actobs; actobs[,1] = 1
actsub = obs$actsub
actdet = obs$actdet

# Save observations
write.table(actobs, "Data/acthmm.txt")


########## V. Run HMM model ##########

# Load model script
source("Code/model_hmm.R")

# Choose run details
nadapt = 1000
nburn = 10000
nchains = 3
niter = 110000
nthin = 10
inc = 10000
par = T

# Run JAGS! 
hmm.jags = autojags(data, inits = inits, params, "Code/fishhmm_dn.jags", iter.increment = inc,
				n.adapt = nadapt, n.burnin = nburn, n.chains = nchains, max.iter = niter, n.thin = nthin,
				parallel = par, save.all.iter = T, factories = "base::Finite sampler FALSE")
print(hmm.jags)
save(object = hmm.jags, file = "Output/hmmjags.rda")

# Actual up-migration transition matrix
trnup

# Estimated up-migration transition matrix
round(hmm.jags$mean$trn[1,,], digits = 2)

# Actual up-migration transition matrix
trndn

# Estimated up-migration transition matrix
round(hmm.jags$mean$trn[2,,], digits = 2)

# Actual hidden transition matrix
trnhm

# Estimated hidden transition matrix
round(hmm.jags$mean$hmm, digits = 2)

# Plot actual and estimated trajectories
for(i in 1:nfish) {
	xl = "Day"
	yl = "State"
	ml = paste("Fish", i, " - click anywhere to advance to the next fish")
	plot(statseq$stat[i,], type = 'l', lwd = 2, ylim = c(1,nstates), xlab = xl, ylab = yl, main = ml)
	lines(seq(ndays), hmm.jags$q50$s[i,], lwd = 2, lty = 2)
	polygon(x = c(seq(ndays),rev(seq(ndays))), y = c(hmm.jags$q2.5$s[i,],rev(hmm.jags$q97.5$s[i,])), border = NA, col = rgb(0,0,0,0.5))
	points(seq(ndays), actobs[i,], pch = 9, lwd = 2)
	legend("topright", bty = 'n', c("Actual", "Median Estimated", "Observations", "95% CI"), lty = c(1,2,NA,1), pch = c(NA,NA,9,NA), lwd = c(2,2,2,10), col = c(1,1,1,rgb(0,0,0,0.5)))
	locator(1)
}


########## VI. Run SMM ##########

# Load model script
source("Code/model_smm.R")

# Choose run details
nadapt = 1000
nburn = 10000
nchains = 3
niter = 110000
nthin = 10
inc = 10000
par = T

# Run JAGS! 
hmmsim.jags = autojags(data, inits = inits, params, "Code/fishsmm_dn.jags", iter.increment = inc,
				n.adapt = nadapt, n.burnin = nburn, n.chains = nchains, max.iter = niter, n.thin = nthin,
				parallel = par, save.all.iter = T, factories = "base::Finite sampler FALSE")
print(hmmsim.jags)
save(object = hmmsim.jags, file = "Output/hmmsimjags.rda")

# Actual up-migration transition matrix
trnup

# Actual up-migration transition matrix
trndn

# Estimated up-migration transition matrix
round(hmmsim.jags$mean$trn, digits = 2)

# Plot actual and estimated trajectories
for(i in 1:nfish) {
	xl = "Day"
	yl = "State"
	ml = paste("Fish", i, " - click anywhere to advance to the next fish")
	plot(statseq$stat[i,], type = 'l', lwd = 2, ylim = c(1,nstates), xlab = xl, ylab = yl, main = ml)
	lines(seq(ndays), hmmsim.jags$q50$s[i,], lwd = 2, lty = 2)
	polygon(x = c(seq(ndays),rev(seq(ndays))), y = c(hmmsim.jags$q2.5$s[i,],rev(hmmsim.jags$q97.5$s[i,])), border = NA, col = rgb(0,0,0,0.5))
	points(seq(ndays), actobs[i,], pch = 9, lwd = 2)
	legend("topright", bty = 'n', c("Actual", "Median Estimated", "Observations", "95% CI"), lty = c(1,2,NA,1), pch = c(NA,NA,9,NA), lwd = c(2,2,2,10), col = c(1,1,1,rgb(0,0,0,0.5)))
	locator(1)
}