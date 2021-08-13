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
#	V. Run SMM
#	VI. Run HMM


########## I. Make some fish and let them swim in your river ##########

# Load simulation and parameter scripts
source("Code/simulation_smm.R")
source("Code/parameters_smm.R")

# Set seed to reproduce results from Duskey et al 2021
set.seed(8235)
#set.seed(499)
#set.seed(2344)
#set.seed(3033)
#set.seed(9845)
#set.seed(294)
#set.seed(8627)
#set.seed(5493)
#set.seed(3621)
#set.seed(1546)

# Get a sequence of locations for your fish
statseq = move(init, trn, nfish, ndays)
write.table(statseq, "Data/statseq_smm.txt")

# Plot sequence of locations for your fish
# for(i in 1:nfish) {
	# xl = "Day"
	# yl = paste("Fish", i)
	# ml = paste("Fish", i, " - click anywhere to advance to the next fish")
	# plot(statseq[i,], type = 'l', lwd = 2, ylim = c(1,nstates), xlab = xl, ylab = yl, main = ml)
	# locator(1)
# }


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
substat = where(statseq, occ)


########## IV. Time to look for your fish ########## 

# Would your fish have been detected by your equipment?  If not, cut those observations out
statdet = see(statseq, p)
hypobs = cut(statseq, statdet)	

# Where did you search with your equipment?  Can't find fish where you don't look
obs = search(rg, hypobs, substat, statdet)

# Store the final results in separate matrices
actobs = obs$actobs; actobs[,1] = 1
actsub = obs$actsub
actdet = obs$actdet

# Save observations
write.table(actobs, "Data/actmat.txt")


########## V. Run SMM model ##########

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
sim.jags = autojags(data, inits = inits, params, "Code/fishsmm_dn.jags", iter.increment = inc,
				n.adapt = nadapt, n.burnin = nburn, n.chains = nchains, max.iter = niter, n.thin = nthin,
				parallel = par, save.all.iter = T, factories = "base::Finite sampler FALSE")
print(sim.jags)
save(object = sim.jags, file = "Output/simjags.rda")

# Actual up-migration transition matrix
trn

# Estimated up-migration transition matrix
round(sim.jags$mean$trn, digits = 2)

# Plot actual and estimated trajectories
for(i in 1:nfish) {
	xl = "Day"
	yl = "State"
	ml = paste("Fish", i, " - click anywhere to advance to the next fish")
	plot(statseq[i,], type = 'l', lwd = 2, ylim = c(1,nstates), xlab = xl, ylab = yl, main = ml)
	lines(seq(ndays), sim.jags$q50$s[i,], lwd = 2, lty = 2)
	polygon(x = c(seq(ndays),rev(seq(ndays))), y = c(sim.jags$q2.5$s[i,],rev(sim.jags$q97.5$s[i,])), border = NA, col = rgb(0,0,0,0.5))
	points(seq(ndays), actobs[i,], pch = 9, lwd = 2)
	legend("topright", bty = 'n', c("Actual", "Median Estimated", "Observations", "95% CI"), lty = c(1,2,NA,1), pch = c(NA,NA,9,NA), lwd = c(2,2,2,10), col = c(1,1,1,rgb(0,0,0,0.5)))
	locator(1)
}


########## V. Run HMM ##########

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
ndhmm.jags = autojags(data, inits = inits, params, "Code/fishhmm_dn.jags", iter.increment = inc,
				n.adapt = nadapt, n.burnin = nburn, n.chains = nchains, max.iter = niter, n.thin = nthin,
				parallel = par, save.all.iter = T, factories = "base::Finite sampler FALSE")
print(ndhmm.jags)
save(object = ndhmm.jags, file = "Output/ndhmmjags.rda")

# Actual transition matrix
trn

# Estimated up-migration transition matrix
round(ndhmm.jags$mean$trn[1,,], digits = 2)

# Estimated up-migration transition matrix
round(ndhmm.jags$mean$trn[2,,], digits = 2)

# Estimated hidden transition matrix
round(ndhmm.jags$mean$hmm, digits = 2)

# Plot actual and estimated trajectories
for(i in 1:nfish) {
	xl = "Day"
	yl = "State"
	ml = paste("Fish", i, " - click anywhere to advance to the next fish")
	plot(statseq[i,], type = 'l', lwd = 2, ylim = c(1,nstates), xlab = xl, ylab = yl, main = ml)
	lines(seq(ndays), ndhmm.jags$q50$s[i,], lwd = 2, lty = 2)
	polygon(x = c(seq(ndays),rev(seq(ndays))), y = c(ndhmm.jags$q2.5$s[i,],rev(ndhmm.jags$q97.5$s[i,])), border = NA, col = rgb(0,0,0,0.5))
	points(seq(ndays), actobs[i,], pch = 9, lwd = 2)
	legend("topright", bty = 'n', c("Actual", "Median Estimated", "Observations", "95% CI"), lty = c(1,2,NA,1), pch = c(NA,NA,9,NA), lwd = c(2,2,2,10), col = c(1,1,1,rgb(0,0,0,0.5)))
	locator(1)
}
