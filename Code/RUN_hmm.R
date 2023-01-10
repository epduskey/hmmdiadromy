# Run the simulation and the model

# Created: May 28, 2021
# Last modified: October 10, 2022 by EPD

# Load packages
library(lattice)
library(jagsUI)
library(rje)
library(coda)
library(mgcv)
library(CARBayesST)
library(mclogit)

# Set working directory
setwd(paste(mypath, "hmmdiadromy-main", sep = ""))

# Contents (ctrl-f):
#	I. Make some fish and let swim in your river
#	II. Make some habitat for your fish
#	III. Which sub-state were your fish in?
#	IV. Time to look for your fish
#	V. Run HMM
#	VI. Run SMM
#	VII. Run GAM
#	VIII. Run STCAR


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
statseq = move(init, trn.list, trnhm, nfish, ndays)

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
statseq$sub = substat

# Save all generated data
write.table(sand, "Data/sand_hmm.txt")
write.table(statseq$stat, "Data/statseq_hmm.txt")
write.table(statseq$hmm, "Data/hidseq_hmm.txt")
write.table(statseq$sub, "Data/subseq_hmm.txt")


########## IV. Time to look for your fish ########## 

# Would your fish have been detected by your equipment?  If not, cut those observations out
statdet = see(statseq$stat, p)
hypobs = cut(statseq$stat, statdet)	

# Where did you search with your equipment?  Can't find fish where you don't look
obs = search(rg, hypobs, substat, statdet)

# Store the final results in separate matrices
allobs = obs$actobs; actobs = allobs; actobs[,1] = 1
actsub = obs$actsub
actdet = obs$actdet

# Treat sub-states as contiguous pieces, and give labels of 1 - nstates*pcs
allsub = allobs * pcs - (pcs - actsub)

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
	lines(seq(ndays), ceiling(hmm.jags$q50$out[i,]/pcs), lwd = 2, lty = 2)
	polygon(x = c(seq(ndays),rev(seq(ndays))), y = ceiling(c(hmm.jags$q2.5$out[i,]/pcs,rev(hmm.jags$q97.5$out[i,]/pcs))), border = NA, col = rgb(0,0,0,0.5))
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
	lines(seq(ndays), ceiling(hmmsim.jags$q50$out[i,]/pcs), lwd = 2, lty = 2)
	polygon(x = c(seq(ndays),rev(seq(ndays))), y = ceiling(c(hmmsim.jags$q2.5$out[i,]/pcs,rev(hmmsim.jags$q97.5$out[i,]/pcs))), border = NA, col = rgb(0,0,0,0.5))
	points(seq(ndays), actobs[i,], pch = 9, lwd = 2)
	legend("topright", bty = 'n', c("Actual", "Median Estimated", "Observations", "95% CI"), lty = c(1,2,NA,1), pch = c(NA,NA,9,NA), lwd = c(2,2,2,10), col = c(1,1,1,rgb(0,0,0,0.5)))
	locator(1)
}


########## VII. Run GAM ##########

# Load model scripts
source("Code/kchoose.R")
source("Code/wcalc.R")

# Organize data into data frame (note multinom GAM requires first category to be '0')
hmmdf = data.frame(fish = factor(rep(seq(nfish),ndays)), 
	day = rep(seq(ndays),each=nfish), 
	state = c(allobs),
	sub = c(actsub))
write.table(hmmdf, "Data/hmmdf.txt")

# Run model and choose number of knots
modlist = kchoose(df = hmmdf, obs = allobs, sub = actsub, sand = sand, tol = 0.01)

# Plot edf against the number of knots
y = do.call(rbind, modlist$edf)
x = matrix(rep(seq(modlist$knots+1-nrow(y)+1,modlist$knots+1),nstates-1), nrow = nrow(y), ncol = ncol(y))
par(mfrow = c(1,1))
matplot(x, y, type = "l", xlab = "Knots", ylab = "EDF")
legend("bottomright", bty = "n", sapply(seq(9), toString), lty = seq(9), col = seq(9), title = "State")

# Run GAM model
hmm.gam = modlist

# Check the output
summary(hmm.gam$model)
plot(hmm.gam$model, pages = 1)
gam.check(hmm.gam$model)

# Convert observations to data frame with substates as factor
indobs = lapply(sapply(c(allsub),data.frame), factor, levels = seq(nstates*pcs))

# Count the number of fish observed in each substate on each day
ctind = lapply(indobs, table)

# Create data frame to house counts per substate
mcldf = data.frame(count = unlist(ctind,use.names=F), state = rep(seq(nstates),each=pcs), sub = seq(pcs))

# Add covariates
mcldf$day = rep(rep(seq(ndays),each=nfish), each = nstates*pcs)
mcldf$sand = sand[matrix(c(mcldf$sub,mcldf$state), ncol = 2)]
mcldf$ind = rep(seq(nfish), each = nstates*pcs)
mcldf$weight = c(sapply(c(allobs), wcalc, nstates = nstates, pcs = pcs))
mcldf$state.day.ind = interaction(mcldf$state, mcldf$day, mcldf$ind)
write.table(mcldf, "Data/hmmmcldf.txt")

# Run occupancy model
hmm.mcl = mclogit(count|state.day.ind ~ sand, weights = weight, random = ~0+sand|ind, data = mcldf)

# Check output
summary(hmm.mcl)

# Save models
hmm.gm = list(hmm.gam = hmm.gam, hmm.mcl = hmm.mcl)
save(object = hmm.gm, file = "Output/hmmgam.rda")


########## VIII. Run STCAR ##########

# Convert observations to data frame with states as factor
facobs = lapply(data.frame(allsub), factor, levels = seq(nstates*pcs))

# Count the number of fish observed in each state on each day
ctobs = lapply(facobs, table)

# Store fish counts per state in KxN table with K substates and N days
carobs = matrix(unlist(ctobs,use.names=F), nrow = nstates*pcs, ncol = ndays)

# Get the number of fish observed on each day and convert to KxN matrix
nfobs = matrix(colSums(carobs), nrow = nstates*pcs, ncol = ndays, byrow = T)

# Put sand covariate in appropriate KxN format
sandmat = matrix(c(sand), nrow = nstates*pcs, ncol = ndays)

# Choose distance over which adjacent states are neighbors
stepsize = abs(apply(actobs, 1, diff))
nbh = median(stepsize[stepsize > 0], na.rm = T)

# Create a data frame of centroids and states for river polygons
xt = data.frame(x = rep(0,nstates*pcs), y = seq(0,nstates*pcs-1), state = rep(seq(10),each=3))

# Create a simple neighborhood matrix
distance = as.matrix(dist(xt[,c(1,3)]))
W = array(0, c(nstates*pcs,nstates*pcs))
W[distance %in% seq(0,nbh)] = 1
diag(W) = 0

# Use AR(1) STCAR from CARBayesST; fit 3 chains
chain1 = ST.CARar(c(carobs) ~ c(sandmat), family = "binomial", W = W, trials = c(nfobs), 
	burnin = 100000, n.sample = 1600000, thin = 100, AR = 1)
chain2 = ST.CARar(c(carobs) ~ c(sandmat), family = "binomial", W = W, trials = c(nfobs), 
	burnin = 100000, n.sample = 1600000, thin = 100, AR = 1)
chain3 = ST.CARar(c(carobs) ~ c(sandmat), family = "binomial", W = W, trials = c(nfobs), 
	burnin = 100000, n.sample = 1600000, thin = 100, AR = 1)

# Combine chains and save
hmm.car = list(chain1, chain2, chain3)
save(object = hmm.car, file = "Output/hmmcar.rda")

# Gather all parameters for diagnostics
beta.samples = mcmc.list(chain1$samples$beta, chain2$samples$beta, chain3$samples$beta)
phi.samples = mcmc.list(chain1$samples$phi, chain2$samples$phi, chain3$samples$phi)
rho.samples = mcmc.list(chain1$samples$rho, chain2$samples$rho, chain3$samples$rho)
tau2.samples = mcmc.list(chain1$samples$tau2, chain2$samples$tau2, chain3$samples$tau2)
fitted.samples = mcmc.list(chain1$samples$fitted, chain2$samples$fitted, chain3$samples$fitted)

# Check Gelman-Rubin
rhbeta = gelman.diag(beta.samples, autoburnin = F, multivariate = F)
rhphi = gelman.diag(phi.samples, autoburnin = F, multivariate = F)
rhrho = gelman.diag(rho.samples, autoburnin = F, multivariate = F)
rhtau2 = gelman.diag(tau2.samples, autoburnin = F, multivariate = F)
rhfitted = gelman.diag(phi.samples, autoburnin = F, multivariate = F)
