# Standard Markov model

cat("model{
	
	# Priors for transition matrix
	for(i in 1:nstates) {
		for(j in 1:nstates) {
			a[i,j] ~ dgamma(1,1)	
		}
		trn[i,1:nstates] ~ ddirich(a[i,1:nstates] + 0.1)
	}
	
	# Estimate states from transition matrix
	for(i in 1:nfish) {
		for(j in 2:ndays) {
			s[i,j] ~ dcat(trn[s[i,j-1], ])
		}
	}
	
	# Estimate detection probability
	for(i in 1:nfish) {
		for(j in 1:ndays) {
			search[i,j] <- ifelse(rg[s[i,j],j] == 1, 1, 0)
			det[i,j] ~ dbin(p*search[i,j], 1)
		}
	}
	
	# Detection prior
	p ~ dunif(0,1)
	
	# Estimate occupancy at the substate scale
	for(i in 1:nstates) {
		for(j in 1:pcs) {
			# Log link for occupancy probability as a function of sand
			log(w[i,j]) <- (a0 + a1*sand[i,j])	
		}
		# Draw from dirichlet
		pocc[i,1:pcs] ~ ddirich(w[i,1:pcs] + 0.1)	
	}
	for(i in 1:nfish) {
		for(j in 1:ndays) {
			# Grid scale occupancy drawn from categorical
			grid[i,j] ~ dcat(pocc[s[i,j],])
		}
	}
	
	# Logistic regression priors - occupancy
	a0 ~ dnorm(0, 0.001)
	a1 ~ dnorm(0, 0.001)
	
}", file = "Code/fishsmm_dn.jags")

# Set up data, initial values, and parameters for jags
data = list(nfish = nfish,
			ndays = ndays,
			nstates = nstates,
			det = actdet,
			grid = actsub,
			sand = t(sand),
			s = actobs,
			rg = rgmat,
			pcs = pcs)

# Load function tmat to calculate via MLE
source("Code/tmat.R")	
	
# Apply tmat to all observations
trnest = tmat(actobs, 1, nstates, est = T)		

# Initial values
inits = function() {		
	list(p = 0.5,
		a0 = 0, a1 = 0, 
		trn = trnest, 
		pocc = matrix(1/pcs, nrow = nstates, ncol = pcs, byrow = T))
}

# Tracked parameters
params = c("s", 
		"a0",
		"a1",
		"trn",
		"p",
		"pocc")