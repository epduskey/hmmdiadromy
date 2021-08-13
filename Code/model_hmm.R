# HMM model

cat("model{
	
	# Priors for transition matrices
	for(i in 1:nstates) {
		for(j in 1:nstates) {
			aup[i,j] ~ dgamma(1,1)
			adn[i,j] ~ dgamma(1,1)
		}
		trn[1,i,1:nstates] ~ ddirich(aup[i,] + 0.1)
		trn[2,i,1:nstates] ~ ddirich(adn[i,] + 0.1)
	}
	
	# Hidden Markov portion of the model
	for(i in 1:2) {
		for(j in 1:2) {
			ahmm[i,j] ~ dgamma(1,1)		
		}
		hmm[i,1:2] ~ ddirich(ahmm[i,1:2] + 0.1)
	}
	
	# Estimate states from transition matrix
	for(i in 1:nfish) {
		m[i,1] <- 1
		for(k in 2:ndays) {
			m[i,k] ~ dcat(hmm[m[i,k-1],])
			s[i,k] ~ dcat(trn[m[i,k],s[i,k-1],])
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
			# Logit link for occupancy probability as a function of sand
			logit(pocc[i,j]) <- a0 + a1*sand[i,j]
		}
	}
	for(i in 1:nfish) {
		for(j in 1:ndays) {
			grid[i,j] ~ dcat(pocc[s[i,j],])
		}
	}
	
	# Logistic regression priors - occupancy
	a0 ~ dnorm(0, 0.001)
	a1 ~ dnorm(0, 0.001)
	
}", file = "Code/fishhmm_dn.jags")

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
	
# Apply tmat to first half of observations for up-migration
trnup.est = tmat(actobs[,1:(ndays/2)], 1, nstates, est = T)	

# Apply tmat to second half of observations for down-migration
trndn.est = tmat(actobs[,(ndays/2+1):ndays], 1, nstates, est = T)	

# Create array to house MLE for both up-migration and down-migration matrices
trnest = array(dim = c(2,10,10))			
trnest[1,,] = trnup.est
trnest[2,,] = trndn.est

# Initial values
inits = function() {		
	list(p = 0.5,
			a0 = 0, a1 = 0, 
			trn = trnest, 
			hmm = matrix(c(0.9,0.1,0.1,0.9), nrow = 2, ncol = 2, byrow = T))
}

# Tracked parameters
params = c("s",
		"a0",
		"a1",
		"trn",
		"p",
		"pocc",
		"hmm")