# Simulation functions

# Functions (ctrl-f):
#	I. Movement
#	II. Habitat
#	III. Detection

########## I. Movement ##########

# Function to move the fish according to tprob matrix
move = function(init, trn, trnhm, nfish, ndays) {				
	
	# Normalize the number of days
	days.norm = seq(0, ndays)/ndays								
	sinit = init[[1]]
	hminit = init[[2]]
	
	# Return matrix of river mile states
	sret = matrix(NA, nrow = nfish, ncol = ndays )
	
	# Return matrix of hidden states				
	hmret = matrix(NA, nrow = nfish, ncol = ndays )	
	
	# Name rows after fish			
	rownames(sret) = paste("Fish", seq(nfish), sep = "")	
	
	# Name columns after the day	
	colnames(sret) = paste("Day", seq(ndays), sep = "")	
	
	# Name rows after fish		
	rownames(hmret) = paste("Fish", seq(nfish), sep = "")
	
	# Name columns after the day		
	colnames(hmret) = paste("Day", seq(ndays), sep = "")		
	
	# Move each fish
	for(i in 1:nfish) {	
		
		# Specify inital location via inital dist										
		sret[i,1] = which(rmultinom(1, 1, sinit) == 1)	
		
		# Specify inital hidden state via inital dist		
		hmret[i,1] = which(rmultinom(1, 1, hminit) == 1)
		
		# Move the fish once each day		
		for(j in 2:(ndays)) {				
			
			# Hidden state probabilities based on previous day					
			hmprob = trnhm[hmret[i,j-1], ]						
			
			# Determine the hmm state
			hmret[i,j] = which(rmultinom(1, 1, hmprob) == 1)	
			
			# Transition matrix is row of previous state	
			sprob = trn[[hmret[i,j]]][sret[i,j-1], ]			
			
			# Move according to sprob
			sret[i,j] = which(rmultinom(1, 1, sprob) == 1)		
		}
	}
	return(list("stat" = sret, "hmm" = hmret))	
}


########## II. Habitat ##########

# Function to divide the river into pieces
hack = function(pcs, nstates) {										
	
	# Return matrix
	ret = matrix(NA, nrow = pcs, ncol = nstates)						
	
	# Call rows after pieces
	rownames(ret) = paste("Pc", seq(pcs), sep = "")					
	
	# Call columns after river miles
	colnames(ret) = paste("Rm", seq(nstates), sep = "")				
	
	# All row indices
	inds = seq(nstates)												
	
	# Maximum proportion for sand (starts at 1)
	throwmax = 1													
	
	# For each piece
	for(i in seq(pcs - 1)) {										
		
		# Throw some sand into the river
		throw = runif(nstates, 0, throwmax)						
		
		# Set first piece to result of throwing sand	
		ret[i, ] = throw											
		
		# Reset throwmax to what's left
		throwmax = 1 - throw										
		
		# If there's only one piece left
		if((i+1) %% pcs == 0) {										
			
			# Set proportion to what's left
			ret[i+1, ] = throwmax									
		}
	}
	return(ret)
}	

# Occupancy probability as a function of sand
pocc = function(a0, a1, sand) {								
	w = expit(a0 + sand*a1)
	p = w/sum(w)
}

# Simulate the location of the fish in the subsections
where = function(statseq, occ) {
	
	# Columns are days
	nd = ncol(statseq)										
	
	# Rows are fish
	nc = nrow(statseq)										
	
	# Return matrix
	ret = matrix(NA, nc, nd)									
	colnames(ret) = paste("Day", seq(nd), sep = "")
	rownames(ret) = paste("Fish", seq(nc), sep = "")	

	# For each day
	for(i in 1:(nc)) {									
		
		# For each fish	
		for(j in 1:nd) {									
			
			# Probabilities from statseq
			prob = occ[, statseq[i,j]]						
			ret[i,j] = which(rmultinom(1, 1, prob) == 1)
		}
	}	
	return(ret)
}


########## III. Detection ##########

# Function for detection
see = function(statseq, p) {								
	
	# Columns are days
	nd = ncol(statseq)									
	
	# Rows are fish
	nc = nrow(statseq)									
	
	# Return matrix
	ret = matrix(NA, nc, ndays)		
						
	# Column names (days)
	colnames(ret) = paste("Day", seq(nd), sep = "")		
	
	# Row names (fish)
	rownames(ret) = paste("Fish", seq(nc), sep = "")		

	# For each fish
	for(i in 1:nc) {										
		
		# For each day
		for(j in 1:nd) {									
			
			# Binomial detection
			ret[i,j] = rbinom(1, 1, p)						
		}
	}
	return(ret)
}

# Function to remove states of the unobserved fish
cut = function(statseq, statdet) {						
	
	# Columns are days
	ndays = ncol(statseq)								
	
	# Rows are fish
	nfish = nrow(statseq)							
	
	# Return matrix	
	ret = matrix(NA, nrow = nfish, ncol = ndays)		
	
	# Column names (days)
	colnames(ret) = paste("Day", seq(ndays), sep = "")	
	
	# Row names (fish)
	rownames(ret) = paste("Fish", seq(nfish), sep = "")	
	
	# Sets observed states equal to actual states
	for(i in 1:nrow(ret)) {								
		for(j in 1:ncol(ret)) {
			
			# If fish was observed, set value equal to actual
			if(statdet[i,j] == 1) {						
				ret[i,j] = statseq[i,j]
			}
		}
	}
	return(ret)	
}

# Function to cut out unobersved states
search = function(rg, hypobs, substat, statdet) {					
	
	# Columns are days
	nd = ncol(hypobs)												
	
	# Rows are fish
	nc = nrow(hypobs)												
	
	# Return matrix for states
	retstat = matrix(NA, nrow = nc, ncol = nd)						
	
	# Return matrix for substates
	retsub = matrix(NA, nrow = nc, ncol = nd)					
	
	# Return matrix for detection	
	retdet = statdet												
	
	# Column names (days)
	colnames(retstat) = paste("Day", seq(nd), sep = "")				
	
	# Row names (fish)
	rownames(retstat) = paste("Fish", seq(nc), sep = "")			
	
	# Column names (days)	
	colnames(retsub) = paste("Day", seq(nd), sep = "")				
	
	# Row names (fish)
	rownames(retsub) = paste("Fish", seq(nc), sep = "")				

	# Checks that length of list range is correct
	if(length(rg) == nd) {										
		
		# For each day	
		for(i in 1:nd) {											
			
			# Indices of fish available for search
			inds = which(hypobs[,i] %in% which(rg[[i]] == 1))					
			
			# If searched, set equal to actual state
			retstat[inds, i] = hypobs[inds, i]						
			
			# If searched, set equal to actual sub-state
			retsub[inds, i] = substat[inds, i]						
			
			# Indices of fish not available for search
			indsdet = which(!(hypobs[,i] %in% which(rg[[i]] == 1)))				
			
			# Set indices of non-searchable fish to 0
			retdet[indsdet, i] = 0									
		}
		return(list(actobs = retstat, 
					actsub = retsub, 
					actdet = retdet))
	
	# If start and end not correct, return error
	} else {														
		print("you must have one vector of search areas for each day")
	}
}