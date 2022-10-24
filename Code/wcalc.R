# A function to calculate mclogit regression weight for each state
#	x: a vector of length nstates*pcs with 0's and one 1 indicating which substate an individual was occupying on a given day
#	nstates: number of spatial states
#	pcs: number of substates in each spatial state
#	returns a vector of length nstates*pcs with 1's indicating which spatial state an individual was occupying on a given day 
wcalc = function(x, nstates, pcs) {
	if(!is.na(x)) {
		y = rep(0, nstates*pcs)
		y[(x*pcs-pcs+1):(x*pcs)] = 1
		return(y)
	} else {
		return(rep(0,nstates*pcs))
	}
}