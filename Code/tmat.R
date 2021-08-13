# Function that estimates the transition matrix from frequency in the data - MLE
tmat = function(X, minlev, maxlev, prob = T, est = F) {
	ret = table(factor(c(X[,-ncol(X)]),
					  levels = c(minlev:maxlev)), 
			   factor(c(X[,-1]),
			   		  levels = c(minlev:maxlev)))
	if(prob) {
		if(est) {
			ret = ret + 0.01; ret = ret / rowSums(ret)
		} else{ret = ret / rowSums(ret)}
	}
	return(ret)
}