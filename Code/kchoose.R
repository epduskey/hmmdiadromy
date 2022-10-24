# A function to choose a number of knots for the multinomial GAM
#	df: data frame with each fish's state on each day
#	obs: all observed state values
#	sub: all observed substate values
#	sand: the observed sand substrate
#	tol: k will no longer increase once edf is not increasing beyong 1 + tol
kchoose = function(df, obs, sub, sand, tol) {

	# Get unique state values and create recoded values so there are no missing states
	stunq = unique(df$state[!is.na(df$state)])
	lookup = data.frame(v1 = stunq[order(stunq)])
	lookup$v2 = seq(0,nrow(lookup)-1)

	# Recode states to renamed states with contiguous integer values
	df$rcstate = lookup$v2[match(as.character(df$state), as.character(lookup$v1))]

	# Add proportion sand in each substate
	df$sand = sand[matrix(c(sub,obs),ncol=2)]

	# Create vector to store edf
	edf = list()
	
	# Initialize knots k
	k = 4

	# Create a formula for each category
	form = as.list(lapply(c(paste("rcstate ~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), 
		rep(paste("~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), nrow(lookup)-2)), as.formula))

	# Run model
	mod = gam(form, data = df, family = multinom(K=nrow(lookup)-1), method = "REML")
	
	# Get indices of spline components
	comp = c("s\\(day\\)", paste("s.", seq(nrow(lookup)-2), "\\(day\\)", sep = ""))
	inds = lapply(comp, grepl, names(mod$edf))
	
	# Sum edf of components to get total edf
	edf[[1]] = unlist(lapply(lapply(inds, subset, x = mod$edf), sum))
	
	# Repeat steps above to create comparable edf vector for k + 1
	k = k + 1
	form = as.list(lapply(c(paste("rcstate ~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), 
		rep(paste("~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), nrow(lookup)-2)), as.formula))
	mod = gam(form, data = df, family = multinom(K=nrow(lookup)-1), method = "REML")
	inds = lapply(comp, grepl, names(mod$edf))
	edf[[2]] = unlist(lapply(lapply(inds, subset, x = mod$edf), sum))
	
	# Use while loop to increase k until edf is no longer increasing by a proportion of 1 + tol
	i = 2
	while(any(edf[[i]] / edf[[i-1]] > 1 + tol)) {
		k = k + 1
		form = as.list(lapply(c(paste("rcstate ~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), 
			rep(paste("~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), nrow(lookup)-2)), as.formula))
		mod = gam(form, data = df, family = multinom(K=nrow(lookup)-1), method = "REML")
		inds = lapply(comp, grepl, names(mod$edf))
		edf[[i+1]] = unlist(lapply(lapply(inds, subset, x = mod$edf), sum))
		i = i + 1
	}
	
	# Set k
	k = k - 1
	
	# Run final model
	form = as.list(lapply(c(paste("rcstate ~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), 
		rep(paste("~ s(day, k = ", k, ") + s(fish, bs = 're')", sep = ""), nrow(lookup)-2)), as.formula))
	mod = gam(form, data = df, family = multinom(K=nrow(lookup)-1), method = "REML")
	
	# Return final model and k value
	return(list(model = mod, knots = k, edf = edf))
	
}