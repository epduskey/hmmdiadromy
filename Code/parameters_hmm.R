# Model parameters

# Parameters (ctrl-f):
#	I. Movement parameters
#	II. Habitat parameters
#	III. Detection parameters
#	IV. Search parameters

########## I. Movement parameters ##########

# Number of states - state 1 = beginning of migration, state 10\highest state = migration destination
nstates = 10

# Number of fish		
nfish = 50	

# Number of days (sampling days - note we are defining movement in reference to observation)	
ndays = 100		

# Up-migration transition probabilities: probability of moving FROM (row) any state TO (column) any other state
trnup = matrix(c(0.3, 0.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		# Up-migration transition probabilities
			     0.0, 0.1, 0.6, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.3, 0.2, 0.2, 0.2, 0.1, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.1, 0.5, 0.3, 0.1, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.4, 0.2, 0.2, 0.2, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6, 0.3, 0.1,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7, 0.3,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
			     nrow = nstates, ncol = nstates, byrow = T)

# Down-migration transition probabilities: probability of moving FROM (row) any state TO (column) any other state
trndn = matrix(c(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,		# Down-migration transition probabilities
			     0.6, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.4, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.3, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.2, 0.4, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.2, 0.2, 0.4, 0.2, 0.0, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.1, 0.3, 0.3, 0.3, 0.0, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.4, 0.2, 0.0, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.3, 0.3, 0.0,
			     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.6, 0.2),
			     nrow = nstates, ncol = nstates, byrow = T)

# Hidden markov transition matrix: probability of switching FROM (row) up-migration movement behavior TO (column) down-migration movement behavior
trnhm = matrix(c(0.98, 0.02,							
				 0.02, 0.98),
				 nrow = 2, ncol = 2, byrow = T)

# State initial distribution (all fish start out in state 1)
sinit = c(1, rep(0, nstates - 1))

# Hidden Markov inital distribution (all start with up-migration)		
hminit = c(1, rep(0, 1))

# Package all initial distributions in one list			
init = list(sinit, hminit)		

# Package all transition matrices in one list
trn = list(trnup, trndn)


########## II. Habitat parameters ##########		

# Specify the number of sub-states in each state (3 here, for simplicity)
pcs = 3	

# Choose intercept of sand preference: lower means a stronger dislike of sandless habitats
a0 = -2						

# Choose slope of sand preference: higher means a stronger selection for sandy habitats
a1 = 2		


########## III. Detection parameters ##########

# Detection probability
p = 0.3	

########## IV. Search parameters ##########

# Make a sampling procedure below i.e. a matrix where rows are states and columns are days: did we search that state on that day?

# Assume we have array that covers entire area (in this case the search function does nothing)
schedule = rep(1, nstates)			

# Matrix of search values - states in rows, days in columns
rgmat = matrix(schedule, nrow = nstates, ncol = ndays)
rg = split(rgmat, col(rgmat))

# Here's an alternative search procedure - UNCOMMENT (remove '#') THE STUFF BELOW AND COMMENT OUT (add '#') THE STUFF ABOVE

# # Five-day rotating schedule; we search 20% of the area every day on a rotation - THIS ONLY WORKS WHEN nstates = 10
# schedule = c(c(1,1,0,0,0,0,0,0,0,0),
				# c(0,0,1,1,0,0,0,0,0,0),
				# c(0,0,0,0,1,1,0,0,0,0),
				# c(0,0,0,0,0,0,1,1,0,0),
				# c(0,0,0,0,0,0,0,0,1,1))

# # Matrix of search values - states in rows, days in columns
# rgmat = matrix(schedule, nrow = nstates, ncol = ndays)
# rg = split(rgmat, col(rgmat))