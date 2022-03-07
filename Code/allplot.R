# Loading and plotting from final simulation models

# Created: January 25, 2019
# Last modified: June 22, 2021

# Load packages
library(jagsUI)
library(HDInterval)
library(vioplot)
library(lattice)

# Set working directory
setwd(paste(mypath, "HMM", sep = ""))

# Contents (ctrl-f):
#	0a. Common values
#	0b. Load all models
#	0c. Load simulated data
#	I. Figure 1
#	II. Figure 2
#	III. Figure 3
#	IV. Figure 4
#	V. Figure 5
#	VI. Figure 6
#	VII. Figure 7
#	VIII. Figure 8
#	IX. Figure 9
#	X. Figure 10
#	XI. Not plotted


########## 0a. Common values ##########

# Number of fish
nfish = 50

# Number of days
ndays = 100

# Number of spatial states
nstates = 10


########## 0b. Load all models ##########

# Load the model files
load("Output/simjags.rda")
load("Output/ndhmmjags.rda")
load("Output/hmmjags.rda")
load("Output/hmmsimjags.rda")


########## 0c. Load simulated data ##########

# Load standard Markov trajectories
statseq_smm = as.matrix(read.table("Data/statseq_smm.txt", header = T))

# Load hidden Markov trajectories
statseq_hmm = as.matrix(read.table("Data/statseq_hmm.txt", header = T))


########## I. Figure 1 ##########

jpeg("Plots/figure_1.jpeg", width = 8900, height = 4450, units = 'px', res = 600)
par(mfrow = c(1,2))
par(mar = c(5,5,3,1))
boxplot.matrix(statseq_smm, outline = F, axes = F, xlab = "", ylab = "", cex.lab = 2)
axis(1, cex.axis = 1.5, at = c(1,20,40,60,80,100))
axis(2, cex.axis = 1.5, at = seq(10))
box()
par(mar = c(5,5,3,1))
boxplot.matrix(statseq_hmm, outline = F, axes = F, xlab = "", ylab = "", cex.lab = 2)
axis(1, cex.axis = 1.5, at = c(1,20,40,60,80,100))
axis(2, cex.axis = 1.5, at = seq(10))
box()
mtext("Day", 1, cex = 2, line = -2, outer = T)
mtext("State", 2, cex = 2, line = -2, outer = T)
dev.off()


########## II. Figure 2 ##########

# This chart was created using Word, and subsequently converted to PDF, then to JPEG and TIFF
# Contact lead author if you would like the original chart file


########## III. Figure 3 ##########

# Actual transition matrix
trn = matrix(c(0.5, 0.2, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			   0.2, 0.5, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			   0.1, 0.2, 0.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.2, 0.5, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.2, 0.5, 0.2, 0.1, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.2, 0.5, 0.2, 0.1, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.5, 0.2, 0.1, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.5, 0.1, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.5, 0.1,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.2, 0.5),
			   nrow = nstates, ncol = nstates, byrow = T)

# Plot transition matrix for variation 1
alltrn = data.frame(y = numeric(), x = character(), z = factor())
z.lev = character()
for(j in 1:10) {
	for(i in 1:10) {
		out = data.frame(y = sim.jags$sims.list$trn[,i,j])
		out$x = rep("1", dim(out)[1])
		out$z = rep(paste(i,",",j,sep=''), dim(out)[1])
		alltrn = rbind(alltrn,out)
		z.lev = c(z.lev, paste(i,",",j,sep=''))
	}
}
alltrn$z = factor(alltrn$z, levels= c(t(apply(matrix(z.lev,nrow=10,ncol=10), 2, rev))))

y01 = c(t(apply(trn, 2, rev)))
dl = c(t(apply(ifelse(diag(10) == 1, "black", NA), 2, rev)))
pt = c(t(apply(sim.jags$q50$trn, 2, rev)))

jpeg("Plots/figure_3.jpeg", width = 6324, height = 6324, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15,15), c(1,15,15))
par(mar = c(0,0,0,0))
frame()
text(seq(0,0.93,length.out = 10), rep(0.5,10), seq(10), font = 2, cex = 2)
frame()
text(rep(0.5,10), seq(1,0.09,length.out=10), seq(10), font = 2, cex = 2)
out = bwplot(y~x|z, data = alltrn, 
		panel = function(x, y, ...) {
		panel.violin(x, y, ...)
		panel.points(1, pt[panel.number()], lwd = 2, pch = 21, col = "black", fill = "white", cex = 0.75)
		panel.segments(0.65,y01[panel.number()],1.35,y01[panel.number()],lwd = 2)
		panel.segments(0.4,1,1.6,0,col=dl[panel.number()])}, 
		strip = F, 
		col = "gray80", lwd = 2, 
		ylim = c(0,1),
		scales = list(y = list(at = NULL, labels = NULL), x = list(at = NULL, labels = "")),
		xlab = "",
		ylab = "")
print(out, newpage = F)
par(mfrow = c(1,1))
dev.off()


########## IV. Figure 4 ##########

# Generate credible intervals for trajectories from standard simulation
actmat = as.matrix(read.table("Data/actmat.txt", header = T))
simmat1.med = sim.jags$q50$s
simmat1.hdi = apply(sim.jags$sims.list$s, c(2,3), hdi)
simmat1.25 = simmat1.hdi[1,,]
simmat1.975 = simmat1.hdi[2,,]
simmat2.med = ndhmm.jags$q50$s
simmat2.hdi = apply(ndhmm.jags$sims.list$s, c(2,3), hdi)
simmat2.25 = simmat2.hdi[1,,]
simmat2.975 = simmat2.hdi[2,,]

# Get RMSE for SMM (best, mid, worst trajectories will change with different simulations)
rmse = function(mod, act, obs) {	
	out = (rowSums(colSums(sweep(mod$sims.list$s, 2:3, act)^2, dims = 1))/mod$mcmc.info$n.samples)/rowSums(is.na(obs))
	return(out)
}
simvar1.rmse = rmse(sim.jags, statseq_smm, actmat)
simvar2.rmse = rmse(ndhmm.jags, statseq_smm, actmat)

# Store best/mid/worst indices
sv1.idx = c(which.min(simvar1.rmse), order(simvar1.rmse)[25], which.max(simvar1.rmse))
sv2.idx = c(which.min(simvar2.rmse), order(simvar2.rmse)[25], which.max(simvar2.rmse))

# For variation 1 applied to data set 1
jpeg("Plots/figure_4a.jpeg", width = 9000, height = 3000, units = 'px', res = 600)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0,0))
plot(simmat1.med[sv1.idx[1],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat1.25[sv1.idx[1],], rev(simmat1.975[sv1.idx[1],])), border = F, col = "gray")
lines(seq(ndays), simmat1.med[sv1.idx[1],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv1.idx[1],i]), NA, 9)}
points(seq(ndays), simmat1.med[sv1.idx[1],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv1.idx[1],], lwd = 2)

par(mar = c(5,5,3,1))
plot(simmat1.med[sv1.idx[2],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat1.25[sv1.idx[2],], rev(simmat1.975[sv1.idx[2],])), border = F, col = "gray")
lines(seq(ndays), simmat1.med[sv1.idx[2],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv1.idx[2],i]), NA, 9)}
points(seq(ndays), simmat1.med[sv1.idx[2],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv1.idx[2],], lwd = 2)

par(mar = c(5,5,3,1))
plot(simmat1.med[sv1.idx[3],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat1.25[sv1.idx[3],], rev(simmat1.975[sv1.idx[3],])), border = F, col = "gray")
lines(seq(ndays), simmat1.med[sv1.idx[3],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv1.idx[3],i]), NA, 9)}
points(seq(ndays), simmat1.med[sv1.idx[3],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv1.idx[3],], lwd = 2)
mtext("State", 2, cex = 2, line = -1, outer = T)
mtext("Day", 1, cex = 2, line = -1, outer = T)
mtext("SM", 3, cex = 2, line = -2, outer = T)
dev.off()

# For variation 2 applied to data set 1
jpeg("Plots/figure_4b.jpeg", width = 9000, height = 3000, units = 'px', res = 600)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0,0))
plot(simmat2.med[sv2.idx[1],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat2.25[sv2.idx[1],], rev(simmat2.975[sv2.idx[1],])), border = F, col = "gray")
lines(seq(ndays), simmat2.med[sv2.idx[1],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv2.idx[1],i]), NA, 9)}
points(seq(ndays), simmat2.med[sv2.idx[1],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv2.idx[1],], lwd = 2)

par(mar = c(5,5,3,1))
plot(simmat2.med[sv2.idx[2],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat2.25[sv2.idx[2],], rev(simmat2.975[sv2.idx[2],])), border = F, col = "gray")
lines(seq(ndays), simmat2.med[sv2.idx[2],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv2.idx[2],i]), NA, 9)}
points(seq(ndays), simmat2.med[sv2.idx[2],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv2.idx[2],], lwd = 2)

par(mar = c(5,5,3,1))
plot(simmat2.med[sv2.idx[3],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat2.25[sv2.idx[3],], rev(simmat2.975[sv2.idx[3],])), border = F, col = "gray")
lines(seq(ndays), simmat2.med[sv2.idx[3],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv2.idx[3],i]), NA, 9)}
points(seq(ndays), simmat2.med[sv2.idx[3],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv2.idx[3],], lwd = 2)
mtext("State", 2, cex = 2, line = -1, outer = T)
mtext("Day", 1, cex = 2, line = -1, outer = T)
mtext("HMM", 3, cex = 2, line = -2, outer = T)
dev.off()


########## V. Figure 5 ##########

# Generate credible intervals for trajectories from non-dynamic hmm simulation
acthmm = as.matrix(read.table("Data/acthmm.txt", header = T))
hmmmat1.med = hmm.jags$q50$s
hmmmat1.hdi = apply(hmm.jags$sims.list$s, c(2,3), hdi)
hmmmat1.25 = hmmmat1.hdi[1,,]
hmmmat1.975 = hmmmat1.hdi[2,,]
hmmmat2.med = hmmsim.jags$q50$s
hmmmat2.hdi = apply(hmmsim.jags$sims.list$s, c(2,3), hdi)
hmmmat2.25 = hmmmat2.hdi[1,,]
hmmmat2.975 = hmmmat2.hdi[2,,]

# Get RMSE for HMM (best, mid, worst trajectories will change with different simulations)
hmmvar1.rmse = rmse(hmmsim.jags, statseq_hmm, acthmm)
hmmvar2.rmse = rmse(hmm.jags, statseq_hmm, acthmm)

jpeg("Plots/figure_5.jpeg", width = 4500, height = 3300, units = 'px', res = 600)
par(mfrow = c(1,1), mar = c(5,5,3,1))
plot(simvar1.rmse, simvar2.rmse, pch = 21, cex = 2, bg = "gray20", xlab = expression("Standard Markov RMSE"["s"["k"]]), ylab = expression("HMM RMSE"["s"["k"]]), cex.lab = 1.5, xlim = c(0,5), ylim = c(0,5), cex.axis = 1.5)
points(hmmvar1.rmse, hmmvar2.rmse, pch = 22, cex = 2, bg = "gray80")
abline(0,1)
legend("topleft", bty = 'n', c("D-SM", "D_HMM"), pch = c(21,22), cex = 1.5, pt.cex = 2, pt.bg = c("gray20","gray80"))
dev.off()


########## VI. Figure 6 ##########

# Up-migration transition probabilities: probability of moving FROM (row) any state TO (column) any other state
trnup = matrix(c(0.3, 0.5, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
trndn = matrix(c(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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

# Plot transition matrix for variation 1
alltrn = data.frame(y = numeric(), x = character(), z = factor())
z.lev = character()
for(j in 1:10) {
	for(i in 1:10) {
		out = data.frame(y = hmmsim.jags$sims.list$trn[,i,j])
		out$x = rep("1", dim(out)[1])
		out$z = rep(paste(i,",",j,sep=''), dim(out)[1])
		alltrn = rbind(alltrn,out)
		z.lev = c(z.lev, paste(i,",",j,sep=''))
	}
}
alltrn$z = factor(alltrn$z, levels= c(t(apply(matrix(z.lev,nrow=10,ncol=10), 2, rev))))

hv1.trn.q50 = hmmsim.jags$q50$trn
pt = c(t(apply(hv1.trn.q50, 2, rev)))
actup = c(t(apply(trnup, 2, rev)))
actdn = c(t(apply(trndn, 2, rev)))

jpeg("Plots/figure_6.jpeg", width = 6324, height = 6324, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15,15), c(1,15,15))
par(mar = c(0,0,0,0))
frame()
text(seq(0,0.93,length.out = 10), rep(0.5,10), seq(10), font = 2, cex = 2)
frame()
text(rep(0.5,10), seq(1,0.09,length.out=10), seq(10), font = 2, cex = 2)
out = bwplot(y~x|z, data = alltrn, 
		panel = function(x, y, ...) {
		panel.violin(x, y, ...)
		panel.points(1, pt[panel.number()], lwd = 2, pch = 21, col = "black", fill = "white", cex = 0.75) 
		panel.points(1, actdn[panel.number()], lwd = 1.5, pch = 25, col = "black", fill = "white", cex = 1)
		panel.points(1, actup[panel.number()], lwd = 1.5, pch = 24, col = "black", fill = "white", cex = 1)
		panel.segments(0.4,1,1.6,0,col=dl[panel.number()])}, 
		strip = F, 
		col = "gray80", lwd = 2, 
		ylim = c(0,1),
		scales = list(y = list(at = NULL, labels = NULL), x = list(at = NULL, labels = "")),
		xlab = "",
		ylab = "")
print(out, newpage = F)
par(mfrow = c(1,1))
dev.off()


########## VII. Figure 7 ##########

hv1.idx = c(which.min(hmmvar1.rmse), order(hmmvar1.rmse)[25], which.max(hmmvar1.rmse))
hv2.idx = c(which.min(hmmvar2.rmse), order(hmmvar2.rmse)[25], which.max(hmmvar2.rmse))

# For variation 1 applied to data set 2
jpeg("Plots/figure_7a.jpeg", width = 9000, height = 3000, units = 'px', res = 600)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0,0))
plot(hmmmat1.med[hv1.idx[1],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat1.25[hv1.idx[1],], rev(hmmmat1.975[hv1.idx[1],])), border = F, col = "gray")
lines(seq(ndays), hmmmat1.med[hv1.idx[1],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv1.idx[1],i]), NA, 9)}
points(seq(ndays), hmmmat1.med[hv1.idx[1],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv1.idx[1],], lwd = 2)

par(mar = c(5,5,3,1))
plot(hmmmat1.med[hv1.idx[2],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat1.25[hv1.idx[2],], rev(hmmmat1.975[hv1.idx[2],])), border = F, col = "gray")
lines(seq(ndays), hmmmat1.med[hv1.idx[2],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv1.idx[2],i]), NA, 9)}
points(seq(ndays), hmmmat1.med[hv1.idx[2],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv1.idx[2],], lwd = 2)

par(mar = c(5,5,3,1))
plot(hmmmat1.med[hv1.idx[3],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat1.25[hv1.idx[3],], rev(hmmmat1.975[hv1.idx[3],])), border = F, col = "gray")
lines(seq(ndays), hmmmat1.med[hv1.idx[3],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv1.idx[3],i]), NA, 9)}
points(seq(ndays), hmmmat1.med[hv1.idx[3],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv1.idx[3],], lwd = 2)
mtext("State", 2, cex = 2, line = -1, outer = T)
mtext("Day", 1, cex = 2, line = -1, outer = T)
mtext("SM", 3, cex = 2, line = -2, outer = T)
dev.off()

# For variation 2 applied to data set 2
jpeg("Plots/figure_7b.jpeg", width = 9000, height = 3000, units = 'px', res = 600)
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0,0))
plot(hmmmat2.med[hv2.idx[1],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[hv2.idx[1],], rev(hmmmat2.975[hv2.idx[1],])), border = F, col = "gray")
lines(seq(ndays), hmmmat2.med[hv2.idx[1],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv2.idx[1],i]), NA, 9)}
points(seq(ndays), hmmmat2.med[hv2.idx[1],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv2.idx[1],], lwd = 2)

par(mar = c(5,5,3,1))
plot(hmmmat2.med[hv2.idx[2],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[hv2.idx[2],], rev(hmmmat2.975[hv2.idx[2],])), border = F, col = "gray")
lines(seq(ndays), hmmmat2.med[hv2.idx[2],], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv2.idx[2],i]), NA, 9)}
points(seq(ndays), hmmmat2.med[hv2.idx[2],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv2.idx[2],], lwd = 2)

par(mar = c(5,5,3,1))
plot(hmmmat1.med[26,], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 2)
axis(2, at = seq(10), cex.axis = 2)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[26,], rev(hmmmat2.975[26,])), border = F, col = "gray")
lines(seq(ndays), hmmmat2.med[26,], lwd = 2, lty = 3)
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[26,i]), NA, 9)}
points(seq(ndays), hmmmat2.med[26,], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[26,], lwd = 2)
mtext("State", 2, cex = 2, line = -1, outer = T)
mtext("Day", 1, cex = 2, line = -1, outer = T)
mtext("HMM", 3, cex = 2, line = -2, outer = T)
dev.off()


########## VIII. Figure 8 ##########

# Plot transition matrix for variation 2 - upriver
alltrn = data.frame(y = numeric(), x = character(), z = factor())
z.lev = character()
for(j in 1:10) {
	for(i in 1:10) {
		out = data.frame(y = hmm.jags$sims.list$trn[,1,i,j])
		out$x = rep("1", dim(out)[1])
		out$z = rep(paste(i,",",j,sep=''), dim(out)[1])
		alltrn = rbind(alltrn,out)
		z.lev = c(z.lev, paste(i,",",j,sep=''))
	}
}
alltrn$z = factor(alltrn$z, levels= c(t(apply(matrix(z.lev,nrow=10,ncol=10), 2, rev))))

y01 = c(t(apply(trnup, 2, rev)))
hv2.trn.q50 = hmm.jags$q50$trn[1,,]
pt = c(t(apply(hv2.trn.q50, 2, rev)))

jpeg("Plots/figure_8a.jpeg", width = 6324, height = 6324, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15,15), c(1,15,15))
par(mar = c(0,0,0,0))
frame()
text(seq(0,0.93,length.out = 10), rep(0.5,10), seq(10), font = 2, cex = 2)
frame()
text(rep(0.5,10), seq(1,0.09,length.out=10), seq(10), font = 2, cex = 2)
out = bwplot(y~x|z, data = alltrn, 
		panel = function(x, y, ...) {
		panel.violin(x, y, ...)
		panel.points(1, pt[panel.number()], lwd = 2, pch = 21, col = "black", fill = "white", cex = 0.75)
		panel.segments(0.65,y01[panel.number()],1.35,y01[panel.number()],lwd = 2)
		panel.segments(0.4,1,1.6,0,col=dl[panel.number()])}, 
		strip = F, 
		col = "gray80", lwd = 2, 
		ylim = c(0,1),
		scales = list(y = list(at = NULL, labels = NULL), x = list(at = NULL, labels = "")),
		xlab = "",
		ylab = "")
print(out, newpage = F)
par(mfrow = c(1,1))
dev.off()

# Plot transition matrix for variation 2 - downriver
alltrn = data.frame(y = numeric(), x = character(), z = factor())
z.lev = character()
for(j in 1:10) {
	for(i in 1:10) {
		out = data.frame(y = hmm.jags$sims.list$trn[,2,i,j])
		out$x = rep("1", dim(out)[1])
		out$z = rep(paste(i,",",j,sep=''), dim(out)[1])
		alltrn = rbind(alltrn,out)
		z.lev = c(z.lev, paste(i,",",j,sep=''))
	}
}
alltrn$z = factor(alltrn$z, levels= c(t(apply(matrix(z.lev,nrow=10,ncol=10), 2, rev))))

y01 = c(t(apply(trndn, 2, rev)))
hv2.trn.q50 = hmm.jags$q50$trn[2,,]
pt = c(t(apply(hv2.trn.q50, 2, rev)))

jpeg("Plots/figure_8b.jpeg", width = 6324, height = 6324, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15,15), c(1,15,15))
par(mar = c(0,0,0,0))
frame()
text(seq(0,0.93,length.out = 10), rep(0.5,10), seq(10), font = 2, cex = 2)
frame()
text(rep(0.5,10), seq(1,0.09,length.out=10), seq(10), font = 2, cex = 2)
out = bwplot(y~x|z, data = alltrn, 
		panel = function(x, y, ...) {
		panel.violin(x, y, ...)
		panel.points(1, pt[panel.number()], lwd = 2, pch = 21, col = "black", fill = "white", cex = 0.75)
		panel.segments(0.65,y01[panel.number()],1.35,y01[panel.number()],lwd = 2)
		panel.segments(0.4,1,1.6,0,col=dl[panel.number()])}, 
		strip = F, 
		col = "gray80", lwd = 2, 
		ylim = c(0,1),
		scales = list(y = list(at = NULL, labels = NULL), x = list(at = NULL, labels = "")),
		xlab = "",
		ylab = "")
print(out, newpage = F)
par(mfrow = c(1,1))
dev.off()


########## IX. Figure 9 ##########

# Actual for variation 2 run on nd hmm simulation
trnhm = matrix(c(0.98, 0.02,
			0.02, 0.98),
			nrow = 2, ncol = 2, byrow = T)

hmmvar2_11.hdi = hdi(hmm.jags$sims.list$hmm[,1,1])
hmmvar2_12.hdi = hdi(hmm.jags$sims.list$hmm[,1,2])
hmmvar2_21.hdi = hdi(hmm.jags$sims.list$hmm[,2,1])
hmmvar2_22.hdi = hdi(hmm.jags$sims.list$hmm[,2,2])

jpeg("Plots/figure_9.jpeg", width = 5000, height = 5000, units = 'px', res = 600)
matlay = layout(matrix(c(1,2,3,4),2,2,T), c(1,1,1,1), c(1,1,1,1))
layout.show(matlay)
par(mar = c(2,5,5,0))
plot(density(hmm.jags$sims.list$hmm[,1,1]), lwd = 2, axes = F, main = "", xlab = "", xlim = c(0.970,0.995), ylim = c(0,400), ylab = "")
polygon(density(hmm.jags$sims.list$hmm[,1,1]), lwd = 2, col = "gray80")
title(main = "Upriver", ylab = "Upriver", cex.main = 2, font.main = 2, cex.lab = 2, font.lab = 2)
lines(rep(0.98,2), c(0,0.7*400), lwd = 2, col = "gray40")
lines(rep(hmm.jags$q50$hmm[1,1],2), c(0,0.7*400), lty = 3, lwd = 2, col = "gray40")
arrows(hmmvar2_11.hdi[1], 0.7*400, hmmvar2_11.hdi[2], 0.7*400, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(hmmvar2_11.hdi), 0.8*400, "95% HDI", cex = 1.5)
axis(1, at = c(0.980, 0.985, 0.990), cex.axis = 1.5)
box()
par(mar = c(2,0,5,5))
plot(density(hmm.jags$sims.list$hmm[,1,2]), lwd = 2, axes = F, main = "", xlab = "", xlim = c(0.005,0.030), ylim = c(0,400), ylab = "")
polygon(density(hmm.jags$sims.list$hmm[,1,2]), lwd = 2, col = "gray80")
title(main = "Downriver", cex.main = 2, font.main = 2)
lines(rep(0.02,2), c(0,0.7*400), lwd = 2, col = "gray40")
lines(rep(hmm.jags$q50$hmm[1,2],2), c(0,0.7*400), lty = 3, lwd = 2, col = "gray40")
arrows(hmmvar2_12.hdi[1], 0.7*400, hmmvar2_12.hdi[2], 0.7*400, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(hmmvar2_12.hdi), 0.8*400, "95% HDI", cex = 1.5)
axis(1, at = c(0.010, 0.015, 0.020), cex.axis = 1.5)
box()
par(mar = c(5,5,2,0))
plot(density(hmm.jags$sims.list$hmm[,2,1]), lwd = 2, axes = F, main = "", xlab = "", xlim = c(0.005,0.030), ylim = c(0,400), ylab = "")
polygon(density(hmm.jags$sims.list$hmm[,2,1]), lwd = 2, col = "gray80")
title(ylab = "Downriver", cex.lab = 2, font.lab = 2)
lines(rep(0.02,2), c(0,0.7*400), lwd = 2, col = "gray40")
lines(rep(hmm.jags$q50$hmm[2,1],2), c(0,0.7*400), lty = 3, lwd = 2, col = "gray40")
arrows(hmmvar2_21.hdi[1], 0.7*400, hmmvar2_21.hdi[2], 0.7*400, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(hmmvar2_21.hdi), 0.8*400, "95% HDI", cex = 1.5)
axis(1, at = c(0.010, 0.015, 0.020), cex.axis = 1.5)
box()
par(mar = c(5,0,2,5))
plot(density(hmm.jags$sims.list$hmm[,2,2]), lwd = 2, axes = F, main = "", xlab = "", xlim = c(0.97,0.995), ylim = c(0,400), ylab = "")
polygon(density(hmm.jags$sims.list$hmm[,2,2]), lwd = 2, col = "gray80")
lines(rep(0.98,2), c(0,0.7*400), lwd = 2, col = "gray40")
lines(rep(hmm.jags$q50$hmm[2,2],2), c(0,0.7*400), lty = 3, lwd = 2, col = "gray40")
arrows(hmmvar2_22.hdi[1], 0.7*400, hmmvar2_22.hdi[2], 0.7*400, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(hmmvar2_22.hdi), 0.8*400, "95% HDI", cex = 1.5)
axis(1, at = c(0.980, 0.985, 0.990), cex.axis = 1.5)
box()
par(mfrow = c(1,1))
dev.off()


########## X. Figure 10 ##########

# Calculate hpdi of a's and b's
a0.mark1 = hdi(sim.jags$sims.list$a0)
a0.hmm1 = hdi(ndhmm.jags$sims.list$a0)
a0.mark2 = hdi(hmmsim.jags$sims.list$a0)
a0.hmm2 = hdi(hmm.jags$sims.list$a0)
a1.mark1 = hdi(sim.jags$sims.list$a1)
a1.hmm1 = hdi(ndhmm.jags$sims.list$a1)
a1.mark2 = hdi(hmmsim.jags$sims.list$a1)
a1.hmm2 = hdi(hmm.jags$sims.list$a1)

# Actual a and b values
a0 = -5
a1 = 10

# Initialize plot
jpeg("Plots/figure_10.jpeg", width = 8000, height = 5000, units = 'px', res = 600)
alay = layout(matrix(c(1,1,1,1,2,3,4,5,6,6,6,6,7,8,9,10),4,4,T), rep(1,4), c(1,15,1,15)) 
layout.show(alay)

# Add title to first row
par(mar = c(0,0,0,0))
frame()
text(0.5,0.5,"a.",font=2,cex=2)

# Plot density for variation 1 run
par(mar = c(5,5,3,1))
plot(density(sim.jags$sims.list$a0), xlab = expression(bold(a[0] ~ "- Standard Markov D-SM*")), ylab = "",
		xlim = c(-7,-4), ylim = c(0,2.0), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
title(ylab = "Density", cex.lab = 2)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.0,0.5))
box()
polygon(density(sim.jags$sims.list$a0), lwd = 2, col = "gray70")
arrows(a0.mark1[1], 0.72*2.0, a0.mark1[2], 0.72*2.0, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.mark1), 0.8*2.0, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.0), lwd = 2, col = "gray40")
lines(c(sim.jags$q50$a0,sim.jags$q50$a0), c(0,0.72*2.0), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,5,3,1))
plot(density(ndhmm.jags$sims.list$a0), xlab = expression(a[0] ~ "- HMM D-SM"), ylab = "",
		xlim = c(-7,-4), ylim = c(0,2.0), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.0,0.5))
box()
polygon(density(ndhmm.jags$sims.list$a0), lwd = 2, col = "gray70")
arrows(a0.hmm1[1], 0.72*2.0, a0.hmm1[2], 0.72*2.0, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.hmm1), 0.8*2.0, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.0), lwd = 2, col = "gray40")
lines(c(ndhmm.jags$q50$a0,ndhmm.jags$q50$a0), c(0,0.72*2.0), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,5,3,1))
plot(density(hmmsim.jags$sims.list$a0), xlab = expression(a[0] ~ "- Standard Markov D-HMM"), ylab = "",
		xlim = c(-7,-4), ylim = c(0,2.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.5,0.5))
box()
polygon(density(hmmsim.jags$sims.list$a0), lwd = 2, col = "gray70")
arrows(a0.mark2[1], 0.72*2.5, a0.mark2[2], 0.72*2.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.mark2), 0.8*2.5, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.5), lwd = 2, col = "gray40")
lines(c(hmmsim.jags$q50$a0,hmmsim.jags$q50$a0), c(0,0.72*2.5), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,5,3,1))
plot(density(hmm.jags$sims.list$a0), xlab = expression(bold(a[0] ~ "- HMM D-HMM*")), ylab = "", font.lab = 2, 
		xlim = c(-7,-4), ylim = c(0,2.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.5,0.5))
box()
polygon(density(hmm.jags$sims.list$a0), lwd = 2, col = "gray70")
arrows(a0.hmm2[1], 0.72*2.5, a0.hmm2[2], 0.72*2.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.hmm2), 0.8*2.5, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.5), lwd = 2, col = "gray40")
lines(c(hmm.jags$q50$a0,hmm.jags$q50$a0), c(0,0.72*2.5), lwd = 2, lty = 3, col = "gray40")

# Add title to second row
par(mar = c(0,0,0,0))
frame()
text(0.5,0.5,"b.",font=2,cex=2)

# Plot density for variation 1 run
par(mar = c(5,5,3,1))
plot(density(sim.jags$sims.list$a1), xlab = expression(bold(a[1] ~ "- Standard Markov D-SM*")), ylab = "",
		xlim = c(6,16), ylim = c(0,0.8), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
title(ylab = "Density", cex.lab = 2)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,0.8,0.2))
box()
polygon(density(sim.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.mark1[1], 0.72*0.8, a1.mark1[2], 0.72*0.8, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.mark1), 0.8*0.8, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*0.8), lwd = 2, col = "gray40")
lines(c(sim.jags$q50$a1,sim.jags$q50$a1), c(0,0.72*0.8), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,5,3,1))
plot(density(ndhmm.jags$sims.list$a1), xlab = expression(a[1] ~ "- HMM D-SM"), ylab = "",
		xlim = c(6,16), ylim = c(0,0.8), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,0.8,0.2))
box()
polygon(density(ndhmm.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.hmm1[1], 0.72*0.8, a1.hmm1[2], 0.72*0.8, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.hmm1), 0.8*0.8, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*0.8), lwd = 2, col = "gray40")
lines(c(ndhmm.jags$q50$a1,ndhmm.jags$q50$a1), c(0,0.72*0.8), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,5,3,1))
plot(density(hmmsim.jags$sims.list$a1), xlab = expression(a[1] ~ "- Standard Markov D-HMM"), ylab = "",
		xlim = c(6,16), ylim = c(0,0.8), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,0.8,0.2))
box()
polygon(density(hmmsim.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.mark2[1], 0.72*0.8, a1.mark2[2], 0.72*0.8, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.mark2), 0.8*0.8, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*0.8), lwd = 2, col = "gray40")
lines(c(hmmsim.jags$q50$a1,hmmsim.jags$q50$a1), c(0,0.72*0.8), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,5,3,1))
plot(density(hmm.jags$sims.list$a1), xlab = expression(bold(a[1] ~ "- HMM D-HMM*")), ylab = "", font.lab = 2, 
		xlim = c(6,16), ylim = c(0,0.8), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,0.8,0.2))
box()
polygon(density(hmm.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.hmm2[1], 0.72*0.8, a1.hmm2[2], 0.72*0.8, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.hmm2), 0.8*0.8, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*0.8), lwd = 2, col = "gray40")
lines(c(hmm.jags$q50$a1,hmm.jags$q50$a1), c(0,0.72*0.8), lwd = 2, lty = 3, col = "gray40")

# Reset and turn off device
par(mfrow = c(1,1))
dev.off()


########## XI. Not plotted ##########

# Plot transition matrix for variation 2 - upriver
alltrn = data.frame(y = numeric(), x = character(), z = factor())
z.lev = character()
for(j in 1:10) {
	for(i in 1:10) {
		out = data.frame(y = ndhmm.jags$sims.list$trn[,1,i,j])
		out$x = rep("1", dim(out)[1])
		out$z = rep(paste(i,",",j,sep=''), dim(out)[1])
		alltrn = rbind(alltrn,out)
		z.lev = c(z.lev, paste(i,",",j,sep=''))
	}
}
alltrn$z = factor(alltrn$z, levels= c(t(apply(matrix(z.lev,nrow=10,ncol=10), 2, rev))))

y01 = c(t(apply(trn, 2, rev)))
pt = c(t(apply(ndhmm.jags$q50$trn[1,,], 2, rev)))

jpeg("Plots/trn_density_simvar2_up.jpeg", width = 6324, height = 6324, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15,15), c(1,15,15))
par(mar = c(0,0,0,0))
frame()
text(seq(0,0.93,length.out = 10), rep(0.5,10), seq(10), font = 2, cex = 2)
frame()
text(rep(0.5,10), seq(1,0.09,length.out=10), seq(10), font = 2, cex = 2)
out = bwplot(y~x|z, data = alltrn, 
		panel = function(x, y, ...) {
		panel.violin(x, y, ...)
		panel.points(1, pt[panel.number()], lwd = 2, pch = 21, col = "black", fill = "white", cex = 0.75) 
		panel.segments(0.65,y01[panel.number()],1.35,y01[panel.number()],lwd = 2)
		panel.segments(0.4,1,1.6,0,col=dl[panel.number()])},
		strip = F, 
		col = "gray80", lwd = 2, 
		ylim = c(0,1),
		scales = list(y = list(at = NULL, labels = NULL), x = list(at = NULL, labels = "")),
		xlab = "",
		ylab = "")
print(out, newpage = F)
par(mfrow = c(1,1))
dev.off()

# Plot transition matrix for variation 2 - downriver
alltrn = data.frame(y = numeric(), x = character(), z = factor())
z.lev = character()
for(j in 1:10) {
	for(i in 1:10) {
		out = data.frame(y = ndhmm.jags$sims.list$trn[,2,i,j])
		out$x = rep("1", dim(out)[1])
		out$z = rep(paste(i,",",j,sep=''), dim(out)[1])
		alltrn = rbind(alltrn,out)
		z.lev = c(z.lev, paste(i,",",j,sep=''))
	}
}
alltrn$z = factor(alltrn$z, levels= c(t(apply(matrix(z.lev,nrow=10,ncol=10), 2, rev))))

pt = c(t(apply(ndhmm.jags$q50$trn[2,,], 2, rev)))

jpeg("Plots/trn_density_simvar2_dn.jpeg", width = 6324, height = 6324, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15,15), c(1,15,15))
par(mar = c(0,0,0,0))
frame()
text(seq(0,0.93,length.out = 10), rep(0.5,10), seq(10), font = 2, cex = 2)
frame()
text(rep(0.5,10), seq(1,0.09,length.out=10), seq(10), font = 2, cex = 2)
out = bwplot(y~x|z, data = alltrn, 
		panel = function(x, y, ...) {
		panel.violin(x, y, ...)
		panel.points(1, pt[panel.number()], lwd = 2, pch = 21, col = "black", fill = "white", cex = 0.75)
		panel.segments(0.4,1,1.6,0,col=dl[panel.number()])}, 
		strip = F, 
		col = "gray80", lwd = 2, 
		ylim = c(0,1),
		scales = list(y = list(at = NULL, labels = NULL), x = list(at = NULL, labels = "")),
		xlab = "",
		ylab = "")
print(out, newpage = F)
par(mfrow = c(1,1))
dev.off()