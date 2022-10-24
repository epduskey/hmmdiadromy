# Loading and plotting from final simulation models

# Created: January 25, 2019
# Last modified: October 10, 2022

# Load packages
library(jagsUI)
library(HDInterval)
library(vioplot)
library(lattice)
library(mgcv)
library(mclogit)
library(rje)

# Set working directory
setwd(paste(mypath, "HMM", sep = ""))

# Contents (ctrl-f):
#	0a. Common values and functions
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
#	IX. Figure S1
#	X. Figure S2
#	XI. Figure S3
#	XII. Figure S4
#	XIII. Figure S5
#	XIV. Figure S6
#	XV. Figure S7


########## 0a. Common values and functions ##########

# Source SMM data generating values
source("Code/parameters_smm.R")

# Source HMM data generating values
source("Code/parameters_hmm.R")

# RMSE for Markov models
rmse_markov = function(mod, act, obs, type = "state") {
	if(type == "state") {	
		out = (rowSums(colSums(sweep(ceiling(mod$sims.list$out/pcs), 2:3, act)^2, dims = 1))/mod$mcmc.info$n.samples)/rowSums(is.na(obs))
	} else if(type == "sub") {
		out = (rowSums(colSums(sweep(mod$sims.list$out, 2:3, act)^2, dims = 1))/mod$mcmc.info$n.samples)/rowSums(is.na(obs))
	} else {
		stop("argument type must by either 'state' or 'sub'")
	}
	return(out)
}

# RMSE for GAM and CAR models
rmse_gc = function(mod, act, obs) {	
	out = (rowSums(colSums(sweep(mod, 2:3, act)^2, dims = 1))/dim(mod)[1])/rowSums(is.na(obs))
	return(out)
}


########## 0b. Load all models ##########

# Load SMM model files
load("Output/simjags.rda")
load("Output/ndhmmjags.rda")

# Load HMM model files
load("Output/hmmjags.rda")
load("Output/hmmsimjags.rda")

# Load GAM model files
load("Output/simgam.rda")
load("Output/hmmgam.rda")

# Load STCAR models
load("Output/simcar.rda")
load("Output/hmmcar.rda")


########## 0c. Load simulated data ##########

# Load standard Markov trajectories
statseq_smm = as.matrix(read.table("Data/statseq_smm.txt", header = T))

# Load standard Markov substates
subseq_smm = as.matrix(read.table("Data/subseq_smm.txt", header = T))

# Load hidden Markov trajectories
statseq_hmm = as.matrix(read.table("Data/statseq_hmm.txt", header = T))

# Load hidden Markov substates
subseq_hmm = as.matrix(read.table("Data/subseq_hmm.txt", header = T))

# Read SMM, HMM, and GAM output data frames
simdf = read.table("Data/simdf.txt", header = T)
hmmdf = read.table("Data/hmmdf.txt", header = T)
mcldf.smm = read.table("Data/smmmcldf.txt", header = T)
mcldf.hmm = read.table("Data/hmmmcldf.txt", header = T)

# Load standard Markov observations
actmat = as.matrix(read.table("Data/actmat.txt", header = T))

# Load hidden Markov observations
acthmm = as.matrix(read.table("Data/acthmm.txt", header = T))

# Load sand observations
sand_hmm = as.matrix(read.table("Data/sand_hmm.txt"))
sand_smm = as.matrix(read.table("Data/sand_smm.txt"))


########## I. Figure 1 ##########

jpeg("Plots/figure_1.jpeg", width = 4450, height = 4450, units = 'px', res = 600)
par(mfrow = c(1,1), mar = c(5,5,3,1))
boxplot.matrix(statseq_hmm, outline = F, axes = F, xlab = "Day", ylab = "State", cex.lab = 2)
axis(1, cex.axis = 1.5, at = c(1,20,40,60,80,100))
axis(2, cex.axis = 1.5, at = seq(10))
box()
dev.off()


########## II. Figure 2 ##########

# This chart was created using Word, and subsequently converted to PDF, then to JPEG and TIFF
# Contact corresponding author if you would like the original chart file


########## III. Figure 3 ##########

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
dl = c(t(apply(ifelse(diag(10) == 1, "black", NA), 2, rev)))
pt = c(t(apply(hv2.trn.q50, 2, rev)))

jpeg("Plots/figure_3a.jpeg", width = 4400, height = 4400, units = 'px', res = 600)
matlay = layout(matrix(c(0,1,2,3),2,2,T), c(1,15), c(1,15))
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
dl = c(t(apply(ifelse(diag(10) == 1, "black", NA), 2, rev)))
pt = c(t(apply(hv2.trn.q50, 2, rev)))

jpeg("Plots/figure_3b.jpeg", width = 4400, height = 4400, units = 'px', res = 600)
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

# Calculate hpdi of hidden transition matrix
hmmvar2_11.hdi = hdi(hmm.jags$sims.list$hmm[,1,1])
hmmvar2_12.hdi = hdi(hmm.jags$sims.list$hmm[,1,2])
hmmvar2_21.hdi = hdi(hmm.jags$sims.list$hmm[,2,1])
hmmvar2_22.hdi = hdi(hmm.jags$sims.list$hmm[,2,2])

jpeg("Plots/figure_4.jpeg", width = 5000, height = 5000, units = 'px', res = 600)
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


########## V. Figure 5 ##########

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
dl = c(t(apply(ifelse(diag(10) == 1, "black", NA), 2, rev)))
actup = c(t(apply(trnup, 2, rev)))
actdn = c(t(apply(trndn, 2, rev)))

jpeg("Plots/figure_5.jpeg", width = 4400, height = 4400, units = 'px', res = 600)
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


########## VI. Figure 6 ##########

# Get actual proportions in each state on each day
hmm.table = apply(statseq_hmm, 2, function(x) table(factor(x, levels = seq(10))))
hmm.prop = unname(hmm.table / colSums(hmm.table))

# Extract fitted CAR values
hmmcar.fitted = rbind(hmm.car[[1]]$samples$fitted, hmm.car[[2]]$samples$fitted, hmm.car[[3]]$samples$fitted) 
hmmcar.mean = array(colMeans(hmmcar.fitted), dim = c(pcs,nstates,ndays))

# Get GAM means for each area on each day
hmmgam.preds = predict(hmm.gm$hmm.gam$model, newdata = hmmdf, type = "response")

# Get probability in each state on each day
linds = seq(1, nfish*ndays, by = nfish)
hinds = seq(nfish, nfish*ndays, by = nfish)
hmmgam.mean = matrix(NA, nrow = ndays, ncol = nstates)
hmmcar.prob = array(dim = c(ndays,pcs,nstates))
for(i in 1:ndays) {
	hmmgam.mean[i,] = colMeans(hmmgam.preds[linds[i]:hinds[i], ])
	hmmcar.prob[i,,] = hmmcar.mean[,,i] / sum(hmmcar.mean[,,i])
}

# Initialize plot
jpeg("Plots/figure_6.jpeg", width = 2700, height = 5100, units = 'px', res = 600)

# Gray colors for bars
col.bp = hcl.colors(n = nstates, palette = "Grays")

# Plot barplot of state occupancy probability across all days, including observations
par(mfrow = c(3,1), mar = c(3,3,1,1), oma = c(3,3,0,4))
hmmobs.bp = barplot(hmm.prop, space = 0, col = col.bp, cex.axis = 1.2, main = "Actual Proportions", cex.main = 1.5)
axis(1, at = hmmobs.bp[c(1,seq(20,100,by=20))], labels = c(1,seq(20,100,by=20)), cex.axis = 1.2)
box()
hmmgam.bp = barplot(t(hmmgam.mean), space = 0, col = col.bp, cex.axis = 1.2, main = "GAM", cex.main = 1.5)
axis(1, at = hmmgam.bp[c(1,seq(20,100,by=20))], labels = c(1,seq(20,100,by=20)), cex.axis = 1.2)
box()
hmmcar.bp = barplot(t(apply(hmmcar.prob, c(1,3), sum)), space = 0, col = col.bp, cex.axis = 1.2, main = "CAR", cex.main = 1.5)
axis(1, at = hmmcar.bp[c(1,seq(20,100,by=20))], labels = c(1,seq(20,100,by=20)), cex.axis = 1.2)
box()
mtext("Day", side = 1, cex = 1.5, outer = T)
mtext("Cumulative Occupancy Probability", side = 2, cex = 1.5, outer = T)

# Legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', sapply(nstates:1, toString), col = rev(col.bp), pch = 15, pt.cex = 3, xpd = TRUE, cex = 1.2, title = "State")

# Reset and turn off device
par(mfrow = c(1,1))
dev.off()


########## VII. Figure 7 ##########

# Plotting colors
color = palette.colors(n = 4, palette = "Okabe-Ito")
col.traj = apply(col2rgb(color), 2, function(x) rgb(x[1], x[2], x[3], 30, maxColorValue = 255))

# Generate credible intervals for HM trajectories from SMM and HMM 
hmmmat1.med = ceiling(hmm.jags$q50$out/pcs)
hmmmat1.hdi = apply(ceiling(hmm.jags$sims.list$out/pcs), c(2,3), hdi)
hmmmat1.25 = hmmmat1.hdi[1,,]
hmmmat1.975 = hmmmat1.hdi[2,,]
hmmmat2.med = ceiling(hmmsim.jags$q50$out/pcs)
hmmmat2.hdi = apply(ceiling(hmmsim.jags$sims.list$out/pcs), c(2,3), hdi)
hmmmat2.25 = hmmmat2.hdi[1,,]
hmmmat2.975 = hmmmat2.hdi[2,,]

# Get RMSE for HMM
hmmvar2.rmse = rmse_markov(hmm.jags, statseq_hmm, acthmm)

# Store best/mid/worst indices for native models
hv2.idx = c(which.min(hmmvar2.rmse), order(hmmvar2.rmse)[25], which.max(hmmvar2.rmse))

# Generate confidence intervals for HM from GAM
set.seed(349)
hmmgam.prob = apply(predict(hmm.gm$hmm.gam$model, newdata = list(day = rep(seq(ndays),nfish), fish = rep(seq(nfish),each=ndays)), type = "response"), 1, rmultinom, n = 10000, size = 1)
hmmgam.samples = aperm(array(apply(hmmgam.prob, 2, function(x) ifelse(which(x == 1) %% nstates == 0, nstates, which(x == 1) %% nstates)), dim = c(10000,ndays,nfish)), c(1,3,2))
hmmgam.all = apply(hmmgam.samples, c(2:3), quantile, prob = c(0.025,0.5,0.975))
hmmgam.traj = list()
hmmgam.traj$best = hmmgam.all[,hv2.idx[1],]
hmmgam.traj$mid = hmmgam.all[,hv2.idx[2],]
hmmgam.traj$worst = hmmgam.all[,hv2.idx[3],]

# Generate confidence intervals for HM from STCAR
set.seed(501)
hmmcar.all = apply(hmmcar.prob, 1, rmultinom, n = 10000, size = 1)
hmmcar.samples = array(apply(hmmcar.all, 2, function(x) ifelse(which(x == 1) %% (pcs*nstates) == 0, pcs*nstates, which(x == 1) %% (pcs*nstates))), dim = c(10000,ndays))
hmmcar.traj = apply(hmmcar.samples, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Initialize plot
jpeg("Plots/figure_7.jpeg", width = 7500, height = 2500, units = 'px', res = 600)

# For all models applied to HMM data set
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0,8))
plot(hmmmat2.med[hv2.idx[1],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 1.5)
axis(2, at = seq(10), cex.axis = 1.5)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[hv2.idx[1],], rev(hmmmat2.975[hv2.idx[1],])), border = F, col = col.traj[1])
lines(seq(ndays), hmmmat2.med[hv2.idx[1],], lwd = 2, lty = 3, col = color[1])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat1.25[hv2.idx[1],], rev(hmmmat1.975[hv2.idx[1],])), border = F, col = col.traj[2])
lines(seq(ndays), hmmmat1.med[hv2.idx[1],], lwd = 2, lty = 4, col = color[2])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmgam.traj$best[1,], rev(hmmgam.traj$best[3,])), border = F, col = col.traj[3])
lines(seq(ndays), hmmgam.traj$best[2,], lwd = 2, lty = 5, col = color[3])
polygon(c(seq(ndays), rev(seq(ndays))), c(ceiling(hmmcar.traj[1,]/pcs), rev(ceiling(hmmcar.traj[3,]/pcs))), border = F, col = col.traj[4])
lines(seq(ndays), ceiling(hmmcar.traj[2,]/pcs), lwd = 2, lty = 6, col = color[4])
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv2.idx[1],i]), NA, 9)}
points(seq(ndays), hmmmat2.med[hv2.idx[1],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv2.idx[1],], lwd = 2)

plot(hmmmat2.med[hv2.idx[2],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, cex.main = 2, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 1.5)
axis(2, at = seq(10), cex.axis = 1.5)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[hv2.idx[2],], rev(hmmmat2.975[hv2.idx[2],])), border = F, col = col.traj[1])
lines(seq(ndays), hmmmat2.med[hv2.idx[2],], lwd = 2, lty = 3, col = color[1])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat1.25[hv2.idx[2],], rev(hmmmat1.975[hv2.idx[2],])), border = F, col = col.traj[2])
lines(seq(ndays), hmmmat1.med[hv2.idx[2],], lwd = 2, lty = 4, col = color[2])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmgam.traj$mid[1,], rev(hmmgam.traj$mid[3,])), border = F, col = col.traj[3])
lines(seq(ndays), hmmgam.traj$mid[2,], lwd = 2, lty = 5, col = color[3])
polygon(c(seq(ndays), rev(seq(ndays))), c(ceiling(hmmcar.traj[1,]/pcs), rev(ceiling(hmmcar.traj[3,]/pcs))), border = F, col = col.traj[4])
lines(seq(ndays), ceiling(hmmcar.traj[2,]/pcs), lwd = 2, lty = 6, col = color[4])
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[hv2.idx[2],i]), NA, 9)}
points(seq(ndays), hmmmat2.med[hv2.idx[2],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[hv2.idx[2],], lwd = 2)

plot(hmmmat1.med[26,], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 1.5)
axis(2, at = seq(10), cex.axis = 1.5)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[hv2.idx[3],], rev(hmmmat2.975[hv2.idx[3],])), border = F, col = col.traj[1])
lines(seq(ndays), hmmmat2.med[hv2.idx[3],], lwd = 2, lty = 3, col = color[1])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat1.25[hv2.idx[3],], rev(hmmmat1.975[hv2.idx[3],])), border = F, col = col.traj[2])
lines(seq(ndays), hmmmat1.med[hv2.idx[3],], lwd = 2, lty = 4, col = color[2])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmgam.traj$worst[1,], rev(hmmgam.traj$worst[3,])), border = F, col = col.traj[3])
lines(seq(ndays), hmmgam.traj$worst[2,], lwd = 2, lty = 5, col = color[3])
polygon(c(seq(ndays), rev(seq(ndays))), c(ceiling(hmmcar.traj[1,]/pcs), rev(ceiling(hmmcar.traj[3,]/pcs))), border = F, col = col.traj[4])
lines(seq(ndays), ceiling(hmmcar.traj[2,]/pcs), lwd = 2, lty = 6, col = color[4])
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(acthmm[26,i]), NA, 9)}
points(seq(ndays), hmmmat2.med[26,], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_hmm[26,], lwd = 2)

# Add outer margin labels
mtext("State", 2, cex = 2, line = -1, outer = T)
mtext("Day", 1, cex = 2, line = -1, outer = T)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("SMM", "HMM", "GAM", "STCAR", "Data"), col = c(col.traj,color[1]), lwd = 2, pch = c(15,15,15,15,18), lty = c(0,0,0,0,0), pt.cex = c(4,4,4,4,1.5), xpd = TRUE, cex = 1.5)
legend("right", bty = 'n', c("SMM", "HMM", "GAM", "STCAR", "Data"), col = c(color,color[1]), pch = c(26,26,26,26,26), lty = c(3,4,5,6,1), lwd = 2, xpd = TRUE, cex = 1.5)

# Finish plots
dev.off()


########## VIII. Figure 8 ##########

# Generate substates for the GAM model
set.seed(582)
gamsub.allp = array(predict(hmm.gm$hmm.mcl, newdata = mcldf.hmm, type = "response"), dim = c(pcs,nstates,nfish,ndays))
gamsub.prob = array(gamsub.allp[cbind(rep(seq(pcs),nfish*ndays),rep(c(hmmgam.all[2,,]),each=pcs),rep(rep(seq(nfish),ndays),each=pcs),rep(seq(ndays),each=nfish*pcs))], dim = c(pcs,nfish,ndays))
gamsub.subs = apply(gamsub.prob, c(2,3), rmultinom, n = 10000, size = 1)
gamsub.samples = apply(gamsub.subs, c(2,3), function(x) ifelse(which(x == 1) %% pcs == 0, pcs, which(x == 1) %% pcs))
gamsub.all = hmmgam.samples*pcs - (pcs-gamsub.samples)

# Calculate HMM sub-state probabilities from each model
hmmact.subprob = expit(a0 + a1*sand_hmm); hmmact.subprob = sweep(hmmact.subprob, 2, colSums(hmmact.subprob), "/")
hmmv1.subprob = expit(mean(hmmsim.jags$sims.list$a0) + mean(hmmsim.jags$sims.list$a1)*sand_hmm); hmmv1.subprob = sweep(hmmv1.subprob, 2, colSums(hmmv1.subprob), "/")
hmmv2.subprob = expit(mean(hmm.jags$sims.list$a0) + mean(hmm.jags$sims.list$a1)*sand_hmm); hmmv2.subprob = sweep(hmmv2.subprob, 2, colSums(hmmv2.subprob), "/")
hmmv3.subprob = sweep(apply(gamsub.allp, c(1,2), sum), 2, colSums(apply(gamsub.allp, c(1,2), sum)), "/")
hmmv4.subprob = sweep(rowSums(hmmcar.mean, dims = 2), 2, colSums(rowSums(hmmcar.mean, dims = 2)), "/")

# Create data frame to house probabilities and model type
hmm.subprob = data.frame(prob = c(hmmact.subprob,hmmv1.subprob,hmmv2.subprob,hmmv3.subprob,hmmv4.subprob), state = rep(seq(nstates),each=pcs), pcs = seq(3))
hmm.subprob$type = factor(rep(c("Actual","HMM","SMM","GAM","STCAR"),each=pcs*nstates), levels = c("Actual","HMM","SMM","STCAR","GAM"))

# Plot all probability estimates together
jpeg("Plots/figure_8.jpeg", width = 6000, height = 3200, units = 'px', res = 600)
colorkey = list(at = seq(0,1,length.out=100), labels = list(at=seq(0,1,0.2),cex = 1.2))
levelplot(prob ~ pcs*state | type, data = hmm.subprob,
	xlab = list("Sub-state",cex=1.5), ylab = list("State",cex=1.5), main = list("Probability of occupying each sub-state",cex=1.5), scales = list(cex.axis=1.2,cex=1.2,x=list(at=seq(pcs)),y=list(at=seq(nstates))),
	colorkey = colorkey, col.regions = hcl.colors(100, "Light Grays", rev = T), 
	layout = c(5,1), par.strip.text = list(font=2), par.settings = list(strip.background=list(col="white")))
dev.off()


########## IX. Figure 9 ##########

# Get RMSE for HMM states
hmmvar1.lnest = log(rmse_markov(hmmsim.jags, statseq_hmm, acthmm) + 1)
hmmvar2.lnest = log(rmse_markov(hmm.jags, statseq_hmm, acthmm) + 1)
hmmvar3.lnest = log(rmse_gc(hmmgam.samples, statseq_hmm, acthmm) + 1)
hmmvar4.lnest = log(rmse_gc(aperm(array(rep(ceiling(hmmcar.samples/pcs), nfish), dim = c(10000,ndays,nfish)), c(1,3,2)), statseq_hmm, acthmm) + 1)

# Get RMSE for HMM substates
hmmvar1.lnesb = log(rmse_markov(hmmsim.jags, statseq_hmm*pcs - (pcs-subseq_hmm), acthmm, "sub") + 1)
hmmvar2.lnesb = log(rmse_markov(hmm.jags, statseq_hmm*pcs - (pcs-subseq_hmm), acthmm, "sub") + 1)
hmmvar3.lnesb = log(rmse_gc(gamsub.all, statseq_hmm*pcs - (pcs-subseq_hmm), acthmm) + 1)
hmmvar4.lnesb = log(rmse_gc(aperm(array(rep(hmmcar.samples, nfish), dim = c(10000,ndays,nfish)), c(1,3,2)), statseq_hmm*pcs - (pcs-subseq_hmm), acthmm) + 1)

jpeg("Plots/figure_9.jpeg", width = 3500, height = 4800, units = 'px', res = 600)
par(mfrow = c(2,1), mar = c(3,5,2,1), oma = c(2,1,0,5))
plot(hmmvar2.lnest, hmmvar1.lnest, pch = 22, cex = 2, bg = "gray20", xlab = "", ylab = "", main = "State", cex.main = 1.5, xlim = c(0,2), ylim = c(0,5), cex.axis = 1.5)
points(hmmvar2.lnest, hmmvar3.lnest, pch = 22, cex = 2, bg = "gray50")
points(hmmvar2.lnest, hmmvar4.lnest, pch = 22, cex = 2, bg = "gray80")
abline(0,1)
plot(hmmvar2.lnesb, hmmvar1.lnesb, pch = 22, cex = 2, bg = "gray20", xlab = "", ylab = "", main = "Sub-state", cex.main = 1.5, xlim = c(1,4), ylim = c(0,7), axes = F); axis(1, at = seq(4), cex.axis = 1.5); axis(2, cex.axis = 1.5); box()
points(hmmvar2.lnesb, hmmvar3.lnesb, pch = 22, cex = 2, bg = "gray50")
points(hmmvar2.lnesb, hmmvar4.lnesb, pch = 22, cex = 2, bg = "gray80")
abline(0,1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
mtext(expression("Native Model RMSE"["s"["k"]]), side = 1, line = -1, cex = 1.5, outer = T)
mtext(expression("Alternative Model RMSE"["s"["k"]]), side = 2, line = -2, cex = 1.5, outer = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("SMM","GAM","STCAR"), pt.bg = c("gray20","gray50","gray80"), pch = 22, pt.cex = 1.5, xpd = TRUE, cex = 1.2)
dev.off()


########## X. Figure S1 ##########

jpeg("Plots/figure_S1.jpeg", width = 4450, height = 4450, units = 'px', res = 600)
par(mfrow = c(1,1), mar = c(5,5,3,1))
boxplot.matrix(statseq_smm, outline = F, axes = F, xlab = "Day", ylab = "State", cex.lab = 2)
axis(1, cex.axis = 1.5, at = c(1,20,40,60,80,100))
axis(2, cex.axis = 1.5, at = seq(10))
box()
dev.off()


########## XI. Figure S2 ##########

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

jpeg("Plots/figure_S2.jpeg", width = 4400, height = 4400, units = 'px', res = 600)
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


########## XII. Figure S3 ##########

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

jpeg("Plots/figure_S3a.jpeg", width = 4400, height = 4400, units = 'px', res = 600)
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

jpeg("Plots/figure_S3b.jpeg", width = 4400, height = 4400, units = 'px', res = 600)
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


########## XIII. Figure S4 ##########

# Get actual proportions in each state on each day
smm.table = apply(statseq_smm, 2, function(x) table(factor(x, levels = seq(10))))
smm.prop = unname(smm.table / colSums(smm.table))

# Extract fitted CAR values
simcar.fitted = rbind(sim.car[[1]]$samples$fitted, sim.car[[2]]$samples$fitted, sim.car[[3]]$samples$fitted) 
simcar.mean = array(colMeans(simcar.fitted), dim = c(pcs,nstates,ndays))

# Get GAM means for each area on each day
simgam.preds = predict(sim.gm$sim.gam$model, newdata = simdf, type = "response")

# Get daily state means across individuals
linds = seq(1, nfish*ndays, by = nfish)
hinds = seq(nfish, nfish*ndays, by = nfish)
simgam.mean = matrix(NA, nrow = ndays, ncol = nstates)
simcar.prob = array(dim = c(ndays,pcs,nstates))
for(i in 1:ndays) {
	simgam.mean[i,] = colMeans(simgam.preds[linds[i]:hinds[i], ])
	simcar.prob[i,,] = simcar.mean[,,i] / sum(simcar.mean[,,i])
}

# Initialize plot
jpeg("Plots/figure_S4.jpeg", width = 2700, height = 5100, units = 'px', res = 600)

# Gray colors for bars
col.bp = hcl.colors(n = nstates, palette = "Grays")

# Plot barplot of state occupancy probability across all days
par(mfrow = c(3,1), mar = c(3,3,1,1), oma = c(3,3,0,4))
smmobs.bp = barplot(smm.prop, space = 0, col = col.bp, cex.axis = 1.2, main = "Actual Proportions", cex.main = 1.5)
axis(1, at = smmobs.bp[c(1,seq(20,100,by=20))], labels = c(1,seq(20,100,by=20)), cex.axis = 1.2)
box()
simgam.bp = barplot(t(simgam.mean), space = 0, col = col.bp, cex.axis = 1.2, main = "GAM", cex.main = 1.5)
axis(1, at = simgam.bp[c(1,seq(20,100,by=20))], labels = c(1,seq(20,100,by=20)), cex.axis = 1.2)
box()
simcar.bp = barplot(t(apply(simcar.prob, c(1,3), sum)), space = 0, col = col.bp, cex.axis = 1.2, main = "CAR", cex.main = 1.5)
axis(1, at = simcar.bp[c(1,seq(20,100,by=20))], labels = c(1,seq(20,100,by=20)), cex.axis = 1.2)
box()
mtext("Day", side = 1, cex = 1.5, outer = T)
mtext("Cumulative Occupancy Probability", side = 2, cex = 1.5, outer = T)

# Legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', sapply(nstates:1, toString), col = rev(col.bp), pch = 15, pt.cex = 3, xpd = TRUE, cex = 1.2, title = "State")

# Reset and turn off device
par(mfrow = c(1,1))
dev.off()


########## XIV. Figure S5 ##########

# Plotting colors
color = palette.colors(n = 4, palette = "Okabe-Ito")
col.traj = apply(col2rgb(color), 2, function(x) rgb(x[1], x[2], x[3], 30, maxColorValue = 255))

# Generate credible intervals for SM trajectories from SMM and HMM
simmat1.med = ceiling(sim.jags$q50$out/pcs)
simmat1.hdi = apply(ceiling(sim.jags$sims.list$out/pcs), c(2,3), hdi)
simmat1.25 = simmat1.hdi[1,,]
simmat1.975 = simmat1.hdi[2,,]
simmat2.med = ceiling(ndhmm.jags$q50$out/pcs)
simmat2.hdi = apply(ceiling(ndhmm.jags$sims.list$out/pcs), c(2,3), hdi)
simmat2.25 = simmat2.hdi[1,,]
simmat2.975 = simmat2.hdi[2,,]

# Get RMSE for HMM
simvar1.rmse = rmse_markov(sim.jags, statseq_smm, actmat)

# Store best/mid/worst indices for native models
sv1.idx = c(which.min(simvar1.rmse), order(simvar1.rmse)[25], which.max(simvar1.rmse))

# Generate confidence intervals for SM from GAM
set.seed(21)
simgam.prob = apply(predict(sim.gm$sim.gam$model, newdata = list(day = rep(seq(ndays),nfish), fish = rep(seq(nfish),each=ndays)), type = "response"), 1, rmultinom, n = 10000, size = 1)
simgam.samples = aperm(array(apply(simgam.prob, 2, function(x) ifelse(which(x == 1) %% nstates == 0, nstates, which(x == 1) %% nstates)), dim = c(10000,ndays,nfish)), c(1,3,2))
simgam.all = apply(simgam.samples, c(2:3), quantile, prob = c(0.025,0.5,0.975))
simgam.traj = list()
simgam.traj$best = simgam.all[,sv1.idx[1],]
simgam.traj$mid = simgam.all[,sv1.idx[2],]
simgam.traj$worst = simgam.all[,sv1.idx[3],]

# Generate confidence intervals for HM from STCAR
set.seed(59)
simcar.all = apply(simcar.prob, 1, rmultinom, n = 10000, size = 1)
simcar.samples = array(apply(simcar.all, 2, function(x) ifelse(which(x == 1) %% (pcs*nstates) == 0, pcs*nstates, which(x == 1) %% (pcs*nstates))), dim = c(10000,ndays))
simcar.traj = apply(simcar.samples, 2, quantile, prob = c(0.025, 0.5, 0.975))

# Initialize plot
jpeg("Plots/figure_S5.jpeg", width = 7500, height = 2500, units = 'px', res = 600)

# For all models applied to SMM data set
par(mfrow = c(1,3))
par(mar = c(5,5,3,1), oma = c(2,2,0,8))
plot(simmat1.med[sv1.idx[1],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 1.5)
axis(2, at = seq(10), cex.axis = 1.5)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat1.25[sv1.idx[1],], rev(simmat1.975[sv1.idx[1],])), border = F, col = col.traj[1])
lines(seq(ndays), simmat2.med[sv1.idx[1],], lwd = 2, lty = 3, col = color[1])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[sv1.idx[1],], rev(simmat2.975[sv1.idx[1],])), border = F, col = col.traj[2])
lines(seq(ndays), simmat2.med[sv1.idx[1],], lwd = 2, lty = 4, col = color[2])
polygon(c(seq(ndays), rev(seq(ndays))), c(simgam.traj$best[1,], rev(simgam.traj$best[3,])), border = F, col = col.traj[3])
lines(seq(ndays), simgam.traj$best[2,], lwd = 2, lty = 5, col = color[3])
polygon(c(seq(ndays), rev(seq(ndays))), c(ceiling(simcar.traj[1,]/pcs), rev(ceiling(simcar.traj[3,]/pcs))), border = F, col = col.traj[4])
lines(seq(ndays), ceiling(simcar.traj[2,]/pcs), lwd = 2, lty = 6, col = color[4])
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv1.idx[1],i]), NA, 9)}
points(seq(ndays), simmat1.med[sv1.idx[1],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv1.idx[1],], lwd = 2)

plot(simmat1.med[sv1.idx[2],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, cex.main = 2, axes = F, xlab = "", ylab = "", main = "SM")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 1.5)
axis(2, at = seq(10), cex.axis = 1.5)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat1.25[sv1.idx[2],], rev(simmat1.975[sv1.idx[2],])), border = F, col = col.traj[1])
lines(seq(ndays), simmat2.med[sv1.idx[2],], lwd = 2, lty = 3, col = color[1])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[sv1.idx[2],], rev(simmat2.975[sv1.idx[2],])), border = F, col = col.traj[2])
lines(seq(ndays), simmat2.med[sv1.idx[2],], lwd = 2, lty = 4, col = color[2])
polygon(c(seq(ndays), rev(seq(ndays))), c(simgam.traj$mid[1,], rev(simgam.traj$mid[3,])), border = F, col = col.traj[3])
lines(seq(ndays), simgam.traj$mid[2,], lwd = 2, lty = 5, col = color[3])
polygon(c(seq(ndays), rev(seq(ndays))), c(ceiling(simcar.traj[1,]/pcs), rev(ceiling(simcar.traj[3,]/pcs))), border = F, col = col.traj[4])
lines(seq(ndays), ceiling(simcar.traj[2,]/pcs), lwd = 2, lty = 6, col = color[4])
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv1.idx[2],i]), NA, 9)}
points(seq(ndays), simmat1.med[sv1.idx[2],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv1.idx[2],], lwd = 2)

plot(simmat1.med[sv1.idx[3],], type = 'n', ylim = c(1,10), xlim = c(1,100), cex.lab = 2.5, axes = F, xlab = "", ylab = "")
axis(1, at = c(1,seq(20,100,by=20)), cex.axis = 1.5)
axis(2, at = seq(10), cex.axis = 1.5)
box()
polygon(c(seq(ndays), rev(seq(ndays))), c(simmat1.25[sv1.idx[3],], rev(simmat1.975[sv1.idx[3],])), border = F, col = col.traj[1])
lines(seq(ndays), simmat2.med[sv1.idx[3],], lwd = 2, lty = 3, col = color[1])
polygon(c(seq(ndays), rev(seq(ndays))), c(hmmmat2.25[sv1.idx[3],], rev(simmat2.975[sv1.idx[3],])), border = F, col = col.traj[2])
lines(seq(ndays), simmat2.med[sv1.idx[3],], lwd = 2, lty = 4, col = color[2])
polygon(c(seq(ndays), rev(seq(ndays))), c(simgam.traj$worst[1,], rev(simgam.traj$worst[3,])), border = F, col = col.traj[3])
lines(seq(ndays), simgam.traj$worst[2,], lwd = 2, lty = 5, col = color[3])
polygon(c(seq(ndays), rev(seq(ndays))), c(ceiling(simcar.traj[1,]/pcs), rev(ceiling(simcar.traj[3,]/pcs))), border = F, col = col.traj[4])
lines(seq(ndays), ceiling(simcar.traj[2,]/pcs), lwd = 2, lty = 6, col = color[4])
ppoint = vector(length = ndays)
for(i in 1:ndays) {ppoint[i] = ifelse(is.na(actmat[sv1.idx[3],i]), NA, 9)}
points(seq(ndays), simmat1.med[sv1.idx[3],], pch = ppoint, lwd = 2)
lines(seq(ndays), statseq_smm[sv1.idx[3],], lwd = 2)

# Add outer margin labels
mtext("State", 2, cex = 2, line = -1, outer = T)
mtext("Day", 1, cex = 2, line = -1, outer = T)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("SMM", "HMM", "GAM", "STCAR", "Data"), col = c(col.traj,color[1]), lwd = 2, pch = c(15,15,15,15,18), lty = c(0,0,0,0,0), pt.cex = c(4,4,4,4,1.5), xpd = TRUE, cex = 1.5)
legend("right", bty = 'n', c("SMM", "HMM", "GAM", "STCAR", "Data"), col = c(color,color[1]), pch = c(26,26,26,26,26), lty = c(3,4,5,6,1), lwd = 2, xpd = TRUE, cex = 1.5)

# Finish plots
dev.off()


########## XV. Figure S6 ##########

# Generate substates for the GAM model
set.seed(582)
gamsim.allp = array(predict(sim.gm$sim.mcl, newdata = mcldf.smm, type = "response"), dim = c(pcs,nstates,nfish,ndays))
gamsim.prob = array(gamsim.allp[cbind(rep(seq(pcs),nfish*ndays),rep(c(simgam.all[2,,]),each=pcs),rep(rep(seq(nfish),ndays),each=pcs),rep(seq(ndays),each=nfish*pcs))], dim = c(pcs,nfish,ndays))
gamsim.subs = apply(gamsim.prob, c(2,3), rmultinom, n = 10000, size = 1)
gamsim.samples = apply(gamsim.subs, c(2,3), function(x) ifelse(which(x == 1) %% pcs == 0, pcs, which(x == 1) %% pcs))
gamsim.all = simgam.samples*pcs - (pcs-gamsim.samples)

# Calculate HMM sub-state probabilities from each model
smmact.subprob = expit(a0 + a1*sand_smm); smmact.subprob = sweep(smmact.subprob, 2, colSums(smmact.subprob), "/")
smmv1.subprob = expit(mean(sim.jags$sims.list$a0) + mean(sim.jags$sims.list$a1)*sand_smm); smmv1.subprob = sweep(smmv1.subprob, 2, colSums(smmv1.subprob), "/")
smmv2.subprob = expit(mean(ndhmm.jags$sims.list$a0) + mean(ndhmm.jags$sims.list$a1)*sand_smm); smmv2.subprob = sweep(smmv2.subprob, 2, colSums(smmv2.subprob), "/")
smmv3.subprob = sweep(apply(gamsim.allp, c(1,2), sum), 2, colSums(apply(gamsim.allp, c(1,2), sum)), "/")
smmv4.subprob = sweep(rowSums(simcar.mean, dims = 2), 2, colSums(rowSums(simcar.mean, dims = 2)), "/")

# Create data frame to house probabilities and model type
smm.subprob = data.frame(prob = c(smmact.subprob,smmv1.subprob,smmv2.subprob,smmv3.subprob,smmv4.subprob), state = rep(seq(nstates),each=pcs), pcs = seq(3))
smm.subprob$type = factor(rep(c("Actual","HMM","SMM","GAM","STCAR"),each=pcs*nstates), levels = c("Actual","HMM","SMM","STCAR","GAM"))

# Plot all probability estimates together
jpeg("Plots/figure_S6.jpeg", width = 6000, height = 3200, units = 'px', res = 600)
colorkey = list(at = seq(0,1,length.out=100), labels = list(at=seq(0,1,0.2),cex = 1.2))
levelplot(prob ~ pcs*state | type, data = smm.subprob,
	xlab = list("Sub-state",cex=1.5), ylab = list("State",cex=1.5), main = list("Probability of occupying each sub-state",cex=1.5), scales = list(cex.axis=1.2,cex=1.2,x=list(at=seq(pcs)),y=list(at=seq(nstates))),
	colorkey = colorkey, col.regions = hcl.colors(100, "Light Grays", rev = T), 
	layout = c(5,1), par.strip.text = list(font=2), par.settings = list(strip.background=list(col="white")))
dev.off()


########## XVI. Figure S7 ##########

# Get RMSE for SMM amd alternative models
simvar1.lnest = log(rmse_markov(sim.jags, statseq_smm, actmat) + 1)
simvar2.lnest = log(rmse_markov(ndhmm.jags, statseq_smm, actmat) + 1)
simvar3.lnest = log(rmse_gc(simgam.samples, statseq_smm, actmat) + 1)
simvar4.lnest = log(rmse_gc(aperm(array(rep(ceiling(simcar.samples/pcs), nfish), dim = c(10000,ndays,nfish)), c(1,3,2)), statseq_smm, actmat) + 1)

# Get RMSE for SMM substates
simvar1.lnesb = log(rmse_markov(sim.jags, statseq_smm*pcs - (pcs-subseq_smm), actmat, "sub") + 1)
simvar2.lnesb = log(rmse_markov(ndhmm.jags, statseq_smm*pcs - (pcs-subseq_smm), actmat, "sub") + 1)
simvar3.lnesb = log(rmse_gc(gamsim.all, statseq_smm*pcs - (pcs-subseq_smm), actmat) + 1)
simvar4.lnesb = log(rmse_gc(aperm(array(rep(simcar.samples, nfish), dim = c(10000,ndays,nfish)), c(1,3,2)), statseq_smm*pcs - (pcs-subseq_smm), actmat) + 1)

jpeg("Plots/figure_S7.jpeg", width = 3500, height = 4800, units = 'px', res = 600)
par(mfrow = c(2,1), mar = c(3,5,2,1), oma = c(2,1,0,5))
plot(simvar1.lnest, simvar2.lnest, pch = 22, cex = 2, bg = "gray20", xlab = "", ylab = "", main = "State", cex.main = 1.5, xlim = c(0.75,1.75), ylim = c(0,3.5), axes = F); axis(1, cex.axis = 1.5); axis(2, at = seq(0,3), cex.axis = 1.5); box()
points(simvar1.lnest, simvar3.lnest, pch = 22, cex = 2, bg = "gray50")
points(simvar1.lnest, simvar4.lnest, pch = 22, cex = 2, bg = "gray80")
abline(0,1)
plot(simvar2.lnesb, simvar1.lnesb, pch = 22, cex = 2, bg = "gray20", xlab = "", ylab = "", main = "Sub-state", cex.main = 1.5, xlim = c(2,4), ylim = c(0,6), cex.axis = 1.5)
points(simvar2.lnesb, simvar3.lnesb, pch = 22, cex = 2, bg = "gray50")
points(simvar2.lnesb, simvar4.lnesb, pch = 22, cex = 2, bg = "gray80")
abline(0,1)
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
mtext(expression("Native Model RMSE"["s"["k"]]), side = 1, line = -1, cex = 1.5, outer = T)
mtext(expression("Alternative Model RMSE"["s"["k"]]), side = 2, line = -2, cex = 1.5, outer = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("HMM","GAM","STCAR"), pt.bg = c("gray20","gray50","gray80"), pch = 22, pt.cex = 1.5, xpd = TRUE, cex = 1.2)
dev.off()


########## XVII. Figure S8 ##########

# Calculate hpdi of a's and b's
a0.mark2 = hdi(hmmsim.jags$sims.list$a0)
a0.hmm2 = hdi(hmm.jags$sims.list$a0)
a0.car2 = hdi(c(hmm.car[[1]]$samples$beta[,1],hmm.car[[2]]$samples$beta[,1],hmm.car[[3]]$samples$beta[,1]))
a1.mark2 = hdi(hmmsim.jags$sims.list$a1)
a1.hmm2 = hdi(hmm.jags$sims.list$a1)
a1.car2 = hdi(c(hmm.car[[1]]$samples$beta[,2],hmm.car[[2]]$samples$beta[,2],hmm.car[[3]]$samples$beta[,2]))

# Initialize plot
jpeg("Plots/figure_S8.jpeg", width = 8000, height = 5000, units = 'px', res = 600)
alay = layout(matrix(c(1,1,1,2,3,4,5,5,5,6,7,8),4,3,T), rep(1,2), c(1,15,1,15)) 
layout.show(alay)

# Add title to first row
par(mar = c(0,0,0,0))
frame()
text(0.5,0.5,"a.",font=2,cex=2)

# Plot density for a0
par(mar = c(5,3,3,1))
plot(density(hmm.jags$sims.list$a0), xlab = expression(bold(a[0] ~ "- HMM*")), ylab = "", font.lab = 2, 
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

par(mar = c(5,3,3,1))
plot(density(hmmsim.jags$sims.list$a0), xlab = expression(a[0] ~ "- SMM"), ylab = "",
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

par(mar = c(5,3,3,1))
plot(density(c(hmm.car[[1]]$samples$beta[,1],hmm.car[[2]]$samples$beta[,1],hmm.car[[3]]$samples$beta[,1])), xlab = expression(a[0] ~ "- CAR"), ylab = "",
		xlim = c(-8,-5), ylim = c(0,2.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.5,0.5))
box()
polygon(density(c(hmm.car[[1]]$samples$beta[,1],hmm.car[[2]]$samples$beta[,1],hmm.car[[3]]$samples$beta[,1])), lwd = 2, col = "gray70")
arrows(a0.car2[1], 0.72*2.5, a0.car2[2], 0.72*2.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.car2), 0.8*2.5, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.5), lwd = 2, col = "gray40")
lines(rep(median(c(hmm.car[[1]]$samples$beta[,1],hmm.car[[2]]$samples$beta[,1],hmm.car[[3]]$samples$beta[,1])),2), c(0,0.72*2.5), lwd = 2, lty = 3, col = "gray40")

# Add title to first row
par(mar = c(0,0,0,0))
frame()
text(0.5,0.5,"b.",font=2,cex=2)

# Plot density for a1
par(mar = c(5,3,3,1))
plot(density(hmm.jags$sims.list$a1), xlab = expression(bold(a[1] ~ "- HMM*")), ylab = "", font.lab = 2, 
		xlim = c(6,16), ylim = c(0,1.0), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,1.0,0.2))
box()
polygon(density(hmm.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.hmm2[1], 0.72*1.0, a1.hmm2[2], 0.72*1.0, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.hmm2), 0.8*1.0, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*1.0), lwd = 2, col = "gray40")
lines(c(hmm.jags$q50$a1,hmm.jags$q50$a1), c(0,0.72*0.8), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,3,3,1))
plot(density(hmmsim.jags$sims.list$a1), xlab = expression(a[1] ~ "- SMM"), ylab = "",
		xlim = c(6,16), ylim = c(0,1.0), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,1.0,0.2))
box()
polygon(density(hmmsim.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.mark2[1], 0.72*1.0, a1.mark2[2], 0.72*1.0, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.mark2), 0.8*1.0, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*1.0), lwd = 2, col = "gray40")
lines(c(hmmsim.jags$q50$a1,hmmsim.jags$q50$a1), c(0,0.72*1.0), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,3,3,1))
plot(density(c(hmm.car[[1]]$samples$beta[,2],hmm.car[[2]]$samples$beta[,2],hmm.car[[3]]$samples$beta[,2])), xlab = expression(a[1] ~ "- CAR"), ylab = "",
		xlim = c(0,10), ylim = c(0,1.0), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,1.0,0.2))
box()
polygon(density(c(hmm.car[[1]]$samples$beta[,2],hmm.car[[2]]$samples$beta[,2],hmm.car[[3]]$samples$beta[,2])), lwd = 2, col = "gray70")
arrows(a1.car2[1], 0.72*1.0, a1.car2[2], 0.72*1.0, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.car2), 0.8*1.0, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*1.0), lwd = 2, col = "gray40")
lines(rep(median(c(hmm.car[[1]]$samples$beta[,2],hmm.car[[2]]$samples$beta[,2],hmm.car[[3]]$samples$beta[,2])),2), c(0,0.72*0.8), lwd = 2, lty = 3, col = "gray40")

# Reset and turn off device
par(mfrow = c(1,1))
dev.off()


########## XVIII. Figure S9 ##########

# Calculate hpdi of a's and b's
a0.mark1 = hdi(sim.jags$sims.list$a0)
a0.hmm1 = hdi(ndhmm.jags$sims.list$a0)
a0.car1 = hdi(c(sim.car[[1]]$samples$beta[,1],sim.car[[2]]$samples$beta[,1],sim.car[[3]]$samples$beta[,1]))
a1.mark1 = hdi(sim.jags$sims.list$a1)
a1.hmm1 = hdi(ndhmm.jags$sims.list$a1)
a1.car1 = hdi(c(sim.car[[1]]$samples$beta[,2],sim.car[[2]]$samples$beta[,2],sim.car[[3]]$samples$beta[,2]))

# Initialize plot
jpeg("Plots/figure_S8.jpeg", width = 8000, height = 5000, units = 'px', res = 600)
alay = layout(matrix(c(1,1,1,2,3,4,5,5,5,6,7,8),4,3,T), rep(1,2), c(1,15,1,15)) 
layout.show(alay)

# Add title to first row
par(mar = c(0,0,0,0))
frame()
text(0.5,0.5,"a.",font=2,cex=2)

# Plot density for a0
par(mar = c(5,3,3,1))
plot(density(ndhmm.jags$sims.list$a0), xlab = expression(a[0] ~ "- HMM"), ylab = "",
		xlim = c(-7,-3), ylim = c(0,2.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.5,0.5))
box()
polygon(density(ndhmm.jags$sims.list$a0), lwd = 2, col = "gray70")
arrows(a0.hmm1[1], 0.72*2.5, a0.hmm1[2], 0.72*2.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.hmm1), 0.8*2.5, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.5), lwd = 2, col = "gray40")
lines(c(ndhmm.jags$q50$a0,ndhmm.jags$q50$a0), c(0,0.72*2.0), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,3,3,1))
plot(density(sim.jags$sims.list$a0), xlab = expression(bold(a[0] ~ "- SMM*")), ylab = "",
		xlim = c(-7,-4), ylim = c(0,2.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.0,0.5))
box()
polygon(density(sim.jags$sims.list$a0), lwd = 2, col = "gray70")
arrows(a0.mark1[1], 0.72*2.5, a0.mark1[2], 0.72*2.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.mark1), 0.8*2.5, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.5), lwd = 2, col = "gray40")
lines(c(sim.jags$q50$a0,sim.jags$q50$a0), c(0,0.72*2.0), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,3,3,1))
plot(density(c(sim.car[[1]]$samples$beta[,1],sim.car[[2]]$samples$beta[,1],sim.car[[3]]$samples$beta[,1])), xlab = expression(a[0] ~ "- STCAR"), ylab = "",
		xlim = c(-8,-5), ylim = c(0,2.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,2.5,0.5))
box()
polygon(density(c(sim.car[[1]]$samples$beta[,1],sim.car[[2]]$samples$beta[,1],sim.car[[3]]$samples$beta[,1])), lwd = 2, col = "gray70")
arrows(a0.car1[1], 0.72*2.5, a0.car1[2], 0.72*2.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a0.car1), 0.8*2.5, "95% HDI", cex = 1.8)
lines(c(a0,a0), c(0,0.72*2.5), lwd = 2, col = "gray40")
lines(rep(median(c(sim.car[[1]]$samples$beta[,1],sim.car[[2]]$samples$beta[,1],sim.car[[3]]$samples$beta[,1])),2), c(0,0.72*2.5), lwd = 2, lty = 3, col = "gray40")

# Add title to second row
par(mar = c(0,0,0,0))
frame()
text(0.5,0.5,"b.",font=2,cex=2)

# Plot density for a1
par(mar = c(5,3,3,1))
plot(density(ndhmm.jags$sims.list$a1), xlab = expression(a[1] ~ "- HMM"), ylab = "",
		xlim = c(6,16), ylim = c(0,1.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,1.5,0.3))
box()
polygon(density(ndhmm.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.hmm1[1], 0.72*1.5, a1.hmm1[2], 0.72*1.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.hmm1), 0.8*1.5, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*1.5), lwd = 2, col = "gray40")
lines(c(ndhmm.jags$q50$a1,ndhmm.jags$q50$a1), c(0,0.72*1.5), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,3,3,1))
plot(density(sim.jags$sims.list$a1), xlab = expression(bold(a[1] ~ "- SMM*")), ylab = "",
		xlim = c(6,16), ylim = c(0,1.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,1.5,0.3))
box()
polygon(density(sim.jags$sims.list$a1), lwd = 2, col = "gray70")
arrows(a1.mark1[1], 0.72*1.5, a1.mark1[2], 0.72*1.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.mark1), 0.8*1.5, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*1.5), lwd = 2, col = "gray40")
lines(c(sim.jags$q50$a1,sim.jags$q50$a1), c(0,0.72*1.5), lwd = 2, lty = 3, col = "gray40")

par(mar = c(5,3,3,1))
plot(density(c(sim.car[[1]]$samples$beta[,2],sim.car[[2]]$samples$beta[,2],sim.car[[3]]$samples$beta[,2])), xlab = expression(a[1] ~ "- STCAR"), ylab = "",
		xlim = c(2,10), ylim = c(0,1.5), 
		lwd = 2, main = "", cex.lab = 1.8, axes = F)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, at = seq(0,1.5,0.3))
box()
polygon(density(c(sim.car[[1]]$samples$beta[,2],sim.car[[2]]$samples$beta[,2],sim.car[[3]]$samples$beta[,2])), lwd = 2, col = "gray70")
arrows(a1.car1[1], 0.72*1.5, a1.car1[2], 0.72*1.5, lwd = 2, angle = 90, length = 0.1, code = 3)
text(mean(a1.car2), 0.8*1.5, "95% HDI", cex = 1.8)
lines(c(a1,a1), c(0,0.72*1.5), lwd = 2, col = "gray40")
lines(rep(median(c(sim.car[[1]]$samples$beta[,2],sim.car[[2]]$samples$beta[,2],sim.car[[3]]$samples$beta[,2])),2), c(0,0.72*1.5), lwd = 2, lty = 3, col = "gray40")

# Reset and turn off device
par(mfrow = c(1,1))
dev.off()