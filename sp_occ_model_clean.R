#############################################################################################
# Spatially Explicit Dynamic Occupancy Model For the Florida Bog Frog on Eglin Air Force Base

# libraries
library(jagsUI)
library(wiqid)
library(tidyverse)
library(geosphere)

# load in detection data and site data
dat = read.csv("bf_surveys.csv")
covs = read.csv("bf_covs.csv")

# calculate distances between sites
dmat = distm(covs[,6:7],fun = distHaversine)

# store first and last years
year1 = min(dat$Year)
yearn = max(dat$Year)
n.years = length(unique(dat$Year))

# create column for number of surveys within years
dat = dat %>%
  group_by(Year, SiteID) %>%
  mutate(survey = 1:n())

# spread columns by year
dat = spread(dat, key = Year, value = RANOKA)

# extract the first 3 surveys
dat = dat %>%
  filter(survey < 4)

# spread by surveys within years
dat = pivot_wider(
  dat, 
  id_cols = 'SiteID', 
  names_from = 'survey', 
  values_from = starts_with("20"), 
  names_glue = '{.value}.{survey}'
)

# Extract the detection histories
DH <- as.matrix(dat[, 2:52])
head(DH)

# convert to sites x occasions x year array
Y <- array(DH, dim=c(nrow(DH), 3, n.years))
head(Y)

# Aggregate detection histories across occasions
y <- apply(Y, c(1, 3), sum, na.rm=TRUE)  # sites by years, number of detections
n <- apply(!is.na(Y), c(1, 3), sum)      # ... and no. of surveys

# specify known z's 
z <- (y > 0)*1
# specifiy unknown occupancy for 0 detections
z[z == 0] <- NA


# standardize covariates
burnS <- wiqid::standardize(covs$hist_fire)
sinS <- wiqid::standardize(covs$sinu)
upS <- wiqid::standardize(covs$upland)
scrubS <- wiqid::standardize(covs$scrub)
wetS <- wiqid::standardize(covs$wetland)
diffS <- wiqid::standardize(covs$current_fire)

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y, method = "spearman")) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

pairs(covs[,2:5], upper.panel = panel.cor)

############
# indicator for distance (prevent self colonization)
indica = dmat
indica[which(indica > 0)] = 1

#############
# JAGS model
#############

sink("bf_sp_occ.txt")
cat("

model {

 # annual random effect of detection

for (k in 1:nYears) {

  logit(p[k]) <- beta.p0 + beta.p1[k]

  beta.p1[k] ~ dnorm(0, p_tau) # Random intercepts

}      

  for (k in 2:nYears) {
  
  logit(eps[k-1]) <- beta.eps0 + beta.eps1[k-1] # extinction prob
  
  beta.eps1[k-1] ~ dnorm(0, eps_tau) # Random intercepts

  } 
 
# Priors
  beta.eps0 ~ dlogis(0, 1) # extinction prior
  beta.psi0 ~ dlogis(0, 1) # extinction prior
  beta.p0 ~ dlogis(0,1)    # detection prior
  beta.rho0 ~ dgamma(1,1)  # colonization prior

  beta.psi1 ~ dnorm(0,0.1)
  beta.psi2 ~ dnorm(0,0.1)
  beta.psi3 ~ dnorm(0,0.1)
  beta.psi4 ~ dnorm(0,0.1)
  beta.psi5 ~ dnorm(0,0.1)

  beta.rho1 ~ dnorm(0, 0.1)


  p_tau <- 1/(p_sigma*p_sigma)
  p_sigma ~ dgamma(1,1)
  eps_tau <- 1/(eps_sigma*eps_sigma)
  eps_sigma ~ dgamma(1,1)


  # Initial State and Site Connectivities
  for (i in 1:nSites) {
    
    logit(psi[i]) <- beta.psi0 + beta.psi1*burn.f[i] + beta.psi2*sin[i] + beta.psi3*scrub[i] + beta.psi4*wetland[i] + beta.psi5*upland[i]
    z[i, 1] ~ dbern(psi[i])

  # baseline colonization rate
    theta[i] <- beta.rho0 + beta.rho1*burn.c[i]



    for (j in 1:nSites) {
      con[i, j] <- exp(-theta[i] * dmat[i, j]) * indica[i,j]
    }
  }

  # Ecological submodel: Define state conditional on parameters
  for (k in 2:nYears) {
    for (i in 1:nSites) {
      for (j in 1:nSites) {
        disp[i, j, k-1] <- z[j, k-1] * con[i, j]
      } #j

      s[i, k-1] <- sum(disp[i, , k-1])
      gamma[i, k-1] <- 1 - exp(-s[i, k-1])


       # Probability of being occupied at time k
       Ez[i,k-1] <- gamma[i,k-1]*(1-z[i,k-1])+(1-(eps[k-1]*1-gamma[i,k]))*z[i,k-1]
       
       z[i,k] ~ dbern(Ez[i,k-1])
      
     } #k
   } #i


  # Observation model

  for (i in 1:nSites) {
      for (k in 1:nYears) {
      y[i, k] ~ dbin(p[k] * z[i, k], n[i, k])
      } #k
  } #i
  
  # Derived variable
  for(k in 1:nYears) {
    N[k] <- sum(z[,k]) # no. sites occupied for each year
  }

} # model

",fill = TRUE)
sink()


# bundle data for jags
jdata <- list(nSites = nrow(y), nYears = ncol(y), y = y,
              n = n, z = z, dmat = (dmat/1000), 
              burn.f = burnS, burn.c = diffS,
              sin = sinS, wetland = wetS, upland = upS, scrub = scrubS,
              indica = indica
              )

# parameters to monitor
wanted <- c("beta.psi0", "beta.psi1","beta.psi2","beta.psi3","beta.psi4","beta.psi5",
            "beta.p0", "p_sigma", "p",
            "beta.eps0", "eps_sigma", "eps",
            "beta.rho0", "beta.rho1",
            "N")

# run specs
nc = 3
ni = 11000
na = 1000
nb = 1000
nt = 50


#Function to create a matrix with initial values for unknown latent states, z[i,t]
inits.state.dynocc <- function(obsarray){
  state <- apply(obsarray,c(1,3),function(x) ifelse(sum(is.na(x))==length(x),NA,max(x,na.rm=T)))
  #Initial value of 1 whenever occupancy state is unknown
  state[state==1] <- 2
  state[is.na(state) | state==0] <- 1
  state[state==2] <- NA
  return(state)
}  

# intial values
inits <- function (){
  list (beta.psi0 = runif(1, -1, -1),beta.p0 = runif(1, -1, 1),beta.eps0 = runif(1, -1, 1),beta.rho0 = runif(1, 1, 5), 
        beta.psi1 = runif(1, -1, 1),beta.psi2 = runif(1, -1, 1),beta.psi3 = runif(1, -1, 1),beta.psi4 = runif(1, -1, 1),beta.psi5 = runif(1, -1, 1), 
        beta.rho1 = runif(1, -1, 1),
        eps_sigma = runif(1, 0.1, 0.2), p_sigma = runif(1, 0.1, 0.2),
        z = inits.state.dynocc(Y))
}


# Run the model
jagsOut <- jags(data = jdata, inits = inits, parameters.to.save = wanted, model.file = "bf_sp_occ.txt",
                n.chains=nc, n.iter=ni, n.adapt=na, n.burnin = nb, n.thin = nt, DIC=FALSE, parallel=TRUE)

traceplot(jagsOut)

###############
# OUTPUT

#occupied sites
n = c("N[1]", "N[2]", "N[3]", "N[4]", "N[5]", "N[6]", "N[7]", "N[8]", "N[9]", "N[10]", "N[11]", "N[12]", "N[13]",  "N[14]", "N[15]", "N[16]", "N[17]")

pdf("bf_occ_N.pdf", width = 7, height = 5)

op <- par(mfrow=c(1,1), mar=c(4,4,1,1)+0.1, pch=4, lwd=2, las=1)

# Number of Occupied Sites
toplot <- jagsOut$summary[n, c(1,3,7)]

plotsegraph <- function(loc, lower, upper, wiskwidth, color = "grey", linewidth = 2) {
  
  w <- wiskwidth/2
  segments(x0 = loc, x1 = loc, y0 = lower, y1 = upper, col = color, 
           lwd = linewidth)
  segments(x0 = loc - w, x1 = loc + w, y0 = upper, y1 = upper, 
           col = color, lwd = linewidth)  # upper whiskers
  segments(x0 = loc - w, x1 = loc + w, y0 = lower, y1 = lower, 
           col = color, lwd = linewidth)  # lower whiskers
}

op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(year1:yearn, toplot[,1], col = "black", pch = 21, bg = "black", cex = 1.5,
     xlim = c(year1, yearn), ylim = c(20,70), ylab = "", xlab = "", axes = FALSE)
lines(year1:yearn, toplot[,1], lwd = 2, type = "c")
plotsegraph(year1:yearn, toplot[,2], toplot[,3], 0.2, color = "black")
axis(1, year1:yearn)
axis(2) 

par(las = 0)
mtext("Year", side = 1, line = 2.5, cex = 1.5)
mtext("Number of Occupied Sites", side = 2, line = 3.7, cex = 1.5)

dev.off()

################
# detection

n = c("p[1]", "p[2]", "p[3]", "p[4]", "p[5]", "p[6]", "p[7]", "p[8]", "p[9]", "p[10]", "p[11]", "p[12]", "p[13]",  "p[14]", "p[15]", "p[16]", "p[17]")

pdf("bf_occ_p.pdf", width = 7, height = 5)
op <- par(mfrow=c(1,1), mar=c(4,4,1,1)+0.1, pch=4, lwd=2, las=1)

# Number of Occupied Sites
toplot <- jagsOut$summary[n, c(1,3,7)]

plotsegraph <- function(loc, lower, upper, wiskwidth, color = "grey", linewidth = 2) {
  
  w <- wiskwidth/2
  segments(x0 = loc, x1 = loc, y0 = lower, y1 = upper, col = color, 
           lwd = linewidth)
  segments(x0 = loc - w, x1 = loc + w, y0 = upper, y1 = upper, 
           col = color, lwd = linewidth)  # upper whiskers
  segments(x0 = loc - w, x1 = loc + w, y0 = lower, y1 = lower, 
           col = color, lwd = linewidth)  # lower whiskers
}

op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(year1:yearn, toplot[,1], col = "black", pch = 21, bg = "black", cex = 1.5,
     xlim = c(year1, yearn), ylim = c(0,1), ylab = "", xlab = "", axes = FALSE)

plotsegraph(year1:yearn, toplot[,2], toplot[,3], 0.2, color = "black")
axis(1, year1:yearn)
axis(2) 

par(las = 0)
mtext("Year", side = 1, line = 2.5, cex = 1.5)
mtext("Detection Probability", side = 2, line = 3.7, cex = 1.5)

dev.off()

#extinction
n = c("eps[1]", "eps[2]", "eps[3]", "eps[4]", "eps[5]", "eps[6]", "eps[7]", "eps[8]", "eps[9]", "eps[10]", "eps[11]", "eps[12]", "eps[13]",  "eps[14]", "eps[15]", "eps[16]")

pdf("bf_occ_eps.pdf", width = 7, height = 5)

op <- par(mfrow=c(1,1), mar=c(4,4,1,1)+0.1, pch=4, lwd=2, las=1)

# Number of Occupied Sites
toplot <- jagsOut$summary[n, c(1,3,7)]

plotsegraph <- function(loc, lower, upper, wiskwidth, color = "grey", linewidth = 2) {
  
  w <- wiskwidth/2
  segments(x0 = loc, x1 = loc, y0 = lower, y1 = upper, col = color, 
           lwd = linewidth)
  segments(x0 = loc - w, x1 = loc + w, y0 = upper, y1 = upper, 
           col = color, lwd = linewidth)  # upper whiskers
  segments(x0 = loc - w, x1 = loc + w, y0 = lower, y1 = lower, 
           col = color, lwd = linewidth)  # lower whiskers
}

op <- par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot((year1:(yearn-1))+0.5, toplot[,1], col = "black", pch = 21, bg = "black", cex = 1.5,
     xlim = c(year1, yearn), ylim = c(0,1), ylab = "", xlab = "", axes = FALSE)

plotsegraph((year1:(yearn-1))+0.5, toplot[,2], toplot[,3], 0.2, color = "black")
axis(1, year1:yearn)
axis(2) 

par(las = 0)
mtext("Year", side = 1, line = 2.5, cex = 1.5)
mtext("Extinction Probability", side = 2, line = 3.7, cex = 1.5)

dev.off()

####################
# Initial Occupancy
####################

####################################
# burn history on initial occupancy
range(covs$burn_interval, na.rm = T)
bh <- seq(2, 10, length=101)
bhS <- wiqid::standardize2match(bh, covs$hist_fire)
bhp1 = c("beta.psi0", "beta.psi1")

toplot <- jagsOut$summary[bhp1, c(1,4,6)]

library(boot)
lpsix <- inv.logit(toplot[1,1] + toplot[2,1] * bhS)
lpsixL <- inv.logit(toplot[1,3] + toplot[2,2] * bhS)
lpsixU <- inv.logit(toplot[1,2] + toplot[2,3] * bhS)

psix <- cbind(plogis(lpsix), plogis(lpsixL), plogis(lpsixU))

pred.dat <- data.frame(bh, preds = lpsix, predsL = lpsixL, predsU = lpsixU)

g1 = ggplot(pred.dat, aes(x = bh, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Historic fire return interval (years)") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()


####################################
# sinuosity on initial occupancy
range(covs$sinu, na.rm = T)
bh <- seq(1, 1.8, length=101)
bhS <- wiqid::standardize2match(bh, covs$sinu)
bhp1 = c("beta.psi0", "beta.psi2")

toplot <- jagsOut$summary[bhp1, c(1,3,7)]

lpsix <- inv.logit(toplot[1,1] + toplot[2,1] * bhS)
lpsixL <- inv.logit(toplot[1,2] + toplot[2,3] * bhS)
lpsixU <- inv.logit(toplot[1,3] + toplot[2,2] * bhS)

psix <- cbind(plogis(lpsix), plogis(lpsixL), plogis(lpsixU))

pred.dat <- data.frame(bh, preds = lpsix, predsL = lpsixL, predsU = lpsixU)

g2 = ggplot(pred.dat, aes(x = bh, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Stream sinuosity") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

####################################
# wetlands on initial occupancy
range(covs$wetland, na.rm = T)
bh <- seq(0, 1, length=101)
bhS <- wiqid::standardize2match(bh, covs$wetland)
bhp1 = c("beta.psi0", "beta.psi4")

toplot <- jagsOut$summary[bhp1, c(1,3,7)]

lpsix <- inv.logit(toplot[1,1] + toplot[2,1] * bhS)
lpsixL <- inv.logit(toplot[1,2] + toplot[2,3] * bhS)
lpsixU <- inv.logit(toplot[1,3] + toplot[2,2] * bhS)

psix <- cbind(plogis(lpsix), plogis(lpsixL), plogis(lpsixU))

pred.dat <- data.frame(bh, preds = lpsix, predsL = lpsixL, predsU = lpsixU)

g3 = ggplot(pred.dat, aes(x = bh, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Proportion wetland") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

####################################
# upland on initial occupancy
range(covs$upland, na.rm = T)
bh <- seq(0, 1, length=101)
bhS <- wiqid::standardize2match(bh, covs$upland)
bhp1 = c("beta.psi0", "beta.psi5")

toplot <- jagsOut$summary[bhp1, c(1,3,7)]

lpsix <- inv.logit(toplot[1,1] + toplot[2,1] * bhS)
lpsixL <- inv.logit(toplot[1,2] + toplot[2,3] * bhS)
lpsixU <- inv.logit(toplot[1,3] + toplot[2,2] * bhS)

psix <- cbind(plogis(lpsix), plogis(lpsixL), plogis(lpsixU))

pred.dat <- data.frame(bh, preds = lpsix, predsL = lpsixL, predsU = lpsixU)

g4 = ggplot(pred.dat, aes(x = bh, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Proportion upland") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

####################################
# flatwoods on initial occupancy
range(covs$scrub, na.rm = T)
bh <- seq(0, 0.6, length=101)
bhS <- wiqid::standardize2match(bh, covs$scrub)
bhp1 = c("beta.psi0", "beta.psi3")

toplot <- jagsOut$summary[bhp1, c(1,3,7)]

lpsix <- inv.logit(toplot[1,1] + toplot[2,1] * bhS)
lpsixL <- inv.logit(toplot[1,2] + toplot[2,3] * bhS)
lpsixU <- inv.logit(toplot[1,3] + toplot[2,2] * bhS)

psix <- cbind(plogis(lpsix), plogis(lpsixL), plogis(lpsixU))

pred.dat <- data.frame(bh, preds = lpsix, predsL = lpsixL, predsU = lpsixU)

g5 = ggplot(pred.dat, aes(x = bh, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Proportion flatwoods") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

library(ggpubr)
gg1 = ggarrange(g2, g3, g4, g5, nrow = 2, ncol = 2)

ggsave("psiburn_preds.pdf", g1, width = 4, height = 4, units = "in", dpi = 300)
ggsave("psiInit_preds.pdf", gg1, width = 5, height = 5, units = "in", dpi = 300)


##############################
# colonization probability
#############################

#################################################################################
# plot predicted relationship between distance and colonization probability

bh <- seq(0, 8, length=101)
bhS <- wiqid::standardize2match(bh, covs$current_fire)
d = 0.5 # scale to 500 m

bhp1 = c("beta.rho0", "beta.rho1")
toplot <- jagsOut$summary[bhp1, c(1,3,7)]

lpsix <- exp(-(toplot[1,1] + toplot[2,1] * bhS) * d)
lpsixL <- exp(-(toplot[1,2] + toplot[2,2] * bhS) * d)
lpsixU <- exp(-(toplot[1,3] + toplot[2,3] * bhS) * d)

gam.mu = 1-exp(-lpsix)
gam.L = 1-exp(-lpsixL)
gam.U = 1-exp(-lpsixU)

pred.dat <- data.frame(bh = bh, preds = gam.mu, predsL = gam.L, predsU = gam.U)

g1 = ggplot(pred.dat, aes(x = bh, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Number of burns during study period (2006-2022)") + 
  ylim(0,0.5) + ylab("Colonization probability") +
  theme_classic()


ggsave("col_fire_pred.pdf", g1, width = 5, height = 4, units = "in", dpi = 300)


#########################################
# plot impact of distance on colonization

bhp1 = c("beta.rho0")

toplot <- jagsOut$summary[bhp1, c(1,3,7)]

dist.pred = seq(0, 2, length.out = 100)

lpsix <- exp(-toplot[1]*dist.pred)
lpsixL <- exp(-toplot[2]*dist.pred)
lpsixU <- exp(-toplot[3]*dist.pred)

gam.mu = 1-exp(-lpsix)
gam.L = 1-exp(-lpsixL)
gam.U = 1-exp(-lpsixU)

pred.dat <- data.frame(d = dist.pred, preds = gam.mu, predsL = gam.L, predsU = gam.U)

g1 = ggplot(pred.dat, aes(x = d, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Distance (km)") + 
  ylim(0,1) + ylab("Colonization probability") +
  theme_classic()

ggsave("col_dist_pred.pdf", g1, width = 5, height = 4, units = "in", dpi = 300)









load("bf_occ_r3.RData")
mcmc = jagsOut$sims.list

library(broom.mixed)
library(dplyr)
library(ggplot2)
library(boot)
library(ggpubr)


####################
# Initial Occupancy
####################

####################################
# burn history on initial occupancy
bh <- seq(2, 10, length=101)
bhS <- wiqid::standardize2match(bh, covs$hist_fire)

## Calculate the fitted values
newdata = data.frame(x = bhS)
Xmat = model.matrix(~x, newdata)
coefs = cbind(mcmc[["beta.psi0"]], mcmc[["beta.psi1"]])
fit = inv.logit(coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), bh)

g1 = ggplot(newdata, aes(x = bh, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Historical fire return interval (years)") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()


####################################
# sinuosity on initial occupancy
bh <- seq(1, 1.8, length=101)
bhS <- wiqid::standardize2match(bh, covs$sinu)

## Calculate the fitted values
newdata = data.frame(x = bhS)
Xmat = model.matrix(~x, newdata)
coefs = cbind(mcmc[["beta.psi0"]], mcmc[["beta.psi2"]])
fit = inv.logit(coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), bh)


g2 = ggplot(newdata, aes(x = bh, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Stream sinuosity") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

####################################
# wetlands on initial occupancy
bh <- seq(0, 1, length=101)
bhS <- wiqid::standardize2match(bh, covs$wetland)

## Calculate the fitted values
newdata = data.frame(x = bhS)
Xmat = model.matrix(~x, newdata)
coefs = cbind(mcmc[["beta.psi0"]], mcmc[["beta.psi4"]])
fit = inv.logit(coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), bh)


g3 = ggplot(newdata, aes(x = bh, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Proportion wetland") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

####################################
# upland on initial occupancy
bh <- seq(0, 1, length=101)
bhS <- wiqid::standardize2match(bh, covs$upland)

## Calculate the fitted values
newdata = data.frame(x = bhS)
Xmat = model.matrix(~x, newdata)
coefs = cbind(mcmc[["beta.psi0"]], mcmc[["beta.psi5"]])
fit = inv.logit(coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), bh)

g4 = ggplot(newdata, aes(x = bh, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Proportion upland") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

####################################
# flatwoods on initial occupancy
bh <- seq(0, 0.6, length=101)
bhS <- wiqid::standardize2match(bh, covs$scrub)

## Calculate the fitted values
newdata = data.frame(x = bhS)
Xmat = model.matrix(~x, newdata)
coefs = cbind(mcmc[["beta.psi0"]], mcmc[["beta.psi3"]])
fit = inv.logit(coefs %*% t(Xmat))
newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), bh)

g5 = ggplot(newdata, aes(x = bh, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Proportion flatwoods") + 
  ylim(0,1) + ylab("Initial occupancy probability") +
  theme_classic()

gg1 = ggarrange(g2, g3, g4, g5, nrow = 2, ncol = 2)

ggsave("psiburn_preds.pdf", g1, width = 4, height = 4, units = "in", dpi = 300, device = cairo_pdf)
ggsave("psiInit_preds.pdf", gg1, width = 5, height = 5, units = "in", dpi = 300, device = cairo_pdf)


##############################
# colonization probability
#############################

#################################################################################
# plot predicted relationship between distance and colonization probability

bh <- seq(0, 8, length=101)
bhS <- wiqid::standardize2match(bh, covs$current_fire)
d = 0.5 # scale to 500 m

## Calculate the fitted values
newdata = data.frame(x = bhS)
Xmat = model.matrix(~x, newdata)
coefs = cbind(mcmc[["beta.rho0"]], mcmc[["beta.rho1"]])
fit = exp(-coefs %*% t(Xmat) * d)
fit = 1-exp(-fit)

newdata = newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = "quantile"), bh)


g1 = ggplot(newdata, aes(x = bh, y = estimate))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  xlab("Number of burns during study period (2006-2022)") + 
  ylim(0,0.5) + ylab("Colonization probability") +
  theme_classic()


ggsave("col_fire_pred.pdf", g1, width = 5, height = 4, units = "in", dpi = 300, device = cairo_pdf)


#########################################
# plot impact of distance on colonization

bhp1 = c("beta.rho0")

toplot <- jagsOut$summary[bhp1, c(1,3,7)]

dist.pred = seq(0, 2, length.out = 100)

lpsix <- exp(-toplot[1]*dist.pred)
lpsixL <- exp(-toplot[2]*dist.pred)
lpsixU <- exp(-toplot[3]*dist.pred)

gam.mu = 1-exp(-lpsix)
gam.L = 1-exp(-lpsixL)
gam.U = 1-exp(-lpsixU)

pred.dat <- data.frame(d = dist.pred, preds = gam.mu, predsL = gam.L, predsU = gam.U)

g1 = ggplot(pred.dat, aes(x = d, y = preds))+
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = predsL, ymax = predsU), alpha = 0.3) +
  xlab("Distance (km)") + 
  ylim(0,1) + ylab("Colonization probability") +
  theme_classic()

ggsave("col_dist_pred.pdf", g1, width = 5, height = 4, units = "in", dpi = 300)

