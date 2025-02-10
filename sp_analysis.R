library(AICcmodavg)
library(cvAUC)
library(ISLR)
library(boot)
library(ncf)
library(spdep)
library(codep)
library(dismo)
library(raster)
library(gbm)
library(ggplot2)
library(wiqid)
library(tidyverse)
library(geosphere)
library(dismo)
library(fields)

#------------------
# SAC analyses
#------------------
# load in detection data and site data
dat <- read.csv("bf_surveys.csv")
covs = read.csv("bf_site_dat.csv")

# create column for number of surveys within years
dat = dat %>%
  group_by(Year, SiteID) %>%
  mutate(survey = 1:n())

# spread columns by year
dat = spread(dat, key = Year, value = RANOKA)

# spread by surveys within years
dat = pivot_wider(
  dat, 
  id_cols = 'SiteID', 
  names_from = 'survey', 
  values_from = starts_with("20"), 
  names_glue = '{.value}.{survey}'
)


# Extract the detection histories
DH <- as.matrix(dat[, 2:358])
df = data.frame(SiteID = dat$SiteID, y = 0)
df[which(rowSums(DH, na.rm = T) > 0),"y"] = 1

df = cbind(df, covs)

#------------
# TSA
#------------
# get polynomials 
xy = as.matrix(cbind(df$Dec_long,df$Dec_lat))

poly3 = poly(xy, degree = 3, raw = T)
colnames(poly3) = c("X","X2","X3","Y1","XY","X2Y","Y2","XY2","Y3")

hab = cbind(df, poly3)
# TSA
poly1fit = glm(y ~ X + Y1, data = hab, family = binomial)
AICc(poly1fit)

# leave one out cross validation
ci.cvAUC(poly1fit$fitted.values, poly1fit$y, label.ordering = NULL, folds = NULL, confidence = 0.95)
# using residuals did not change predictions

#------------
# Correlogram
#------------
# popdists matrix
popdists <- as.matrix(rdist.earth(cbind(df$Dec_long, df$Dec_lat), miles = F, R = NULL))
diag(popdists) <- 0

# autocorr function
autocorr <- function(w,x,dist=0.9){
  aa <- ceiling(max(w)/dist)
  dists <- seq(0,aa*dist,dist)
  cors <- NULL
  for(i in 1:aa){
    w1 <- ifelse(w > dists[i] & w <= dists[i+1], 1, 0) 
    w2 <- w1
    for(j in 1:dim(w1)[1]){
      nu <- sum(w1[j,])
      if(nu>0){
        w2[j,] <- w1[j,]/nu
      }  
    }
    lag <- w2 %*% x
    cors <- c(cors,cor(x,lag))
  }
  return(cors)
}

ac1 <- autocorr(w=popdists,x=df$y,dist=2)
  
  
  it <- 1000
  mc <- matrix(NA,nrow=it,ncol=length(ac1))
  for(j in 1:it){
    df$rand <- sample(df$y,length(df$y),replace=F)
    mc[j,] <- autocorr(w=popdists,x=df$rand,dist=2)
  }
  
  aa = ceiling(max(popdists)/2)
  dists <- seq(2,aa*2,2)
  
  ac1 <- data.frame(cbind(ac1,dists))
  ac1 <- cbind(ac1,t(apply(mc,2,quantile, probs = c(0.025,0.975), na.rm = T)))
  names(ac1) <- c("ac","dist","lci","uci")
  
g1 = ggplot(ac1, aes(dist, ac)) + 
    geom_point(colour = "darkblue", size = 3) + 
    geom_line(colour = "red") +
    scale_x_continuous('Distance (km)',limits=c(0,30)) + 
    scale_y_continuous('Autocorrelation',limits=c(-1,1)) +
    theme_bw() + 
    geom_hline(yintercept=0) +   
    geom_smooth(aes(ymin = lci, ymax = uci), stat="identity",fill="blue",colour="darkblue") +
    theme(axis.title.x = element_text(face="bold", size=20),
        axis.text.x  = element_text(size=15),
        axis.title.y = element_text(face="bold", size=20),
        axis.text.y  = element_text(size=15))

ggsave(filename = "FigS3.pdf", plot = g1, device = cairo_pdf, width = 5, height = 4, units = "in", dpi =300)
#----------------------------
# VARIANCE PARTITIONING
#----------------------------

# extract latlons
xy = as.matrix(cbind(df$Dec_long, df$Dec_lat))

# create neighbourhood matrix <10km based on GPS locations
nbnear10 = dnearneigh(xy, 0, 7, row.names = df$SiteID, longlat = TRUE) 

par(mfrow = c(1,1))
plot(nbnear10, xy, col = "red", pch = 20, cex = 1)
title(main = "neighbours if 0<d<10km")

###We now use the multiscale co-dependence analysis "codep" package to perform the more general Moran Eigenvector Maps (MEM) analysis
mem = eigenmap(xy, weighting = wf.Drayf1, boundaries = c(0, 7))#performing MEM with the 1-dij/max(dij) weighting function and threshold = 4km
names(mem)#Of immediate interest are lambda and U
par(mfrow = c(1, 1), mar = c(4, 4, 4, 4))
barplot(mem$lambda, ylab = "Lambda")#Note differences in the number of positive eigenvectors compared to the PCNM

par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
plot(mem)#A scree plot of eigenvalues and spatial pattern (map) plots of spatial eigenvectors (i.e., the MEMs)



library(maptools)
library(ggmap)

windows()
# all US counties for plotting
usa <- maps::map("county", fill=TRUE)
gcounty <- map_data("county")
gcounty <- mutate(gcounty, polyname = paste(region, subregion, sep = ","))


# all US states for plotting
IDs <- usa$names
albers.proj <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
usa <- map2SpatialPolygons(usa, IDs=IDs, proj4string=albers.proj)
states <- maps::map("state", fill=TRUE, plot=FALSE)
states.df2 <- fortify(states, region='ID')

# just FL
eg.map <- gcounty[which(gcounty$region == "florida"),]
eg.map = cbind(eg.map, eggs = rep(0, length(eg.map[,1])))

eg.map[which(eg.map$subregion == "okaloosa"),8] = 1
eg.map[which(eg.map$subregion == "santa rosa"),8] = 1


fl.map <- ggplot() +
  geom_polygon(data=eg.map,
               aes(x=long+0.1, y=lat-0.1, group=group), fill='gray80') +
  geom_polygon(data=eg.map, 
               aes(x=long+0.05, y=lat-0.05, group=group, fill=factor(eggs)),
               color="black", size=0.25) +
  coord_equal() +
  scale_fill_manual(name="", values = c("white", 'lightblue'), guide = F) +
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-87.5, -85.9),
                     breaks = c(-87:-86), 
                     labels=c(paste(c(87:86), '°W', sep=''))) +
  scale_y_continuous(limits=c(30,31.1),
                     breaks=c(seq(30, 31, 1)),
                     labels=c(paste(seq(30, 31, 1), '°N', sep=''))) +
  labs(x='', y='')+  
  geom_point(data=df[,2:11], aes(x=Dec_long, y=Dec_lat, color=factor(y)))



##Now take your large positive eigenvectors and work with them as you :please
myEigenvecs = mem$U[, 1:3]


library(vegan)
vapa = varpart(Y = df$y, myEigenvecs, ~ sinu + upland + scrub + wetland + burn_interval, data = df)
par(mfrow = c(1, 1))
plot(vapa, cutoff = 0, cex = 1.5, bg = 2:5)

#-----------------------
# SAC test & Correlogram
#-----------------------
##Compute pairwise Euclidean distances with the function nbdist and calculate spatial weight matrix 
#The Euclidean distances are reported in units of decimal degrees 
dist_ecase = nbdists(nbnear10, xy)
#Now, define weight as a function of distance (e.g. 1-dij/max(dij))
#You could use other definitions here but this should work for starters
fdist = lapply(dist_ecase, function(x) 1 - x / max(dist(xy)))
#Now define your weight matrix, with 0 for non-neighbors and inverse-weighted distance for neighbors 
listw_ecase = nb2listw(nbnear10, glist = fdist, style = "W")
#Creating a neighborhood matrix from the nblistw file
mynbmatrix = listw2mat(listw_ecase)


#---------------------------------------------------------
# neighbourhood matrices and Moran's I for 1st eigenvector
#---------------------------------------------------------
n = length(myEigenvecs[, 1])
W = sum(mynbmatrix)
I_obs = (1 / W * t(myEigenvecs[, 1] - mean(myEigenvecs[, 1])) %*% mynbmatrix %*% (myEigenvecs[, 1] - mean(myEigenvecs[, 1])))/(1 / n * sum((myEigenvecs[, 1] - mean(myEigenvecs[, 1])) ^ 2))#Calculating observed value of Moran's I

###Randomization test for first eigenvector

nreps = 1000         # Note: 10,000 replications
I_ran = vector()

for (i in 1:nreps) {
  randScore <- sample(myEigenvecs[, 1], n, replace = FALSE)
  I_ran[i] = (1 / W * t(randScore - mean(randScore)) %*% mynbmatrix %*% (randScore - mean(randScore)))/(1 / n * sum((randScore - mean(randScore)) ^ 2))#Calculating observed value of Moran's I
}

##We can plot the distribution of Morans and test it for normality
hist(I_ran)
qqnorm(I_ran)
qqline(I_ran)
shapiro.test(I_ran)
###We know that the expected value of Moran's I is -1/(n-1) = -0.2 in this case, but we can't proceed with a Z-test due to lack of normality###
##Our other option is a permutation test based on the empirical CDF and the observed value of Moran's I
#Let's plot the hostogram and shade the tail of values > I_obs
par(mfrow = c(1, 1))
histMorans = hist(I_ran, breaks = 50, plot = FALSE)
plot(histMorans, col = ifelse(histMorans$breaks > I_obs[, 1], "red", "grey50"), main = "Randomization Test for 1st Eigenvector", xlab = "Moran's I")
text(0.03, 40, "Moran's I = 0.99", cex = 1.5)
#The Test
perc.rank = ecdf(I_ran)#This establishes percentile values for the empirical CDF
perc.rank(I_obs)
p = 1 - (perc.rank(I_obs))
p
###So with a p-value of 0, there is spatial autocorrelation in our data

#---------------------------------
# repeat with 2nd eigenvector
#----------------------------------
##Moran's I and Randomization test for 2nd positive eigenvector
I_obs = (1 / W * t(myEigenvecs[, 2] - mean(myEigenvecs[, 2])) %*% mynbmatrix %*% (myEigenvecs[, 2] - mean(myEigenvecs[, 2]))) / (1 / n * sum((myEigenvecs[, 2] - mean(myEigenvecs[, 2])) ^ 2))#Calculating observed value of Moran's I

nreps = 1000         # Note: 10,000 replications
I_ran = vector()

for (i in 1:nreps) {
  randScore <- sample(myEigenvecs[, 2], n, replace = FALSE)
  I_ran[i] = (1 / W * t(randScore - mean(randScore)) %*% mynbmatrix %*% (randScore - mean(randScore))) / (1 / n * sum((randScore - mean(randScore)) ^ 2))#Calculating observed value of Moran's I
}

##We can plot the distribution of Morans and test it for normality
hist(I_ran)
qqnorm(I_ran)
qqline(I_ran)
shapiro.test(I_ran)
###We know that the expected value of Moran's I is -1/(n-1) = -0.2 in this case, but we can't proceed with a Z-test due to lack of normality###
##Our other option is a permutation test based on the empirical CDF and the observed value of Moran's I
#Let's plot the hostogram and shade the tail of values > I_obs
histMorans = hist(I_ran, breaks = 50, plot = FALSE)
plot(histMorans, col = ifelse(histMorans$breaks > I_obs[,1], "red", "grey50"), main = "Randomization Test for 2th Eigenvector", xlab = "Moran's I")
text(0.02, 40, "Moran's I = 0.94", cex = 1.5)
#The Test
perc.rank = ecdf(I_ran)#This establishes percentile values for the empirical CDF
perc.rank(I_obs)
p = 1 - (perc.rank(I_obs))
p
###So with a p-value of 0, there is spatial autocorrelation in our 2nd eigenvector



#-------------------------
# BRT without eigenvectors
#-------------------------
sub.eig = df[, -c(3,8,9,11)]

sub.eig.boost = gbm.step(data = sub.eig, gbm.x = 3:7, gbm.y = 2,  family = "bernoulli", tree.complexity = 2, max.trees = 100000, learning.rate = 0.0001, bag.fraction = 0.85)
#Model has some good outcomes for both training and CV. You can tweak some of the tuning parameters to see how they affect the model performance
#Variable importance & plot
par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))
summary(sub.eig.boost, las = 1)

#Partial dependence plots
gbm.plot(sub.eig.boost, write.title = FALSE, plot.layout = c(2,3))

#-----------------------
# BRT with eigenvectors
#-----------------------
sub.eig = cbind(df[, -c(3,8,9,11)], myEigenvecs)

sub.eig.boost = gbm.step(data = sub.eig,gbm.x = 3:10, gbm.y = 2,  family = "bernoulli", 
                         tree.complexity = 3, learning.rate = 0.0001, bag.fraction = 0.85, max.trees = 1000000)
#Model has some good outcomes for both training and CV. You can tweak some of the tuning parameters to see how they affect the model performance
#Variable importance & plot
par(mfrow = c(1, 1), mar = c(5, 8, 4, 2))
summary(sub.eig.boost, las = 1)

#Partial dependence plots
par(mfrow = c(2,4), mar = c(5, 2, 2, 2))
gbm.plot(sub.eig.boost, write.title = FALSE, plot.layout = c(2, 4))




