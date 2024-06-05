setwd("C:/Users/Lab/OneDrive/Academia/Rothier et al/Locomotion in fluids/5.Functional Ecology/Review 1/analyses_review")

dat <- read.csv("data_2analyses.csv", sep=",", header=T, na.strings = "NA")
rownames(dat) <- dat$Species
dat <- dat[2:31]

library(ape); library(geiger)
library(evobiR); library(phytools)
library(ggplot2)
library(dispRity)

tree <- read.tree("mytree.nex")
tree <- rotate(tree, 846)
plot(tree, show.tip.label = FALSE); axisPhylo()

to_remove1 <- tree$tip.label[!tree$tip.label%in%rownames(dat)]
t1 <- drop.tip(tree, to_remove1)
name.check(t1, dat) 
dat1 <- ReorderData(t1, dat, taxa.names="row names")

dat1$Transf_Mass <- log10(dat1$Body_Mass^(1/3)) # transform bm in linear scale
str(dat1)




##==========================================================================================================##
##                                                                                                          ##
##   Subset 2                                                                                               ##
##                                                                                                          ##
##==========================================================================================================##




### ------ subset 2: dataset preparation ------ ###

dat2 <- na.omit(dat1[,c(1:28,30,31)])
str(dat2) # 837 species

name.check(tree, dat2)

to_remove2 <- tree$tip.label[!tree$tip.label%in%rownames(dat2)]
t2 <- drop.tip(tree, to_remove2)
name.check(t2, dat2) 

library(psych)
lgm2 <- log10(apply(dat2[,c(9:28,30)], 1, geometric.mean))
limb2 <- as.matrix(log10(dat2[,9:28]))
traits2 <- cbind(lgm2,limb2)
head(limb2)

media <- dat2$Media
names(media) <- dat2$Species  

subset2_dat <- list(tree=t2, log_limb=limb2, log_gm=lgm2, media=media)
save(subset2_dat, file="subset2_dat.RData")



##==========================================================================================================##
##                                                                                                          ##
####                         2.    Linear models                                                          ####
##                                                                                                          ##
##==========================================================================================================##

library(mvMORPH)

#### ------ Humerus ------####


    ## 1) With allometry:

#fit_BM_hum <- mvgls(limb2[,1:5]~1, tree=t2, model="BM", penalty="LASSO")
#fit_OU_hum <- mvgls(limb2[,1:5]~1, tree=t2, model="OU", penalty="LASSO")
#fit_EB_hum <- mvgls(limb2[,1:5]~1, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_hum); GIC(fit_OU_hum); GIC(fit_EB_hum) # OU


    ## 2) Allometry removed:

#fit_BM_hum2 <- mvgls(limb2[,1:5]~lgm2, tree=t2, model="BM", penalty="LASSO")
#fit_OU_hum2 <- mvgls(limb2[,1:5]~lgm2, tree=t2, model="OU", penalty="LASSO")
#fit_EB_hum2 <- mvgls(limb2[,1:5]~lgm2, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_hum2); GIC(fit_OU_hum2); GIC(fit_EB_hum2) # OU



##### ------ Radius ------ ####


    ## 1) With allometry

#fit_BM_rad <- mvgls(limb2[,6:10]~1, tree=t2, model="BM", penalty="LASSO")
#fit_OU_rad <- mvgls(limb2[,6:10]~1, tree=t2, model="OU", penalty="LASSO")
#fit_EB_rad <- mvgls(limb2[,6:10]~1, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_rad); GIC(fit_OU_rad); GIC(fit_EB_rad) # OU



    ## 2) Allometry removed:

#fit_BM_rad2 <- mvgls(limb2[,6:10]~lgm2, tree=t2, model="BM", penalty="LASSO")
#fit_OU_rad2 <- mvgls(limb2[,6:10]~lgm2, tree=t2, model="OU", penalty="LASSO")
#fit_EB_rad2 <- mvgls(limb2[,6:10]~lgm2, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_rad2); GIC(fit_OU_rad2); GIC(fit_EB_rad2) # OU




#### ------ Metacarpal ------ ####


    ## 1) With allometry

#fit_BM_met <- mvgls(limb2[,11:15]~1, tree=t2, model="BM", penalty="LASSO")
#fit_OU_met <- mvgls(limb2[,11:15]~1, tree=t2, model="OU", penalty="LASSO")
#fit_EB_met <- mvgls(limb2[,11:15]~1, tree=t2, model="EB", penalty="LASSO")

    #GIC(fit_BM_met); GIC(fit_OU_met); GIC(fit_EB_met) # OU



    ## 2) Allometry removed:

#fit_BM_met2 <- mvgls(limb2[,11:15]~lgm2, tree=t2, model="BM", penalty="LASSO")
#fit_OU_met2 <- mvgls(limb2[,11:15]~lgm2, tree=t2, model="OU", penalty="LASSO")
#fit_EB_met2 <- mvgls(limb2[,11:15]~lgm2, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_met2); GIC(fit_OU_met2); GIC(fit_EB_met2) # OU




#### ------ Phalanx ------ ####


    ## 1) With allometry

#fit_BM_phal <- mvgls(limb2[,16:20]~1, tree=t2, model="BM", penalty="LASSO")
#fit_OU_phal <- mvgls(limb2[,16:20]~1, tree=t2, model="OU", penalty="LASSO")
#fit_EB_phal <- mvgls(limb2[,16:20]~1, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_phal); GIC(fit_OU_phal); GIC(fit_EB_phal) # OU



    ## 2) Allometry removed:


#fit_BM_phal2 <- mvgls(limb2[,16:20]~lgm2, tree=t2, model="BM", penalty="LASSO")
#fit_OU_phal2 <- mvgls(limb2[,16:20]~lgm2, tree=t2, model="OU", penalty="LASSO")
#fit_EB_phal2 <- mvgls(limb2[,16:20]~lgm2, tree=t2, model="EB", penalty="LASSO")

#GIC(fit_BM_phal2); GIC(fit_OU_phal2); GIC(fit_EB_phal2) # OU


#bones_models <- list(fit_BM_hum=fit_BM_hum, fit_OU_hum=fit_OU_hum, fit_EB_hum=fit_EB_hum, 
#                    fit_BM_rad=fit_BM_rad,fit_OU_rad=fit_OU_rad, fit_EB_rad=fit_EB_rad, 
 #                   fit_BM_met=fit_BM_met, fit_OU_met=fit_OU_met, fit_EB_met=fit_EB_met, 
  #                  fit_BM_phal=fit_BM_phal, fit_OU_phal=fit_OU_phal, fit_EB_phal=fit_EB_phal, 
   #                 fit_BM_hum2=fit_BM_hum2,fit_OU_hum2=fit_OU_hum2,fit_EB_hum2=fit_EB_hum2,
    #                fit_BM_rad2=fit_BM_rad2,fit_OU_rad2=fit_OU_rad2,fit_EB_rad2=fit_EB_rad2,
     #               fit_BM_met2=fit_BM_met2,fit_OU_met2=fit_OU_met2,fit_EB_met2=fit_EB_met2,
      #              fit_BM_phal2=fit_BM_phal2,fit_OU_phal2=fit_OU_phal2,fit_EB_phal2=fit_EB_phal2)


#save(bones_models, file="bones_models.RData")
load("bones_models.RData")

lapply(bones_models, function(x) sum(diag(x$sigma$Pinv)))


propvar <- function(x) {
  x/sum(x)} # function to calculate variance explained by each pc


media <- dat2$Media
names(media) <- row.names(dat2)





##==========================================================================================================##
##                                                                                                          ##
####                                    3. manova gls                                                     ####
##                                                                                                          ##
##==========================================================================================================##


data.hum = list(hum=limb2[,1:5], tree=t2, lgm=lgm2, media=media)
data.rad = list(rad = limb2[,6:10], tree=t2, lgm=lgm2, media=media)
data.met = list(met = limb2[,11:15], tree=t2, lgm=lgm2, media=media)
data.phal = list( phal = limb2[,16:20], tree=t2, lgm=lgm2, media=media)

pgls.hum <-  mvgls(hum~lgm2*media, tree=t2, data=data.hum, model="OU", method="LL")
pgls.rad <-  mvgls(rad~lgm2*media, tree=t2, data=data.rad, model="OU", method="LL")
pgls.met <-  mvgls(met~lgm2*media, tree=t2, data=data.met, model="OU", method="LL")
pgls.phal <-  mvgls(phal~lgm2*media, tree=t2, data=data.phal, model="OU", method="LL")

pairwise.glh(pgls.hum, test="Pillai")
(multivariate_test <- manova.gls(pgls.hum, nperm=999, type="II", test="Pillai"))

pairwise.glh(pgls.rad, test="Pillai")
(multivariate_test <- manova.gls(pgls.rad, nperm=999, type="II", test="Pillai"))

pairwise.glh(pgls.met, test="Pillai")
(multivariate_met <- manova.gls(pgls.met, nperm=999, type="II", test="Pillai"))

pairwise.glh(pgls.phal, test="Pillai")
(multivariate_test <- manova.gls(pgls.phal, nperm=999, type="II", test="Pillai"))







##==========================================================================================================##
##                                                                                                          ##
####                                    4.   Analyses HUMERUS                                             ####
##                                                                                                          ##
##==========================================================================================================##



#### ------ SHAPE without allometry (traits ~ gm) ------ ####

fit_OU_hum2 <- bones_models$fit_OU_hum2
pca_hum_gm <- mvgls.pca(fit_OU_hum2)
row.names(pca_hum_gm$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_hum_gm <- as.data.frame(pca_hum_gm$scores)
sc_pca_hum_gm$Group <-  dat2$Group
sc_pca_hum_gm$Media <- dat2$Media
sc_pca_hum_gm$MonophyleticGroup <- dat2$major_clade

data.frame(pca_hum_gm$values)->eigen.df_hum2


propvar_hum2 <- apply(eigen.df_hum2, 2, propvar)
cumvar_hum2 <- cumsum(propvar_hum2)


#### ------ regular PCA ------ ####

opca_hum_shape <- prcomp(fit_OU_hum2$residuals)
summary(opca_hum_shape)

opca_hum_shape$sdev^2 


#### ------ RAW with allometry (traits ~ 1) ------ ####



## 1) phylo PCA

fit_OU_hum <- bones_models$fit_OU_hum
pca_hum_raw <- mvgls.pca(fit_OU_hum)
row.names(pca_hum_raw$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_hum_raw <- as.data.frame(pca_hum_raw$scores)
sc_pca_hum_raw$Group <-  dat2$Group
sc_pca_hum_raw$Media <- dat2$Media
sc_pca_hum_raw$MonophyleticGroup <- dat2$major_clade

data.frame(pca_hum_raw$values)->eigen.df_hum



propvar_hum <- apply(eigen.df_hum, 2, propvar)
cumvar_hum <- cumsum(propvar_hum)



#### ------ regular PCA ------ ####

opca_hum_size <- prcomp(fit_OU_hum$residuals)
summary(opca_hum_size)

opca_hum_size$sdev^2 



##==========================================================================================================##
##                                                                                                          ##
####                                    5.   Analyses RADIUS                                              ####
##                                                                                                          ##
##==========================================================================================================##




#### ------ SHAPE without allometry (traits ~ gm) ------ ####

fit_OU_rad2 <- bones_models$fit_OU_rad2
pca_rad_gm <- mvgls.pca(fit_OU_rad2)
row.names(pca_rad_gm$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_rad_gm <- as.data.frame(pca_rad_gm$scores)
sc_pca_rad_gm$Group <-  dat2$Group
sc_pca_rad_gm$Media <- dat2$Media
sc_pca_rad_gm$MonophyleticGroup <- dat2$major_clade

data.frame(pca_rad_gm$values)->eigen.df_rad2


propvar_rad2 <- apply(eigen.df_rad2, 2, propvar)
cumvar_rad2 <- cumsum(propvar_rad2)



#### ------ regular PCA ------ ####

opca_rad_shape <- prcomp(fit_OU_rad2$residuals)
summary(opca_rad_shape)

opca_rad_shape$sdev^2 



# ===============================================================


#### ------ with allometry (traits ~ 1) ------ ####

fit_OU_rad <- bones_models$fit_OU_rad

pca_rad_raw <- mvgls.pca(fit_OU_rad)
row.names(pca_rad_raw$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_rad_raw <- as.data.frame(pca_rad_raw$scores)
sc_pca_rad_raw$Group <-  dat2$Group
sc_pca_rad_raw$Media <- dat2$Media
sc_pca_rad_raw$MonophyleticGroup <- dat2$major_clade

data.frame(pca_rad_raw$values)->eigen.df_rad



propvar_rad <- apply(eigen.df_rad, 2, propvar)
cumvar_rad <- cumsum(propvar_rad)



#### ------ regular PCA ------ ####

opca_rad_size <- prcomp(fit_OU_rad$residuals)
summary(opca_rad_size)

opca_rad_size$sdev^2 




##==========================================================================================================##
##                                                                                                          ##
####                                    6.   Analyses METACARPAL                                          ####
##                                                                                                          ##
##==========================================================================================================##




#### ------ SHAPE without allometry (traits ~ gm) ------ ####

fit_OU_met2 <- bones_models$fit_OU_met2
pca_met_gm <- mvgls.pca(fit_OU_met2)
row.names(pca_met_gm$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_met_gm <- as.data.frame(pca_met_gm$scores)
sc_pca_met_gm$Group <-  dat2$Group
sc_pca_met_gm$Media <- dat2$Media
sc_pca_met_gm$MonophyleticGroup <- dat2$major_clade

data.frame(pca_met_gm$values)->eigen.df_met2


propvar_met2 <- apply(eigen.df_met2, 2, propvar)
cumvar_met2 <- cumsum(propvar_met2)


#### ------ regular PCA ------ ####

opca_met_shape <- prcomp(fit_OU_met2$residuals)
summary(opca_met_shape)

opca_met_shape$sdev^2 



# ===============================================================


#### ------ with allometry (traits ~ 1) ------ ####


## 1) phylo PCA

fit_OU_met <- bones_models$fit_OU_met
pca_met_raw <- mvgls.pca(fit_OU_met)
row.names(pca_met_raw$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_met_raw <- as.data.frame(pca_met_raw$scores)
sc_pca_met_raw$Group <-  dat2$Group
sc_pca_met_raw$Media <- dat2$Media
sc_pca_met_raw$MonophyleticGroup <- dat2$major_clade

data.frame(pca_met_raw$values)->eigen.df_met



propvar_met <- apply(eigen.df_met, 2, propvar)
cumvar_met <- cumsum(propvar_met)


#### ------ regular PCA ------ ####

opca_met_size <- prcomp(fit_OU_met$residuals)
summary(opca_met_size)

opca_met_size$sdev^2 


##==========================================================================================================##
##                                                                                                          ##
####                                    7.   Analyses PHALANX                                             ####
##                                                                                                          ##
##==========================================================================================================##




#### ------ without allometry (traits ~ gm) ------ ####

fit_OU_phal2 <- bones_models$fit_OU_phal2
pca_phal_gm <- mvgls.pca(fit_OU_phal2)
row.names(pca_phal_gm$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_phal_gm <- as.data.frame(pca_phal_gm$scores)
sc_pca_phal_gm$Group <-  dat2$Group
sc_pca_phal_gm$Media <- dat2$Media
sc_pca_phal_gm$MonophyleticGroup <- dat2$major_clade

data.frame(pca_phal_gm$values)->eigen.df_phal2


propvar_phal2 <- apply(eigen.df_phal2, 2, propvar)
cumvar_phal2 <- cumsum(propvar_phal2)



#### ------ regular PCA ------ ####

opca_phal_shape <- prcomp(fit_OU_phal2$residuals)
summary(opca_phal_shape)

opca_phal_shape$sdev^2 


# ===============================================================

#### ------ with allometry (traits ~ 1) ------ ####




## 1) phylo PCA

fit_OU_phal <- bones_models$fit_OU_phal
pca_phal_raw <- mvgls.pca(fit_OU_phal)
row.names(pca_phal_raw$vectors) <- c("Length", "PWidth", "MWidth", "DWidth", "Height")

sc_pca_phal_raw <- as.data.frame(pca_phal_raw$scores)
sc_pca_phal_raw$Group <-  dat2$Group
sc_pca_phal_raw$Media <- dat2$Media
sc_pca_phal_raw$MonophyleticGroup <- dat2$major_clade

data.frame(pca_phal_raw$values)->eigen.df_phal



propvar_phal <- apply(eigen.df_phal, 2, propvar)
cumvar_phal <- cumsum(propvar_phal)


#### ------ regular PCA ------ ####

opca_phal_size <- prcomp(fit_OU_phal$residuals)
summary(opca_phal_size)

opca_phal_size$sdev^2 



#---------------------------------------------------------------------


bones_pca <- list(sc_pca_hum_raw=sc_pca_hum_raw,
                  sc_pca_hum_gm=sc_pca_hum_gm, 
                  sc_pca_rad_raw=sc_pca_rad_raw,
                  sc_pca_rad_gm=sc_pca_rad_gm, 
                  sc_pca_met_raw=sc_pca_met_raw,
                  sc_pca_met_gm=sc_pca_met_gm, 
                  sc_pca_phal_raw=sc_pca_phal_raw,
                  sc_pca_phal_gm=sc_pca_phal_gm)

save(bones_pca, file="bones_pca.Rdata")

bones_opca <- list(sc_opca_hum_size=opca_hum_size$x,
                  sc_opca_hum_shape=opca_hum_shape$x, 
                  sc_opca_rad_size=opca_rad_size$x,
                  sc_opca_rad_shape=opca_rad_shape$x, 
                  sc_opca_met_size=opca_met_size$x,
                  sc_opca_met_shape=opca_met_shape$x, 
                  sc_opca_phal_size=opca_phal_size$x,
                  sc_opca_phal_shape=opca_phal_shape$x)

save(bones_opca, file="bones_opca.Rdata")

##==========================================================================================================##
##                                                                                                          ##
####                                    8. Disparity                                                      ####
##                                                                                                          ##
##==========================================================================================================##



##=========================================================================##
#### Humerus #####
##=========================================================================##

m2 <- as.character(dat2$Media)
names(m2) <- rownames(dat2)


mg <- as.character(dat2$major_clade)
names(mg) <- rownames(dat2)
unique(mg)

#### ------ Environment ------ ####

##  SHAPE


ord_matrix_hum_shape <- fit_OU_hum2$residuals 

env_subsets_hum_shape <- custom.subsets(ord_matrix_hum_shape, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


#env_bootstrapped_hum_shape <- boot.matrix(env_subsets_hum_shape, bootstraps = 1000,
#                                rarefaction = TRUE, verbose=TRUE) 

#env_disparity_hum_shape <- dispRity(env_bootstrapped_hum_shape, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_hum_shape, file="env_disparity_hum_shape.RData")
load("env_disparity_hum_shape.RData")

summary(env_disparity_hum_shape)
plot.dispRity(env_disparity_hum_shape, rarefaction=TRUE)
plot(env_disparity_hum_shape, xlab="Media of locomotion", ylab="Disparity - Humerus Shape", cex.axis=0.8)  

disparity_hum_shape <- env_disparity_hum_shape$disparity

means_hum_shape <- lapply(disparity_hum_shape, function(x) mean(x[[2]]))
sd_hum_shape <- lapply(disparity_hum_shape, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_hum_shape, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_hum_shape,bhatt.coeff, "pairwise")





##  RAW

ord_matrix_hum_size <- fit_OU_hum$residuals 

env_subsets_hum_size <- custom.subsets(ord_matrix_hum_size, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_hum_size <- boot.matrix(env_subsets_hum_size, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_hum_size <- dispRity(env_bootstrapped_hum_size, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_hum_size, file="env_disparity_hum_size.RData")
load("env_disparity_hum_size.RData")

summary(env_disparity_hum_size)
plot.dispRity(env_disparity_hum_size, rarefaction=TRUE)
plot(env_disparity_hum_size, xlab="Media of locomotion", ylab="Disparity - Humerus size", cex.axis=0.8)  

disparity_hum_size <- env_disparity_hum_size$disparity

means_hum_size <- lapply(disparity_hum_size, function(x) mean(x[[2]]))
sd_hum_sizee <- lapply(disparity_hum_size, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_hum_size, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_hum_size,bhatt.coeff, "pairwise")




#### ------ Group ------ ####


## Shape


## Creating the table that contain the elements and their attributes
group_subsets_hum_shape <- custom.subsets(ord_matrix_hum_shape, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_hum_shape <- boot.matrix(group_subsets_hum_shape, bootstraps = 1000,
                                  rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

#group_disparity_hum_shape <- dispRity(group_bootstrapped_hum_shape, metric = c(sum, variances), verbose=TRUE)
#summary(group_disparity_hum_shape)
#save(group_disparity_hum_shape, file="group_disparity_hum_shape.RData")
load("group_disparity_hum_shape.RData")

plot.dispRity(group_disparity_hum_shape, rarefaction=TRUE)
plot(group_disparity_hum_shape, xlab="Group", ylab="Disparity - Humerus Shape", cex.axis=0.8)  

gr_disparity_hum_shape <- group_disparity_hum_shape$disparity
gr_means_hum_shape <- lapply(gr_disparity_hum_shape, function(x) mean(x[[2]]))
gr_sd_hum_shape <- lapply(gr_disparity_hum_shape, function(x) sd(x[[2]]))


test.dispRity(group_disparity_hum_shape, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_hum_shape,bhatt.coeff, "pairwise")



## Raw


## Creating the table that contain the elements and their attributes
group_subsets_hum_size <- custom.subsets(ord_matrix_hum_size, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_hum_size <- boot.matrix(group_subsets_hum_size, bootstraps = 1000,
                                            rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_hum_size <- dispRity(group_bootstrapped_hum_size, metric = c(sum, variances), verbose=TRUE)
save(group_disparity_hum_size, file="group_disparity_hum_size.RData")
#load("group_disparity_hum_size.RData")

plot.dispRity(group_disparity_hum_size, rarefaction=TRUE)
plot(group_disparity_hum_size, xlab="Group", ylab="Disparity - Humerus size", cex.axis=0.8)  

gr_disparity_hum_size <- group_disparity_hum_size$disparity
gr_means_hum_size <- lapply(gr_disparity_hum_size, function(x) mean(x[[2]]))
gr_sd_hum_size <- lapply(gr_disparity_hum_size, function(x) sd(x[[2]]))


test.dispRity(group_disparity_hum_size, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_hum_size,bhatt.coeff, "pairwise")




##=========================================================================##
#### Radius #####
##=========================================================================##

#### ------ Environment ------ ####

##  SHAPE

ord_matrix_rad_shape <- fit_OU_rad2$residuals 

env_subsets_rad_shape <- custom.subsets(ord_matrix_rad_shape, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_rad_shape <- boot.matrix(env_subsets_rad_shape, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_rad_shape <- dispRity(env_bootstrapped_rad_shape, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_rad_shape, file="env_disparity_rad_shape.RData")
#load("env_disparity_rad_shape.RData")

plot.dispRity(env_disparity_rad_shape, rarefaction=TRUE)
plot(env_disparity_rad_shape, xlab="Media of locomotion", ylab="Disparity - radius Shape", cex.axis=0.8)  

disparity_rad_shape <- env_disparity_rad_shape$disparity

means_rad_shape <- lapply(disparity_rad_shape, function(x) mean(x[[2]]))
sd_rad_shape <- lapply(disparity_rad_shape, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_rad_shape, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_rad_shape,bhatt.coeff, "pairwise")





##  RAW

ord_matrix_rad_size <- fit_OU_rad$residuals 

env_subsets_rad_size <- custom.subsets(ord_matrix_rad_size, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_rad_size <- boot.matrix(env_subsets_rad_size, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_rad_size <- dispRity(env_bootstrapped_rad_size, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_rad_size, file="env_disparity_rad_size.RData")
#load("env_disparity_rad_size.RData")

plot.dispRity(env_disparity_rad_size, rarefaction=TRUE)
plot(env_disparity_rad_size, xlab="Media of locomotion", ylab="Disparity - radius size", cex.axis=0.8)  

disparity_rad_size <- env_disparity_rad_size$disparity

means_rad_size <- lapply(disparity_rad_size, function(x) mean(x[[2]]))
sd_rad_sizee <- lapply(disparity_rad_size, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_rad_size, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_rad_size,bhatt.coeff, "pairwise")




#### ------ Group ------ ####


## Shape


## Creating the table that contain the elements and their attributes
group_subsets_rad_shape <- custom.subsets(ord_matrix_rad_shape, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_rad_shape <- boot.matrix(group_subsets_rad_shape, bootstraps = 1000,
                                            rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_rad_shape <- dispRity(group_bootstrapped_rad_shape, metric = c(sum, variances), verbose=TRUE)
save(group_disparity_rad_shape, file="group_disparity_rad_shape.RData")
#load("group_disparity_rad_shape.RData")

plot.dispRity(group_disparity_rad_shape, rarefaction=TRUE)
plot(group_disparity_rad_shape, xlab="Group", ylab="Disparity - radius Shape", cex.axis=0.8)  

gr_disparity_rad_shape <- group_disparity_rad_shape$disparity
gr_means_rad_shape <- lapply(gr_disparity_rad_shape, function(x) mean(x[[2]]))
gr_sd_rad_shape <- lapply(gr_disparity_rad_shape, function(x) sd(x[[2]]))


test.dispRity(group_disparity_rad_shape, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_rad_shape,bhatt.coeff, "pairwise")



## Raw


## Creating the table that contain the elements and their attributes
group_subsets_rad_size <- custom.subsets(ord_matrix_rad_size, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_rad_size <- boot.matrix(group_subsets_rad_size, bootstraps = 1000,
                                           rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_rad_size <- dispRity(group_bootstrapped_rad_size, metric = c(sum, variances), verbose=TRUE)
save(group_disparity_rad_size, file="group_disparity_rad_size.RData")
#load("group_disparity_rad_size.RData")

plot.dispRity(group_disparity_rad_size, rarefaction=TRUE)
plot(group_disparity_rad_size, xlab="Group", ylab="Disparity - radius size", cex.axis=0.8)  

gr_disparity_rad_size <- group_disparity_rad_size$disparity
gr_means_rad_size <- lapply(gr_disparity_rad_size, function(x) mean(x[[2]]))
gr_sd_rad_size <- lapply(gr_disparity_rad_size, function(x) sd(x[[2]]))


test.dispRity(group_disparity_rad_size, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_rad_size,bhatt.coeff, "pairwise")




##=========================================================================##
# Metacarpal #
##=========================================================================##

#### ------ Environment ------ ####

##  SHAPE

ord_matrix_met_shape <- fit_OU_met2$residuals 

env_subsets_met_shape <- custom.subsets(ord_matrix_met_shape, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_met_shape <- boot.matrix(env_subsets_met_shape, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_met_shape <- dispRity(env_bootstrapped_met_shape, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_met_shape, file="env_disparity_met_shape.RData")
#load("env_disparity_met_shape.RData")

plot.dispRity(env_disparity_met_shape, rarefaction=TRUE)
plot(env_disparity_met_shape, xlab="Media of locomotion", ylab="Disparity - metacarpal Shape", cex.axis=0.8)  

disparity_met_shape <- env_disparity_met_shape$disparity

means_met_shape <- lapply(disparity_met_shape, function(x) mean(x[[2]]))
sd_met_shape <- lapply(disparity_met_shape, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_met_shape, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_met_shape,bhatt.coeff, "pairwise")





##  RAW

ord_matrix_met_size <- fit_OU_met$residuals 

env_subsets_met_size <- custom.subsets(ord_matrix_met_size, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_met_size <- boot.matrix(env_subsets_met_size, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_met_size <- dispRity(env_bootstrapped_met_size, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_met_size, file="env_disparity_met_size.RData")
#load("env_disparity_met_size.RData")

plot.dispRity(env_disparity_met_size, rarefaction=TRUE)
plot(env_disparity_met_size, xlab="Media of locomotion", ylab="Disparity - metacarpal size", cex.axis=0.8)  

disparity_met_size <- env_disparity_met_size$disparity

means_met_size <- lapply(disparity_met_size, function(x) mean(x[[2]]))
sd_met_sizee <- lapply(disparity_met_size, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_met_size, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_met_size,bhatt.coeff, "pairwise")




#### ------ Group ------ ####


## Shape


## Creating the table that contain the elements and their attributes
group_subsets_met_shape <- custom.subsets(ord_matrix_met_shape, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_met_shape <- boot.matrix(group_subsets_met_shape, bootstraps = 1000,
                                            rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_met_shape <- dispRity(group_bootstrapped_met_shape, metric = c(sum, variances), verbose=TRUE)
save(group_disparity_met_shape, file="group_disparity_met_shape.RData")
#load("group_disparity_met_shape.RData")

plot.dispRity(group_disparity_met_shape, rarefaction=TRUE)
plot(group_disparity_met_shape, xlab="Group", ylab="Disparity - metacarpal Shape", cex.axis=0.8)  

gr_disparity_met_shape <- group_disparity_met_shape$disparity
gr_means_met_shape <- lapply(gr_disparity_met_shape, function(x) mean(x[[2]]))
gr_sd_met_shape <- lapply(gr_disparity_met_shape, function(x) sd(x[[2]]))


test.dispRity(group_disparity_met_shape, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_met_shape,bhatt.coeff, "pairwise")



## Raw


## Creating the table that contain the elements and their attributes
group_subsets_met_size <- custom.subsets(ord_matrix_met_size, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_met_size <- boot.matrix(group_subsets_met_size, bootstraps = 1000,
                                           rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_met_size <- dispRity(group_bootstrapped_met_size, metric = c(sum, variances), verbose=TRUE)
save(group_disparity_met_size, file="group_disparity_met_size.RData")
#load("group_disparity_met_size.RData")


plot.dispRity(group_disparity_met_size, rarefaction=TRUE)
plot(group_disparity_met_size, xlab="Group", ylab="Disparity - metacarpal size", cex.axis=0.8)  

gr_disparity_met_size <- group_disparity_met_size$disparity
gr_means_met_size <- lapply(gr_disparity_met_size, function(x) mean(x[[2]]))
gr_sd_met_size <- lapply(gr_disparity_met_size, function(x) sd(x[[2]]))


test.dispRity(group_disparity_met_size, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_met_size,bhatt.coeff, "pairwise")





##=========================================================================##
# Phalanx #
##=========================================================================##

#### ------ Environment ------ ####

##  SHAPE

ord_matrix_phal_shape <- fit_OU_phal2$residuals 

env_subsets_phal_shape <- custom.subsets(ord_matrix_phal_shape, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_phal_shape <- boot.matrix(env_subsets_phal_shape, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_phal_shape <- dispRity(env_bootstrapped_phal_shape, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_phal_shape, file="env_disparity_phal_shape.RData")
#load("env_disparity_phal_shape.RData")

plot.dispRity(env_disparity_phal_shape, rarefaction=TRUE)
plot(env_disparity_phal_shape, xlab="Media of locomotion", ylab="Disparity - phalanx Shape", cex.axis=0.8)  

disparity_phal_shape <- env_disparity_phal_shape$disparity

means_phal_shape <- lapply(disparity_phal_shape, function(x) mean(x[[2]]))
sd_phal_shape <- lapply(disparity_phal_shape, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_phal_shape, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_phal_shape,bhatt.coeff, "pairwise")





##  RAW

ord_matrix_phal_size <- fit_OU_phal$residuals 

env_subsets_phal_size <- custom.subsets(ord_matrix_phal_size, tree=t2, group = list(
  "terrestrial" = which(m2 == "terrestrial"),
  "aquatic" = which(m2 == "aquatic"),
  "semi-aquatic" = which(m2 == "semi-aquatic"),
  "air" = which(m2 == "air"),
  "semi-aerial" = which(m2 == "semi-aerial")))


env_bootstrapped_phal_size <- boot.matrix(env_subsets_phal_size, bootstraps = 1000,
                                rarefaction = TRUE, verbose=TRUE) 

env_disparity_phal_size <- dispRity(env_bootstrapped_phal_size, metric = c(sum, variances), verbose=TRUE)
save(env_disparity_phal_size, file="env_disparity_phal_size.RData")
#load("env_disparity_phal_size.RData")


plot.dispRity(env_disparity_phal_size, rarefaction=TRUE)
plot(env_disparity_phal_size, xlab="Media of locomotion", ylab="Disparity - phalanx size", cex.axis=0.8)  

disparity_phal_size <- env_disparity_phal_size$disparity

means_phal_size <- lapply(disparity_phal_size, function(x) mean(x[[2]]))
sd_phal_sizee <- lapply(disparity_phal_size, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_phal_size, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_phal_size,bhatt.coeff, "pairwise")




#### ------ Group ------ ####


## Shape


## Creating the table that contain the elements and their attributes
group_subsets_phal_shape <- custom.subsets(ord_matrix_phal_shape, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_phal_shape <- boot.matrix(group_subsets_phal_shape, bootstraps = 1000,
                                            rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_phal_shape <- dispRity(group_bootstrapped_phal_shape, metric = c(sum, variances), verbose=TRUE)

save(group_disparity_phal_shape, file="group_disparity_phal_shape.RData")
#load("group_disparity_phal_shape.RData")

plot.dispRity(group_disparity_phal_shape, rarefaction=TRUE)
plot(group_disparity_phal_shape, xlab="Group", ylab="Disparity - phalanx Shape", cex.axis=0.8)  

gr_disparity_phal_shape <- group_disparity_phal_shape$disparity
gr_means_phal_shape <- lapply(gr_disparity_phal_shape, function(x) mean(x[[2]]))
gr_sd_phal_shape <- lapply(gr_disparity_phal_shape, function(x) sd(x[[2]]))


test.dispRity(group_disparity_phal_shape, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_phal_shape,bhatt.coeff, "pairwise")



## Raw


## Creating the table that contain the elements and their attributes
group_subsets_phal_size <- custom.subsets(ord_matrix_phal_size, tree=t2, group = list(
  "others" = which(mg == "others"),
  "Cetacea" = which(mg == "Cetacea"),
  "Pinnipedia" = which(mg == "Pinnipedia"),
  "Chiroptera" = which(mg == "Chiroptera")))


## Bootstrapping the data
group_bootstrapped_phal_size <- boot.matrix(group_subsets_phal_size, bootstraps = 1000,
                                           rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

group_disparity_phal_size <- dispRity(group_bootstrapped_phal_size, metric = c(sum, variances), verbose=TRUE)
save(group_disparity_phal_size, file="group_disparity_phal_size.RData")
#load("group_disparity_phal_size.RData")


plot.dispRity(group_disparity_phal_size, rarefaction=TRUE)
plot(group_disparity_phal_size, xlab="Group", ylab="Disparity - phalanx size", cex.axis=0.8)  

gr_disparity_phal_size <- group_disparity_phal_size$disparity
gr_means_phal_size <- lapply(gr_disparity_phal_size, function(x) mean(x[[2]]))
gr_sd_phal_size <- lapply(gr_disparity_phal_size, function(x) sd(x[[2]]))


test.dispRity(group_disparity_phal_size, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_phal_size,bhatt.coeff, "pairwise")







##==========================================================================================================##
##                                                                                                          ##
####                                   9. Phyl ANOVA                                                      ####
##                                                                                                          ##
##==========================================================================================================##


bones_pca <- list(sc_pca_hum_raw=sc_pca_hum_raw,
                  sc_pca_hum_gm=sc_pca_hum_gm, 
                  sc_pca_rad_raw=sc_pca_rad_raw,
                  sc_pca_rad_gm=sc_pca_rad_gm, 
                  sc_pca_met_raw=sc_pca_met_raw,
                  sc_pca_met_gm=sc_pca_met_gm, 
                  sc_pca_phal_raw=sc_pca_phal_raw,
                  sc_pca_phal_gm=sc_pca_phal_gm)


phylaov_hum.shape <- list()
for(i in 1:4){
  phylaov_hum.shape[[i]] <-  phylANOVA(t2, media, sc_pca_hum_gm[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}
phylaov_hum.size <- list()
for(i in 1:4){
  phylaov_hum.size[[i]] <-  phylANOVA(t2, media, sc_pca_hum_raw[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}


phylaov_rad.shape <- list()
for(i in 1:4){
  phylaov_rad.shape[[i]] <-  phylANOVA(t2, media, sc_pca_rad_gm[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}
phylaov_rad.size <- list()
for(i in 1:4){
  phylaov_rad.size[[i]] <-  phylANOVA(t2, media, sc_pca_rad_raw[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}


phylaov_met.shape <- list()
for(i in 1:4){
  phylaov_met.shape[[i]] <-  phylANOVA(t2, media, sc_pca_met_gm[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}
phylaov_met.size <- list()
for(i in 1:4){
  phylaov_met.size[[i]] <-  phylANOVA(t2, media, sc_pca_met_raw[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}


phylaov_phal.shape <- list()
for(i in 1:4){
  phylaov_phal.shape[[i]] <-  phylANOVA(t2, media, sc_pca_phal_gm[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}
phylaov_phal.size <- list()
for(i in 1:4){
  phylaov_phal.size[[i]] <-  phylANOVA(t2, media, sc_pca_phal_raw[,i], posthoc = TRUE, p.adj="holm", nsim=10000)
}



##==========================================================================================================##
##                                                                                                          ##
####                                   10. saving data                                                    ####
##                                                                                                          ##
##==========================================================================================================##


## plot Humerus ------------

df_s_hum_shape<- (disparity_hum_shape$terrestrial[[2]])
df_aq_hum_shape <- (disparity_hum_shape$aquatic[[2]])
df_saq_hum_shape <- (disparity_hum_shape$`semi-aquatic`[[2]])
df_a_hum_shape <- (disparity_hum_shape$air[[2]])
df_sa_hum_shape <- (disparity_hum_shape$`semi-aerial`[[2]])

df_hum_shape_env <- t(rbind(df_s_hum_shape,df_aq_hum_shape,df_saq_hum_shape, df_a_hum_shape, df_sa_hum_shape))
colnames(df_hum_shape_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )

 
df_s_hum_size<- (disparity_hum_size$terrestrial[[2]])
df_aq_hum_size <- (disparity_hum_size$aquatic[[2]])
df_saq_hum_size <- (disparity_hum_size$`semi-aquatic`[[2]])
df_a_hum_size <- (disparity_hum_size$air[[2]])
df_sa_hum_size <- (disparity_hum_size$`semi-aerial`[[2]])

df_hum_size_env <- t(rbind(df_s_hum_size,df_aq_hum_size,df_saq_hum_size, df_a_hum_size, df_sa_hum_size))
colnames(df_hum_size_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )



#--

df_oth_hum_shape <- (gr_disparity_hum_shape$others[[2]])
df_cet_hum_shape <- (gr_disparity_hum_shape$Cetacea[[2]])
df_pinn_hum_shape <- (gr_disparity_hum_shape$Pinnipedia[[2]])
df_chir_hum_shape <- (gr_disparity_hum_shape$Chiroptera[[2]])

df_hum_shape_group <- t(rbind(df_oth_hum_shape, df_cet_hum_shape, df_pinn_hum_shape, df_chir_hum_shape))
colnames(df_hum_shape_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")


df_oth_hum_size <- (gr_disparity_hum_size$others[[2]])
df_cet_hum_size <- (gr_disparity_hum_size$Cetacea[[2]])
df_pinn_hum_size <- (gr_disparity_hum_size$Pinnipedia[[2]])
df_chir_hum_size <- (gr_disparity_hum_size$Chiroptera[[2]])

df_hum_size_group <- t(rbind(df_oth_hum_size, df_cet_hum_size, df_pinn_hum_size, df_chir_hum_size))
colnames(df_hum_size_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")

library(reshape2)


disparity_hum_shape_env <- melt(df_hum_shape_env)
disparity_hum_size_env <- melt(df_hum_size_env)
disparity_hum_shape_group <- melt(df_hum_shape_group)
disparity_hum_size_group <- melt(df_hum_size_group)






## plot Radius ------------

df_s_rad_shape<- (disparity_rad_shape$terrestrial[[2]])
df_aq_rad_shape <- (disparity_rad_shape$aquatic[[2]])
df_saq_rad_shape <- (disparity_rad_shape$`semi-aquatic`[[2]])
df_a_rad_shape <- (disparity_rad_shape$air[[2]])
df_sa_rad_shape <- (disparity_rad_shape$`semi-aerial`[[2]])

df_rad_shape_env <- t(rbind(df_s_rad_shape,df_aq_rad_shape,df_saq_rad_shape, df_a_rad_shape, df_sa_rad_shape))
colnames(df_rad_shape_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )


df_s_rad_size<- (disparity_rad_size$terrestrial[[2]])
df_aq_rad_size <- (disparity_rad_size$aquatic[[2]])
df_saq_rad_size <- (disparity_rad_size$`semi-aquatic`[[2]])
df_a_rad_size <- (disparity_rad_size$air[[2]])
df_sa_rad_size <- (disparity_rad_size$`semi-aerial`[[2]])

df_rad_size_env <- t(rbind(df_s_rad_size,df_aq_rad_size,df_saq_rad_size, df_a_rad_size, df_sa_rad_size))
colnames(df_rad_size_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )



#--

df_oth_rad_shape <- (gr_disparity_rad_shape$others[[2]])
df_cet_rad_shape <- (gr_disparity_rad_shape$Cetacea[[2]])
df_pinn_rad_shape <- (gr_disparity_rad_shape$Pinnipedia[[2]])
df_chir_rad_shape <- (gr_disparity_rad_shape$Chiroptera[[2]])

df_rad_shape_group <- t(rbind(df_oth_rad_shape, df_cet_rad_shape, df_pinn_rad_shape, df_chir_rad_shape))
colnames(df_rad_shape_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")


df_oth_rad_size <- (gr_disparity_rad_size$others[[2]])
df_cet_rad_size <- (gr_disparity_rad_size$Cetacea[[2]])
df_pinn_rad_size <- (gr_disparity_rad_size$Pinnipedia[[2]])
df_chir_rad_size <- (gr_disparity_rad_size$Chiroptera[[2]])

df_rad_size_group <- t(rbind(df_oth_rad_size, df_cet_rad_size, df_pinn_rad_size, df_chir_rad_size))
colnames(df_rad_size_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")

library(reshape2)


disparity_rad_shape_env <- melt(df_rad_shape_env)
disparity_rad_size_env <- melt(df_rad_size_env)
disparity_rad_shape_group <- melt(df_rad_shape_group)
disparity_rad_size_group <- melt(df_rad_size_group)



## plot Metacarpal ------------


df_s_met_shape<- (disparity_met_shape$terrestrial[[2]])
df_aq_met_shape <- (disparity_met_shape$aquatic[[2]])
df_saq_met_shape <- (disparity_met_shape$`semi-aquatic`[[2]])
df_a_met_shape <- (disparity_met_shape$air[[2]])
df_sa_met_shape <- (disparity_met_shape$`semi-aerial`[[2]])

df_met_shape_env <- t(rbind(df_s_met_shape,df_aq_met_shape,df_saq_met_shape, df_a_met_shape, df_sa_met_shape))
colnames(df_met_shape_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )


df_s_met_size<- (disparity_met_size$terrestrial[[2]])
df_aq_met_size <- (disparity_met_size$aquatic[[2]])
df_saq_met_size <- (disparity_met_size$`semi-aquatic`[[2]])
df_a_met_size <- (disparity_met_size$air[[2]])
df_sa_met_size <- (disparity_met_size$`semi-aerial`[[2]])

df_met_size_env <- t(rbind(df_s_met_size,df_aq_met_size,df_saq_met_size, df_a_met_size, df_sa_met_size))
colnames(df_met_size_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )



#--

df_oth_met_shape <- (gr_disparity_met_shape$others[[2]])
df_cet_met_shape <- (gr_disparity_met_shape$Cetacea[[2]])
df_pinn_met_shape <- (gr_disparity_met_shape$Pinnipedia[[2]])
df_chir_met_shape <- (gr_disparity_met_shape$Chiroptera[[2]])

df_met_shape_group <- t(rbind(df_oth_met_shape, df_cet_met_shape, df_pinn_met_shape, df_chir_met_shape))
colnames(df_met_shape_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")


df_oth_met_size <- (gr_disparity_met_size$others[[2]])
df_cet_met_size <- (gr_disparity_met_size$Cetacea[[2]])
df_pinn_met_size <- (gr_disparity_met_size$Pinnipedia[[2]])
df_chir_met_size <- (gr_disparity_met_size$Chiroptera[[2]])

df_met_size_group <- t(rbind(df_oth_met_size, df_cet_met_size, df_pinn_met_size, df_chir_met_size))
colnames(df_met_size_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")

library(reshape2)


disparity_met_shape_env <- melt(df_met_shape_env)
disparity_met_size_env <- melt(df_met_size_env)
disparity_met_shape_group <- melt(df_met_shape_group)
disparity_met_size_group <- melt(df_met_size_group)



## plot phalanx ------------


df_s_phal_shape<- (disparity_phal_shape$terrestrial[[2]])
df_aq_phal_shape <- (disparity_phal_shape$aquatic[[2]])
df_saq_phal_shape <- (disparity_phal_shape$`semi-aquatic`[[2]])
df_a_phal_shape <- (disparity_phal_shape$air[[2]])
df_sa_phal_shape <- (disparity_phal_shape$`semi-aerial`[[2]])

df_phal_shape_env <- t(rbind(df_s_phal_shape,df_aq_phal_shape,df_saq_phal_shape, df_a_phal_shape, df_sa_phal_shape))
colnames(df_phal_shape_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )


df_s_phal_size<- (disparity_phal_size$terrestrial[[2]])
df_aq_phal_size <- (disparity_phal_size$aquatic[[2]])
df_saq_phal_size <- (disparity_phal_size$`semi-aquatic`[[2]])
df_a_phal_size <- (disparity_phal_size$air[[2]])
df_sa_phal_size <- (disparity_phal_size$`semi-aerial`[[2]])

df_phal_size_env <- t(rbind(df_s_phal_size,df_aq_phal_size,df_saq_phal_size, df_a_phal_size, df_sa_phal_size))
colnames(df_phal_size_env) <- c("terrestrial","aquatic","semi-aquatic","air","semi-aerial" )



#--

df_oth_phal_shape <- (gr_disparity_phal_shape$others[[2]])
df_cet_phal_shape <- (gr_disparity_phal_shape$Cetacea[[2]])
df_pinn_phal_shape <- (gr_disparity_phal_shape$Pinnipedia[[2]])
df_chir_phal_shape <- (gr_disparity_phal_shape$Chiroptera[[2]])

df_phal_shape_group <- t(rbind(df_oth_phal_shape, df_cet_phal_shape, df_pinn_phal_shape, df_chir_phal_shape))
colnames(df_phal_shape_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")


df_oth_phal_size <- (gr_disparity_phal_size$others[[2]])
df_cet_phal_size <- (gr_disparity_phal_size$Cetacea[[2]])
df_pinn_phal_size <- (gr_disparity_phal_size$Pinnipedia[[2]])
df_chir_phal_size <- (gr_disparity_phal_size$Chiroptera[[2]])

df_phal_size_group <- t(rbind(df_oth_phal_size, df_cet_phal_size, df_pinn_phal_size, df_chir_phal_size))
colnames(df_phal_size_group) <- c("others","Cetacea","Pinnipedia","Chiroptera")

library(reshape2)


disparity_phal_shape_env <- melt(df_phal_shape_env)
disparity_phal_size_env <- melt(df_phal_size_env)
disparity_phal_shape_group <- melt(df_phal_shape_group)
disparity_phal_size_group <- melt(df_phal_size_group)




#### saving R data ---------------


disparity_bones_env <- list(disparity_hum_shape_env=disparity_hum_shape_env, disparity_hum_size_env=disparity_hum_size_env,
                            disparity_rad_shape_env=disparity_rad_shape_env, disparity_rad_size_env=disparity_rad_size_env,
                            disparity_met_shape_env=disparity_met_shape_env, disparity_met_size_env=disparity_met_size_env,
                            disparity_phal_shape_env=disparity_phal_shape_env, disparity_phal_size_env=disparity_phal_size_env)
save(disparity_bones_env, file="disparity_bones_env.RData")


disparity_bones_group <- list(disparity_hum_shape_group=disparity_hum_shape_group, disparity_hum_size_group=disparity_hum_size_group,
                            disparity_rad_shape_group=disparity_rad_shape_group, disparity_rad_size_group=disparity_rad_size_group,
                            disparity_met_shape_group=disparity_met_shape_group, disparity_met_size_group=disparity_met_size_group,
                            disparity_phal_shape_group=disparity_phal_shape_group, disparity_phal_size_group=disparity_phal_size_group)
save(disparity_bones_group, file="disparity_bones_group.RData")



