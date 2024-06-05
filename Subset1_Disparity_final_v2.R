##==========================================================================================================##
##                                                                                                          ##
##   Subset 1                                                                                               ##
##                                                                                                          ##
##==========================================================================================================##


# two subsets
# subs1 = species without na (remove those without digit length)
      # Calculating disparity and convergence of overall dataset, with and without allometry

# subs2 = species without ph length (golden moles), don't consider digit length


setwd("C:/Users/Lab/OneDrive/Academia/Rothier et al/Locomotion in fluids/5.Functional Ecology/Review 1/analyses_review")


dat <- read.csv("data_2analyses.csv", sep=",", header=T, na.strings = "NA")
rownames(dat) <- dat$Species
dat <- dat[2:31]

library(ape); library(geiger)
library(evobiR); library(phytools); library(reshape2)
library(ggplot2)


tree <- read.tree("mytree.nex")
plot(tree, show.tip.label = FALSE); axisPhylo()

dat <- ReorderData(tree, dat, taxa.names="rownames")

dat$Transf_Mass <- log10(dat$Body_Mass^(1/3)) # transform bm in linear scale
str(dat)




##==========================================================================================================##
##                                                                                                          ##
####                                   1.      Dataset preparation                                        ####
##                                                                                                          ##
##==========================================================================================================##


str(dat)
dat1 <- na.omit(dat)
str(dat1) # 797 tips
name.check(tree, dat1)

to_remove1 <- tree$tip.label[!tree$tip.label%in%rownames(dat1)]
t1 <- drop.tip(tree, to_remove1)
name.check(t1, dat1) 


library(psych)
lgm <- log10(apply(dat1[,c(9:29,31)], 1, geometric.mean))
limb <- as.matrix(log10(dat1[,9:29]))
traits1 <- cbind(lgm,limb)
head(traits1)



##==========================================================================================================##
##                                                                                                          ##
####                                      2.    Linear models                                             ####
##                                                                                                          ##
##==========================================================================================================##

library(mvMORPH)




#### ------ allometry removed ------ ####

#fit_BM2 <- mvgls(limb~lgm, tree=t1, model="BM", penalty="LASSO")
#fit_OU2 <- mvgls(limb~lgm, tree=t1, model="OU", penalty="LASSO")
#fit_EB2 <- mvgls(limb~lgm, tree=t1, model="EB", penalty="LASSO")
#GIC(fit_BM2); GIC(fit_OU2); GIC(fit_EB2) # OU

#d1_BM_removed <- list(fit_BM2, fit_OU2, fit_EB2)
#names(d1_BM_removed) <- c("fit_BM2", "fit_OU2", "fit_EB2")



#### ------ with allometry ------ ####

#fit_BM1 <- mvgls(limb~1, tree=t1, model="BM", penalty="LASSO")
#fit_OU1 <- mvgls(limb~1, tree=t1, model="OU", penalty="LASSO")
#fit_EB1 <- mvgls(limb~1, tree=t1, model="EB", penalty="LASSO")
#GIC(fit_BM1); GIC(fit_OU1); GIC(fit_EB1)
#OU has the best fit

#d1_with_BM <- list(fit_BM1, fit_OU1, fit_EB1)
#names(d1_with_BM) <- c("fit_BM1", "fit_OU1", "fit_EB1")
#lapply(d1_with_BM, function(x) sum(diag(x$sigma$Pinv)))



#subset1 <- list(d1_with_BM, d1_BM_removed)
#names(subset1) <- c("d1_with_BM", "d1_BM_removed")
#save(subset1, file="subset1_all_measurements.RData") 




## Lasso penalization serves to:
# 1) the estimation of the correlations is improved because some shrinkage has been applied to elements of the 
#covariance matrix that have some uncertainty
# 2) the lasso penalization is introducing zeros in the (partial) correlations/covariances. This means that it's
#automatically performing some kind of model selection by setting to zero the values that are statistically 
#insignificant.
# The OU model implemented in "mvgls" is simpler than the one in mvOU. But for large datasets, this function is better 
# (mvOU cannot handle big dimensions, say when #traits*#species > 2000-3000, this starts to be extremely time consuming)





#### ------ load models ------ ####

load("subset1_all_measurements.Rdata")
d1_with_BM <- subset1$d1_with_BM
fit_OU_raw <- d1_with_BM$fit_OU1
d1_without_BM <- subset1$d1_BM_removed
fit_OU_gm <- d1_without_BM$fit_OU2





##==========================================================================================================##
##                                                                                                          ##
####                                     3.    Manova per group                                           ####
##                                                                                                          ##
##==========================================================================================================##


media <- as.factor(dat1$Media)
names(media) <- rownames(dat1)
data.gls = list(limb=limb, tree=t1, lgm=lgm, media=media)


pgls.resid <- mvgls(limb~lgm*media, tree=t1, data=data.gls, model="OU", method="LL")
pairwise.glh(pgls.resid, test="Pillai")
(multivariate_test <- manova.gls(pgls.resid, nperm=999, type="II", test="Pillai"))

pgls.raw<- mvgls(limb~media, tree=t1, data=data.gls, model="OU", method="LL")
pairwise.glh(pgls.raw, test="Pillai")
(multivariate_test <- manova.gls(pgls.raw, nperm=999, type="II", test="Pillai"))


##==========================================================================================================##
##                                                                                                          ##
####                      4.      Shape linear models (traits ~ gm)                                       ####
##                                                                                                          ##
##==========================================================================================================##



#### ------ phylo PCA ------ ####


pca2 <- mvgls.pca(fit_OU_gm)
row.names(pca2$vectors) <- colnames(limb)

eigen.df2 <- as.data.frame(pca2$values)

propvar <- function(x) {
  x/sum(x)}

propvar_pca2 <- apply(eigen.df2, 2, propvar)
cumvar_pca2 <- cumsum(propvar_pca2)

sc_pca2 <- as.data.frame(pca2$scores)
sc_pca2$Group <-  dat1$Group
sc_pca2$Media <- dat1$Media
sc_pca2$MonophyleticGroup <- dat1$major_clade



pca_plot_shape <- ggplot(sc_pca2, aes(V1, V2, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=as.vector(glasbey(26)))+
  theme_minimal()+
  xlab("PC1 (30.12%)")+
  ylab("PC2 (14.28%)")+
  guides(col=guide_legend(ncol=1))+
 # theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.61), breaks=seq(-0.4, 0.6, by=0.2))

library(plotly)

ggplotly(pca_plot_shape, tooltip=c("text", "x", "y", "group"))



#### ------ regular PCA ------ ####

pca2.r <- prcomp(fit_OU_gm$residuals)
summary(pca2.r)
?prcomp

pca2.r$sdev^2 



##==========================================================================================================##
##                                                                                                          ##
####                      5.      Shape disparity                                                         ####
##                                                                                                          ##
##==========================================================================================================##
library(dispRity)

m1 <- as.character(dat1$Media)
names(m1) <- rownames(dat1$Species)
unique(m1)

ord_matrix_shape <- as.matrix(fit_OU_gm$residuals)

env_subsets_shape <- custom.subsets(ord_matrix_shape, tree=t1, group = list(
  "terrestrial" = which(m1 == "terrestrial"),
  "aquatic" = which(m1 == "aquatic"),
  "semi-aquatic" = which(m1 == "semi-aquatic"),
  "air" = which(m1 == "air"),
  "semi-aerial" = which(m1 == "semi-aerial")))

#env_bootstrapped_shape <- boot.matrix(env_subsets_shape, bootstraps = 1000,
#                                rarefaction = TRUE, verbose=TRUE) 

#env_disparity_shape <- dispRity(env_bootstrapped_shape, metric = c(sum, variances), verbose=TRUE)

#save(env_disparity_shape, file="env_disparity_shape_wholelimb_v2.RData")
load("env_disparity_shape_wholelimb_v2.RData")
str(env_disparity_shape)
summary(env_disparity_shape)
class(env_disparity_shape)
plot.dispRity(env_disparity_shape, rarefaction=TRUE)
plot.dispRity(env_disparity_shape)

disparity_shape <- env_disparity_shape$disparity
str(env_disparity_shape$subsets)

means_shape <- lapply(disparity_shape, function(x) mean(x[[2]]))
sd_shape <- lapply(disparity_shape, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_shape, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_shape,bhatt.coeff, "pairwise")

env_disparity_shape$subsets



#### ------ Disparity considering monophyletic groups------ ####


m2 <- as.character(sc_pca2$MonophyleticGroup)
names(m2) <- rownames(sc_pca2)
unique(m2)

## Creating the table that contain the elements and their attributes
group_subsets <- custom.subsets(ord_matrix_shape, tree=t1, group = list(
  "others" = which(m2 == "others"),
  "Cetacea" = which(m2 == "Cetacea"),
  "Pinnipedia" = which(m2 == "Pinnipedia"),
  "Chiroptera" = which(m2 == "Chiroptera")))


## Bootstrapping the data
#group_bootstrapped_shape <- boot.matrix(group_subsets, bootstraps = 1000,
#                                  rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

#group_disparity_shape <- dispRity(group_bootstrapped_shape, metric = c(sum, variances), verbose=TRUE)
#save(group_disparity_shape, file="group_disparity_shape_wholelimb_v2.RData")
load("group_disparity_shape_wholelimb_v2.RData")

summary(group_disparity_shape)

plot.dispRity(group_disparity_shape, rarefaction=TRUE)
plot.dispRity(group_disparity_shape)

gr_disparity_shape <- group_disparity_shape$disparity
gr_means_shape <- lapply(gr_disparity_shape, function(x) mean(x[[2]]))
gr_sd_shape <- lapply(gr_disparity_shape, function(x) sd(x[[2]]))


test.dispRity(group_disparity_shape, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_shape,bhatt.coeff, "pairwise")





##==========================================================================================================##
##                                                                                                          ##
####                         6.    Linear models with allometry  (traits ~ 1)                             ####
##                                                                                                          ##
##==========================================================================================================##




#### ------ phylo PCA ------ ####


pca1 <- mvgls.pca(fit_OU_raw)
row.names(pca1$vectors) <- colnames(limb)
data.frame(pca1$values)->eigen.df
propvar_pca1 <- apply(eigen.df, 2, propvar)
cumvar_pca1 <- cumsum(propvar_pca1)

sc_pca1 <- as.data.frame(pca1$scores)
sc_pca1$Group <-  dat1$Group
sc_pca1$Media <- dat1$Media
sc_pca1$MonophyleticGroup <- dat1$major_clade
sc_pca1$lgm <- lgm
sc_pca1$bm <- dat1$Transf_Mass

library(nlme)


#### ------ regular PCA ------ ####

pca1.r <- prcomp(fit_OU_raw$residuals)
summary(pca1.r)
pca1.r$sdev^2 

##==========================================================================================================##
##                                                                                                          ##
####                      7.      size disparity                                                          ####
##                                                                                                          ##
##==========================================================================================================##


#### ------ Disparity between ecological groups ------ ####


ord_matrix2 <- pca1$scores


m1 <- as.character(dat1$Media)
names(m1) <- rownames(sc_pca2)
unique(m1)

ord_matrix_size <- as.matrix(fit_OU_raw$residuals)

env_subsets_size <- custom.subsets(ord_matrix_size, tree=t1, group = list(
  "terrestrial" = which(m1 == "terrestrial"),
  "aquatic" = which(m1 == "aquatic"),
  "semi-aquatic" = which(m1 == "semi-aquatic"),
  "air" = which(m1 == "air"),
  "semi-aerial" = which(m1 == "semi-aerial")))

#env_bootstrapped_size <- boot.matrix(env_subsets_size, bootstraps = 1000,
#                                 rarefaction = TRUE, verbose=TRUE) 

#env_disparity_size <- dispRity(env_bootstrapped_size, metric = c(sum, variances), verbose=TRUE)

#save(env_disparity_size, file="env_disparity_size_wholelimb_v2.RData")
load("env_disparity_size_wholelimb_v2.RData")

summary(env_disparity_size)
plot.dispRity(env_disparity_size, rarefaction=TRUE)
plot.dispRity(env_disparity_size)
disparity_size <- env_disparity_size$disparity

means_size <- lapply(disparity_size, function(x) mean(x[[2]]))
sd_size <- lapply(disparity_size, function(x) sd(x[[2]]))

par(bty = "n")

test.dispRity(env_disparity_size, test = t.test, correction="bonferroni")  
test.dispRity(env_disparity_size,bhatt.coeff, "pairwise")




#### ------ Disparity considering monophyletic groups------ ####




## Creating the table that contain the elements and their attributes
group_subsets_size <- custom.subsets(ord_matrix_size, tree=t1, group = list(
  "others" = which(m2 == "others"),
  "Cetacea" = which(m2 == "Cetacea"),
  "Pinnipedia" = which(m2 == "Pinnipedia"),
  "Chiroptera" = which(m2 == "Chiroptera")))


## Bootstrapping the data
#group_bootstrapped_size <- boot.matrix(group_subsets_size, bootstraps = 1000,
#                                  rarefaction = TRUE, verbose=TRUE)

##Calculating disparity

#group_disparity_size <- dispRity(group_bootstrapped_size, metric = c(sum, variances), verbose=TRUE)
#save(group_disparity_size, file="group_disparity_size_wholelimb.RData")
load("group_disparity_size_wholelimb.RData")

summary(group_disparity_size)

plot.dispRity(group_disparity_size, rarefaction=TRUE)
plot.dispRity(group_disparity_size)

gr_disparity_size <- group_disparity_size$disparity
gr_means_size <- lapply(gr_disparity_size, function(x) mean(x[[2]]))
gr_sd_size <- lapply(gr_disparity_size, function(x) sd(x[[2]]))


test.dispRity(group_disparity_size, test = t.test, correction="bonferroni")  
test.dispRity(group_disparity_size,bhatt.coeff, "pairwise")



##==========================================================================================================##
##                                                                                                          ##
####                                    5.  Phyl Anova                                                    ####
##                                                                                                          ##
##==========================================================================================================##



### phyl aov using pcs from raw data PCA


media <- dat1$Media
names(media) <- rownames(dat1)
pc1 <- sc_pca1$V1
names(pc1) <- rownames(sc_pca1)
pc2 <- sc_pca1$V2
names(pc2) <- rownames(sc_pca1)
pc3 <- sc_pca1$V3
names(pc3) <- rownames(sc_pca1)
pc4 <- sc_pca1$V4
names(pc4) <- rownames(sc_pca1)

pc1.gm <- sc_pca2$V1
names(pc1.gm) <- rownames(sc_pca2)
pc2.gm <- sc_pca2$V2
names(pc2.gm) <- rownames(sc_pca2)
pc3.gm <- sc_pca2$V3
names(pc3.gm) <- rownames(sc_pca2)
pc4.gm <- sc_pca2$V4
names(pc4.gm) <- rownames(sc_pca2)


mass <- dat1$Transf_Mass
names(mass) <- rownames(dat1)
aov_bodysize <- phylANOVA(t1, media, mass, posthoc = TRUE, p.adj="holm", nsim=1000)

boxplot(mass~media)
means_bm <- as.data.frame(dat1 %>% group_by(Media) %>% summarise(mean(Body_Mass)))

aov_bodysize <- phylANOVA(t1, media, lgm, posthoc = TRUE, p.adj="holm", nsim=1000) 

aov1 <- phylANOVA(t1, media, pc1, posthoc = TRUE, p.adj="holm", nsim=10000) 
aov2 <- phylANOVA(t1, media, pc2, posthoc = TRUE, p.adj="holm", nsim=10000)
aov3 <- phylANOVA(t1, media, pc3, posthoc = TRUE, p.adj="holm", nsim=10000)
aov4 <- phylANOVA(t1, media, pc3, posthoc = TRUE, p.adj="holm", nsim=10000)
aov5 <- phylANOVA(t1, media, pc1.gm, posthoc = TRUE, p.adj="holm", nsim=10000)
aov6 <- phylANOVA(t1, media, pc2.gm, posthoc = TRUE, p.adj="holm", nsim=10000)
aov7 <- phylANOVA(t1, media, pc3.gm, posthoc = TRUE, p.adj="holm", nsim=10000)
aov8 <- phylANOVA(t1, media, pc3.gm, posthoc = TRUE, p.adj="holm", nsim=10000)



?phylANOVA

### phyl aov using pcs from raw data PCA






##==========================================================================================================##
##                                                                                                          ##
####                                        6.  Convergence                                               ####
##                                                                                                          ##
##==========================================================================================================##



# Preparing dataset

names(media) ==t1$tip.label

unique(media)

library(dplyr)
Full_Aquatic <- recode(media, terrestrial="nostate",
                       'semi-aquatic'="nostate",
                       aquatic= "Full_Aquatic",
                       air="nostate", 
                       'semi-aerial'= "nostate")
names(Full_Aquatic) <- names(media)


Semi_Aquatic <- recode(media, terrestrial="nostate",
                       'semi-aquatic'="Semi_Aquatic",
                       aquatic= "nostate",
                       air="nostate", 
                       'semi-aerial'= "nostate")
names(Semi_Aquatic) <- names(media)


Semi_Aerial <- recode(media, terrestrial="nostate",
                       'semi-aquatic'="nostate",
                       aquatic= "nostate",
                       air="nostate", 
                       'semi-aerial'= "Semi_Aerial")
names(Semi_Aerial ) <- names(media)


Aerial <- recode(media, terrestrial="nostate",
                 'semi-aquatic'="nostate",
                 aquatic= "nostate",
                 air="Full_Aerial", 
                 'semi-aerial'= "Semi_Aerial")
names(Aerial ) <- names(media)


Aquatic <- recode(media, terrestrial="nostate",
                       'semi-aquatic'="Semi_Aquatic",
                       aquatic= "Full_Aquatic",
                       air="nostate", 
                       'semi-aerial'= "nostate")
names(Aquatic) <- names(media)




library(RRphylo)


#### ------ size removed ------ ####

shape_resid <- as.matrix(fit_OU_gm$residuals)


## searching convergence within a single state
name.check(t1, shape_resid)

SC.Full_Aquatic_form <- search.conv(tree=t1, y=shape_resid, state=Full_Aquatic,declust=FALSE) # clades clustred, so declust =F
plotConv(SC.Full_Aquatic_form,shape_resid,variable=1,state=Full_Aquatic)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(23,22),bg="gray"),
              legend.args = list(pch=c(23,22),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))


SC.Semi_Aquatic_form <- search.conv(tree=t1, y=shape_resid, state=Semi_Aquatic,declust=TRUE)
plotConv(SC.Semi_Aquatic_form, shape_resid,variable=1,state=Semi_Aquatic)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(23,22),bg="gray"),
              legend.args = list(pch=c(23,22),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))


SC.Semi_Aerial_form <- search.conv(tree=t1, y=shape_resid, state=Semi_Aerial,declust=TRUE,)
plotConv(SC.Semi_Aerial_form,shape_resid,variable=1,state=Semi_Aerial)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(23,22),bg="gray"),
              legend.args = list(pch=c(23,22),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))



## searching convergence within between full states and semi-states

SC_Aerial_form <- search.conv(tree=t1,y=shape_resid, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
plotConv(SC_Aerial_form,shape_resid,variable=1,state=Aerial)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue","red"),lty=1),
              points.args = list(pch=c(23,22,21),bg="gray"),
              legend.args = list(pch=c(23,22,21),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))



SC_Aquatic_form <- search.conv(tree=t1,y=shape_resid, state=Aquatic, nsim=100,clus=2/parallel::detectCores())
plotConv(SC_Aquatic_form,shape_resid,variable=1,state=Aquatic)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue","red"),lty=1),
              points.args = list(pch=c(23,22,21),bg="gray"),
              legend.args = list(pch=c(23,22,21),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))




conv_1st_form<- rbind(SC.Full_Aquatic_form$state.res, SC.Semi_Aquatic_form$state.res, SC.Semi_Aerial_form$state.res)
conv_2st_form <- rbind(SC_Aerial_form$state.res, SC_Aquatic_form$state.res)





#### ------ with size ------ ####


raw_resid <- as.matrix(fit_OU_raw$residuals)



## searching convergence within a single state
name.check(t1, raw_resid)

SC.Full_Aquatic_size <- search.conv(tree=t1, y=raw_resid, state=Full_Aquatic,declust=FALSE) # clades clustred, so declust =F
plotConv(SC.Full_Aquatic_size,raw_resid,variable=1,state=Full_Aquatic)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(23,22),bg="gray"),
              legend.args = list(pch=c(23,22),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))


SC.Semi_Aquatic_size <- search.conv(tree=t1, y=raw_resid, state=Semi_Aquatic,declust=TRUE)
plotConv(SC.Semi_Aquatic_size, raw_resid,variable=1,state=Semi_Aquatic)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(23,22),bg="gray"),
              legend.args = list(pch=c(23,22),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))


SC.Semi_Aerial_size <- search.conv(tree=t1, y=raw_resid, state=Semi_Aerial,declust=TRUE,)
plotConv(SC.Semi_Aerial_size,raw_resid,variable=1,state=Semi_Aerial)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(23,22),bg="gray"),
              legend.args = list(pch=c(23,22),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))



## searching convergence within between full states and semi-states

SC_Aerial_size <- search.conv(tree=t1,y=raw_resid, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
plotConv(SC_Aerial_size,raw_resid,variable=1,state=Aerial)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue","red"),lty=1),
              points.args = list(pch=c(23,22,21),bg="gray"),
              legend.args = list(pch=c(23,22,21),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))



SC_Aquatic_size <- search.conv(tree=t1,y=raw_resid, state=Aquatic, nsim=100,clus=2/parallel::detectCores())
plotConv(SC_Aquatic_size,raw_resid,variable=1,state=Aquatic)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue","red"),lty=1),
              points.args = list(pch=c(23,22,21),bg="gray"),
              legend.args = list(pch=c(23,22,21),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))




conv_1st_size<- rbind(SC.Full_Aquatic_size$state.res, SC.Semi_Aquatic_size$state.res, SC.Semi_Aerial_size$state.res)
conv_2st_size <- rbind(SC_Aerial_size$state.res, SC_Aquatic_size$state.res)







##==========================================================================================================##
##                                                                                                          ##
####             7.    Creating dataframes for plots                                                      ####
##                                                                                                          ##
##==========================================================================================================##



#### pcas ####

#### disparity ####

# shape environment
df_s.shape.e <- (disparity_shape$terrestrial[[2]]) # taking just bootstrap rep from non-rarefied elements
df_aq.shape.e <- (disparity_shape$aquatic[[2]])
df_saq.shape.e <- (disparity_shape$`semi-aquatic`[[2]])
df_a.shape.e <- (disparity_shape$air[[2]])
df_sa.shape.e <- (disparity_shape$`semi-aerial`[[2]])

df.shape.e <- t(rbind(df_s.shape.e,df_aq.shape.e,df_saq.shape.e, df_a.shape.e, df_sa.shape.e))
colnames(df.shape.e) <- c("terrestrial", "aquatic", "semi-aquatic", "air", "semi-aerial")
melt.df.shape.e <- melt(df.shape.e)

# shape group
df_o.shape.g <- (gr_disparity_shape$others[[2]])
df_ce.shape.g <- (gr_disparity_shape$Cetacea[[2]])
df_p.shape.g <- (gr_disparity_shape$Pinnipedia[[2]])
df_ch.shape.g <- (gr_disparity_shape$Chiroptera[[2]])

df.shape.g <- t(rbind(df_o.shape.g, df_ce.shape.g, df_p.shape.g, df_ch.shape.g))
colnames(df.shape.g) <- c("others", "Cetacea", "Pinnipedia", "Chiroptera")


# size environment
df_s.size.e <- (disparity_size$terrestrial[[2]])
df_aq.size.e <- (disparity_size$aquatic[[2]])
df_saq.size.e <- (disparity_size$`semi-aquatic`[[2]])
df_a.size.e <- (disparity_size$air[[2]])
df_sa.size.e <- (disparity_size$`semi-aerial`[[2]])

df.size.e <- t(rbind(df_s.size.e,df_aq.size.e,df_saq.size.e, df_a.size.e, df_sa.size.e))
colnames(df.size.e) <- c("terrestrial", "aquatic", "semi-aquatic", "air", "semi-aerial")


# size group
df_o.size.g <- (gr_disparity_size$others[[2]])
df_ce.size.g <- (gr_disparity_size$Cetacea[[2]])
df_p.size.g <- (gr_disparity_size$Pinnipedia[[2]])
df_ch.size.g <- (gr_disparity_size$Chiroptera[[2]])

df.size.g <- t(rbind(df_o.size.g, df_ce.size.g, df_p.size.g, df_ch.size.g))
colnames(df.size.g) <- c("others", "Cetacea", "Pinnipedia", "Chiroptera")


## saving DF for plots
disparity_shape_env <- melt(df.shape.e)
disparity_size_env <- melt(df.size.e )
disparity_shape_group <- melt(df.shape.g )
disparity_size_group <- melt(df.size.g )


df_plots <- list(dat1=dat1, t1=t1, sc_pca_raw=sc_pca1, sc_pca_shape=sc_pca2, reg_pca_raw=pca1.r, reg_pca_shape=pca2.r,
                 disparity_shape_env=disparity_shape_env, disparity_shape_group=disparity_shape_group,
                 disparity_size_env=disparity_size_env, disparity_size_group=disparity_size_group)

save(df_plots, file="df_plots.RData")
str(df_plots)








