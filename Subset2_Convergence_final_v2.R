
##==========================================================================================================##
##                                                                                                          ##
##                                       CONVERGENCE PER BONE                                               ##
##                                                                                                          ##
##==========================================================================================================##

setwd("C:/Users/Lab/OneDrive/Academia/Rothier et al/Locomotion in fluids/Analyses/analyses_2ndround")


library(mvMORPH); library(geiger)


load("subset2_dat.Rdata")

tree <- subset2_dat$tree
media <- subset2_dat$media

load("bones_models.Rdata")

#loading the best fit models
fit_hum2 <- bones_models$fit_OU_hum2
fit_hum1 <- bones_models$fit_OU_hum

fit_rad2 <- bones_models$fit_OU_rad2
fit_rad1 <- bones_models$fit_OU_rad


fit_met2 <- bones_models$fit_OU_met2
fit_met1 <- bones_models$fit_OU_met


fit_phal2 <- bones_models$fit_OU_phal2
fit_phal1 <- bones_models$fit_OU_phal



library(RRphylo)
# read search.conv debugged version




  ## preparing dataset


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




#### 1. Humerus ---------------------------------------------------------------------- ####

sc_hum1 <- fit_hum1$residuals #size
sc_hum2 <- fit_hum2$residuals #shape

  ## search.conv with size -----
name.check(tree, sc_hum1)

SC.Full_Aquatic_hum1 <- search.conv(tree=tree, y=sc_hum1, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_hum1 <- search.conv(tree=tree, y=sc_hum1, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_hum1 <- search.conv(tree=tree, y=sc_hum1, state=Semi_Aerial,declust=TRUE)


# testing the new plot function:
plotConv(SC.Semi_Aerial_hum1,sc_hum1,variable=1,state=Semi_Aerial)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue"),lty=1),
              points.args = list(pch=c(21,21),bg="gray"),
              legend.args = list(pch=c(21,21),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))

?plotConv
  # searching convergence within between full states and semi-states

SC_Aerial_hum1 <- search.conv(tree=tree,y=sc_hum1, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_hum1 <- search.conv(tree=tree,y=sc_hum1, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_hum1 <- rbind(SC.Full_Aquatic_hum1$state.res, SC.Semi_Aquatic_hum1$state.res, SC.Semi_Aerial_hum1$state.res)
conv_2st_hum1 <- rbind(SC_Aerial_hum1$state.res, SC_Aquatic_hum1$state.res)



## search.conv without size -----

name.check(tree, sc_hum2)

SC.Full_Aquatic_hum2 <- search.conv(tree=tree, y=sc_hum2, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_hum2 <- search.conv(tree=tree, y=sc_hum2, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_hum2 <- search.conv(tree=tree, y=sc_hum2, state=Semi_Aerial,declust=TRUE)


SC_Aerial_hum2 <- search.conv(tree=tree,y=sc_hum2, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_hum2 <- search.conv(tree=tree,y=sc_hum2, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_hum2 <- rbind(SC.Full_Aquatic_hum2$state.res, SC.Semi_Aquatic_hum2$state.res, SC.Semi_Aerial_hum2$state.res)
conv_2st_hum2 <- rbind(SC_Aerial_hum2$state.res, SC_Aquatic_hum2$state.res)




#### 2. Radius ---------------------------------------------------------------------- ####


## pca -----

sc_rad1 <- fit_rad1$residuals #size
sc_rad2 <- fit_rad2$residuals #shape


## search.conv with size -----

SC.Full_Aquatic_rad1 <- search.conv(tree=tree, y=sc_rad1, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_rad1 <- search.conv(tree=tree, y=sc_rad1, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_rad1 <- search.conv(tree=tree, y=sc_rad1, state=Semi_Aerial,declust=TRUE)



# searching convergence within between full states and semi-states

SC_Aerial_rad1 <- search.conv(tree=tree,y=sc_rad1, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_rad1 <- search.conv(tree=tree,y=sc_rad1, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_rad1 <- rbind(SC.Full_Aquatic_rad1$state.res, SC.Semi_Aquatic_rad1$state.res, SC.Semi_Aerial_rad1$state.res)
conv_2st_rad1 <- rbind(SC_Aerial_rad1$state.res, SC_Aquatic_rad1$state.res)



## search.conv without size -----

name.check(tree, sc_rad2)
SC.Full_Aquatic_rad2 <- search.conv(tree=tree, y=sc_rad2, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_rad2 <- search.conv(tree=tree, y=sc_rad2, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_rad2 <- search.conv(tree=tree, y=sc_rad2, state=Semi_Aerial,declust=TRUE)


SC_Aerial_rad2 <- search.conv(tree=tree,y=sc_rad2, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_rad2 <- search.conv(tree=tree,y=sc_rad2, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_rad2 <- rbind(SC.Full_Aquatic_rad2$state.res, SC.Semi_Aquatic_rad2$state.res, SC.Semi_Aerial_rad2$state.res)
conv_2st_rad2 <- rbind(SC_Aerial_rad2$state.res, SC_Aquatic_rad2$state.res)




#### 3. Metacarpal ---------------------------------------------------------------------- ####

## pca -----

sc_met1 <- fit_met1$residuals #size
sc_met2 <- fit_met2$residuals #shape

## search.conv with size -----

name.check(tree, sc_met1)


# searching convergence within between full states and semi-states

SC.Full_Aquatic_met1 <- search.conv(tree=tree, y=sc_met1, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_met1 <- search.conv(tree=tree, y=sc_met1, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_met1 <- search.conv(tree=tree, y=sc_met1, state=Semi_Aerial,declust=TRUE)



# searching convergence within between full states and semi-states

SC_Aerial_met1 <- search.conv(tree=tree,y=sc_met1, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_met1 <- search.conv(tree=tree,y=sc_met1, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_met1 <- rbind(SC.Full_Aquatic_met1$state.res, SC.Semi_Aquatic_met1$state.res, SC.Semi_Aerial_met1$state.res)
conv_2st_met1 <- rbind(SC_Aerial_met1$state.res, SC_Aquatic_met1$state.res)



## search.conv without size -----

name.check(tree, sc_met2)

SC.Full_Aquatic_met2 <- search.conv(tree=tree, y=sc_met2, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_met2 <- search.conv(tree=tree, y=sc_met2, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_met2 <- search.conv(tree=tree, y=sc_met2, state=Semi_Aerial,declust=TRUE)


SC_Aerial_met2 <- search.conv(tree=tree,y=sc_met2, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_met2 <- search.conv(tree=tree,y=sc_met2, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_met2 <- rbind(SC.Full_Aquatic_met2$state.res, SC.Semi_Aquatic_met2$state.res, SC.Semi_Aerial_met2$state.res)
conv_2st_met2 <- rbind(SC_Aerial_met2$state.res, SC_Aquatic_met2$state.res)


plotConv(SC_Aerial_met2,sc_met2,variable=1,state=Aerial)->pc
pc$plotPChull(chull.args = list(border=c("gray70","blue","red"),lty=1),
              points.args = list(pch=c(23,22,21),bg="gray"),
              legend.args = list(pch=c(23,22,21),x="top"))

pc$plotPolar(polar.args = list(clockwise=TRUE,start=0,rad.col="black",grid.col="black"),
             polygon.args = list(line.col="green",poly.col=NA,lwd=2),
             line.args = list(line.col="deeppink",lty=2,lwd=3))


#### 4. Phalanx ---------------------------------------------------------------------- ####


sc_phal1 <- fit_phal1$residuals #size
sc_phal2 <- fit_phal2$residuals #shape


## search.conv with size -----

name.check(tree, sc_phal1)

SC.Full_Aquatic_phal1 <- search.conv(tree=tree, y=sc_phal1, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_phal1 <- search.conv(tree=tree, y=sc_phal1, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_phal1 <- search.conv(tree=tree, y=sc_phal1, state=Semi_Aerial,declust=TRUE)



# searching convergence within between full states and semi-states

SC_Aerial_phal1 <- search.conv(tree=tree,y=sc_phal1, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_phal1 <- search.conv(tree=tree,y=sc_phal1, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_phal1 <- rbind(SC.Full_Aquatic_phal1$state.res, SC.Semi_Aquatic_phal1$state.res, SC.Semi_Aerial_phal1$state.res)
conv_2st_phal1 <- rbind(SC_Aerial_phal1$state.res, SC_Aquatic_phal1$state.res)



## search.conv without size -----

name.check(tree, sc_phal2)

SC.Full_Aquatic_phal2 <- search.conv(tree=tree, y=sc_phal2, state=Full_Aquatic,declust=FALSE) 
SC.Semi_Aquatic_phal2 <- search.conv(tree=tree, y=sc_phal2, state=Semi_Aquatic,declust=TRUE)
SC.Semi_Aerial_phal2 <- search.conv(tree=tree, y=sc_phal2, state=Semi_Aerial,declust=TRUE)


SC_Aerial_phal2 <- search.conv(tree=tree,y=sc_phal2, state=Aerial, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)
SC_Aquatic_phal2 <- search.conv(tree=tree,y=sc_phal2, state=Aquatic, nsim=100,clus=2/parallel::detectCores(), declust=TRUE)


conv_1st_phal2 <- rbind(SC.Full_Aquatic_phal2$state.res, SC.Semi_Aquatic_phal2$state.res, SC.Semi_Aerial_phal2$state.res)
conv_2st_phal2 <- rbind(SC_Aerial_phal2$state.res, SC_Aquatic_phal2$state.res)


