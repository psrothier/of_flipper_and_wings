setwd("C:/Users/Lab/OneDrive/Academia/Rothier et al/Locomotion in fluids/5.Functional Ecology/Review 1/analyses_review")

##==========================================================================================================##
##                                                                                                          ##
####                                            - Plots -                                                 ####
##                                                                                                          ##
##==========================================================================================================##
library(ggplot2)
library(pals)
library(cowplot)
library(phytools)
library(grid)
library(gridExtra)
library(ggeasy)

## loading objects
## ----- Whole limb ----- ##

load("df_plots.RData") ##### include dat2

tree <- df_plots$t1
dat1 <- df_plots$dat1
sc_pca1 <- df_plots$sc_pca_raw
sc_pca2<- df_plots$sc_pca_shape
reg_pca_raw <- df_plots$reg_pca_raw
reg_pca_shape <- df_plots$reg_pca_shape
disparity_shape_env  <- df_plots$disparity_shape_env
disparity_shape_group <- df_plots$disparity_shape_group
disparity_size_env  <- df_plots$disparity_size_env
disparity_size_group <- df_plots$disparity_size_group

## ----- Each bone ----- ##

load("bones_pca.Rdata")

sc_pca_hum_gm <- bones_pca$sc_pca_hum_gm
sc_pca_rad_gm <- bones_pca$sc_pca_rad_gm
sc_pca_met_gm <- bones_pca$sc_pca_met_gm
sc_pca_phal_gm <- bones_pca$sc_pca_phal_gm

sc_pca_hum_raw <- bones_pca$sc_pca_hum_raw
sc_pca_rad_raw <- bones_pca$sc_pca_rad_raw
sc_pca_met_raw <- bones_pca$sc_pca_met_raw
sc_pca_phal_raw <- bones_pca$sc_pca_phal_raw


load("bones_opca.Rdata")

sc_opca_hum_shape <- bones_opca$sc_opca_hum_shape
sc_opca_rad_shape <- bones_opca$sc_opca_rad_shape
sc_opca_met_shape <- bones_opca$sc_opca_met_shape
sc_opca_phal_shape <- bones_opca$sc_opca_phal_shape

sc_opca_hum_size <- bones_opca$sc_opca_hum_size
sc_opca_rad_size <- bones_opca$sc_opca_rad_size
sc_opca_met_size <- bones_opca$sc_opca_met_size
sc_opca_phal_size <- bones_opca$sc_opca_phal_size



load("disparity_bones_env.RData")
str(disparity_bones_env)

disp_hum_shape.e <- disparity_bones_env$disparity_hum_shape_env
disp_rad_shape.e <- disparity_bones_env$disparity_rad_shape
disp_met_shape.e <- disparity_bones_env$disparity_met_shape_env
disp_phal_shape.e <- disparity_bones_env$ disparity_phal_shape_env

disp_hum_size.e <- disparity_bones_env$disparity_hum_size_env
disp_rad_size.e <- disparity_bones_env$disparity_rad_size_env
disp_met_size.e <- disparity_bones_env$disparity_met_size_env
disp_phal_size.e <- disparity_bones_env$disparity_phal_size_env



load("disparity_bones_group.RData")
str(disparity_bones_group)

disp_hum_shape.g <- disparity_bones_group$disparity_hum_shape_group
disp_rad_shape.g <- disparity_bones_group$disparity_rad_shape
disp_met_shape.g <- disparity_bones_group$disparity_met_shape_group
disp_phal_shape.g <- disparity_bones_group$ disparity_phal_shape_group

disp_hum_size.g <- disparity_bones_group$disparity_hum_size_group
disp_rad_size.g <- disparity_bones_group$disparity_rad_size_group
disp_met_size.g <- disparity_bones_group$disparity_met_size_group
disp_phal_size.g <- disparity_bones_group$disparity_phal_size_group





cols2<-setNames(c("#F4B06F","#F27489","seagreen", "skyblue", "gray88"),c("air","semi-aerial", "semi-aquatic", "aquatic", "terrestrial")) 
colsg <- setNames(c("#665665", "skyblue", "skyblue", "#F4B06F"),c("others","Cetacea","Pinnipedia","Chiroptera"))

cols.taxon <- setNames( c("#B3F94E","black","#BFBFBF",
                          "#FCDE06","#FFA300",
                          "#00FFF8","#008CFF","#005CFF","#025091","#0000FF","#4C1CF7","#8436FF",
                          "#838E54" ,"#04B72F","#215B15","#3CAF82",
                          "#FF00DB","#6B3030","#F99DF7","#BA93BC","#F4C1C1","#FF0000"),                           
                        c("Rodentia","Monotremata","Marsupialia",
                          "Cingulata","Pilosa",
                          "Eulipotyphla","Perissodactyla","Cetartiodactyla_Non_Cet","Cetacea","Carnivora","Pholidota","Chiroptera",
                          "Dermoptera","Scandentia","Primates","Lagomorpha",
                          "Sirenia","Proboscidea","Hyracoidea","Tubulidentata","Macroscelidea","Afrosoricida"))


cols.taxon.gray <- setNames( c("white",rep(c("#BFBFBF","black"), 10), "black"),
                        c("Monotremata","Marsupialia",
                          "Pilosa","Cingulata",
                          "Sirenia","Hyracoidea","Proboscidea","Tubulidentata","Macroscelidea","Afrosoricida",
                          "Dermoptera","Scandentia","Primates","Lagomorpha","Rodentia",
                          "Eulipotyphla","Chiroptera","Pholidota","Carnivora","Perissodactyla","Cetacea","Cetartiodactyla_Non_Cet"))

##==========================================================================================================##
##                                                                                                          ##
####  - Figure 1 -----                                                                                    ####
##                                                                                                          ##
##==========================================================================================================##



##---------------------------------------------------------##
## fan                                                     ##
##---------------------------------------------------------##


habitat <- dat1$Media
names(habitat) <- tree$tip.label
mapped_tree <- make.simmap(ladderize(tree), habitat, model="ER", nsim=1)
Clade <- dat1$Group 
names(Clade) <- dat1$Species


# plot the mapped ancestral states
plotSimmap(mapped_tree, cols2,pts=F,lwd=1.5, type="fan", fsize=0.25, ftype="off")
#add.simmap.legend(colors=cols2,vertical=TRUE)
# add a legend


#legend("topleft", c("air","semi-aquatic", "aquatic", "solid"),
#       fill=c("lightcoral","#3AC9BA", "#002953",  "gray"))
par("mar")

par(mar=c(4, 4,0.1,0.1))


obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
n<-Ntip(mapped_tree)
cols.taxon
par(lend=3)
for(i in 1:Ntip(mapped_tree)){
  cc<-if(obj$xx[i]>0) 14 else -14
  th<-atan(obj$yy[i]/obj$xx[i])
  segments(obj$xx[i],obj$yy[i],
           obj$xx[i]+cc*cos(th),
           obj$yy[i]+cc*sin(th),
           lwd=4,
           col=cols.taxon.gray[Clade[mapped_tree$tip.label[i]]])
  }




##==========================================================================================================##
##                                                                                                          ##
####  - Figure 2 -----                                                                                    ####
##                                                                                                          ##
##==========================================================================================================##



##---------------------------------------------------------##
## Phenogram                                               ##
##---------------------------------------------------------##




#phenogram
sc1 <- sc_pca1$V1
names(sc1) <- row.names(sc_pca1)
sc2 <- sc_pca2$V1
names(sc2) <- row.names(sc_pca2)


phenogram(mapped_tree, sc1,colors=cols2, spread.labels=TRUE, fsize=0.1)
phenogram(mapped_tree, sc2,colors=cols2, spread.labels=TRUE, fsize=0.1)





##---------------------------------------------------------##
#### phylogenetic pPC1 vs pPC2                           ####
##---------------------------------------------------------##



## - PCA without size - ##


pc1pc2_g2 <- ggplot(sc_pca2, aes(V1, V2, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.61), breaks=seq(-0.4, 0.6, by=0.2))


pc1pc2_m2 <- ggplot(sc_pca2, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.61), breaks=seq(-0.4, 0.6, by=0.2))



library(plotly)

ggplotly(pc1pc2_g2, tooltip=c("text", "x", "y", "group"))




## - PCA with size - ##


pc1pc2_g <- ggplot(sc_pca1, aes(V1*-1, V2, group=Group, text=rownames(sc_pca1))) + # multiplying scores *-1 to help visualization
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  xlim(-5,7.5)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")


pc1pc2_m <- ggplot(sc_pca1, aes(V1*-1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")


## Get legend

library(cowplot)

leg_g <- ggplot(sc_pca1, aes(V1, V3, group=Group)) +
  geom_point(aes(color=Group), size=2.5, alpha=0.5) +
  scale_color_manual(values=(cols.taxon))+
  theme_classic()+
  xlab("pPC1")+
  ylab("pPC3")+
  guides(col=guide_legend(ncol=1))+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))



leg_m <- ggplot(sc_pca1, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2.5, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_classic()+
  xlab("pPC1")+
  ylab("pPC2")+
  guides(col=guide_legend(ncol=1))+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)


lg <- get_legend(leg_g)
lm <- get_legend(leg_m)




### ------ Get plot ------ ###


p1 <- plot_grid(pc1pc2_m2, pc1pc2_g2, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p1 <- textGrob("pPC2 (14.3%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p1 <- textGrob("pPC1 (30.1%)",  gp=gpar(col="black", fontsize=10))
p1.p <- grid.arrange(arrangeGrob(p1, left = y.lab.p1, bottom = x.lab.p1))

?plot_grid

p2 <- plot_grid(pc1pc2_m, pc1pc2_g, ncol=2, labels=c("Media", "Taxa"),  align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p2 <- textGrob("pPC2 (2.7%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p2 <- textGrob("pPC1 (90.2%)",  gp=gpar(col="black", fontsize=10))
p2.p <- grid.arrange(arrangeGrob(p2, left = y.lab.p2, bottom = x.lab.p2))

p_leg <- plot_grid(lg, lm, ncol=1,rel_heights = c(3,1), align="V")


title.p1 <- ggdraw() + draw_label("A   shape",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

title.p2 <- ggdraw() +  draw_label("B   size",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p1.with.title <- plot_grid(title.p1, p1.p, ncol = 1,rel_heights = c(0.15, 1))
p2.with.title <- plot_grid(title.p2, p2.p, ncol = 1,rel_heights = c(0.15, 1))


ppca <- plot_grid(p1.with.title, p2.with.title, ncol=1)
plot_grid(ppca, p_leg, ncol=2, rel_widths = c(4,1))





  
  







##==========================================================================================================##
##                                                                                                          ##
####  - Figure 3 -----                                                                                    ####
##                                                                                                          ##
##==========================================================================================================##


# pca plot of bone **SHAPE** data 


library(ggpubr)
library(patchwork)



## ----------- Humerus ----------- ##


f2_hum <- ggplot(sc_pca_hum_gm, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_minimal()+
  xlab("shape pPC1 (36.04%)")+
  ylab("shape pPC2 (20.34%)")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits= c(-0.4, 0.42), breaks=seq(-0.4, 0.4, by=0.2))+
  scale_y_continuous(limits= c(-0.68, 0.4), breaks=seq(-0.5, 0.25, by=0.25))


x_f2_hum <- ggplot(sc_pca_hum_gm, aes(V1, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-0.4, 0.42))+
  theme_void() + 
  theme(legend.position = "none")


y_f2_hum <- ggplot(sc_pca_hum_gm, aes(V2, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-0.68, 0.4))+
  theme_void() + 
  theme(legend.position = "none")+
  coord_flip()

pf2_hum <- x_f2_hum + plot_spacer() + f2_hum + y_f2_hum + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))




## ----------- Radius ----------- ##


f2_rad <- ggplot(sc_pca_rad_gm, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_minimal()+
  xlab("shape pPC1 (34.06%)")+
  ylab("shape pPC2 (21.92%)")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits= c(-0.55, 0.5), breaks=seq(-0.4, 0.4, by=0.2)) +
  scale_y_continuous(limits= c(-0.25, 0.53), breaks=seq(-0.2, 0.4, by=0.2))


x_f2_rad <- ggplot(sc_pca_rad_gm, aes(V1, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-0.55, 0.5))+
  theme_void() + 
  theme(legend.position = "none")


y_f2_rad <- ggplot(sc_pca_rad_gm, aes(V2, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-0.25, 0.53))+
  theme_void() + 
  theme(legend.position = "none")+
  coord_flip()

pf2_rad <- x_f2_rad + plot_spacer() + f2_rad + y_f2_rad + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))




## ----------- Metacarpal ----------- ##


f2_met <- ggplot(sc_pca_met_gm, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_minimal()+
  xlab("shape pPC1 (48.19%)")+
  ylab("shape pPC2 (18.27%)")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits= c(-1.0, 1), breaks=seq(-1.0, 1.0, by=0.5))+
  scale_y_continuous(limits= c(-0.42, 0.41), breaks=seq(-0.4, 0.4, by=0.2))



x_f2_met <- ggplot(sc_pca_met_gm, aes(V1, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-1.0, 1.0))+
  theme_void() + 
  theme(legend.position = "none")


y_f2_met <- ggplot(sc_pca_met_gm, aes(V2, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  theme_void() + 
  xlim(c(-0.42, 0.41))+
  theme(legend.position = "none")+
  coord_flip()

pf2_met <- x_f2_met + plot_spacer() + f2_met + y_f2_met + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))



## ----------- Phalanx ----------- ##


f2_phal <- ggplot(sc_pca_phal_gm, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_minimal()+
  xlab("shape pPC1 (52.68%)")+
  ylab("shape pPC2 (19.46%)")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_x_continuous(limits= c(-0.78, 0.9), breaks=seq(-0.8, 0.8, by=0.4))+
  scale_y_continuous(limits= c(-0.4, 0.3), breaks=seq(-0.4, 0.2, by=0.2))


x_f2_phal <- ggplot(sc_pca_phal_gm, aes(V1, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-0.78, 0.9))+
  theme_void() + 
  theme(legend.position = "none")


y_f2_phal <- ggplot(sc_pca_phal_gm, aes(V2, fill = Media)) + 
  geom_density(alpha = 0.4, color = NA) + 
  scale_fill_manual(values=cols2) +
  xlim(c(-0.4, 0.3))+
  theme_void() + 
  theme(legend.position = "none")+
  coord_flip()

pf2_phal <- x_f2_phal + plot_spacer() + f2_phal + y_f2_phal + 
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))




plot_grid(pf2_hum, pf2_rad, 
          pf2_met, pf2_phal, labels=c("A  humerus", "B  radius", "C  metacarpal", "D  phalanx"))


sc_pca_met_gm

##==========================================================================================================##
##                                                                                                          ##
####  - Figure 4 -----                                                                                    ####
##                                                                                                          ##
##==========================================================================================================##



##---------------------------------------------------------##
## Disparity                                               ##
##---------------------------------------------------------##


## 1) whole limb

disp_env_shape <- ggplot(disparity_shape_env, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  ylim(c(0.0,0.35))+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")



disp_env_size <- ggplot(disparity_size_env, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 8), 
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.limb.e <- plot_grid(disp_env_shape, disp_env_size, ncol=2, align = "hv", scale=0.95)
y.lab.d.limb.e<- textGrob("whole limb disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.limb.env<- grid.arrange(arrangeGrob(d.limb.e , left = y.lab.d.limb.e))





## ----------- Humerus ----------- ##



dhum_shape.e <- ggplot(disp_hum_shape.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")



dhum_size.e <- ggplot(disp_hum_size.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.hum.e <- plot_grid(dhum_shape.e, dhum_size.e, ncol=2, align = "hv", scale=0.95)
y.lab.d.hum.e<- textGrob("humerus disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.hum.env<- grid.arrange(arrangeGrob(d.hum.e , left = y.lab.d.hum.e))







## ----------- Radius ----------- ##


drad_shape.e <- ggplot(disp_rad_shape.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")



drad_size.e <- ggplot(disp_rad_size.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.rad.e <- plot_grid(drad_shape.e, drad_size.e, ncol=2, align = "hv", scale=0.95)
y.lab.d.rad.e<- textGrob("radius disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.rad.env<- grid.arrange(arrangeGrob(d.rad.e , left = y.lab.d.rad.e))




## ----------- Metacarpal ----------- ##


dmet_shape.e <- ggplot(disp_met_shape.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")



dmet_size.e <- ggplot(disp_met_size.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.met.e <- plot_grid(dmet_shape.e, dmet_size.e, ncol=2, align = "hv", scale=0.95)
y.lab.d.met.e<- textGrob("metacarpal disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.met.env<- grid.arrange(arrangeGrob(d.met.e , left = y.lab.d.met.e))





## ----------- Phalanx ----------- ##

dphal_shape.e <- ggplot(disp_phal_shape.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 8, angle=45, hjust=0.5))



dphal_size.e <- ggplot(disp_phal_size.e, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=cols2)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=cols2)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 8, angle=45, hjust=0.5))




d.phal.e <- plot_grid(dphal_shape.e, dphal_size.e, ncol=2, align = "hv")
y.lab.d.phal.e<- textGrob("phalanx disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.phal.env<- grid.arrange(arrangeGrob(d.phal.e , left = y.lab.d.phal.e))




## ----------- all plots ----------- ##


plot_grid(p.d.limb.env, p.d.hum.env, p.d.rad.env, p.d.met.env, p.d.phal.env, ncol=1, labels=c("A","B","C","D","E"), align = "hv")


##==========================================================================================================##
##                                                                                                          ##
####  - Figure 5  -----                                                                                   ####
##                                                                                                          ##
##==========================================================================================================##

#---------------------------------------------------------##
## Disparity between groups                                ##
##---------------------------------------------------------##



## 1) shape groups

disp1g <- ggplot(disparity_shape_group, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.4, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=1, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("whole limb shape disparity
(size removed)")+
  labs(fill="Group")+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 8, angle=45, hjust=1))+
  scale_y_continuous(limits= c(-0.02,0.3), breaks=seq(0.0, 0.3, by=0.1))



## 2) size groups

disp2g <- ggplot(disparity_size_group, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.4, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=1, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("whole limb size disparity
(raw)")+
  labs(fill="Group")+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 8, angle=45, hjust=1))+
  scale_y_continuous(limits= c(-0.1,6.5), breaks=seq(0.0, 6.0, by=2),  labels = function(x) round(x,  digits=2))



plot_grid(disp1g, disp2g, nrow=1, labels=c("A","B"))






##==========================================================================================================##
##                                                                                                          ##
####  - Supporting Information                                                                            ####
##                                                                                                          ##
##==========================================================================================================##

## get legend

leg_g <- ggplot(sc_pca_hum_gm, aes(V1, V2, group=Group)) + 
  geom_point(aes(color=Group), size=1.5, alpha=0.35) +
  scale_color_manual(values=as.vector(cols.taxon))+
  theme_minimal()+
  guides(col=guide_legend(ncol=1))
  

leg_m <- ggplot(sc_pca_hum_gm, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=1.5, alpha=0.35) +
  scale_color_manual(values=cols2)+
  theme_minimal()+
  guides(col=guide_legend(ncol=1))+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)


lg <- get_legend(leg_g)
lm <- get_legend(leg_m)






##==========================================================================================================##
##                                                                                                          ##
####   Figure S1 -----                                                                                  ####
##                                                                                                          ##
##==========================================================================================================##


### Panel A: shape
# Panel A.1 : ppc1 x ppc2   ppc3 x pp4
# Panel A.2 : pc1 x pc2     pc3 x p4

### Panel B: size
# Panel B.1 : ppc1 x ppc2   ppc3 x pp4
# Panel B.2 : pc1 x pc2     pc3 x p4


phy.pca_limb_size <- df_plots$sc_pca_raw
phy.pca_limb_shape<- df_plots$sc_pca_shape
opca_limb_size <- df_plots$reg_pca_raw
opca_limb_shape <- df_plots$reg_pca_shape

ord.pca_limb_shape <- as.data.frame(opca_shape$x)
ord.pca_limb_shape$Group <- phy.pca_shape$Group
ord.pca_limb_shape$Media <- phy.pca_shape$Media

ord.pca_limb_size <- as.data.frame(opca_size$x)
ord.pca_limb_size$Group <- phy.pca_size$Group
ord.pca_limb_size$Media <- phy.pca_size$Media


##-----------------------------------------------------------------------------------------------------##
### Panel A: limb shape #### 
##-----------------------------------------------------------------------------------------------------##


##  Panel A phy PCA ------------------------------------------------

## Panel a1 L ppc1 x ppc2

phy.pc1pc2_m_limb_shape <- ggplot(phy.pca_limb_shape, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.61), breaks=seq(-0.4, 0.6, by=0.2))



phy.pc1pc2_g_limb_shape <- ggplot(phy.pca_limb_shape, aes(V1, V2, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.61), breaks=seq(-0.4, 0.6, by=0.2))


p.a1.L_phy_limb_shape_grid <- plot_grid(phy.pc1pc2_m_limb_shape, phy.pc1pc2_g_limb_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.L_phy_limb_shape <- textGrob("pPC2 (14.28%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.L_phy_limb_shape <- textGrob("pPC1 (30.12%)",  gp=gpar(col="black", fontsize=10))
p.a1.L_phy_limb_shape <- grid.arrange(arrangeGrob(p.a1.L_phy_limb_shape_grid, left = y.lab.p.a1.L_phy_limb_shape, bottom = x.lab.p.a1.L_phy_limb_shape))




## Panel a1 R ppc3 x ppc4


phy.pc3pc4_m_limb_shape <- ggplot(phy.pca_limb_shape, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.55, 0.35), breaks=seq(-0.4, 0.2, by=0.2))



phy.pc3pc4_g_limb_shape <- ggplot(phy.pca_limb_shape, aes(V3, V4, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.55, 0.35), breaks=seq(-0.4, 0.2, by=0.2))


p.a1_phy_limb_shape_R_grid <- plot_grid(phy.pc3pc4_m_limb_shape, phy.pc3pc4_g_limb_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.R_phy_limb_shape<- textGrob("pPC4 (6.64%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.R_phy_limb_shape<- textGrob("pPC3 (7.94%)",  gp=gpar(col="black", fontsize=10))
p.a1.R_phy_limb_shape<- grid.arrange(arrangeGrob(p.a1_phy_limb_shape_R_grid, left = y.lab.p.a1.R_phy_limb_shape, bottom = x.lab.p.a1.R_phy_limb_shape))



title.a1_phy_limb_shape <- ggdraw() + draw_label("     A.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a1_phy_limb_shape <- plot_grid(p.a1.L_phy_limb_shape, p.a1.R_phy_limb_shape, align="H", ncol=2)
a1.with.title_phy_limb_shape <- plot_grid(title.a1_phy_limb_shape,p.a1_phy_limb_shape , ncol = 1,rel_heights = c(0.15, 1))



## Panel A ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_limb_shape <- ggplot(ord.pca_limb_shape, aes(PC1*-1, PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.52, 0.83), breaks=seq(-0.4, 0.8, by=0.4))



ord.pc1pc2_g_limb_shape <- ggplot(ord.pca_limb_shape, aes(PC1*-1, PC2*-1, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.52, 0.83), breaks=seq(-0.4, 0.8, by=0.4))



p.a2.L_ord_limb_shape_grid <- plot_grid(ord.pc1pc2_m_limb_shape, ord.pc1pc2_g_limb_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.L_ord_limb_shape <- textGrob("PC2 (9.28%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.L_ord_limb_shape <- textGrob("PC1 (67.26%)",  gp=gpar(col="black", fontsize=10))
p.a2.L_ord_limb_shape <- grid.arrange(arrangeGrob(p.a2.L_ord_limb_shape_grid, left = y.lab.p.a2.L_ord_limb_shape, bottom = x.lab.p.a2.L_ord_limb_shape))




## Panel R pc3 x pc4


ord.pc3pc4_m_limb_shape <- ggplot(ord.pca_limb_shape, aes(x=PC3*-1, y=PC4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.6, 0.58), breaks=seq(-0.6, 0.3, by=0.3))



ord.pc3pc4_g_limb_shape <- ggplot(ord.pca_limb_shape, aes(x=PC3*-1, y=PC4, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.6, 0.58), breaks=seq(-0.6, 0.3, by=0.3))


p.a2_ord_limb_shape_R_grid <- plot_grid(ord.pc3pc4_m_limb_shape, ord.pc3pc4_g_limb_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.R_ord_limb_shape<- textGrob("PC4 (3.73%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.R_ord_limb_shape<- textGrob("PC3 (5.48%)",  gp=gpar(col="black", fontsize=10))
p.a2.R_ord_limb_shape<- grid.arrange(arrangeGrob(p.a2_ord_limb_shape_R_grid, left = y.lab.p.a2.R_ord_limb_shape, bottom = x.lab.p.a2.R_ord_limb_shape))



title.a2_ord_limb_shape <- ggdraw() + draw_label("     A.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a2_ord_limb_shape <- plot_grid(p.a2.L_ord_limb_shape, p.a2.R_ord_limb_shape, align="H", ncol=2)
a2.with.title_ord_limb_shape <- plot_grid(title.a2_ord_limb_shape,p.a2_ord_limb_shape , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot A ------ ###


title.a_limb <- ggdraw() + draw_label("A.  WHOLE LIMB SHAPE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.limb_shape <- plot_grid(a1.with.title_phy_limb_shape, a2.with.title_ord_limb_shape, ncol=1)

p.limb_shape_with.title <- plot_grid(title.a_limb, p.limb_shape, ncol = 1,rel_heights = c(0.15, 1))




##-----------------------------------------------------------------------------------------------------##
### Panel B: limb size #### 
##-----------------------------------------------------------------------------------------------------##

## Panel B phy PCA ------------------------------------------------

## Panel b1 L ppc1 x ppc2

phy.pc1pc2_m_limb_size <- ggplot(phy.pca_limb_size, aes(V1*-1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-1.3, 1.4), breaks=seq(-1.0, 1.0, by=0.5))



phy.pc1pc2_g_limb_size <- ggplot(phy.pca_limb_size, aes(V1*-1, V2, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-1.3, 1.4), breaks=seq(-1.0, 1.0, by=0.5))


p.b1.L_phy_limb_size_grid <- plot_grid(phy.pc1pc2_m_limb_size, phy.pc1pc2_g_limb_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.L_phy_limb_size <- textGrob("pPC2 (2.72%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.L_phy_limb_size <- textGrob("pPC1 (90.25%)",  gp=gpar(col="black", fontsize=10))
p.b1.L_phy_limb_size <- grid.arrange(arrangeGrob(p.b1.L_phy_limb_size_grid, left = y.lab.p.b1.L_phy_limb_size, bottom = x.lab.p.b1.L_phy_limb_size))




## Panel b1 R ppc3 x ppc4


phy.pc3pc4_m_limb_size <- ggplot(phy.pca_limb_size, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.45, 0.52), breaks=seq(-0.4, 0.4, by=0.2))+
  scale_x_continuous(limits= c(-0.50, 0.63), breaks=seq(-0.3, 0.6, by=0.3))



phy.pc3pc4_g_limb_size <- ggplot(phy.pca_limb_size, aes(V3, V4, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.45, 0.52), breaks=seq(-0.4, 0.4, by=0.2))+
  scale_x_continuous(limits= c(-0.50, 0.63), breaks=seq(-0.3, 0.6, by=0.3))


p.b1_phy_limb_size_R_grid <- plot_grid(phy.pc3pc4_m_limb_size, phy.pc3pc4_g_limb_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.R_phy_limb_size<- textGrob("pPC4 (0.75%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.R_phy_limb_size<- textGrob("pPC3 (1.36%)",  gp=gpar(col="black", fontsize=10))
p.b1.R_phy_limb_size<- grid.arrange(arrangeGrob(p.b1_phy_limb_size_R_grid, left = y.lab.p.b1.R_phy_limb_size, bottom = x.lab.p.b1.R_phy_limb_size))



title.b1_phy_limb_size <- ggdraw() + draw_label("     B.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b1_phy_limb_size <- plot_grid(p.b1.L_phy_limb_size, p.b1.R_phy_limb_size, align="H", ncol=2)
b1.with.title_phy_limb_size <- plot_grid(title.b1_phy_limb_size,p.b1_phy_limb_size , ncol = 1,rel_heights = c(0.15, 1))



## Panel B ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_limb_size <- ggplot(ord.pca_limb_size, aes(x=PC1*-1, y=PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-1.5, 1.2), breaks=seq(-1.5, 1.0, by=0.5),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_limb_size <- ggplot(ord.pca_limb_size, aes(PC1*-1, PC2*-1, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-1.5, 1.2), breaks=seq(-1.5, 1.0, by=0.5),  labels = function(x) round(x,  digits=1))



p.b2.L_ord_limb_size_grid <- plot_grid(ord.pc1pc2_m_limb_size, ord.pc1pc2_g_limb_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.L_ord_limb_size <- textGrob("PC2 (4.68%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.L_ord_limb_size <- textGrob("PC1 (93.22%)",  gp=gpar(col="black", fontsize=10))
p.b2.L_ord_limb_size <- grid.arrange(arrangeGrob(p.b2.L_ord_limb_size_grid, left = y.lab.p.b2.L_ord_limb_size, bottom = x.lab.p.b2.L_ord_limb_size))




## Panel R pc3 x pc4


ord.pc3pc4_m_limb_size <- ggplot(ord.pca_limb_size, aes(x=PC3*-1, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.7, 0.38), breaks=seq(-0.6, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_limb_size <- ggplot(ord.pca_limb_size, aes(x=PC3*-1, y=PC4*-1, group=Group, text=rownames(sc_pca2))) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.7, 0.38), breaks=seq(-0.6, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



p.b2_ord_limb_size_R_grid <- plot_grid(ord.pc3pc4_m_limb_size, ord.pc3pc4_g_limb_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.R_ord_limb_size<- textGrob("PC4 (0.38%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.R_ord_limb_size<- textGrob("PC3 (0.54%)",  gp=gpar(col="black", fontsize=10))
p.b2.R_ord_limb_size<- grid.arrange(arrangeGrob(p.b2_ord_limb_size_R_grid, left = y.lab.p.b2.R_ord_limb_size, bottom = x.lab.p.b2.R_ord_limb_size))



title.b2_ord_limb_size <- ggdraw() + draw_label("     B.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b2_ord_limb_size <- plot_grid(p.b2.L_ord_limb_size, p.b2.R_ord_limb_size, align="H", ncol=2)
b2.with.title_ord_limb_size <- plot_grid(title.b2_ord_limb_size,p.b2_ord_limb_size , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot B ------ ###


title.b_limb <- ggdraw() + draw_label("B.   WHOLE LIMB SIZE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.limb_size <- plot_grid(b1.with.title_phy_limb_size, b2.with.title_ord_limb_size, ncol=1)

p.limb_size_with.title <- plot_grid(title.b_limb, p.limb_size, ncol = 1,rel_heights = c(0.15, 1))



### ------ Get Figure S2.1 ------ ###

plot_grid(p.limb_shape_with.title, p.limb_size_with.title, ncol = 1, align="h")





##==========================================================================================================##
##                                                                                                          ##
####  - Figure S2 -----                                                                                 ####
##                                                                                                          ##
##==========================================================================================================##


### Panel A: shape
# Panel A.1 : ppc1 x ppc2   ppc3 x pp4
# Panel A.2 : pc1 x pc2     pc3 x p4

### Panel B: size
# Panel B.1 : ppc1 x ppc2   ppc3 x pp4
# Panel B.2 : pc1 x pc2     pc3 x p4
bones_pca

phy.pca_hum_size <-bones_pca$sc_pca_hum_raw
phy.pca_hum_shape<- bones_pca$sc_pca_hum_gm
opca_hum_size <- bones_opca$sc_opca_hum_size
opca_hum_shape <- bones_opca$sc_opca_hum_shape

ord.pca_hum_shape <- as.data.frame(opca_hum_shape)
ord.pca_hum_shape$Group <- phy.pca_hum_size$Group
ord.pca_hum_shape$Media <- phy.pca_hum_size$Media

ord.pca_hum_size <- as.data.frame(opca_hum_size)
ord.pca_hum_size$Group <- phy.pca_hum_size$Group
ord.pca_hum_size$Media <- phy.pca_hum_size$Media


##-----------------------------------------------------------------------------------------------------##
### Panel A: humerus shape #### 
##-----------------------------------------------------------------------------------------------------##


## Panel A phy PCA ------------------------------------------------

## Panel a1 L ppc1 x ppc2

phy.pc1pc2_m_hum_shape <- ggplot(phy.pca_hum_shape, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.7,0.4), breaks=seq(-0.6, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_hum_shape <- ggplot(phy.pca_hum_shape, aes(V1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.7,0.4), breaks=seq(-0.6, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))


p.a1.L_phy_hum_shape_grid <- plot_grid(phy.pc1pc2_m_hum_shape, phy.pc1pc2_g_hum_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.L_phy_hum_shape <- textGrob("pPC2 (20.34%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.L_phy_hum_shape <- textGrob("pPC1 (36.04%)",  gp=gpar(col="black", fontsize=10))
p.a1.L_phy_hum_shape <- grid.arrange(arrangeGrob(p.a1.L_phy_hum_shape_grid, left = y.lab.p.a1.L_phy_hum_shape, bottom = x.lab.p.a1.L_phy_hum_shape))




## Panel a1 R ppc3 x ppc4


phy.pc3pc4_m_hum_shape <- ggplot(phy.pca_hum_shape, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.3, 0.2), breaks=seq(-0.2, 0.2, by=0.2))



phy.pc3pc4_g_hum_shape <- ggplot(phy.pca_hum_shape, aes(V3, V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.3, 0.2), breaks=seq(-0.2, 0.2, by=0.2))


p.a1_phy_hum_shape_R_grid <- plot_grid(phy.pc3pc4_m_hum_shape, phy.pc3pc4_g_hum_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.R_phy_hum_shape<- textGrob("pPC4 (14.15%%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.R_phy_hum_shape<- textGrob("pPC3 (18.32%%)",  gp=gpar(col="black", fontsize=10))
p.a1.R_phy_hum_shape<- grid.arrange(arrangeGrob(p.a1_phy_hum_shape_R_grid, left = y.lab.p.a1.R_phy_hum_shape, bottom = x.lab.p.a1.R_phy_hum_shape))



title.a1_phy_hum_shape <- ggdraw() + draw_label("     A.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a1_phy_hum_shape <- plot_grid(p.a1.L_phy_hum_shape, p.a1.R_phy_hum_shape, align="H", ncol=2)
a1.with.title_phy_hum_shape <- plot_grid(title.a1_phy_hum_shape,p.a1_phy_hum_shape , ncol = 1,rel_heights = c(0.15, 1))



## Panel A ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_hum_shape <- ggplot(ord.pca_hum_shape, aes(PC1, PC2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.73,0.35), breaks=seq(-0.6, 0.2, by=0.2),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_hum_shape <- ggplot(ord.pca_hum_shape, aes(PC1, PC2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.73,0.35), breaks=seq(-0.6, 0.2, by=0.2),  labels = function(x) round(x,  digits=1))



p.a2.L_ord_hum_shape_grid <- plot_grid(ord.pc1pc2_m_hum_shape, ord.pc1pc2_g_hum_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.L_ord_hum_shape <- textGrob("PC2 (30.04%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.L_ord_hum_shape <- textGrob("PC1 (51.65%)",  gp=gpar(col="black", fontsize=10))
p.a2.L_ord_hum_shape <- grid.arrange(arrangeGrob(p.a2.L_ord_hum_shape_grid, left = y.lab.p.a2.L_ord_hum_shape, bottom = x.lab.p.a2.L_ord_hum_shape))




## Panel R pc3 x pc4


ord.pc3pc4_m_hum_shape <- ggplot(ord.pca_hum_shape, aes(x=PC3, y=PC4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.32, 0.22), breaks=seq(-0.3, 0.2, by=0.1), labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_hum_shape <- ggplot(ord.pca_hum_shape, aes(x=PC3, y=PC4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.32, 0.22), breaks=seq(-0.3, 0.2, by=0.1), labels = function(x) round(x,  digits=1))


p.a2_ord_hum_shape_R_grid <- plot_grid(ord.pc3pc4_m_hum_shape, ord.pc3pc4_g_hum_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.R_ord_hum_shape<- textGrob("PC4 (5.36%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.R_ord_hum_shape<- textGrob("PC3 (8.43%)",  gp=gpar(col="black", fontsize=10))
p.a2.R_ord_hum_shape<- grid.arrange(arrangeGrob(p.a2_ord_hum_shape_R_grid, left = y.lab.p.a2.R_ord_hum_shape, bottom = x.lab.p.a2.R_ord_hum_shape))



title.a2_ord_hum_shape <- ggdraw() + draw_label("     A.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a2_ord_hum_shape <- plot_grid(p.a2.L_ord_hum_shape, p.a2.R_ord_hum_shape, align="H", ncol=2)
a2.with.title_ord_hum_shape <- plot_grid(title.a2_ord_hum_shape,p.a2_ord_hum_shape , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot A ------ ###


title.a_hum <- ggdraw() + draw_label("A.  HUMERUS SHAPE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.hum_shape <- plot_grid(a1.with.title_phy_hum_shape, a2.with.title_ord_hum_shape, ncol=1)

p.hum_shape_with.title <- plot_grid(title.a_hum, p.hum_shape, ncol = 1,rel_heights = c(0.15, 1))




##-----------------------------------------------------------------------------------------------------##
### Panel B: humerus size #### 
##-----------------------------------------------------------------------------------------------------##

## Panel B phy PCA ------------------------------------------------

## Panel b1 L ppc1 x ppc2

phy.pc1pc2_m_hum_size <- ggplot(phy.pca_hum_size, aes(V1*-1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.45), breaks=seq(-0.5, 0.25, by=0.25))



phy.pc1pc2_g_hum_size <- ggplot(phy.pca_hum_size, aes(V1*-1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.5, 0.45), breaks=seq(-0.5, 0.25, by=0.25))


p.b1.L_phy_hum_size_grid <- plot_grid(phy.pc1pc2_m_hum_size, phy.pc1pc2_g_hum_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.L_phy_hum_size <- textGrob("pPC2 (1.82%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.L_phy_hum_size <- textGrob("pPC1 (94.99%)",  gp=gpar(col="black", fontsize=10))
p.b1.L_phy_hum_size <- grid.arrange(arrangeGrob(p.b1.L_phy_hum_size_grid, left = y.lab.p.b1.L_phy_hum_size, bottom = x.lab.p.b1.L_phy_hum_size))




## Panel b1 R ppc3 x ppc4


phy.pc3pc4_m_hum_size <- ggplot(phy.pca_hum_size, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.28, 0.17), breaks=seq(-0.3, 0.1, by=0.1), labels = function(x) round(x,  digits=1))


phy.pc3pc4_g_hum_size <- ggplot(phy.pca_hum_size, aes(V3, V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.28, 0.17), breaks=seq(-0.3, 0.1, by=0.1), labels = function(x) round(x,  digits=1))


p.b1_phy_hum_size_R_grid <- plot_grid(phy.pc3pc4_m_hum_size, phy.pc3pc4_g_hum_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.R_phy_hum_size<- textGrob("pPC4 (1.04%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.R_phy_hum_size<- textGrob("pPC3 (1.30%)",  gp=gpar(col="black", fontsize=10))
p.b1.R_phy_hum_size<- grid.arrange(arrangeGrob(p.b1_phy_hum_size_R_grid, left = y.lab.p.b1.R_phy_hum_size, bottom = x.lab.p.b1.R_phy_hum_size))



title.b1_phy_hum_size <- ggdraw() + draw_label("     B.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b1_phy_hum_size <- plot_grid(p.b1.L_phy_hum_size, p.b1.R_phy_hum_size, align="H", ncol=2)
b1.with.title_phy_hum_size <- plot_grid(title.b1_phy_hum_size,p.b1_phy_hum_size , ncol = 1,rel_heights = c(0.15, 1))



## Panel B ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_hum_size <- ggplot(ord.pca_hum_size, aes(x=PC1*-1, y=PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.55,0.42), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_hum_size <- ggplot(ord.pca_hum_size, aes(PC1*-1, PC2*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.55,0.42), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



p.b2.L_ord_hum_size_grid <- plot_grid(ord.pc1pc2_m_hum_size, ord.pc1pc2_g_hum_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.L_ord_hum_size <- textGrob("PC2 (1.62%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.L_ord_hum_size <- textGrob("PC1 (97.39%)",  gp=gpar(col="black", fontsize=10))
p.b2.L_ord_hum_size <- grid.arrange(arrangeGrob(p.b2.L_ord_hum_size_grid, left = y.lab.p.b2.L_ord_hum_size, bottom = x.lab.p.b2.L_ord_hum_size))




## Panel R pc3 x pc4


ord.pc3pc4_m_hum_size <- ggplot(ord.pca_hum_size, aes(x=PC3*-1, y=PC4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.24, 0.17), breaks=seq(-0.2, 0.1, by=0.1),  labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_hum_size <- ggplot(ord.pca_hum_size, aes(x=PC3*-1, y=PC4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.24, 0.17), breaks=seq(-0.2, 0.1, by=0.1),  labels = function(x) round(x,  digits=1))




p.b2_ord_hum_size_R_grid <- plot_grid(ord.pc3pc4_m_hum_size, ord.pc3pc4_g_hum_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.R_ord_hum_size<- textGrob("PC4 (0.37%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.R_ord_hum_size<- textGrob("PC3 (0.42%)",  gp=gpar(col="black", fontsize=10))
p.b2.R_ord_hum_size<- grid.arrange(arrangeGrob(p.b2_ord_hum_size_R_grid, left = y.lab.p.b2.R_ord_hum_size, bottom = x.lab.p.b2.R_ord_hum_size))



title.b2_ord_hum_size <- ggdraw() + draw_label("     B.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b2_ord_hum_size <- plot_grid(p.b2.L_ord_hum_size, p.b2.R_ord_hum_size, align="H", ncol=2)
b2.with.title_ord_hum_size <- plot_grid(title.b2_ord_hum_size,p.b2_ord_hum_size , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot B ------ ###


title.b_hum <- ggdraw() + draw_label("B.   HUMERUS SIZE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.hum_size <- plot_grid(b1.with.title_phy_hum_size, b2.with.title_ord_hum_size, ncol=1)

p.hum_size_with.title <- plot_grid(title.b_hum, p.hum_size, ncol = 1,rel_heights = c(0.15, 1))



### ------ Get Figure S2.1 ------ ###

plot_grid(p.hum_shape_with.title, p.hum_size_with.title, ncol = 1, align="h")







##==========================================================================================================##
##                                                                                                          ##
####  - Figure S3 -----                                                                                 ####
##                                                                                                          ##
##==========================================================================================================##


### Panel A: shape
# Panel A.1 : ppc1 x ppc2   ppc3 x pp4
# Panel A.2 : pc1 x pc2     pc3 x p4

### Panel B: size
# Panel B.1 : ppc1 x ppc2   ppc3 x pp4
# Panel B.2 : pc1 x pc2     pc3 x p4

phy.pca_rad_size <-bones_pca$sc_pca_rad_raw
phy.pca_rad_shape<- bones_pca$sc_pca_rad_gm
opca_rad_size <- bones_opca$sc_opca_rad_size
opca_rad_shape <- bones_opca$sc_opca_rad_shape

ord.pca_rad_shape <- as.data.frame(opca_rad_shape)
ord.pca_rad_shape$Group <- phy.pca_rad_size$Group
ord.pca_rad_shape$Media <- phy.pca_rad_size$Media

ord.pca_rad_size <- as.data.frame(opca_rad_size)
ord.pca_rad_size$Group <- phy.pca_rad_size$Group
ord.pca_rad_size$Media <- phy.pca_rad_size$Media


##-----------------------------------------------------------------------------------------------------##
### Panel A: radius shape #### 
##-----------------------------------------------------------------------------------------------------##


## Panel A phy PCA ------------------------------------------------

## Panel a1 L ppc1 x ppc2

phy.pc1pc2_m_rad_shape <- ggplot(phy.pca_rad_shape, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.24,0.53), breaks=seq(-0.2, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_rad_shape <- ggplot(phy.pca_rad_shape, aes(V1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.24,0.53), breaks=seq(-0.2, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))


p.a1.L_phy_rad_shape_grid <- plot_grid(phy.pc1pc2_m_rad_shape, phy.pc1pc2_g_rad_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.L_phy_rad_shape <- textGrob("pPC2 (21.92%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.L_phy_rad_shape <- textGrob("pPC1 (34.06%)",  gp=gpar(col="black", fontsize=10))
p.a1.L_phy_rad_shape <- grid.arrange(arrangeGrob(p.a1.L_phy_rad_shape_grid, left = y.lab.p.a1.L_phy_rad_shape, bottom = x.lab.p.a1.L_phy_rad_shape))




## Panel a1 R ppc3 x ppc4


phy.pc3pc4_m_rad_shape <- ggplot(phy.pca_rad_shape, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.25, 0.17), breaks=seq(-0.2, 0.1, by=0.1))



phy.pc3pc4_g_rad_shape <- ggplot(phy.pca_rad_shape, aes(V3, V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.25, 0.17), breaks=seq(-0.2, 0.1, by=0.1))


p.a1_phy_rad_shape_R_grid <- plot_grid(phy.pc3pc4_m_rad_shape, phy.pc3pc4_g_rad_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.R_phy_rad_shape<- textGrob("pPC4 (12.04%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.R_phy_rad_shape<- textGrob("pPC3 (20.49%)",  gp=gpar(col="black", fontsize=10))
p.a1.R_phy_rad_shape<- grid.arrange(arrangeGrob(p.a1_phy_rad_shape_R_grid, left = y.lab.p.a1.R_phy_rad_shape, bottom = x.lab.p.a1.R_phy_rad_shape))



title.a1_phy_rad_shape <- ggdraw() + draw_label("     A.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a1_phy_rad_shape <- plot_grid(p.a1.L_phy_rad_shape, p.a1.R_phy_rad_shape, align="H", ncol=2)
a1.with.title_phy_rad_shape <- plot_grid(title.a1_phy_rad_shape,p.a1_phy_rad_shape , ncol = 1,rel_heights = c(0.15, 1))



## Panel A ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_rad_shape <- ggplot(ord.pca_rad_shape, aes(PC1, PC2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.4,0.4), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_rad_shape <- ggplot(ord.pca_rad_shape, aes(PC1, PC2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.4,0.4), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



p.a2.L_ord_rad_shape_grid <- plot_grid(ord.pc1pc2_m_rad_shape, ord.pc1pc2_g_rad_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.L_ord_rad_shape <- textGrob("PC2 (18.90%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.L_ord_rad_shape <- textGrob("PC1 (56.89%)",  gp=gpar(col="black", fontsize=10))
p.a2.L_ord_rad_shape <- grid.arrange(arrangeGrob(p.a2.L_ord_rad_shape_grid, left = y.lab.p.a2.L_ord_rad_shape, bottom = x.lab.p.a2.L_ord_rad_shape))




## Panel R pc3 x pc4


ord.pc3pc4_m_rad_shape <- ggplot(ord.pca_rad_shape, aes(x=PC3*-1, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.2, 0.21), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_rad_shape <- ggplot(ord.pca_rad_shape, aes(x=PC3*-1, y=PC4*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.2, 0.21), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))


p.a2_ord_rad_shape_R_grid <- plot_grid(ord.pc3pc4_m_rad_shape, ord.pc3pc4_g_rad_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.R_ord_rad_shape<- textGrob("PC4 (7.99%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.R_ord_rad_shape<- textGrob("PC3 (9.80%)",  gp=gpar(col="black", fontsize=10))
p.a2.R_ord_rad_shape<- grid.arrange(arrangeGrob(p.a2_ord_rad_shape_R_grid, left = y.lab.p.a2.R_ord_rad_shape, bottom = x.lab.p.a2.R_ord_rad_shape))



title.a2_ord_rad_shape <- ggdraw() + draw_label("     A.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a2_ord_rad_shape <- plot_grid(p.a2.L_ord_rad_shape, p.a2.R_ord_rad_shape, align="H", ncol=2)
a2.with.title_ord_rad_shape <- plot_grid(title.a2_ord_rad_shape,p.a2_ord_rad_shape , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot A ------ ###


title.a_rad <- ggdraw() + draw_label("A.  RADIUS SHAPE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.rad_shape <- plot_grid(a1.with.title_phy_rad_shape, a2.with.title_ord_rad_shape, ncol=1)

p.rad_shape_with.title <- plot_grid(title.a_rad, p.rad_shape, ncol = 1,rel_heights = c(0.15, 1))




##-----------------------------------------------------------------------------------------------------##
### Panel B: radius size #### 
##-----------------------------------------------------------------------------------------------------##

## Panel B phy PCA ------------------------------------------------

## Panel b1 L ppc1 x ppc2

phy.pc1pc2_m_rad_size <- ggplot(phy.pca_rad_size, aes(V1*-1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.6, 0.4), breaks=seq(-0.6, 0.4, by=0.2), labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_rad_size <- ggplot(phy.pca_rad_size, aes(V1*-1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.6, 0.4), breaks=seq(-0.6, 0.4, by=0.2), labels = function(x) round(x,  digits=1))


p.b1.L_phy_rad_size_grid <- plot_grid(phy.pc1pc2_m_rad_size, phy.pc1pc2_g_rad_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.L_phy_rad_size <- textGrob("pPC2 (2.46%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.L_phy_rad_size <- textGrob("pPC1 (93.73%)",  gp=gpar(col="black", fontsize=10))
p.b1.L_phy_rad_size <- grid.arrange(arrangeGrob(p.b1.L_phy_rad_size_grid, left = y.lab.p.b1.L_phy_rad_size, bottom = x.lab.p.b1.L_phy_rad_size))




## Panel b1 R ppc3 x ppc4


phy.pc3pc4_m_rad_size <- ggplot(phy.pca_rad_size, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.2, 0.23), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))


phy.pc3pc4_g_rad_size <- ggplot(phy.pca_rad_size, aes(V3, V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.2, 0.23), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))



p.b1_phy_rad_size_R_grid <- plot_grid(phy.pc3pc4_m_rad_size, phy.pc3pc4_g_rad_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.R_phy_rad_size<- textGrob("pPC4 (1.23%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.R_phy_rad_size<- textGrob("pPC3 (1.73%)",  gp=gpar(col="black", fontsize=10))
p.b1.R_phy_rad_size<- grid.arrange(arrangeGrob(p.b1_phy_rad_size_R_grid, left = y.lab.p.b1.R_phy_rad_size, bottom = x.lab.p.b1.R_phy_rad_size))



title.b1_phy_rad_size <- ggdraw() + draw_label("     B.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b1_phy_rad_size <- plot_grid(p.b1.L_phy_rad_size, p.b1.R_phy_rad_size, align="H", ncol=2)
b1.with.title_phy_rad_size <- plot_grid(title.b1_phy_rad_size,p.b1_phy_rad_size , ncol = 1,rel_heights = c(0.15, 1))



## Panel B ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_rad_size <- ggplot(ord.pca_rad_size, aes(x=PC1*-1, y=PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.65,0.48), breaks=seq(-0.6, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_rad_size <- ggplot(ord.pca_rad_size, aes(PC1*-1, PC2*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.65,0.48), breaks=seq(-0.6, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



p.b2.L_ord_rad_size_grid <- plot_grid(ord.pc1pc2_m_rad_size, ord.pc1pc2_g_rad_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.L_ord_rad_size <- textGrob("PC2 (1.74%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.L_ord_rad_size <- textGrob("PC1 (97.38%)",  gp=gpar(col="black", fontsize=10))
p.b2.L_ord_rad_size <- grid.arrange(arrangeGrob(p.b2.L_ord_rad_size_grid, left = y.lab.p.b2.L_ord_rad_size, bottom = x.lab.p.b2.L_ord_rad_size))




## Panel R pc3 x pc4


ord.pc3pc4_m_rad_size <- ggplot(ord.pca_rad_size, aes(x=PC3*-1, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.24, 0.2), breaks=seq(-0.2, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_rad_size <- ggplot(ord.pca_rad_size, aes(x=PC3*-1, y=PC4*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.24, 0.2), breaks=seq(-0.2, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))




p.b2_ord_rad_size_R_grid <- plot_grid(ord.pc3pc4_m_rad_size, ord.pc3pc4_g_rad_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.R_ord_rad_size<- textGrob("PC4 (0.28%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.R_ord_rad_size<- textGrob("PC3 (0.35%)",  gp=gpar(col="black", fontsize=10))
p.b2.R_ord_rad_size<- grid.arrange(arrangeGrob(p.b2_ord_rad_size_R_grid, left = y.lab.p.b2.R_ord_rad_size, bottom = x.lab.p.b2.R_ord_rad_size))



title.b2_ord_rad_size <- ggdraw() + draw_label("     B.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b2_ord_rad_size <- plot_grid(p.b2.L_ord_rad_size, p.b2.R_ord_rad_size, align="H", ncol=2)
b2.with.title_ord_rad_size <- plot_grid(title.b2_ord_rad_size,p.b2_ord_rad_size , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot B ------ ###


title.b_rad <- ggdraw() + draw_label("B.   RADIUS SIZE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.rad_size <- plot_grid(b1.with.title_phy_rad_size, b2.with.title_ord_rad_size, ncol=1)

p.rad_size_with.title <- plot_grid(title.b_rad, p.rad_size, ncol = 1,rel_heights = c(0.15, 1))



### ------ Get Figure S2.1 ------ ###

plot_grid(p.rad_shape_with.title, p.rad_size_with.title, ncol = 1, align="h")





##==========================================================================================================##
##                                                                                                          ##
####  - Figure S4 -----                                                                                 ####
##                                                                                                          ##
##==========================================================================================================##


### Panel A: shape
# Panel A.1 : ppc1 x ppc2   ppc3 x pp4
# Panel A.2 : pc1 x pc2     pc3 x p4

### Panel B: size
# Panel B.1 : ppc1 x ppc2   ppc3 x pp4
# Panel B.2 : pc1 x pc2     pc3 x p4

phy.pca_met_size <-bones_pca$sc_pca_met_raw
phy.pca_met_shape<- bones_pca$sc_pca_met_gm
opca_met_size <- bones_opca$sc_opca_met_size
opca_met_shape <- bones_opca$sc_opca_met_shape

ord.pca_met_shape <- as.data.frame(opca_met_shape)
ord.pca_met_shape$Group <- phy.pca_met_size$Group
ord.pca_met_shape$Media <- phy.pca_met_size$Media

ord.pca_met_size <- as.data.frame(opca_met_size)
ord.pca_met_size$Group <- phy.pca_met_size$Group
ord.pca_met_size$Media <- phy.pca_met_size$Media


##-----------------------------------------------------------------------------------------------------##
### Panel A: metacarpal shape #### 
##-----------------------------------------------------------------------------------------------------##


## Panel A phy PCA ------------------------------------------------

## Panel a1 L ppc1 x ppc2

phy.pc1pc2_m_met_shape <- ggplot(phy.pca_met_shape, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.42,0.4), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_met_shape <- ggplot(phy.pca_met_shape, aes(V1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.42,0.4), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))


p.a1.L_phy_met_shape_grid <- plot_grid(phy.pc1pc2_m_met_shape, phy.pc1pc2_g_met_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.L_phy_met_shape <- textGrob("pPC2 (18.27%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.L_phy_met_shape <- textGrob("pPC1 (48.19%)",  gp=gpar(col="black", fontsize=10))
p.a1.L_phy_met_shape <- grid.arrange(arrangeGrob(p.a1.L_phy_met_shape_grid, left = y.lab.p.a1.L_phy_met_shape, bottom = x.lab.p.a1.L_phy_met_shape))




## Panel a1 R ppc3 x ppc4


phy.pc3pc4_m_met_shape <- ggplot(phy.pca_met_shape, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.3, 0.25), breaks=seq(-0.3, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))



phy.pc3pc4_g_met_shape <- ggplot(phy.pca_met_shape, aes(x=V3, y=V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.3, 0.25), breaks=seq(-0.3, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))


p.a1_phy_met_shape_R_grid <- plot_grid(phy.pc3pc4_m_met_shape, phy.pc3pc4_g_met_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.R_phy_met_shape<- textGrob("pPC4 (11.62%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.R_phy_met_shape<- textGrob("pPC3 (14.42%)",  gp=gpar(col="black", fontsize=10))
p.a1.R_phy_met_shape<- grid.arrange(arrangeGrob(p.a1_phy_met_shape_R_grid, left = y.lab.p.a1.R_phy_met_shape, bottom = x.lab.p.a1.R_phy_met_shape))



title.a1_phy_met_shape <- ggdraw() + draw_label("     A.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a1_phy_met_shape <- plot_grid(p.a1.L_phy_met_shape, p.a1.R_phy_met_shape, align="H", ncol=2)
a1.with.title_phy_met_shape <- plot_grid(title.a1_phy_met_shape,p.a1_phy_met_shape , ncol = 1,rel_heights = c(0.15, 1))



## Panel A ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_met_shape <- ggplot(ord.pca_met_shape, aes(PC1, PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.49,0.47), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_met_shape <- ggplot(ord.pca_met_shape, aes(PC1, PC2*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.49,0.47), breaks=seq(-0.4, 0.4, by=0.2),  labels = function(x) round(x,  digits=1))



p.a2.L_ord_met_shape_grid <- plot_grid(ord.pc1pc2_m_met_shape, ord.pc1pc2_g_met_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.L_ord_met_shape <- textGrob("PC2 (9.15%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.L_ord_met_shape <- textGrob("PC1 (82.13%)",  gp=gpar(col="black", fontsize=10))
p.a2.L_ord_met_shape <- grid.arrange(arrangeGrob(p.a2.L_ord_met_shape_grid, left = y.lab.p.a2.L_ord_met_shape, bottom = x.lab.p.a2.L_ord_met_shape))




## Panel R pc3 x pc4


ord.pc3pc4_m_met_shape <- ggplot(ord.pca_met_shape, aes(x=PC3, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.22, 0.2), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_met_shape <- ggplot(ord.pca_met_shape, aes(x=PC3, y=PC4*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.22, 0.2), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))


p.a2_ord_met_shape_R_grid <- plot_grid(ord.pc3pc4_m_met_shape, ord.pc3pc4_g_met_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.R_ord_met_shape<- textGrob("PC4 (2.63%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.R_ord_met_shape<- textGrob("PC3 (4.67%)",  gp=gpar(col="black", fontsize=10))
p.a2.R_ord_met_shape<- grid.arrange(arrangeGrob(p.a2_ord_met_shape_R_grid, left = y.lab.p.a2.R_ord_met_shape, bottom = x.lab.p.a2.R_ord_met_shape))



title.a2_ord_met_shape <- ggdraw() + draw_label("     A.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a2_ord_met_shape <- plot_grid(p.a2.L_ord_met_shape, p.a2.R_ord_met_shape, align="H", ncol=2)
a2.with.title_ord_met_shape <- plot_grid(title.a2_ord_met_shape,p.a2_ord_met_shape , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot A ------ ###


title.a_met <- ggdraw() + draw_label("A.  METACARPAL SHAPE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.met_shape <- plot_grid(a1.with.title_phy_met_shape, a2.with.title_ord_met_shape, ncol=1)

p.met_shape_with.title <- plot_grid(title.a_met, p.met_shape, ncol = 1,rel_heights = c(0.15, 1))




##-----------------------------------------------------------------------------------------------------##
### Panel B: metacarpal size #### 
##-----------------------------------------------------------------------------------------------------##

## Panel B phy PCA ------------------------------------------------

## Panel b1 L ppc1 x ppc2

phy.pc1pc2_m_met_size <- ggplot(phy.pca_met_size, aes(V1*-1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.95, 0.95), breaks=seq(-0.8, 0.8, by=0.4), labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_met_size <- ggplot(phy.pca_met_size, aes(V1*-1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.95, 0.95), breaks=seq(-0.8, 0.8, by=0.4), labels = function(x) round(x,  digits=1))


p.b1.L_phy_met_size_grid <- plot_grid(phy.pc1pc2_m_met_size, phy.pc1pc2_g_met_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.L_phy_met_size <- textGrob("pPC2 (5.24%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.L_phy_met_size <- textGrob("pPC1 (90.37%)",  gp=gpar(col="black", fontsize=10))
p.b1.L_phy_met_size <- grid.arrange(arrangeGrob(p.b1.L_phy_met_size_grid, left = y.lab.p.b1.L_phy_met_size, bottom = x.lab.p.b1.L_phy_met_size))




## Panel b1 R ppc3 x ppc4


phy.pc3pc4_m_met_size <- ggplot(phy.pca_met_size, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.3, 0.28), breaks=seq(-0.3, 0.3, by=0.1), labels = function(x) round(x,  digits=1))


phy.pc3pc4_g_met_size <- ggplot(phy.pca_met_size, aes(V3, V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.3, 0.28), breaks=seq(-0.3, 0.3, by=0.1), labels = function(x) round(x,  digits=1))




p.b1_phy_met_size_R_grid <- plot_grid(phy.pc3pc4_m_met_size, phy.pc3pc4_g_met_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.R_phy_met_size<- textGrob("pPC4 (1.56%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.R_phy_met_size<- textGrob("pPC3 (1.88%)",  gp=gpar(col="black", fontsize=10))
p.b1.R_phy_met_size<- grid.arrange(arrangeGrob(p.b1_phy_met_size_R_grid, left = y.lab.p.b1.R_phy_met_size, bottom = x.lab.p.b1.R_phy_met_size))



title.b1_phy_met_size <- ggdraw() + draw_label("     B.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b1_phy_met_size <- plot_grid(p.b1.L_phy_met_size, p.b1.R_phy_met_size, align="H", ncol=2)
b1.with.title_phy_met_size <- plot_grid(title.b1_phy_met_size,p.b1_phy_met_size , ncol = 1,rel_heights = c(0.15, 1))



## Panel B ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_met_size <- ggplot(ord.pca_met_size, aes(x=PC1, y=PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-1.1,0.8), breaks=seq(-1.0, 0.5, by=0.5),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_met_size <- ggplot(ord.pca_met_size, aes(PC1, PC2*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-1.1,0.8), breaks=seq(-1.0, 0.5, by=0.5),  labels = function(x) round(x,  digits=1))


p.b2.L_ord_met_size_grid <- plot_grid(ord.pc1pc2_m_met_size, ord.pc1pc2_g_met_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.L_ord_met_size <- textGrob("PC2 (7.00%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.L_ord_met_size <- textGrob("PC1 (92.16%)",  gp=gpar(col="black", fontsize=10))
p.b2.L_ord_met_size <- grid.arrange(arrangeGrob(p.b2.L_ord_met_size_grid, left = y.lab.p.b2.L_ord_met_size, bottom = x.lab.p.b2.L_ord_met_size))




## Panel R pc3 x pc4


ord.pc3pc4_m_met_size <- ggplot(ord.pca_met_size, aes(x=PC3*-1, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.22, 0.23), breaks=seq(-0.2, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_met_size <- ggplot(ord.pca_met_size, aes(x=PC3*-1, y=PC4*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
scale_y_continuous(limits= c(-0.22, 0.23), breaks=seq(-0.2, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))


p.b2_ord_met_size_R_grid <- plot_grid(ord.pc3pc4_m_met_size, ord.pc3pc4_g_met_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.R_ord_met_size<- textGrob("PC4 (0.29%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.R_ord_met_size<- textGrob("PC3 (0.42%)",  gp=gpar(col="black", fontsize=10))
p.b2.R_ord_met_size<- grid.arrange(arrangeGrob(p.b2_ord_met_size_R_grid, left = y.lab.p.b2.R_ord_met_size, bottom = x.lab.p.b2.R_ord_met_size))



title.b2_ord_met_size <- ggdraw() + draw_label("     B.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b2_ord_met_size <- plot_grid(p.b2.L_ord_met_size, p.b2.R_ord_met_size, align="H", ncol=2)
b2.with.title_ord_met_size <- plot_grid(title.b2_ord_met_size,p.b2_ord_met_size , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot B ------ ###


title.b_met <- ggdraw() + draw_label("B.   METACARPAL SIZE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.met_size <- plot_grid(b1.with.title_phy_met_size, b2.with.title_ord_met_size, ncol=1)

p.met_size_with.title <- plot_grid(title.b_met, p.met_size, ncol = 1,rel_heights = c(0.15, 1))



### ------ Get Figure S2.3 ------ ###

plot_grid(p.met_shape_with.title, p.met_size_with.title, ncol = 1, align="h")

##==========================================================================================================##
##                                                                                                          ##
####  - Figure S5 -----                                                                                 ####
##                                                                                                          ##
##==========================================================================================================##


### Panel A: shape
# Panel A.1 : ppc1 x ppc2   ppc3 x pp4
# Panel A.2 : pc1 x pc2     pc3 x p4

### Panel B: size
# Panel B.1 : ppc1 x ppc2   ppc3 x pp4
# Panel B.2 : pc1 x pc2     pc3 x p4

phy.pca_phal_size <-bones_pca$sc_pca_phal_raw
phy.pca_phal_shape<- bones_pca$sc_pca_phal_gm
opca_phal_size <- bones_opca$sc_opca_phal_size
opca_phal_shape <- bones_opca$sc_opca_phal_shape

ord.pca_phal_shape <- as.data.frame(opca_phal_shape)
ord.pca_phal_shape$Group <- phy.pca_phal_size$Group
ord.pca_phal_shape$Media <- phy.pca_phal_size$Media

ord.pca_phal_size <- as.data.frame(opca_phal_size)
ord.pca_phal_size$Group <- phy.pca_phal_size$Group
ord.pca_phal_size$Media <- phy.pca_phal_size$Media


##-----------------------------------------------------------------------------------------------------##
### Panel A: phalanx shape #### 
##-----------------------------------------------------------------------------------------------------##


## Panel A phy PCA ------------------------------------------------

## Panel a1 L ppc1 x ppc2

phy.pc1pc2_m_phal_shape <- ggplot(phy.pca_phal_shape, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.35,0.29), breaks=seq(-0.2, 0.2, by=0.2),  labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_phal_shape <- ggplot(phy.pca_phal_shape, aes(V1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.35,0.29), breaks=seq(-0.2, 0.2, by=0.2),  labels = function(x) round(x,  digits=1))


p.a1.L_phy_phal_shape_grid <- plot_grid(phy.pc1pc2_m_phal_shape, phy.pc1pc2_g_phal_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.L_phy_phal_shape <- textGrob("pPC2 (19.46%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.L_phy_phal_shape <- textGrob("pPC1 (52.68%)",  gp=gpar(col="black", fontsize=10))
p.a1.L_phy_phal_shape <- grid.arrange(arrangeGrob(p.a1.L_phy_phal_shape_grid, left = y.lab.p.a1.L_phy_phal_shape, bottom = x.lab.p.a1.L_phy_phal_shape))




## Panel a1 R ppc3 x ppc4


phy.pc3pc4_m_phal_shape <- ggplot(phy.pca_phal_shape, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.23, 0.20), breaks=seq(-0.2, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))



phy.pc3pc4_g_phal_shape <- ggplot(phy.pca_phal_shape, aes(x=V3, y=V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.23, 0.20), breaks=seq(-0.2, 0.2, by=0.1),  labels = function(x) round(x,  digits=1))


p.a1_phy_phal_shape_R_grid <- plot_grid(phy.pc3pc4_m_phal_shape, phy.pc3pc4_g_phal_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a1.R_phy_phal_shape<- textGrob("pPC4 (7.08%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a1.R_phy_phal_shape<- textGrob("pPC3 (15.49%)",  gp=gpar(col="black", fontsize=10))
p.a1.R_phy_phal_shape<- grid.arrange(arrangeGrob(p.a1_phy_phal_shape_R_grid, left = y.lab.p.a1.R_phy_phal_shape, bottom = x.lab.p.a1.R_phy_phal_shape))



title.a1_phy_phal_shape <- ggdraw() + draw_label("     A.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a1_phy_phal_shape <- plot_grid(p.a1.L_phy_phal_shape, p.a1.R_phy_phal_shape, align="H", ncol=2)
a1.with.title_phy_phal_shape <- plot_grid(title.a1_phy_phal_shape,p.a1_phy_phal_shape , ncol = 1,rel_heights = c(0.15, 1))



## Panel A ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_phal_shape <- ggplot(ord.pca_phal_shape, aes(PC1, PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.47,0.3), breaks=seq(-0.4, 0.2, by=0.2),  labels = function(x) round(x,  digits=1))


ord.pc1pc2_g_phal_shape <- ggplot(ord.pca_phal_shape, aes(PC1, PC2*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.47,0.3), breaks=seq(-0.4, 0.2, by=0.2),  labels = function(x) round(x,  digits=1))



p.a2.L_ord_phal_shape_grid <- plot_grid(ord.pc1pc2_m_phal_shape, ord.pc1pc2_g_phal_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.L_ord_phal_shape <- textGrob("PC2 (7.71%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.L_ord_phal_shape <- textGrob("PC1 (82.14%)",  gp=gpar(col="black", fontsize=10))
p.a2.L_ord_phal_shape <- grid.arrange(arrangeGrob(p.a2.L_ord_phal_shape_grid, left = y.lab.p.a2.L_ord_phal_shape, bottom = x.lab.p.a2.L_ord_phal_shape))




## Panel R pc3 x pc4


ord.pc3pc4_m_phal_shape <- ggplot(ord.pca_phal_shape, aes(x=PC3*-1, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.22, 0.14), breaks=seq(-0.2, 0.1, by=0.1), labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_phal_shape <- ggplot(ord.pca_phal_shape, aes(x=PC3*-1, y=PC4*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.22, 0.14), breaks=seq(-0.2, 0.1, by=0.1), labels = function(x) round(x,  digits=1))


p.a2_ord_phal_shape_R_grid <- plot_grid(ord.pc3pc4_m_phal_shape, ord.pc3pc4_g_phal_shape, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.a2.R_ord_phal_shape<- textGrob("PC4 (1.93%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.a2.R_ord_phal_shape<- textGrob("PC3 (6.49%)",  gp=gpar(col="black", fontsize=10))
p.a2.R_ord_phal_shape<- grid.arrange(arrangeGrob(p.a2_ord_phal_shape_R_grid, left = y.lab.p.a2.R_ord_phal_shape, bottom = x.lab.p.a2.R_ord_phal_shape))



title.a2_ord_phal_shape <- ggdraw() + draw_label("     A.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.a2_ord_phal_shape <- plot_grid(p.a2.L_ord_phal_shape, p.a2.R_ord_phal_shape, align="H", ncol=2)
a2.with.title_ord_phal_shape <- plot_grid(title.a2_ord_phal_shape,p.a2_ord_phal_shape , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot A ------ ###


title.a_phal <- ggdraw() + draw_label("A.  PHALANX SHAPE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.phal_shape <- plot_grid(a1.with.title_phy_phal_shape, a2.with.title_ord_phal_shape, ncol=1)

p.phal_shape_with.title <- plot_grid(title.a_phal, p.phal_shape, ncol = 1,rel_heights = c(0.15, 1))




##-----------------------------------------------------------------------------------------------------##
### Panel B: phalanx size #### 
##-----------------------------------------------------------------------------------------------------##

## Panel B phy PCA ------------------------------------------------

## Panel b1 L ppc1 x ppc2

phy.pc1pc2_m_phal_size <- ggplot(phy.pca_phal_size, aes(V1, V2, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.75, 0.95), breaks=seq(-0.8, 0.8, by=0.4), labels = function(x) round(x,  digits=1))



phy.pc1pc2_g_phal_size <- ggplot(phy.pca_phal_size, aes(V1, V2, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.75, 0.95), breaks=seq(-0.8, 0.8, by=0.4), labels = function(x) round(x,  digits=1))


p.b1.L_phy_phal_size_grid <- plot_grid(phy.pc1pc2_m_phal_size, phy.pc1pc2_g_phal_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.L_phy_phal_size <- textGrob("pPC2 (5.33%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.L_phy_phal_size <- textGrob("pPC1 (91.16%)",  gp=gpar(col="black", fontsize=10))
p.b1.L_phy_phal_size <- grid.arrange(arrangeGrob(p.b1.L_phy_phal_size_grid, left = y.lab.p.b1.L_phy_phal_size, bottom = x.lab.p.b1.L_phy_phal_size))




## Panel b1 R ppc3 x ppc4


phy.pc3pc4_m_phal_size <- ggplot(phy.pca_phal_size, aes(V3, V4, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.23, 0.19), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))


phy.pc3pc4_g_phal_size <- ggplot(phy.pca_phal_size, aes(V3, V4, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.23, 0.19), breaks=seq(-0.2, 0.2, by=0.1), labels = function(x) round(x,  digits=1))




p.b1_phy_phal_size_R_grid <- plot_grid(phy.pc3pc4_m_phal_size, phy.pc3pc4_g_phal_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b1.R_phy_phal_size<- textGrob("pPC4 (0.94%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b1.R_phy_phal_size<- textGrob("pPC3 (1.94%)",  gp=gpar(col="black", fontsize=10))
p.b1.R_phy_phal_size<- grid.arrange(arrangeGrob(p.b1_phy_phal_size_R_grid, left = y.lab.p.b1.R_phy_phal_size, bottom = x.lab.p.b1.R_phy_phal_size))



title.b1_phy_phal_size <- ggdraw() + draw_label("     B.1  phylogenetic PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b1_phy_phal_size <- plot_grid(p.b1.L_phy_phal_size, p.b1.R_phy_phal_size, align="H", ncol=2)
b1.with.title_phy_phal_size <- plot_grid(title.b1_phy_phal_size,p.b1_phy_phal_size , ncol = 1,rel_heights = c(0.15, 1))



## Panel B ordinary PCA -----------------------------------------------

## Panel L pc1 x pc2

ord.pc1pc2_m_phal_size <- ggplot(ord.pca_phal_size, aes(x=PC1, y=PC2*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.85,0.8), breaks=seq(-0.5, 0.5, by=0.5),  labels = function(x) round(x,  digits=1))



ord.pc1pc2_g_phal_size <- ggplot(ord.pca_phal_size, aes(PC1, PC2*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.85,0.8), breaks=seq(-0.5, 0.5, by=0.5),  labels = function(x) round(x,  digits=1))



p.b2.L_ord_phal_size_grid <- plot_grid(ord.pc1pc2_m_phal_size, ord.pc1pc2_g_phal_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.L_ord_phal_size <- textGrob("PC2 (5.70%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.L_ord_phal_size <- textGrob("PC1 (93.52%)",  gp=gpar(col="black", fontsize=10))
p.b2.L_ord_phal_size <- grid.arrange(arrangeGrob(p.b2.L_ord_phal_size_grid, left = y.lab.p.b2.L_ord_phal_size, bottom = x.lab.p.b2.L_ord_phal_size))




## Panel R pc3 x pc4


ord.pc3pc4_m_phal_size <- ggplot(ord.pca_phal_size, aes(x=PC3*-1, y=PC4*-1, color=Media))+
  geom_point(aes(color=Media), size=2, alpha=0.8) +
  scale_color_manual(values=cols2)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(plot.title = element_text(size=12, face="bold"), 
        axis.text.y = element_text(size = 10), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  stat_ellipse(geom = "polygon",
               aes(fill = Media),
               alpha = 0.2) + scale_fill_manual(values=cols2)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.26, 0.15), breaks=seq(-0.2, 0.1, by=0.1),  labels = function(x) round(x,  digits=1))



ord.pc3pc4_g_phal_size <- ggplot(ord.pca_phal_size, aes(x=PC3*-1, y=PC4*-1, group=Group)) +
  geom_point(aes(color=Group), size=2, alpha=0.5) +
  scale_color_manual(values=cols.taxon)+
  theme_void()+
  xlab("")+
  ylab("")+
  #guides(col=guide_legend(ncol=1))+
  theme(legend.position="none")+
  theme(axis.text.y = element_blank(), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 10), axis.title.x = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=10),
        panel.grid.major.x = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1),
        panel.grid.major.y = element_line(color = "lightgray",
                                          size = 0.5,
                                          linetype = 1))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  scale_y_continuous(limits= c(-0.26, 0.15), breaks=seq(-0.2, 0.1, by=0.1),  labels = function(x) round(x,  digits=1))



p.b2_ord_phal_size_R_grid <- plot_grid(ord.pc3pc4_m_phal_size, ord.pc3pc4_g_phal_size, ncol=2, labels=c("Media", "Taxa"), align = "hv", label_x = 0.05, label_y = 0.98,label_size = 11, scale=0.95)
y.lab.p.b2.R_ord_phal_size<- textGrob("PC4 (0.13%)",  gp=gpar(col="black", fontsize=10), rot=90)
x.lab.p.b2.R_ord_phal_size<- textGrob("PC3 (0.54%)",  gp=gpar(col="black", fontsize=10))
p.b2.R_ord_phal_size<- grid.arrange(arrangeGrob(p.b2_ord_phal_size_R_grid, left = y.lab.p.b2.R_ord_phal_size, bottom = x.lab.p.b2.R_ord_phal_size))



title.b2_ord_phal_size <- ggdraw() + draw_label("     B.2  ordinary PCA",fontface = 'bold', x = 0, hjust = 0, size=11) +
  theme(plot.margin = margin(0, 0, 0, 7))


p.b2_ord_phal_size <- plot_grid(p.b2.L_ord_phal_size, p.b2.R_ord_phal_size, align="H", ncol=2)
b2.with.title_ord_phal_size <- plot_grid(title.b2_ord_phal_size,p.b2_ord_phal_size , ncol = 1,rel_heights = c(0.15, 1))



### ------ Get plot B ------ ###


title.b_phal <- ggdraw() + draw_label("B.   PHALANX SIZE",fontface = 'bold', x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

p.phal_size <- plot_grid(b1.with.title_phy_phal_size, b2.with.title_ord_phal_size, ncol=1)

p.phal_size_with.title <- plot_grid(title.b_phal, p.phal_size, ncol = 1,rel_heights = c(0.15, 1))



### ------ Get Figure S2.5 ------ ###

plot_grid(p.phal_shape_with.title, p.phal_size_with.title, ncol = 1, align="h")



##==========================================================================================================##
##                                                                                                          ##
####  - Figure Sx -----                                                                                 ####
##                                                                                                          ##
##==========================================================================================================##



## ----------- Humerus ----------- ##


dhum_shape.g <- ggplot(disp_hum_shape.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




dhum_size.g <- ggplot(disp_hum_size.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.hum.g <- plot_grid(dhum_shape.g, dhum_size.g, ncol=2, align = "hv", scale=0.95)
y.lab.d.hum.g<- textGrob("humerus
disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.hum.group<- grid.arrange(arrangeGrob(d.hum.g , left = y.lab.d.hum.g))







## ----------- Radius ----------- ##


drad_shape.g <- ggplot(disp_rad_shape.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")



drad_size.g <- ggplot(disp_rad_size.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.rad.g <- plot_grid(drad_shape.g, drad_size.g, ncol=2, align = "hv", scale=0.95)
y.lab.d.rad.g<- textGrob("radius
disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.rad.group<- grid.arrange(arrangeGrob(d.rad.g , left = y.lab.d.rad.g))




## ----------- Metacarpal ----------- ##


dmet_shape.g <- ggplot(disp_met_shape.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")+
 scale_y_continuous(limits= c(0.0,0.095), breaks=seq(0.0, 0.08, by=0.02))


dmet_size.g <- ggplot(disp_met_size.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  easy_remove_axes("x")




d.met.g <- plot_grid(dmet_shape.g, dmet_size.g, ncol=2, align = "hv", scale=0.95)
y.lab.d.met.g<- textGrob("metacarpal
disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.met.group<- grid.arrange(arrangeGrob(d.met.g , left = y.lab.d.met.g))





## ----------- Phalanx ----------- ##

dphal_shape.g <- ggplot(disp_phal_shape.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 8, angle=45, hjust=0.5))



dphal_size.g <- ggplot(disp_phal_size.g, aes(x=Var2, y=value, fill=Var2, color=Var2)) +
  theme_classic()+
  geom_boxplot(width=0.5, color="gray50", outlier.colour = NA,  alpha=0.2) +
  scale_fill_manual(values=colsg)+
  geom_point(size=0.9, shape=21, position=position_jitter(width=0.13, height=0.01), alpha=0.3) +
  scale_fill_manual(values=colsg)+
  scale_color_manual(values=colsg)+
  xlab(NULL)+
  ylab("")+
  labs(fill="Media")+
  theme(legend.position="none")+
  #ylim(c())+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 8), axis.title.y = element_text(size=12),
        axis.text.x = element_text(size = 8, angle=45, hjust=0.5))




d.phal.g <- plot_grid(dphal_shape.g, dphal_size.g, ncol=2, align = "hv")
y.lab.d.phal.g<- textGrob("phalanx
disparity", gp=gpar(col="black", fontsize=10), rot=90)
p.d.phal.group<- grid.arrange(arrangeGrob(d.phal.g , left = y.lab.d.phal.g))




## ----------- all plots ----------- ##


plot_grid(p.d.hum.group, p.d.rad.group, p.d.met.group, p.d.phal.group, ncol=1, labels=c("A","B","C","D"), align = "hv")




##==========================================================================================================##
##                                                                                                          ##
####  - Figure S6 -----                                                                                   ####
##                                                                                                          ##
##==========================================================================================================##

library(dispRity)

# plot name 

## environment

load("env_disparity_shape_wholelimb_v2.RData")
plot.dispRity(env_disparity_shape, rarefaction=TRUE)


load("env_disparity_size_wholelimb_v2.RData")
plot.dispRity(env_disparity_size, rarefaction=TRUE)


load("env_disparity_hum_shape.RData")
plot.dispRity(env_disparity_hum_shape, rarefaction=TRUE)

load("env_disparity_hum_size.RData")
plot.dispRity(env_disparity_hum_size, rarefaction=TRUE)


load("env_disparity_rad_shape.RData")
plot.dispRity(env_disparity_rad_shape, rarefaction=TRUE)

load("env_disparity_rad_size.RData")
plot.dispRity(env_disparity_rad_size, rarefaction=TRUE)


load("env_disparity_met_shape.RData")
plot.dispRity(env_disparity_met_shape, rarefaction=TRUE)

load("env_disparity_met_size.RData")
plot.dispRity(env_disparity_met_size, rarefaction=TRUE)

load("env_disparity_phal_shape.RData")
plot.dispRity(env_disparity_phal_shape, rarefaction=TRUE)

load("env_disparity_phal_size.RData")
plot.dispRity(env_disparity_phal_size, rarefaction=TRUE)


## group
rar_limb_shape_gr

load("group_disparity_shape_wholelimb_v2.RData")
plot.dispRity(group_disparity_shape, rarefaction=TRUE)

load("group_disparity_size_wholelimb_v2.RData")
plot.dispRity(group_disparity_size, rarefaction=TRUE)

load("group_disparity_hum_shape.RData")
plot.dispRity(group_disparity_hum_shape, rarefaction=TRUE)

#load("group_disparity_hum_size.RData")
#plot.dispRity(group_disparity_hum_size, rarefaction=TRUE)


load("group_disparity_rad_shape.RData")
plot.dispRity(group_disparity_rad_shape, rarefaction=TRUE)

load("group_disparity_rad_size.RData")
plot.dispRity(group_disparity_rad_size, rarefaction=TRUE)


load("group_disparity_met_shape.RData")
plot.dispRity(group_disparity_met_shape, rarefaction=TRUE)

load("group_disparity_met_size.RData")
plot.dispRity(group_disparity_met_size, rarefaction=TRUE)

load("group_disparity_phal_shape.RData")
plot.dispRity(group_disparity_phal_shape, rarefaction=TRUE)

load("group_disparity_phal_size.RData")
plot.dispRity(group_disparity_phal_size, rarefaction=TRUE)

