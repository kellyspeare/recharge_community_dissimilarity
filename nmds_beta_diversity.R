library(cowplot)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(tidyr)
library(sjPlot)
library(DHARMa)
library(car)
library(effects)

# color palette: "#5637C8","#C83761",  "#A9C837", "#37C89E"
# "#EE964B", "#F4D35E","#28AFB0","#19647E"
#  "#D1495B","#EDAE49",  "#00798C", "#003D5B"
# "#EE6055","#FFD97D",  "#17BEBB", "#AAF683"
#"#F8766D", "#7CAE00","#00BFC4","#C77CFF"

#  scale_color_manual(values=c("#e41b76", "#d01a8e", "#bd1aa4", "#a919ba", "#9619d1","#8318e7",
# "#d1f00f","#b2e716","#92de1c","#73d523","#53cc29","#32c230", "#109641",
# "#e1551e","#de751b","#db9418", "#867d79","#a8a19e"))+

cover<-read.csv("data/disturbed_plots_percent_cover.csv", header=TRUE, stringsAsFactors = TRUE)
#reordering factor levels
cover$Herb_Trt<-factor(cover$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))

#subsetting cover for T0
cover0<-subset(cover, cover$Timepoint=="T0")
cover3<-subset(cover, cover$Timepoint=="T3")
cover6<-subset(cover, cover$Timepoint=="T6")
cover9<-subset(cover, cover$Timepoint=="T9")
cover12<-subset(cover, cover$Timepoint=="T12")

# NOTE! when you make your species matrix dataframes make sure that you don't include Sum_Macroalgae or Sum_Coral in the species matrix
# these are sums of other columns in the species matrix
# the code below puts these columns in the 'meta' dataframes

# load the functional group key 
fun_groups<-read.csv("data/functional_group_key.csv", header = T, stringsAsFactors = T)
fun_groups$functional_group<-factor(fun_groups$functional_group, c("branching corals","encrusting corals","other corals","hydrocorals","fungids","soft corals",
                                                                   "calcareous algae","filamentous algae","foliose algae","macrophytes","microalgae","other macroalgae",
                                                                   "turf, bare space, encrusting algae", "anemones","sponges","tridacna","rubble","sand"))

#load timepoints key
timepoints<-read.csv("Time_point_key.csv", header = T, stringsAsFactors = T)
timepoints$Date<-as.Date(timepoints$Date, "%m/%d/%Y")

startdate <- as.Date("2018-07-10","%Y-%m-%d") # setting the start date of the experiment
timepoints$Days_since_start  <- as.numeric(difftime(timepoints$Date, startdate, units="days"))

# labeller for plots
herb_names <- c(
  `Open` = "High",
  `3X3` = "Medium",
  `2X2` = "Low",
  `1X1` = "Very Low"
)

# -------------------------------------------------------------------#
# NMDS Plots --------------------------------------------------------
# -------------------------------------------------------------------#

## NMDS T0 ----------------------------------------------------------

# Make a matrix of any species functional group that is present

sp_matrix0<-as.matrix(cover0[,c(13:63)])
meta0<-cover0[,c(1:12)]

# NMDS 
nmds0<-metaMDS(sp_matrix0, distance="bray", trymax=1000, k=2, autotransform =FALSE) 
stressplot(nmds0)
print(nmds0) #stress is 0.2  

nmds0_df<-data.frame(MDS1=nmds0$points[,1], MDS2=nmds0$points[,2])
nmds0_df_meta<-cbind(meta0, nmds0_df)
vec.sp0<-envfit(nmds0$points, sp_matrix0, perm=1000, choices=c(1,2))
spp.scrs0<-as.data.frame(scores(vec.sp0, display = "vectors"))
spp.scrs0$species<-rownames(spp.scrs0)

#renaming factor levels
spp.scrs0$species<-as.factor(spp.scrs0$species)
spp.scrs0$species<-factor(spp.scrs0$species)
#joining with functional group names
spp.scrs0<-left_join(spp.scrs0, fun_groups, by="species")

# nmds with color for herb treatment
nmds_t0_plot<-ggplot(nmds0_df_meta)+ 
  geom_point(aes(x=MDS1, y=MDS2, color=Herb_Trt, shape=Nutrient_Trt), size=2.5, stroke = 1)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                     name = "Herbivory & Predation", labels = c("High","Medium", "Low", "Very low"))+
  stat_ellipse(type = "t", mapping=aes(x=MDS1, y=MDS2, color=Herb_Trt), linetype = 1)+
  coord_fixed()+
  xlim(c(-1.24,1.24))+ ylim(c(-1.24,1.24))+
  theme_bw()+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)
#ggsave("figures/nmds_t0_plot.pdf", width=6, height=6, units="in")

# Plot with legend
nmds_t0_plot_legend<-ggplot(nmds0_df_meta)+ 
  geom_point(aes(x=MDS1, y=MDS2, color=Herb_Trt, shape=Nutrient_Trt), stroke = 1)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                     name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low"))+
  stat_ellipse(type = "t", mapping=aes(x=MDS1, y=MDS2, color=Herb_Trt), linetype = 1)+
  coord_fixed()+
  xlim(c(-1.24,1.24))+ ylim(c(-1.24,1.24))+
  theme_bw()+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  #theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)

nmds_legend <- get_legend(nmds_t0_plot_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))

nmds_t0_vectors<-ggplot(nmds0_df_meta, aes(x=MDS1, y=MDS2))+ 
  geom_segment(aes(x=0, y=0, xend=MDS1, yend=MDS2, color=basic_group), data=spp.scrs0, arrow=arrow(length = unit(0.0, "cm")), 
               linewidth=0.5)+
  #geom_text(data=spp.scrs0,aes(x=MDS1,y=MDS2,label=functional_group),size=3)+
  scale_shape_manual(values=c(1,2))+
  scale_color_manual(values=c("#06c914", "#a8a19e", "#964b00", "#FFC53F"))+ # algae, colonizable, corals, other
  #stat_ellipse(type = "t", mapping=aes(x=MDS1, y=MDS2, color=Herb_Trt))+
  xlim(c(-1.24,1.24))+ ylim(c(-1.24,1.24))+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))

nmds_t0_vectors_legend<-ggplot(nmds0_df_meta, aes(x=MDS1, y=MDS2))+ 
  geom_segment(aes(x=0, y=0, xend=MDS1, yend=MDS2, color=basic_group), data=spp.scrs0, arrow=arrow(length = unit(0.0, "cm")), 
               linewidth=0.5)+
  #geom_text(data=spp.scrs0,aes(x=MDS1,y=MDS2,label=functional_group),size=3)+
  scale_shape_manual(values=c(1,2))+
  scale_color_manual(values=c("#06c914", "#a8a19e", "#964b00", "#FFC53F"), name="Benthic groups")+
  #stat_ellipse(type = "t", mapping=aes(x=MDS1, y=MDS2, color=Herb_Trt))+
  xlim(c(-1.24,1.24))+ ylim(c(-1.24,1.24))+
  coord_fixed()+
  theme_bw()+
  #theme(legend.position="none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))

vector_legend <- get_legend(nmds_t0_vectors_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))


## --------------------------- NMDS T12 -------------------------------

#Make a matrix of any species functional group that is present

sp_matrix12<-as.matrix(cover12[,c(13:63)])
meta12<-cover12[,c(1:12)]

# NMDS 
nmds12<-metaMDS(sp_matrix12, distance="bray", trymax=1000, k=2, autotransform =FALSE) 
stressplot(nmds12)
print(nmds12) #stress is 0.09  

nmds12_df<-data.frame(MDS1=nmds12$points[,1], MDS2=nmds12$points[,2])
nmds12_df_meta<-cbind(meta12, nmds12_df)
vec.sp12<-envfit(nmds12$points, sp_matrix12, perm=1000, choices=c(1,2))
spp.scrs12<-as.data.frame(scores(vec.sp12, display = "vectors"))
spp.scrs12$species<-rownames(spp.scrs12)

#renaming factor levels
spp.scrs12$species<-as.factor(spp.scrs12$species)
spp.scrs12$species<-factor(spp.scrs12$species)
spp.scrs12<-left_join(spp.scrs12, fun_groups, by="species")

nmds_t12_plot<-ggplot(nmds12_df_meta)+ 
  geom_point(aes(x=MDS1, y=MDS2, color=Herb_Trt, shape=Nutrient_Trt), size=2.5, stroke = 1)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                     name = "Herbivory & Predation", labels = c("High","Medium", "Low", "Very low"))+
  stat_ellipse(type = "t", mapping=aes(x=MDS1, y=MDS2, color=Herb_Trt), linetype = 1)+
  coord_fixed()+
  xlim(c(-1.24,1.24))+ ylim(c(-1.24,1.24))+
  theme_bw()+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)
#ggsave("figures/nmds_t12_plot.pdf", width=6, height=6, units="in")

nmds_t12_vectors<-ggplot(nmds12_df_meta, aes(x=MDS1, y=MDS2))+ 
  geom_segment(aes(x=0, y=0, xend=MDS1, yend=MDS2, color=basic_group), data=spp.scrs12, arrow=arrow(length = unit(0.0, "cm")), 
               linewidth=0.5)+
  #geom_text(data=spp.scrs0,aes(x=MDS1,y=MDS2,label=functional_group),size=3)+
  scale_shape_manual(values=c(1,2))+
  scale_color_manual(values=c("#06c914", "#a8a19e", "#964b00", "#FFC53F"), name="Benthic groups")+
  #stat_ellipse(type = "t", mapping=aes(x=MDS1, y=MDS2, color=Herb_Trt))+
  xlim(c(-1.24,1.24))+ ylim(c(-1.24,1.24))+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))
#ggsave("figures/nmds_t12_vector_plot.pdf", width=6, height=6, units="in")


# -------------------------------------------------------------------#
#  Beta Dispersion---------------------------------------------------#
# -------------------------------------------------------------------#

# data summary function that makes mean points for plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-(sd(x)/sqrt(length(x)))
  ymax <- m+(sd(x)/sqrt(length(x)))
  return(c(y=m,ymin=ymin,ymax=ymax))
}

## Beta Dispersion T0 - herb trt, 4 groups  -----------------------------

sp_matrix0<-as.matrix(cover0[,c(13:63)])
meta0<-cover0[,c(1:12)]

bray_dist0<-vegdist(sp_matrix0, method = "bray")
bray_mod0<-betadisper(bray_dist0, meta0$Herb_Trt, type = "centroid")
anova(bray_mod0) # NS, P= 0.2721
plot(bray_mod0)
boxplot(bray_mod0)
TukeyHSD(bray_mod0)

betadisp0<-data.frame(group=bray_mod0$group, distances=bray_mod0$distances)
betadisp0<-cbind(betadisp0, meta0)

t0_dist_plot<-ggplot(betadisp0, aes(x=Herb_Trt, y=distances, color=Herb_Trt, shape=Nutrient_Trt))+
  geom_point(position=position_jitter(w = 0.3, h = 0.01), aes(shape=Nutrient_Trt), size=2.5, stroke=1)+
  stat_summary(fun.data=data_summary, color="black", size=1, shape=16)+
  scale_color_manual(values=c( "#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Herbivory", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2))+
  scale_x_discrete(name ="Consumer pressure", labels=c("High","Medium", "Low", "Very low"))+
  ylim(c(0.03,0.4))+
  coord_fixed()+
  theme_bw()+
  ylab("Distance to centroid")+
  labs(color="Herbivory & Predation")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)

## Beta Dispersion T0 - herb trt, 4 groups  -----------------------------

sp_matrix12<-as.matrix(cover12[,c(13:63)])
meta12<-cover12[,c(1:12)]

bray_dist12<-vegdist(sp_matrix12, method = "bray")
bray_mod12<-betadisper(bray_dist12, meta12$Herb_Trt, type = "centroid")
anova(bray_mod12)
plot(bray_mod12)
boxplot(bray_mod12)
TukeyHSD(bray_mod12)

betadisp12<-data.frame(group=bray_mod12$group, distances=bray_mod12$distances)
betadisp12<-cbind(betadisp12, meta12)

t12_dist_plot<-ggplot(betadisp12, aes(x=Herb_Trt, y=distances, color=Herb_Trt, shape=Nutrient_Trt))+
  geom_point(position=position_jitter(w = 0.3, h = 0.01), aes(shape=Nutrient_Trt), size=2.5, stroke=1)+
  stat_summary(fun.data=data_summary, color="black", size=1, shape=16)+
  scale_color_manual(values=c( "#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Herbivory", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2))+
  scale_x_discrete(name ="Consumer pressure", labels=c("High","Medium", "Low", "Very low"))+
  ylim(c(0.03,0.4))+
  coord_fixed()+
  theme_bw()+
  ylab("Distance to centroid")+
  labs(color="Herbivory & Predation")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(legend.position = "none")+
  theme(aspect.ratio = 4/4)

# Note to self, presenting the results with the above method. but results w/below method are quantitatively the same.
# ------------- Beta Dispersion T0 - herb and nutrients as one factor, 8 groups ------
# ### other ways to look at the T12 data that consider both treatments in factorial design
# 
# meta0_twofac<-meta0
# meta0_twofac$Herb_Nuts_Trt<-paste(meta0_twofac$Herb_Trt, meta0_twofac$Nutrient_Trt, sep = "_")
# 
# sp_matrix0_twofac<-as.matrix(cover0[,c(13:63)])
# bray_disT0_twofac<-vegdist(sp_matrix0_twofac, method = "bray")
# bray_mod0_twofac<-betadisper(bray_disT0_twofac, meta0_twofac$Herb_Nuts_Trt, type = "centroid")
# 
# anova(bray_mod0_twofac)
# plot(bray_mod0_twofac)
# boxplot(bray_mod0_twofac)
# TukeyHSD(bray_mod0_twofac)
# 
# betadisp0_twofac_df<-data.frame(group=bray_mod0_twofac$group, distances=bray_mod0_twofac$distances)
# 
# betadisp0_twofac_df <- betadisp0_twofac_df %>% separate(group, c("Herb_Trt", "Nutrient_Trt"), "_")
# betadisp0_twofac_df$Herb_Trt<-factor(betadisp0_twofac_df$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))
# betadisp0_twofac_df$Timepoint<-"T0"
# 
## Beta Dispersion T12 - herb and nutrients as one factor, 8 groups ---------------
# 
# meta12_twofac<-meta12
# meta12_twofac$Herb_Nuts_Trt<-paste(meta12_twofac$Herb_Trt, meta12_twofac$Nutrient_Trt, sep = "_")
# 
# sp_matrix12_twofac<-as.matrix(cover12[,c(13:63)])
# bray_dist12_twofac<-vegdist(sp_matrix12_twofac, method = "bray")
# bray_mod12_twofac<-betadisper(bray_dist12_twofac, meta12_twofac$Herb_Nuts_Trt, type = "centroid")
# 
# anova(bray_mod12_twofac)
# plot(bray_mod12_twofac)
# boxplot(bray_mod12_twofac)
# TukeyHSD(bray_mod12_twofac)
# 
# betadisp12_twofac_df<-data.frame(group=bray_mod12_twofac$group, distances=bray_mod12_twofac$distances)
# 
# betadisp12_twofac_df <- betadisp12_twofac_df %>% separate(group, c("Herb_Trt", "Nutrient_Trt"), "_")
# betadisp12_twofac_df$Herb_Trt<-factor(betadisp12_twofac_df$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))
# betadisp12_twofac_df$Timepoint<-"T12"
# 
# t0_dist_plot<-ggplot(betadisp0_twofac_df, aes(x=Herb_Trt, y=distances, color=Herb_Trt, shape=Nutrient_Trt))+
#   geom_point(position=position_jitter(w = 0.3, h = 0.01), aes(shape=Nutrient_Trt), size=2.5, stroke=1)+
#   stat_summary(fun.data=data_summary, color="black", size=1, shape=16)+
#   scale_color_manual(values=c( "#5637C8","#C83761",  "#A9C837", "#37C89E"),
#                      name = "Herbivory", labels = c("High","Medium", "Low", "Very low" ))+
#   #scale_shape_manual(values=c(1,2))+
#   scale_x_discrete(name ="Consumer pressure", labels=c("High","Medium", "Low", "Very low"))+
#   ylim(c(0,0.35))+
#   coord_fixed()+
#   theme_bw()+
#   ylab("Distance to centroid")+
#   labs(color="Herbivory & Predation")+
#   theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
#   theme(legend.position = "none")+
#   theme(aspect.ratio = 4/4)
# #ggsave("figures/T0_dist_plot.pdf", width=4, height=4, units="in")
# 
# t12_dist_plot<-ggplot(betadisp12_twofac_df, aes(x=Herb_Trt, y=distances, color=Herb_Trt, shape=Nutrient_Trt))+
#   geom_point(position=position_jitter(w = 0.3, h = 0.01), aes(shape=Nutrient_Trt), size=2.5, stroke=1)+
#   stat_summary(fun.data=data_summary, color="black", size=1, shape=16)+
#   scale_color_manual(values=c( "#5637C8","#C83761",  "#A9C837", "#37C89E"),
#                      name = "Herbivory", labels = c("High","Medium", "Low", "Very low" ))+
#   scale_shape_manual(values=c(1,2))+
#   scale_x_discrete(name ="Consumer pressure", labels=c("High","Medium", "Low", "Very low"))+
#   ylim(c(0,0.35))+
#   coord_fixed()+
#   theme_bw()+
#   ylab("Distance to centroid")+
#   labs(color="Herbivory & Predation")+
#   theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
#   theme(legend.position = "none")+
#   theme(aspect.ratio = 4/4)
# #ggsave("figures/t12_dist_plot.pdf", width=4, height=4, units="in")

# Plot Figure 2 ----------------------------------------------------------------

nmds_distance_plot<-cowplot::plot_grid(nmds_t0_plot, nmds_t0_vectors, t0_dist_plot, nmds_t12_plot, nmds_t12_vectors, t12_dist_plot, 
                                       align = c("v","h"), labels = c("(a)","(b)","(c)","(d)","(e)","(f)"))
legend<-cowplot::plot_grid(nmds_legend, vector_legend, ncol=1, nrow = 2,  align="vh")
nmds_distance_plot_legend<-cowplot::plot_grid(nmds_distance_plot, legend, align="vh", rel_widths = c(5, 1))
ggsave("figures/nmds_distance_plot_legend.pdf", width=11, height=6, units="in")



nmds_distance_plot<-cowplot::plot_grid(nmds_t0_plot, nmds_t0_vectors, t0_dist_plot, nmds_t12_plot, nmds_t12_vectors, t12_dist_plot, 
                                       align = c("v","h"), ncol=2, nrow=3, labels = c("(a)","(b)","(c)","(d)","(e)","(f)"))



# PERMANOVA -------------------------------------------------------------------
# does community composition differ by treatment?

## PERMANOVA T0 ----------------------------------------------------------------
### Herbivory ------------------------------------------------------------------

adonis2(sp_matrix0 ~ Herb_Trt, data=meta0, permutations = 999, method = "bray")
# NS, P=0.759, F(3,28)=0.643, R2=0.064
adonis2(sp_matrix0_twofac ~ Herb_Trt, data=meta0_twofac, permutations = 999, method = "bray")
# NS, P=0.759, F(3,28)=0.643
# fyi running this with the two-factor data frame and the regular dataframe return same rsults. duh they should

### Nutrients ------------------------------------------------------------------

adonis2(sp_matrix0 ~ Nutrient_Trt, data=meta0, permutations = 999, method = "bray")
# NS, P=0.958, F(1,30)=0.0382

### two-factor Herbivory and Nutrients -----------------------------------------

adonis2(sp_matrix0_twofac ~ Twofac, data=meta0_twofac, permutations = 999, method = "bray")
adonis2(sp_matrix0 ~ Nutrient_Trt + Herb_Trt + Nutrient_Trt*Herb_Trt, data=meta0, permutations = 999, method = "bray")
####################### Df SumOfSqs      R2       F Pr(>F) 
#Nutrient_Trt           1  0.00048 0.00127 0.0373  0.956
#Herb_Trt               3  0.02426 0.06443 0.6299  0.786
#Nutrient_Trt:Herb_Trt  3  0.04368 0.11600 1.1341  0.326
#Residual              24  0.30815 0.81830 

#3 test of disprsions separately for two factors
#------- Herbivory treatment -----
bray_mod0_herb<-betadisper(bray_dist0, meta0$Herb_Trt, type = "centroid")
anova(bray_mod0_herb)
#          Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     3 0.005148 0.0017161  1.3706 0.2721
#Residuals 28 0.035060 0.0012521   

#------- Nutrient treatment -----
bray_mod0_nutrients<-betadisper(bray_dist0, meta0$Nutrient_Trt, type = "centroid")
anova(bray_mod0_nutrients)
#Df  Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.000792 0.00079205  0.6991 0.4097
#Residuals 30 0.033991 0.00113303  

## PERMANOVA T12 ---------------------------------------------------------------
### Herbivory ------------------------------------------------------------------

adonis2(sp_matrix12 ~ Herb_Trt, data=meta12, permutations = 999, method = "bray")
# *** P=0.001, F(3,28)=9.699, R2=0.510

### Nutrients ------------------------------------------------------------------

adonis2(sp_matrix12 ~ Nutrient_Trt, data=meta12, permutations = 999, method = "bray")
# NS, P=0.153, F(1,30)=0.145, R2=0.052

### two-factor Herbivory and Nutrients -----------------------------------------

adonis2(sp_matrix12 ~ Nutrient_Trt + Herb_Trt + Nutrient_Trt*Herb_Trt, data=meta12, permutations = 999, method = "bray")
####################### Df SumOfSqs      R2       F Pr(>F)    
#Nutrient_Trt           1  0.14536 0.05158  3.2099  0.023 *  
#Herb_Trt               3  1.43600 0.50960 10.5703  0.001 ***
#Nutrient_Trt:Herb_Trt  3  0.14971 0.05313  1.1020  0.354    
#Residual              24  1.08682 0.38569            

#3 test of disprsions separately for two factors
#------- Herbivory treatment -----
bray_mod12_herb<-betadisper(bray_dist12, meta12$Herb_Trt, type = "centroid")
anova(bray_mod12_herb)
#          Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     3 0.083691 0.0278971  10.728 7.253e-05 ***
#Residuals 28 0.072811 0.0026004     
TukeyHSD(bray_mod12_herb)

#------- Nutrient treatment -----
bray_mod12_nutrients<-betadisper(bray_dist12, meta0$Nutrient_Trt, type = "centroid")
anova(bray_mod12_nutrients)
#Df  Sum Sq   Mean Sq F value Pr(>F)
#Groups     1 0.00000 0.0000001       0 0.9975
#Residuals 30 0.24882 0.0082939   

#--------------------------- SIMPER T12 --------------------------
#only doing for t12 because there were not significant differences at t0
sim_herb <- with(meta12_twofac, simper(sp_matrix12, Herb_Trt, permutations = 99))
summary(sim_herb)
sim_herb_df<-summary(sim_herb)
sim_herb_df<-data.frame(unclass(sim_herb_df))
write.csv(sim_herb_df, "model outputs/simper/sim_herb_df.csv")

sim_nutrients <- with(meta12_twofac, simper(sp_matrix12, Nutrient_Trt, permutations = 99))
summary(sim_nutrients)
sim_nutrients_df<-summary(sim_nutrients)
sim_nutrients_df<-data.frame(unclass(sim_nutrients_df))
write.csv(sim_nutrients_df, "model outputs/simper/sim_nutrients_df.csv")

#### KELLY start here
# Beta Dispersion all timepoints - two factor --------------------------

cover_twofac<-cover
# creating a single factor that combines herbivory and nutrients. 8 levels total
cover_twofac$Twofac<-as.factor(paste(cover_twofac$Herb_Trt, cover_twofac$Nutrient_Trt, sep="_"))
cover_twofac$Twofac<-factor(cover_twofac$Twofac)

# subsetting data for each timepoint 
cover0_twofac<-subset(cover_twofac, Timepoint_numeric==0)
cover1_twofac<-subset(cover_twofac, Timepoint_numeric==1)
cover2_twofac<-subset(cover_twofac, Timepoint_numeric==2)
cover3_twofac<-subset(cover_twofac, Timepoint_numeric==3)
cover4_twofac<-subset(cover_twofac, Timepoint_numeric==4)
cover6_twofac<-subset(cover_twofac, Timepoint_numeric==6)
cover8_twofac<-subset(cover_twofac, Timepoint_numeric==8)
cover9_twofac<-subset(cover_twofac, Timepoint_numeric==9)
cover10_twofac<-subset(cover_twofac, Timepoint_numeric==10)
cover11_twofac<-subset(cover_twofac, Timepoint_numeric==11)
cover12_twofac<-subset(cover_twofac, Timepoint_numeric==12)

# making the species matrix for each timepoint
sp_matrix0_twofac<-as.matrix(cover0_twofac[,c(13:63)])
sp_matrix1_twofac<-as.matrix(cover1_twofac[,c(13:63)])
sp_matrix2_twofac<-as.matrix(cover2_twofac[,c(13:63)])
sp_matrix3_twofac<-as.matrix(cover3_twofac[,c(13:63)])
sp_matrix4_twofac<-as.matrix(cover4_twofac[,c(13:63)])
sp_matrix6_twofac<-as.matrix(cover6_twofac[,c(13:63)])
sp_matrix8_twofac<-as.matrix(cover8_twofac[,c(13:63)])
sp_matrix9_twofac<-as.matrix(cover9_twofac[,c(13:63)])
sp_matrix10_twofac<-as.matrix(cover10_twofac[,c(13:63)])
sp_matrix11_twofac<-as.matrix(cover11_twofac[,c(13:63)])
sp_matrix12_twofac<-as.matrix(cover12_twofac[,c(13:63)])

#metadata for each timepoint
meta0_twofac<-cover0_twofac[,c(1:12,64)]
meta1_twofac<-cover1_twofac[,c(1:12,64)]
meta2_twofac<-cover2_twofac[,c(1:12,64)]
meta3_twofac<-cover3_twofac[,c(1:12,64)]
meta4_twofac<-cover4_twofac[,c(1:12,64)]
meta6_twofac<-cover6_twofac[,c(1:12,64)]
meta8_twofac<-cover8_twofac[,c(1:12,64)]
meta9_twofac<-cover9_twofac[,c(1:12,64)]
meta10_twofac<-cover10_twofac[,c(1:12,64)]
meta11_twofac<-cover11_twofac[,c(1:12,64)]
meta12_twofac<-cover12_twofac[,c(1:12,64)]

# Bray-Curtis dissimilarity matrix 
bray_dist0_twofac<-vegdist(sp_matrix0_twofac, method = "bray")
bray_dist1_twofac<-vegdist(sp_matrix1_twofac, method = "bray")
bray_dist2_twofac<-vegdist(sp_matrix2_twofac, method = "bray")
bray_dist3_twofac<-vegdist(sp_matrix3_twofac, method = "bray")
bray_dist4_twofac<-vegdist(sp_matrix4_twofac, method = "bray")
bray_dist6_twofac<-vegdist(sp_matrix6_twofac, method = "bray")
bray_dist8_twofac<-vegdist(sp_matrix8_twofac, method = "bray")
bray_dist9_twofac<-vegdist(sp_matrix9_twofac, method = "bray")
bray_dist10_twofac<-vegdist(sp_matrix10_twofac, method = "bray")
bray_dist11_twofac<-vegdist(sp_matrix11_twofac, method = "bray")
bray_dist12_twofac<-vegdist(sp_matrix12_twofac, method = "bray")

# Betadisper - Testing for differences in dispersion among 8 factors (one-way test)
bray_mod0_twofac<-betadisper(bray_dist0_twofac, meta0_twofac$Twofac, type = "centroid")
bray_mod1_twofac<-betadisper(bray_dist1_twofac, meta1_twofac$Twofac, type = "centroid")
bray_mod2_twofac<-betadisper(bray_dist2_twofac, meta2_twofac$Twofac, type = "centroid")
bray_mod3_twofac<-betadisper(bray_dist3_twofac, meta3_twofac$Twofac, type = "centroid")
bray_mod4_twofac<-betadisper(bray_dist4_twofac, meta4_twofac$Twofac, type = "centroid")
bray_mod6_twofac<-betadisper(bray_dist6_twofac, meta6_twofac$Twofac, type = "centroid")
bray_mod8_twofac<-betadisper(bray_dist8_twofac, meta8_twofac$Twofac, type = "centroid")
bray_mod9_twofac<-betadisper(bray_dist9_twofac, meta9_twofac$Twofac, type = "centroid")
bray_mod10_twofac<-betadisper(bray_dist10_twofac, meta10_twofac$Twofac, type = "centroid")
bray_mod11_twofac<-betadisper(bray_dist11_twofac, meta11_twofac$Twofac, type = "centroid")
bray_mod12_twofac<-betadisper(bray_dist12_twofac, meta12_twofac$Twofac, type = "centroid")

# Anova on Betadisepr models - Testing for differences in dispersion among 8 factors (one-way test)
anova(bray_mod0_twofac)
anova(bray_mod1_twofac)
anova(bray_mod2_twofac)
anova(bray_mod3_twofac)
anova(bray_mod4_twofac)
anova(bray_mod6_twofac)
anova(bray_mod8_twofac)
anova(bray_mod9_twofac)
anova(bray_mod10_twofac)
anova(bray_mod11_twofac)
anova(bray_mod12_twofac)

# saving model outputs
write.csv(anova(bray_mod0_twofac), "model outputs/betadisper/betadisp_anova_bray_mod0_twofac.csv")
write.csv(anova(bray_mod1_twofac), "model outputs/betadisper/betadisp_anova_bray_mod1_twofac.csv")
write.csv(anova(bray_mod2_twofac), "model outputs/betadisper/betadisp_anova_bray_mod2_twofac.csv")
write.csv(anova(bray_mod3_twofac), "model outputs/betadisper/betadisp_anova_bray_mod3_twofac.csv")
write.csv(anova(bray_mod4_twofac), "model outputs/betadisper/betadisp_anova_bray_mod4_twofac.csv")
write.csv(anova(bray_mod6_twofac), "model outputs/betadisper/betadisp_anova_bray_mod6_twofac.csv")
write.csv(anova(bray_mod8_twofac), "model outputs/betadisper/betadisp_anova_bray_mod8_twofac.csv")
write.csv(anova(bray_mod9_twofac), "model outputs/betadisper/betadisp_anova_bray_mod9_twofac.csv")
write.csv(anova(bray_mod10_twofac), "model outputs/betadisper/betadisp_anova_bray_mod10_twofac.csv")
write.csv(anova(bray_mod11_twofac), "model outputs/betadisper/betadisp_anova_bray_mod11_twofac.csv")
write.csv(anova(bray_mod12_twofac), "model outputs/betadisper/betadisp_anova_bray_mod12_twofac.csv")

# post hoc Tukey test - testing for significant differences in dispersion between individual groups (think we don't need these)
tukey_0<-TukeyHSD(bray_mod0_twofac, "group")
tukey_1<-TukeyHSD(bray_mod1_twofac, "group")
tukey_2<-TukeyHSD(bray_mod2_twofac, "group")
tukey_3<-TukeyHSD(bray_mod3_twofac, "group")
tukey_4<-TukeyHSD(bray_mod4_twofac, "group")
tukey_6<-TukeyHSD(bray_mod6_twofac, "group")
tukey_8<-TukeyHSD(bray_mod8_twofac, "group")
tukey_9<-TukeyHSD(bray_mod9_twofac, "group")
tukey_10<-TukeyHSD(bray_mod10_twofac, "group")
tukey_11<-TukeyHSD(bray_mod11_twofac, "group")
tukey_12<-TukeyHSD(bray_mod12_twofac, "group")

# saving results as CSV
# write.csv(as.data.frame(tukey_0$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod0_twofac.csv") 
# write.csv(as.data.frame(tukey_1$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod1_twofac.csv") 
# write.csv(as.data.frame(tukey_2$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod2_twofac.csv") 
# write.csv(as.data.frame(tukey_3$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod3_twofac.csv") 
# write.csv(as.data.frame(tukey_4$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod4_twofac.csv") 
# write.csv(as.data.frame(tukey_6$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod6_twofac.csv") 
# write.csv(as.data.frame(tukey_8$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod8_twofac.csv") 
# write.csv(as.data.frame(tukey_9$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod9_twofac.csv") 
# write.csv(as.data.frame(tukey_10$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod10_twofac.csv") 
# write.csv(as.data.frame(tukey_11$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod11_twofac.csv") 
# write.csv(as.data.frame(tukey_12$group), "model outputs/betadisper_tukey/betadisp_tukey_bray_mod12_twofac.csv") 

# testing for interaction in two-way permanova - testing for differences in location of group centroids
# write.csv(adonis2(sp_matrix0_twofac ~ Herb_Trt*Nutrient_Trt, data=meta0_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod0_twofac.csv")
# write.csv(adonis2(sp_matrix1_twofac ~ Herb_Trt*Nutrient_Trt, data=meta1_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod1_twofac.csv")
# write.csv(adonis2(sp_matrix2_twofac ~ Herb_Trt*Nutrient_Trt, data=meta2_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod2_twofac.csv")
# write.csv(adonis2(sp_matrix3_twofac ~ Herb_Trt*Nutrient_Trt, data=meta3_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod3_twofac.csv")
# write.csv(adonis2(sp_matrix4_twofac ~ Herb_Trt*Nutrient_Trt, data=meta4_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod4_twofac.csv")
# write.csv(adonis2(sp_matrix6_twofac ~ Herb_Trt*Nutrient_Trt, data=meta6_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod6_twofac.csv")
# write.csv(adonis2(sp_matrix8_twofac ~ Herb_Trt*Nutrient_Trt, data=meta8_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod8_twofac.csv")
# write.csv(adonis2(sp_matrix9_twofac ~ Herb_Trt*Nutrient_Trt, data=meta9_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod9_twofac.csv")
# write.csv(adonis2(sp_matrix10_twofac ~ Herb_Trt*Nutrient_Trt, data=meta10_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod10_twofac.csv")
# write.csv(adonis2(sp_matrix11_twofac ~ Herb_Trt*Nutrient_Trt, data=meta11_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod11_twofac.csv")
# write.csv(adonis2(sp_matrix12_twofac ~ Herb_Trt*Nutrient_Trt, data=meta12_twofac, permutations = 999, method = "bray"), "model outputs/adonis/adonis_bray_mod12_twofac.csv")

# interaction is not significant at any timpoint
# therefore we can proceed to test for effects of herbivory and nutrients separately 

#---------------------- pulling distance to centroid data for use in mixed effects models -----------------------

# we calculated distances to centroid for 8 groups. groups represent a single  factor that combines herbivory and nutrients
# since we want to test for the effects of nutrients and herbivory on betadiversity we will use those data 

#pulling distances into a dataframe for each timepoint
betadisp0_twofac_df<-data.frame(group=bray_mod0_twofac$group, distances=bray_mod0_twofac$distances)
betadisp1_twofac_df<-data.frame(group=bray_mod1_twofac$group, distances=bray_mod1_twofac$distances)
betadisp2_twofac_df<-data.frame(group=bray_mod2_twofac$group, distances=bray_mod2_twofac$distances)
betadisp3_twofac_df<-data.frame(group=bray_mod3_twofac$group, distances=bray_mod3_twofac$distances)
betadisp4_twofac_df<-data.frame(group=bray_mod4_twofac$group, distances=bray_mod4_twofac$distances)
betadisp6_twofac_df<-data.frame(group=bray_mod6_twofac$group, distances=bray_mod6_twofac$distances)
betadisp8_twofac_df<-data.frame(group=bray_mod8_twofac$group, distances=bray_mod8_twofac$distances)
betadisp9_twofac_df<-data.frame(group=bray_mod9_twofac$group, distances=bray_mod9_twofac$distances)
betadisp10_twofac_df<-data.frame(group=bray_mod10_twofac$group, distances=bray_mod10_twofac$distances)
betadisp11_twofac_df<-data.frame(group=bray_mod11_twofac$group, distances=bray_mod11_twofac$distances)
betadisp12_twofac_df<-data.frame(group=bray_mod12_twofac$group, distances=bray_mod12_twofac$distances)

# adding a factor column for timepoint
betadisp0_twofac_df$Timepoint<-"T0"
betadisp1_twofac_df$Timepoint<-"T1"
betadisp2_twofac_df$Timepoint<-"T2"
betadisp3_twofac_df$Timepoint<-"T3"
betadisp4_twofac_df$Timepoint<-"T4"
betadisp6_twofac_df$Timepoint<-"T6"
betadisp8_twofac_df$Timepoint<-"T8"
betadisp9_twofac_df$Timepoint<-"T9"
betadisp10_twofac_df$Timepoint<-"T10"
betadisp11_twofac_df$Timepoint<-"T11"
betadisp12_twofac_df$Timepoint<-"T12"

# recombining with the meta dataframe to get plot id's
# vegan spits out output data frame in same order as input. same as with metamds function
betadisp0_twofac_df<-cbind(betadisp0_twofac_df, meta0_twofac[,c(4:8,10)])
betadisp1_twofac_df<-cbind(betadisp1_twofac_df, meta1_twofac[,c(4:8,10)])
betadisp2_twofac_df<-cbind(betadisp2_twofac_df, meta2_twofac[,c(4:8,10)])
betadisp3_twofac_df<-cbind(betadisp3_twofac_df, meta3_twofac[,c(4:8,10)])
betadisp4_twofac_df<-cbind(betadisp4_twofac_df, meta4_twofac[,c(4:8,10)])
betadisp6_twofac_df<-cbind(betadisp6_twofac_df, meta6_twofac[,c(4:8,10)])
betadisp8_twofac_df<-cbind(betadisp8_twofac_df, meta8_twofac[,c(4:8,10)])
betadisp9_twofac_df<-cbind(betadisp9_twofac_df, meta9_twofac[,c(4:8,10)])
betadisp10_twofac_df<-cbind(betadisp10_twofac_df, meta10_twofac[,c(4:8,10)])
betadisp11_twofac_df<-cbind(betadisp11_twofac_df, meta11_twofac[,c(4:8,10)])
betadisp12_twofac_df<-cbind(betadisp12_twofac_df, meta12_twofac[,c(4:8,10)])

# combining all the dataframes from individual timepoints
betadisp_twofac_df<-rbind(betadisp0_twofac_df, betadisp1_twofac_df, betadisp2_twofac_df, betadisp3_twofac_df, betadisp4_twofac_df, 
                          betadisp6_twofac_df, betadisp8_twofac_df, betadisp9_twofac_df, betadisp10_twofac_df, betadisp11_twofac_df, betadisp12_twofac_df)
# making timepoint a factor
betadisp_twofac_df$Timepoint<-as.factor(betadisp_twofac_df$Timepoint) 
betadisp_twofac_df$Timepoint<-factor(betadisp_twofac_df$Timepoint)
# reordering factor levels
betadisp_twofac_df$Timepoint <- ordered(betadisp_twofac_df$Timepoint, levels = 
                                          c("T0","T1","T2","T3","T4","T6","T8","T9","T10","T11","T12"))
# taking the single factor that combined herbivory and nutrients and splitting it into 2 factors
betadisp_twofac_df$Twofac<-betadisp_twofac_df$group
betadisp_twofac_df<-(separate(data = betadisp_twofac_df, col=group, into = c("Herb_Trt", "Nutrient_Trt"), sep = "_"))
# herbivory 
betadisp_twofac_df$Herb_Trt<-as.factor(betadisp_twofac_df$Herb_Trt)
betadisp_twofac_df$Herb_Trt<-factor(betadisp_twofac_df$Herb_Trt)
# nutrients
betadisp_twofac_df$Nutrient_Trt<-as.factor(betadisp_twofac_df$Nutrient_Trt)
betadisp_twofac_df$Nutrient_Trt<-factor(betadisp_twofac_df$Nutrient_Trt)

# merging with timepoint dataframe
betadisp_twofac_df<-merge(betadisp_twofac_df, timepoints, by="Timepoint")
betadisp_twofac_df$Herb_Trt<-factor(betadisp_twofac_df$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))

# save as csv for use in a later script
write.csv(betadisp_twofac_df, "data_summaries/betadisp_twofac_df.csv")


#---------------------- summarizing the data for plotting  -----------------------

## one factor - consumer pressure ----------------------------------------------
# summarizing for Herbivory and Nutrient treatments
betadisp_twofac_df_sum<-ddply(betadisp_twofac_df, .(Timepoint, Nutrient_Trt, Herb_Trt, Twofac), summarize,
                              mean_dist=mean(distances), 
                              se=sd(distances)/sqrt(length(distances)))

# merging with timepoints datafraame so we can plot with dates
betadisp_twofac_df_sum<-merge(betadisp_twofac_df_sum, timepoints, by="Timepoint")


# one factor only
betadisp_df_sum<-ddply(betadisp_twofac_df, .(Timepoint, Herb_Trt), summarize,
                       mean_dist=mean(distances), 
                       se=sd(distances)/sqrt(length(distances)))

# merging with timepoints datafraame so we can plot with dates
betadisp_df_sum<-merge(betadisp_df_sum, timepoints, by="Timepoint")


## two factors - consumer pressure x nutrients ---------------------------------
betadisp_df_twofac_sum<-ddply(betadisp_twofac_df, .(Timepoint, Herb_Trt, Nutrient_Trt), summarize,
                              mean_dist=mean(distances), 
                              se=sd(distances)/sqrt(length(distances)))

# merging with timepoints datafraame so we can plot with dates
betadisp_df_twofac_sum<-merge(betadisp_df_twofac_sum, timepoints, by="Timepoint")
betadisp_df_twofac_sum$Twofac<-as.factor(paste(betadisp_df_twofac_sum$Herb_Trt, betadisp_df_twofac_sum$Nutrient_Trt, sep="_"))


# Mixed effects model ------------------------------------------------------------
# evaluating effects of herbivory and nutrients on community dispersion 
#  Block_Plot_Herb_Trt as a random effect
# factors are Days_since_start, Herb_Trt, Nutrient_trt
# included plot (Block_Plot_Herb_Trt) as a random effect to account for repeated measures through time.

dist.mod.1<-lmer(distances ~ Herb_Trt + Nutrient_Trt + Days_since_start+Herb_Trt*Nutrient_Trt*Days_since_start  + (1|Block_Plot_Herb_Trt), data=betadisp_twofac_df)
car::Anova(dist.mod.1) # no three-way interaction

dist.mod.2<-lmer(distances ~ Herb_Trt + Nutrient_Trt + Days_since_start + Herb_Trt*Nutrient_Trt + Herb_Trt*Days_since_start + Nutrient_Trt*Days_since_start + (1|Block_Plot_Herb_Trt), data=betadisp_twofac_df)
car::Anova(dist.mod.2) #Type II Wald chisquare tests
plot(allEffects(dist.mod.2)) ## this is the plot you want for the supp!
plot_model(dist.mod.2)


qqnorm(ranef(dist.mod.2)[[1]][,1])
qqline(ranef(dist.mod.2)[[1]][,1])

# getting the predicted relationships
dist_ref<-ref_grid(dist.mod.2, at=list(Days_since_start=c(0,200,400, 600, 800, 1000, 1200, 1362, 1500)))
dist_pred<-data.frame(emmip(dist_ref, Herb_Trt + Nutrient_Trt ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

# plotting the predicted relationships on top of the real data 
betadisper_days_plot<-ggplot()+ 
  geom_point(data=betadisp_twofac_df, aes(x=Days_since_start, y=distances, color=Herb_Trt, fill=Herb_Trt, shape=Nutrient_Trt),
             position=position_jitter(w = 25), size=2, stroke = 1, alpha=0.5)+
  geom_ribbon(data=dist_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Herb_Trt, group=tvar),alpha= 0.3) +
  geom_line(data=dist_pred, aes(x=Days_since_start, y = yvar, color=Herb_Trt, group=tvar, linetype=Nutrient_Trt), size = 0.9)+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                     name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_fill_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                    name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  scale_linetype_manual(values=c(1,2), name="Nutrients")+
  coord_fixed()+
  theme_bw()+
  ylab("Distance to centroid")+
  xlab("Days since start of experiment")+
  labs(color="Herbivory")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 4/6)
ggsave("figures/betadisper_days_plot.pdf", width=7, height=4, units="in")

# or plotting the effects of consumers and nutrients separately

# Herb_Trt
dist_herb_ref<-ref_grid(dist.mod.2, at=list(Days_since_start=c(0,200,400, 600, 800, 1000, 1200, 1362, 1500)))
dist_herb_pred<-data.frame(emmip(dist_herb_ref, Herb_Trt  ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

ggplot()+ 
  geom_point(data=betadisp_twofac_df, aes(x=Days_since_start, y=distances, color=Herb_Trt, fill=Herb_Trt, shape=Nutrient_Trt),
             position=position_jitter(w = 25), size=2, stroke = 1, alpha=0.5)+
  geom_ribbon(data=dist_herb_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Herb_Trt, group=Herb_Trt),alpha= 0.3) +
  geom_line(data=dist_herb_pred, aes(x=Days_since_start, y = yvar, color=Herb_Trt, group=Herb_Trt), size = 0.9)+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                     name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_fill_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                    name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  coord_fixed()+
  theme_bw()+
  ylab("Distance to centroid")+
  xlab("Days since start of experiment")+
  labs(color="Herbivory")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 4/6)

# Nutrient_Trt - Don't really like this
dist_nut_ref<-ref_grid(dist.mod.2, at=list(Days_since_start=c(0,200,400, 600, 800, 1000, 1200, 1362, 1500)))
dist_nut_pred<-data.frame(emmip(dist_nut_ref, Nutrient_Trt  ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

ggplot()+ 
  geom_point(data=betadisp_twofac_df, aes(x=Days_since_start, y=distances, color=Herb_Trt, fill=Herb_Trt, shape=Nutrient_Trt),
             position=position_jitter(w = 25), size=2, stroke = 1, alpha=0.5)+
  geom_ribbon(data=dist_nut_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, group=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=dist_nut_pred, aes(x=Days_since_start, y = yvar, group=Nutrient_Trt, linetype=Nutrient_Trt), size = 0.9)+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"), 
                     name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  # scale_fill_manual(values=c("#5637C8","#C83761",  "#A9C837", "#37C89E"), 
  #                   name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  scale_linetype_manual(values=c(1,2), name="Nutrients")+
  coord_fixed()+
  theme_bw()+
  ylab("Distance to centroid")+
  xlab("Days since start of experiment")+
  labs(color="Herbivory")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(aspect.ratio = 4/6)
