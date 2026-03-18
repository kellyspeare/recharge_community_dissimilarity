# ----------------------------------------------------------------------- #
# Reducing consumer pressure increases community dissimilarity leading to 
# coral- or algal- dominance on a coral reef

# 3_rank_abundance_shifts

# ----------------------------------------------------------------------- #
# this script analyzes data on rank abundance shifts (mean rank shifts)
# Plots Figures 5, 6

# Packages --------------------------------------------------------------- #

library(codyn)
library(plyr)
#library(dplyr)
library(tidyverse)
library(lme4)
library(pbkrtest)
library(effects)
library(emmeans)
library(lmerTest)
library(cowplot)
library(ggplot2)
library(ciTools)
library(forcats)
library(car)
library(DHARMa)

# dataframe with all of the cover data
cover<-read.csv("data/data_clean/plots_benthic_cover.csv", header=TRUE, stringsAsFactors = TRUE)

#reordering factor levels
cover$Herb_Trt<-factor(cover$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))

cover_all<-cover
cover_all_l<-cover_all %>% pivot_longer(cols=13:56, names_to="category", values_to="percent") 
cover_all_l$Date<-as.Date(cover_all_l$Date, "%m/%d/%Y")
cover_all_l$category<-as.factor(cover_all_l$category)

# labeller for plots
herb_names <- c(
  `Open` = "High",
  `3X3` = "Medium",
  `2X2` = "Low",
  `1X1` = "Very Low"
)

# ----------------------------------------------------------------------------------#
# ----------------- Calculating Mean Rank Shifts through time -----------------------
# ----------------------------------------------------------------------------------#

rank_cover<-cover_all_l
rank_cover$Timepoint_ordered<-ifelse(rank_cover$Time_point==0, 0,
                                     ifelse(rank_cover$Time_point==1, 1,
                                            ifelse(rank_cover$Time_point==2, 2,
                                                   ifelse(rank_cover$Time_point==3,3,
                                                          ifelse(rank_cover$Time_point==4,4,
                                                                 ifelse(rank_cover$Time_point==6,5,
                                                                        ifelse(rank_cover$Time_point==8,6,
                                                                               ifelse(rank_cover$Time_point==9,7,
                                                                                      ifelse(rank_cover$Time_point==10,8,
                                                                                             ifelse(rank_cover$Time_point==11,9,
                                                                                                    ifelse(rank_cover$Time_point==12,10, NA)))))))))))

ranks<-rank_shift(rank_cover, time.var="Timepoint_ordered", species.var="category", abundance.var="percent", replicate.var="Block_Plot_Herb_Trt")
ranks<-rank_shift(subset(rank_cover, category!="Turf_Bare_Encrusting_Algae" & category!="Sand" &category!="Rubble"), time.var="Timepoint_ordered", species.var="category", abundance.var="percent", replicate.var="Block_Plot_Herb_Trt")
summary(ranks)
ranks$Block_Plot_Herb_Trt<-as.factor(ranks$Block_Plot_Herb_Trt)

# timepoint at the START of the timepoint pair
ranks$Time_point_start<-ifelse(ranks$year_pair=="0-1",0, 
                                 ifelse(ranks$year_pair=="1-2",1, 
                                        ifelse(ranks$year_pair=="2-3",2, 
                                               ifelse(ranks$year_pair=="3-4",3,
                                                      ifelse(ranks$year_pair=="4-5",4,
                                                             ifelse(ranks$year_pair=="5-6",6,
                                                                    ifelse(ranks$year_pair=="6-7",8,
                                                                           ifelse(ranks$year_pair=="7-8",9,
                                                                                  ifelse(ranks$year_pair=="8-9",10,
                                                                                         ifelse(ranks$year_pair=="9-10",11,
                                                                                                ifelse(ranks$year_pair=="10-11",12,
                                                                                                       ifelse(ranks$year_pair=="3-5", 3, NA)))))))))))) # this is for those two plots where we don't have T4 data
# timepoint at the END of the timepoint pair
ranks$Time_point_end<-ifelse(ranks$year_pair=="0-1",1, 
                                ifelse(ranks$year_pair=="1-2",2, 
                                       ifelse(ranks$year_pair=="2-3",3, 
                                              ifelse(ranks$year_pair=="3-4",4,
                                                     ifelse(ranks$year_pair=="4-5",6,
                                                            ifelse(ranks$year_pair=="5-6",8,
                                                                   ifelse(ranks$year_pair=="6-7",9,
                                                                          ifelse(ranks$year_pair=="7-8",10,
                                                                                 ifelse(ranks$year_pair=="8-9",11,
                                                                                        ifelse(ranks$year_pair=="9-10",12,
                                                                                               ifelse(ranks$year_pair=="3-5", 6, NA))))))))))) # this is for those two plots where we don't have T4 data




plots<-cover[,c(1,2,3,7:12,53)] # colum 53 is Turf_Bare_Encrusting_Algae
plots$Date<-as.Date(as.character(plots$Date),format="%m/%d/%Y")
#plots<-data.frame(unique(plots, plots$Block_Plot_Herb_Trt, incomparables = FALSE))

ranks_plots<-left_join(ranks, plots, by=c("Block_Plot_Herb_Trt", "Time_point_end"="Time_point"))

# reorder factor levels
ranks_plots$Herb_Trt<-factor(ranks_plots$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))

#creating a column for ratio of coral to macroalgae
ranks_plots$ratio_coral_macro<-(1+ ranks_plots$Sum_Coral) / (1 + ranks_plots$Sum_Macroalgae)
ranks_plots$log_ratio_coral_macro<-log(ranks_plots$ratio_coral_macro)
#binary column for dominated by corals or macroalgae
ranks_plots$binary_ratio<-as.factor(ifelse(ranks_plots$log_ratio_coral_macro>0, "coral", "macro"))
# binary column 1 for coral
ranks_plots$binary_coral<-ifelse(ranks_plots$log_ratio_coral_macro>0, 1, 0)
# binary column 1 for macroalgae
ranks_plots$binary_macro<-ifelse(ranks_plots$log_ratio_coral_macro>0, 0, 1)


#### how correlated is time (days since start of experiment) and colonizable space?
cor.test(ranks_plots$Days_since_start, ranks_plots$Turf_Bare_Encrusting_Algae, method=c("pearson"))
# cor= -0.545466 
ggplot(ranks_plots, aes(x=Days_since_start, y=Turf_Bare_Encrusting_Algae, color=Herb_Trt))+
  geom_point()
# so it really depends on herbivory treatment

#make date into a numeric
startdate <- as.Date("2018-07-10","%Y-%m-%d")
ranks_plots$Days_since_start  <- as.numeric(difftime(ranks_plots$Date, startdate ,units="days"))

# ----------------------------------------------------------------------#
# -------------------- Colonizable Space Model -------------------------
# ----------------------------------------------------------------------#

col_space<-subset(cover_all_l, category=="Turf_Bare_Encrusting_Algae") # col_space should have 351 rows
startdate <- as.Date("2018-07-10","%Y-%m-%d")
col_space$Days_since_start  <- as.numeric(difftime(col_space$Date, startdate ,units="days"))

space.mod.1<-lmer(percent ~ Days_since_start + Herb_Trt + Nutrient_Trt + Days_since_start*Herb_Trt*Nutrient_Trt + (1|Block_Plot_Herb_Trt), data=col_space)

space.mod.1<-lmer(percent ~ Days_since_start + Herb_Trt + Nutrient_Trt + Days_since_start*Herb_Trt*Nutrient_Trt + (1|Block_Plot_Herb_Trt), data=col_space)
car::Anova(space.mod.1) #NS 3-way interaction, so remove it

space.mod.2<-lmer(percent ~ Days_since_start + Herb_Trt + Nutrient_Trt + Days_since_start*Herb_Trt + Herb_Trt*Nutrient_Trt + Days_since_start*Nutrient_Trt + (1|Block_Plot_Herb_Trt), data=col_space)
car::Anova(space.mod.2)
plot(simulateResiduals(fittedModel=space.mod.2, quantileFunction = qnorm))
#                                 Chisq Df Pr(>Chisq)    
# Days_since_start              507.8720  1  < 2.2e-16 ***
# Herb_Trt                      213.0150  3  < 2.2e-16 ***
# Nutrient_Trt                    7.1503  1   0.007495 ** 
# Days_since_start:Herb_Trt     218.4827  3  < 2.2e-16 ***
# Herb_Trt:Nutrient_Trt           1.6052  3   0.658216    
# Days_since_start:Nutrient_Trt   4.5663  1   0.032606 * 

plot(allEffects(space.mod.2))

qqnorm(ranef(space.mod.2)[[1]][,1])
qqline(ranef(space.mod.2)[[1]][,1])

### predicted relationships
# Herb_Trt are slopes different?
emtrends(space.mod.2, pairwise ~ Herb_Trt, var = "Days_since_start")
cld(emtrends(space.mod.2, pairwise ~ Herb_Trt, var = "Days_since_start"))
# Nutrient treatment are slopes different?
emtrends(space.mod.2, pairwise ~ Nutrient_Trt, var = "Days_since_start")
cld(emtrends(space.mod.2, pairwise ~ Nutrient_Trt, var = "Days_since_start"))

emtrends(space.mod.2, pairwise ~ Nutrient_Trt + Herb_Trt, var = "Days_since_start")

space_herb_ref<-ref_grid(space.mod.2, at=list(Days_since_start=c(0,200,400, 600, 800, 1000, 1200, 1300, 1400, 1472)))
space_herb_pred<-data.frame(emmip(space_herb_ref, Herb_Trt ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

## --------------------- Colonizable Space Herb Plot -------------------------

#Herb plot 
space_herb_plot<-ggplot()+
  geom_point(data=col_space, aes(x=Days_since_start, y=percent, color=Herb_Trt, group=Herb_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 20),  size=2, stroke = 1, alpha=0.5)+
  geom_ribbon(data=space_herb_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Herb_Trt),alpha= 0.3) +
  geom_line(data=space_herb_pred, aes(x=Days_since_start, y = yvar, color=Herb_Trt), size = 0.9)+
  #geom_smooth(method="lm")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_fill_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                    name = "Consumer pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  xlab("Days since start of experiment")+
  ylab("% colonizable space")+
  xlim(c(-50, 1600))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
space_herb_plot


## --------------------- Colonizable Space Nutrient  Plot -------------------------

space_nut_ref<-ref_grid(space.mod.2, at=list(Days_since_start=c(0,200,400, 600, 800, 1000, 1200, 1300, 1400, 1472)))
space_nut_pred<-data.frame(emmip(space_herb_ref, Nutrient_Trt ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

#Nutrient plot 
space_nutrient_plot<-ggplot()+
  geom_point(data=col_space, aes(x=Days_since_start, y=percent, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 20),  size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=space_nut_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=space_nut_pred, aes(x=Days_since_start, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Days since start of experiment")+
  ylab("% colonizable space")+
  xlim(c(-50, 1600))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
space_nutrient_plot

## --------------------- Plot Legends -------------------------
## legends
space_herb_plot_legend<-ggplot()+
  geom_point(data=ranks_plots, aes(x=Days_since_start, y=Turf_Bare_Encrusting_Algae, color=Herb_Trt, group=Herb_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 20),  size=2, stroke = 1, alpha=0.5)+
  geom_ribbon(data=space_herb_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Herb_Trt),alpha= 0.3) +
  geom_line(data=space_herb_pred, aes(x=Days_since_start, y = yvar, color=Herb_Trt), size = 0.9)+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumer\npressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_fill_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                    name = "Consumer\npressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_shape_manual(values=c(1,2), name="Nutrients", guide=NULL)+
  xlab("Days since start of experiment")+
  ylab("% colonizable space")+
  xlim(c(-50, 1600))+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
space_herb_legend<- get_legend(space_herb_plot_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))

space_nutrient_plot_legend<-ggplot()+
  geom_point(data=ranks_plots, aes(x=Days_since_start, y=Turf_Bare_Encrusting_Algae, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 20),  size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=space_nut_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=space_nut_pred, aes(x=Days_since_start, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Days since start of experiment")+
  ylab("% colonizable space (% cover)")+
  xlim(c(-50, 1600))+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
space_nutrient_legend<- get_legend(space_nutrient_plot_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))


# ----------------------------------------------------------------------#
# ------------------ Mean Rank Shifts (MRS) Model -----------------------
# ----------------------------------------------------------------------#

# mod with colonizable space and time
mrs.mod.1<-lmer(MRS ~ Days_since_start + Herb_Trt + Nutrient_Trt + Days_since_start*Herb_Trt*Nutrient_Trt + (1|Block_Plot_Herb_Trt), data=ranks_plots)
car::Anova(mrs.mod.1) #3-way interaction NS, remove

mrs.mod.2<-lmer(MRS ~ Days_since_start + Herb_Trt + Nutrient_Trt + Days_since_start*Herb_Trt + Herb_Trt*Nutrient_Trt + Days_since_start*Nutrient_Trt + (1|Block_Plot_Herb_Trt), data=ranks_plots)
car::Anova(mrs.mod.2)

#                                  Chisq Df Pr(>Chisq)    
# Days_since_start              38.6151  1  5.162e-10 ***
# Herb_Trt                      11.9838  3   0.007439 ** 
# Nutrient_Trt                   0.2792  1   0.597223    
# Days_since_start:Herb_Trt     43.5856  3  1.848e-09 ***
# Herb_Trt:Nutrient_Trt          1.7107  3   0.634553    
# Days_since_start:Nutrient_Trt  0.1450  1   0.703404    

# interaction btw days and herb_trt, no effect of nutrients or interaction 
plot(simulateResiduals(fittedModel=mrs.mod.2, quantileFunction = qnorm))

## --------------------- Mean Rank Shifts Herb Plot ---------------------

### predicted relationships
emtrends(mrs.mod.2, pairwise ~ Herb_Trt, var = "Days_since_start")
# herb trt
cld(emtrends(mrs.mod.2, pairwise ~ Herb_Trt, var = "Days_since_start"))
# slopes of lines are  statistically different. 
#open is different from all the rest. none of the others differ
# nutrient trt
cld(emtrends(mrs.mod.2, pairwise ~ Nutrient_Trt, var = "Days_since_start"))

# Herbivores Plot
mrs_herb_ref<-ref_grid(mrs.mod.2, at=list(Days_since_start=c(80,200,400,600,800,1000,1200,1350,1500)))
mrs_herb_pred<-data.frame(emmip(mrs_herb_ref, Herb_Trt ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

legend<- get_legend(plot_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))

mrs_herb_plot<-ggplot()+
  geom_point(data=ranks_plots, aes(x=Days_since_start, y=MRS, color=Herb_Trt, group=Herb_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 20), size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=mrs_herb_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Herb_Trt),alpha= 0.3) +
  geom_line(data=mrs_herb_pred, aes(x=Days_since_start, y = yvar, color=Herb_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumers pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_fill_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumers pressure", labels = c("High","Medium", "Low", "Very low" ))+
  xlab("Days since start of experiment")+
  ylab("Mean rank shifts")+
  xlim(c(-50, 1600))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
mrs_herb_plot

## --------------------- Mean Rank Shifts Nutrient Plot --------------------- 

mrs_nuts_ref<-ref_grid(mrs.mod.2, at=list(Days_since_start=c(80,200,400, 600, 800, 1000, 1200, 1362)))
mrs_nuts_pred<-data.frame(emmip(mrs_nuts_ref, Nutrient_Trt ~ Days_since_start, cov.reduce=range, plotit=FALSE)) #look at them visually

mrs_nutrient_plot<-ggplot()+
  geom_point(data=ranks_plots, aes(x=Days_since_start, y=MRS, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 20),  size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=mrs_nuts_pred, aes(x=Days_since_start, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=mrs_nuts_pred, aes(x=Days_since_start, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Days since start of experiment")+
  ylab("Mean rank shifts")+
  xlim(c(-50, 1600))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
mrs_nutrient_plot


# ----------------------------------------------------------------------------------#
# --------------------- MRS vs Beta Dispersion  model ------------------------------
# ----------------------------------------------------------------------------------#

#read in beta dispersion dataframe generated in script 2 

#betadisp_twofac_df<-read.csv("data_summaries/betadisp_twofac_df.csv", header=TRUE, stringsAsFactors = TRUE)

betadisp_dat<-betadisp_twofac_df
betadisp_dat$Timepoint<-as.character(betadisp_dat$Timepoint)
betadisp_dat$Timepoint<-as.factor(betadisp_dat$Timepoint)

#use ranks_plots for mrs
ranks_dat<-ranks_plots

# combine MRS data (second of 2 timepoints) with distance data that corresponds to the second timepoint
betadisp_mrs<-left_join(ranks_dat, betadisp_dat[,c(1,4,8,9)], by=c("Timepoint", "Block_Plot_Herb_Trt"))
# dataframe should have 319 rows

### mixed effects model
mrs.dist.mod.1<-lmer(distances ~ MRS + Herb_Trt + Nutrient_Trt + MRS*Herb_Trt*Nutrient_Trt + (1|Block_Plot_Herb_Trt), data=betadisp_mrs)
car::Anova(mrs.dist.mod.1) #no 3-way interaction

mrs.dist.mod.2<-lmer(distances ~ MRS*Herb_Trt + MRS*Nutrient_Trt + Nutrient_Trt*Herb_Trt + (1|Block_Plot_Herb_Trt), data=betadisp_mrs)
car::Anova(mrs.dist.mod.2) #no 2-way interactions

#                           Chisq Df Pr(>Chisq)    
# MRS                   16.2348  1  5.596e-05 ***
# Herb_Trt              40.6449  3  7.777e-09 ***
# Nutrient_Trt           0.0032  1   0.954980    
# MRS:Herb_Trt          14.7523  3   0.002041 ** 
# MRS:Nutrient_Trt       1.5620  1   0.211367    
# Herb_Trt:Nutrient_Trt  0.7565  3   0.859838   

plot(allEffects(mrs.dist.mod.2))
plot(mrs.dist.mod.2)
plot(simulateResiduals(mrs.dist.mod.2, quantileFunction = qnorm))

#herb trt - different slopes?
emtrends(mrs.dist.mod.2, pairwise ~ Herb_Trt, var = "MRS")
cld(emtrends(mrs.dist.mod.2, pairwise ~ Herb_Trt, var = "MRS"))

# nutrient trt - different slopes
emtrends(mrs.dist.mod.2, pairwise ~ Nutrient_Trt, var = "MRS")
cld(emtrends(mrs.dist.mod.2, pairwise ~ Nutrient_Trt, var = "MRS"))

### Consumer Pressure

mrs_dist_herb_ref<-ref_grid(mrs.dist.mod.2, at=list(MRS=c(0,1,2,3,4,5,6,7,7.2)))
mrs_dist_herb_pred<-data.frame(emmip(mrs_dist_herb_ref, Herb_Trt ~ MRS, cov.reduce=range, plotit=FALSE)) #look at them visually

#legend<- get_legend(plot_legend + theme(legend.box.margin = margin(0, 0, 0, 12)))
dist_herb_plot<-ggplot()+
  geom_point(data=betadisp_mrs, aes(x=MRS, y=distances, color=Herb_Trt, group=Herb_Trt, shape=Nutrient_Trt), 
             size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=mrs_dist_herb_pred, aes(x=MRS, ymin=yvar-SE, ymax=yvar+SE, fill=Herb_Trt),alpha= 0.3) +
  geom_line(data=mrs_dist_herb_pred, aes(x=MRS, y = yvar, color=Herb_Trt), size = 0.9)+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumers pressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_fill_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                    name = "Consumers pressure", labels = c("High","Medium", "Low", "Very low" ))+
  xlab("Mean rank shifts")+
  ylab("Beta dispersion")+
  xlim(c(-0.75,7.2))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
dist_herb_plot

### Nutrients
mrs_dist_nut_ref<-ref_grid(mrs.dist.mod.2, at=list(MRS=c(0,1,2,3,4,5,6,7,7.2)))
mrs_dist_nut_pred<-data.frame(emmip(mrs_dist_nut_ref, Nutrient_Trt ~ MRS, cov.reduce=range, plotit=FALSE)) #look at them visually

dist_nutrient_plot<-ggplot()+
  geom_point(data=betadisp_mrs, aes(x=MRS, y=distances, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=mrs_dist_nut_pred, aes(x=MRS, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=mrs_dist_nut_pred, aes(x=MRS, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Mean rank shifts")+
  ylab("Beta dispersion")+
  xlim(c(-0.75, 7.2))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
dist_nutrient_plot

# Plot Figure 5  -------------------------------------------------------------------

six_panel_plot<-cowplot::plot_grid(space_herb_plot, space_nutrient_plot, 
                                   mrs_herb_plot, mrs_nutrient_plot, 
                                   dist_herb_plot, dist_nutrient_plot, 
                                   align="vh", ncol=2, 
                                   labels=c("(a)","(b)","(c)","(d)","(e)","(f)"))

six_panel_legend<-cowplot::plot_grid(space_herb_legend, space_nutrient_legend, align="vh", ncol=2)
six_panel_plot_with_legend<-cowplot::plot_grid(six_panel_plot, six_panel_legend, align="vh", ncol=1, rel_heights = c(5, 0.5))
six_panel_plot_with_legend
ggsave("figures/six_panel_plot_with_legend_median_v2.pdf", width=7.2, height=10, units="in")


# ----------------------------------------------------------------------------------#
# ----------------- Ratio Coral to Macroalgae through time--------------------------
# ----------------------------------------------------------------------------------#

### simple stats to go with this plot
# number of times poisitive or negative
ratio_sums<-ddply(ranks_plots, .(Herb_Trt), summarize,
                  perc_coral=(sum(binary_coral)/length(binary_coral))*100,
                  perc_macro=(sum(binary_macro)/length(binary_macro))*100)

#  note to self need to subset to T11 becuase that was the starting timpeoint. the way you set up this df there is no t12
ratio_sums_t12<-ddply(subset(ranks_plots, Timepoint=="T11"), .(Herb_Trt), summarize,
                      perc_coral=(sum(binary_coral)/length(binary_coral))*100,
                      perc_macro=(sum(binary_macro)/length(binary_macro))*100)

# Plot Figure 6 -----------------------------------------------------------------------
coral_macro_ratio_time<-ggplot(ranks_plots, aes(x=Date, y=log_ratio_coral_macro, group=Block_Plot_Herb_Trt, color=Herb_Trt, shape=Nutrient_Trt))+
  geom_hline(yintercept=0, linetype="dashed", color = "gray")+
  geom_point(aes(color=Herb_Trt), size=2, stroke = 1, alpha=0.8)+
  geom_line(data=ranks_plots, aes(x=Date, y=log_ratio_coral_macro, color=Herb_Trt, group=Block_Plot_Herb_Trt), size=0.5, alpha=0.8)+
  scale_shape_manual(values=c(1,2),  name="Nutrients")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumer Pressure", labels = c("High","Medium", "Low", "Very low" ))+
  facet_wrap(~Herb_Trt, labeller = as_labeller(herb_names), nrow=2)+ 
  xlab("")+
  ylab("log(Coral:Macroalgae)")+
  coord_fixed()+
  theme_minimal()+
  theme(panel.border = element_rect(color = "gray 50", fill = NA))+
  theme(strip.text.x = element_text(size=12, face="bold", color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 4/6)
coral_macro_ratio_time
ggsave("figures/coral_macro_ratio_time_v2.pdf", width=9, height=5, units="in")


