# ----------------------------------------------------------------------- #
# Reducing consumer pressure increases community dissimilarity leading to 
# coral- or algal- dominance on a coral reef
 
# 1_coral_macroalgae_patterns

# ----------------------------------------------------------------------- #
# this script analyzes data on the abundance of corals and algae in plots
# Plots Figures 1b, 1c, 1d, 5

# Packages --------------------------------------------------------------- #

library(cowplot)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(tidyr)
library(sjPlot)
library(DHARMa)
library(car)
library(effects)
library(multcomp)

# ----------------------------------------------------------------------------------#
# -- Load Data ----------------------------------------------------------------------
# ----------------------------------------------------------------------------------#


cover<-read.csv("data/data_clean/plots_benthic_cover.csv", header=TRUE, stringsAsFactors = TRUE)

#reordering factor levels
cover$Herb_Trt<-factor(cover$Herb_Trt, levels=c("Open", "3X3", "2X2", "1X1"))


# ----------------------------------------------------------------------------------#
# -- calculating range of abundance of corals and macroalgae, rate of change -------
# ----------------------------------------------------------------------------------#

cover_all<-cover

t12_coral_macro_sums<- cover_all %>% subset(Timepoint=="T12") %>%
  group_by(Herb_Trt) %>%
  summarize(high_coral = max(Sum_Coral),
            low_coral= min(Sum_Coral),
            mean_coral = mean(Sum_Coral),
            median_coral=median(Sum_Coral),
            high_del_coral=high_coral/4,
            low_del_coral=low_coral/4,
            mean_del_coral=mean_coral/4,
            median_del_coral=median_coral/4,
            high_macro=max(Sum_Macroalgae),
            low_macro=min(Sum_Macroalgae),
            mean_macro=mean(Sum_Macroalgae),
            median_macro=median(Sum_Macroalgae),
            high_del_macro=high_macro/4,
            low_del_macro=low_macro/4,
            mean_del_macro=mean_macro/4,
            median_del_macro=median_macro/4) %>%
  ungroup()


# calculating rates the correct way for individual plots
t12_coral_macro<-subset(cover_all[,c(1:12)], Timepoint=="T12")

t0_coral_macro<-subset(cover_all[,c(1:12)], Timepoint=="T0")

coral_macro_combined<-rbind(t12_coral_macro, t0_coral_macro)
coral_macro_combined_w<-pivot_wider(coral_macro_combined, id_cols=c(Block, Block_Plot_Herb_Trt, Nutrient_Trt, Herb_Trt), names_from =Timepoint, values_from = c(Sum_Coral, Sum_Macroalgae))

#change in coral cover
coral_macro_combined_w$coral_change<-coral_macro_combined_w$Sum_Coral_T12 - coral_macro_combined_w$Sum_Coral_T0
coral_macro_combined_w$coral_change_rate<-coral_macro_combined_w$coral_change/4

#change in macroalgae cover
coral_macro_combined_w$macro_change<-coral_macro_combined_w$Sum_Macroalgae_T12 - coral_macro_combined_w$Sum_Macroalgae_T0
coral_macro_combined_w$macro_change_rate<-coral_macro_combined_w$macro_change/4

# means by herb treatment
coral_macro_change_means_herbs<-ddply(coral_macro_combined_w, .(Herb_Trt), summarize,
                                      mean_coral_change=mean(coral_change),
                                      se_coral_change=sd(coral_change)/sqrt(length(coral_change)),
                                      mean_coral_change_rate=mean(coral_change_rate),
                                      se_coral_change_rate=sd(coral_change_rate)/sqrt(length(coral_change_rate)),
                                      mean_macro_change=mean(macro_change),
                                      se_macro_change=sd(macro_change)/sqrt(length(macro_change)),
                                      mean_macro_change_rate=mean(macro_change_rate),
                                      se_macro_change_rate=sd(macro_change_rate)/sqrt(length(macro_change_rate)))

# means by nutrient treatment
coral_macro_change_means_nuts<-ddply(coral_macro_combined_w, .(Nutrient_Trt), summarize,
                                     mean_coral_change=mean(coral_change),
                                     se_coral_change=sd(coral_change)/sqrt(length(coral_change)),
                                     mean_coral_change_rate=mean(coral_change_rate),
                                     se_coral_change_rate=sd(coral_change_rate)/sqrt(length(coral_change_rate)),
                                     mean_macro_change=mean(macro_change),
                                     se_macro_change=sd(macro_change)/sqrt(length(macro_change)),
                                     mean_macro_change_rate=mean(macro_change_rate),
                                     se_macro_change_rate=sd(macro_change_rate)/sqrt(length(macro_change_rate)))

# ----------------------------------------------------------------------------------#
# ---------- rate of change mixed effects model ------------------------------------
# ----------------------------------------------------------------------------------#

# rate of change of coral cover
coral_change_mod.1<-lmer(coral_change_rate ~ Herb_Trt*Nutrient_Trt + (1|Block), data=coral_macro_combined_w)
car::Anova(coral_change_mod.1)

# rate of change of macroalgae cover
macro_change_mod.1<-lmer(macro_change_rate ~ Herb_Trt*Nutrient_Trt + (1|Block), data=coral_macro_combined_w)
car::Anova(macro_change_mod.1)


# ----------------------------------------------------------------------------------#
# ---------- plot rate of change  ------------------------------------
# ----------------------------------------------------------------------------------#

coral_macro_combined_l<-coral_macro_combined_w %>% 
  select(Block, Block_Plot_Herb_Trt, Herb_Trt, Nutrient_Trt, coral_change_rate, macro_change_rate) %>% 
  pivot_longer(cols=c(coral_change_rate, macro_change_rate), names_to="group", values_to = "change_rate")

# data summary function that makes mean points for plots #
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-(sd(x)/sqrt(length(x)))
  ymax <- m+(sd(x)/sqrt(length(x)))
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# labeller for plots
group_names <- c(
  `coral_change_rate` = "Coral",
  `macro_change_rate` = "Macroalgae"
)
 
# Plot Figure 2b --------------------------------------------------------------
coral_macro_change_plot<-ggplot(coral_macro_combined_l, aes(x=Herb_Trt, y=change_rate))+
  geom_jitter(size=2.5, stroke=1, aes(color=Herb_Trt, shape=Nutrient_Trt), width=0.2, alpha=1)+
  stat_summary(fun.data=data_summary, color="black", size=1, shape=16)+
  scale_color_manual(values=c( "#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumer\npressure", labels = c("High","Medium", "Low", "Very low" ))+
  scale_x_discrete(name ="Consumer pressure", labels=c("High","Medium", "Low", "Very low"))+
  scale_shape_manual(values=c(1,2), name ="Nutrients")+
  facet_wrap(~group, ncol=4, labeller = as_labeller(group_names),)+
  theme_classic()+
  ylab(expression(paste("% change ", year^-1)))+
  xlab("")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "gray 50", fill = NA))+
  theme(strip.text.x = element_text(size=16, face="bold", color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 1/1)
coral_macro_change_plot
#ggsave("figures/coral_macro_change_plot.pdf",width=6.35, height=4,units="in")

# ----------------------------------------------------------------------------------#
# --------------------- biplot of coral and macroalgae ------------------------------
# ----------------------------------------------------------------------------------#

# Plot Figure 2c --------------------------------------------------------------------
coral_macro_biplot<-ggplot(cover, aes(x=Sum_Coral, y=Sum_Macroalgae, color=Date, shape=Nutrient_Trt, group=Block_Plot_Herb_Trt))+
  geom_point(size=2, stroke=1)+
  scale_color_manual(name="Timepoint", values=c("#440154","#482475","#414487","#355f8d", "#2a788e","#21918c","#22a884","#44bf70","#7ad151","#bddf26","#fde725"), 
                     labels=c("Jul 2018", "Nov 2018", "Mar 2019", "Jul 2019", "Nov 2019", "Aug 2020", "May 2021", "Jul 2021", "Nov 2021", "Apr 2022", "Jul 2022"))+
  scale_shape_manual(values=c(1,2),  name="Nutrients", guide="none")+
  facet_wrap(~Herb_Trt, labeller = as_labeller(herb_names), nrow=1)+ 
  xlim(c(0,100))+
  ylim(c(0,100))+
  ylab("Abundance of Macroalgae \n(% cover)")+
  xlab("Abundance of Coral (% cover)")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "gray 50", fill = NA))+
  theme(strip.text.x = element_text(size=14, face="bold", color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 1/1)
coral_macro_biplot
#ggsave("figures/coral_macro_biplot.pdf", width=10, height=3.5, units="in")

# cowplot::plot_grid(coral_macro_change_plot + theme(legend.justification = c(0,1)),
#                    coral_macro_biplot + theme(legend.justification = c(0,1)), 
#                    nrow=2, align = c("vh"), axis = "bt",labels = c("(a)","(b)"))

# ----------------------------------------------------------------------------------#
# ---------------- plots of timeseries of coral and macroalgae ----------------------
# ----------------------------------------------------------------------------------#

cover$Date<-as.Date(cover$Date, "%m/%d/%Y")

cover_l<-cover[,c(1:12)] %>% pivot_longer(cols=11:12, names_to="category", values_to="percent")
cover_l$category<-as.factor(cover_l$category)

#calculate mean cover of coral and macroalgae at each timepoint
cover_means<-ddply(cover_l, .(Herb_Trt, Date, Timepoint, category), summarize,
                   mean=mean(percent),
                   se=sd(percent)/sqrt(length(percent)))

# Plot Figure 2d --------------------------------------------------------------

# labeller for plots
category_names <- c(
  `Sum_Coral` = "Coral",
  `Sum_Macroalgae` = "Macroalgae"
)

coral_algae_timeseries<-ggplot(cover_means, aes(x=Date, y=mean, group=Herb_Trt, color=Herb_Trt))+
  geom_jitter(data=cover_l, aes(x=Date, y=percent, color=Herb_Trt, shape=Nutrient_Trt), size=2, stroke = 1, alpha=0.5)+
  geom_point(aes(color=Herb_Trt), size=2, stroke = 1, alpha=0.8)+
  geom_line(data=cover_means, aes(x=Date, y=mean, color=Herb_Trt, group=Herb_Trt), size=0.5, alpha=0.8)+
  scale_shape_manual(values=c(1,2),  name="Nutrients")+
  scale_color_manual(values=c("#0d0887","#7e03a8",  "#cc4778", "#f89540"),
                     name = "Consumer\npressure", labels = c("High","Medium", "Low", "Very low" ))+
  facet_wrap(~category, nrow=1, labeller = as_labeller(category_names))+ 
  xlab("")+
  ylab("Abundance\n(% cover)")+
  coord_fixed()+
  theme_minimal()+
  theme(panel.border = element_rect(color = "gray 50", fill = NA))+
  theme(strip.text.x = element_text(size=14, face="bold", color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 4/6)
coral_algae_timeseries
ggsave("figures/coral_algae_timeseries.pdf", width=10, height=3.5, units="in")




# simple stats coral and macro abundance at start ------------------------------

start_abund_sums<-ddply(subset(cover, Timepoint=="T0"), .(), summarize,
                        mean_coral=mean(Sum_Coral),
                        min_coral=min(Sum_Coral),
                        max_coral=max(Sum_Coral),
                        mean_macro=mean(Sum_Macroalgae),
                        min_macro=min(Sum_Macroalgae),
                        max_macro=max(Sum_Macroalgae))

# simple stats coral and macro abundance at end ------------------------------

end_abund_sums_herb<-ddply(subset(cover, Timepoint=="T12"), .(Herb_Trt), summarize,
                           mean_coral=mean(Sum_Coral),
                           se_coral=sd(Sum_Coral)/sqrt(length(Sum_Coral)),
                           min_coral=min(Sum_Coral),
                           max_coral=max(Sum_Coral),
                           mean_macro=mean(Sum_Macroalgae),
                           se_macro=sd(Sum_Macroalgae)/sqrt(length(Sum_Macroalgae)),
                           min_macro=min(Sum_Macroalgae),
                           max_macro=max(Sum_Macroalgae))

end_abund_sums_nuts<-ddply(subset(cover, Timepoint=="T12"), .(Nutrient_Trt), summarize,
                           mean_coral=mean(Sum_Coral),
                           se_coral=sd(Sum_Coral)/sqrt(length(Sum_Coral)),
                           min_coral=min(Sum_Coral),
                           max_coral=max(Sum_Coral),
                           mean_macro=mean(Sum_Macroalgae),
                           se_macro=sd(Sum_Macroalgae)/sqrt(length(Sum_Macroalgae)),
                           min_macro=min(Sum_Macroalgae),
                           max_macro=max(Sum_Macroalgae))
