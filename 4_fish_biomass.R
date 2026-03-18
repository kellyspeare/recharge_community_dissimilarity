# ----------------------------------------------------------------------- #

# Reducing consumer pressure increases community dissimilarity leading to 
# coral- or algal- dominance on a coral reef
#  
# 4_fish_biomass

# ----------------------------------------------------------------------- #

# this script analyzes data on fish abundance and corallivory in experimental plots
# Plots Figures S2, S3

# Packages --------------------------------------------------------------- #

library(ggplot2)
library(plyr)
library(cowplot)
library(multcomp)
library(glmmTMB)

# data --------------------------------------------------------------------

fish<-read.csv("data/data_clean/fish_biomass.csv", header=TRUE, stringsAsFactors = TRUE)

#how many unique survey dates?
unique(fish$monthyr) #8

# labeller for plots
nutrient_names <- c(
  `ambient` = "Ambient Nutrients",
  `enriched` = "Enriched Nutrients"
)

# ordering facotr levels
fish$consumers<-factor(fish$consumers, levels=c("High", "Medium", "Low", "Very Low"))

#calculating mean biomass by functional group
fish_means<-ddply(fish, .(consumers, nutrients, functionalgroup), summarize,
                  meanbiomass=mean(biomass),
                  sebiomass=sd(biomass)/sqrt(length(biomass)))
# calculating total biomass for herbivores only
herbivore_totals<-ddply(subset(fish, functionalgroup!="Corallivore"), .(monthyr, blockplot, nutrients, consumers), summarize,
                  totalbiomass=sum(biomass))
# calculating mean biomass for herbivores only
herbivore_means<-ddply(herbivore_totals, .(nutrients, consumers), summarize,
                       meanbiomass=mean(totalbiomass),
                       sebiomass=sd(totalbiomass)/sqrt(length(totalbiomass)))

# model herbivores --------------------------------------------------------------

herb.mod<-lmer(log(totalbiomass+1) ~ consumers*nutrients + (1|blockplot), data=herbivore_totals)
car::Anova(herb.mod)
plot(simulateResiduals(fittedModel=herb.mod, quantileFunction = qnorm))


# plot herbivores by consumer pressure ----------------------------------------- 

herbivore_means_no_nutrients<-ddply(herbivore_totals, .(consumers), summarize,
                       meanbiomass=mean(totalbiomass),
                       sebiomass=sd(totalbiomass)/sqrt(length(totalbiomass)))

fish_means_no_nutrients<-ddply(fish, .(consumers, functionalgroup), summarize,
                  meanbiomass=mean(biomass),
                  sebiomass=sd(biomass)/sqrt(length(biomass)))

# Tukey Post hoc
herb.mod.tukey<-emmeans(herb.mod, type="response", ~consumers)
herb.cld<-multcomp::cld(herb.mod.tukey, Letters=letters)
names(herb.cld)[names(herb.cld)==".group"]<-"group"

herbivore_means<-left_join(herbivore_means, herb.cld[c(1,2,7)], by=c("consumers"))

herbivore_plot<-ggplot()+
  geom_bar(data=subset(fish_means_no_nutrients, functionalgroup!="Corallivore"), aes(x=consumers, y=meanbiomass, fill=functionalgroup), stat="identity")+
  geom_errorbar(data=herbivore_means_no_nutrients, aes(x=consumers, y=meanbiomass, ymin=meanbiomass-sebiomass, ymax=meanbiomass+sebiomass), width = .25)+
  #geom_text(data=herbivore_means, aes(label = group, x=consumers, y= meanbiomass+sebiomass), vjust =-2) +
  scale_fill_manual(values=c("#ffd166","#06d6a0","#118ab2","#073b4c"), name = "Functional group")+
  #facet_wrap(~nutrients, labeller = as_labeller(nutrient_names))+
  xlab("Consumer pressure")+
  ylab("Mean biomass (g/exclosure)")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  theme(strip.text.x = element_text(size=12, color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
herbivore_plot

# model corallivores ------------------------------------------------------------
corallivore_totals<-subset(fish, functionalgroup=="Corallivore")

corallivore.mod<-lmer(biomass ~ consumers*nutrients + (1|blockplot), data=corallivore_totals)
car::Anova(corallivore.mod)
plot(simulateResiduals(fittedModel=corallivore.mod, quantileFunction = qnorm))

corallivore.mod.2 <-lmer(log(biomass+1) ~ consumers*nutrients + (1|blockplot), data=corallivore_totals)
car::Anova(corallivore.mod.2)
plot(allEffects(corallivore.mod.2))
plot(simulateResiduals(fittedModel=corallivore.mod.2, quantileFunction = qnorm))


# plot corallivores by consumer pressure ---------------------------------------- 

corallivore_means_no_nutrients<-ddply(corallivore_totals, .(consumers), summarize,
                                    meanbiomass=mean(biomass),
                                    sebiomass=sd(biomass)/sqrt(length(biomass)))


corallivore.mod.tukey<-emmeans(corallivore.mod, type="response", ~consumers)
corallivore.cld<-multcomp::cld(corallivore.mod.tukey, Letters=letters)
names(corallivore.cld)[names(corallivore.cld)==".group"]<-"group"

corallivore_means<-left_join(corallivore_means_no_nutrients, corallivore.cld[c(1,2,7)], by=c("consumers"))

# plot 
corallivore_plot<-ggplot()+
  geom_bar(data=subset(fish_means_no_nutrients, functionalgroup=="Corallivore"), aes(x=consumers, y=meanbiomass, fill=functionalgroup), stat="identity")+
  geom_errorbar(data=corallivore_means_no_nutrients, aes(x=consumers, y=meanbiomass, ymin=meanbiomass-sebiomass, ymax=meanbiomass+sebiomass), width = .25)+
  #geom_text(data=corallivore_means, aes(label = group, x=consumers, y= meanbiomass+sebiomass), vjust =-2) +
  scale_fill_manual(values=c("#ef476f"), name = "Functional group")+
  #facet_wrap(~nutrients, labeller = as_labeller(nutrient_names))+
  xlab("Consumer pressure")+
  ylab("Mean biomass (g/exclosure)")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  theme(strip.text.x = element_text(size=12, color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
corallivore_plot

# plot Figure S2 ----------------------------------------------------------------- 

## cowplot 
herbivore_corallivore_plot<-cowplot::plot_grid(herbivore_plot, corallivore_plot, nrow=1, ncol=2, align="vh",  labels=c("(a)","(b)"))
herbivore_corallivore_plot
ggsave("figures/herbivore_corallivore_plot.pdf", width=10, height=4, units="in")


# Fish biomass by species -------------------------------------- 

fish_sp<-read.csv("data/data_clean/fish_biomass_species.csv", header = TRUE, stringsAsFactors = TRUE)

#reordering factor levels
fish_sp$consumers<-factor(fish_sp$consumers, levels=c("High", "Medium", "Low", "Very Low"))

fish_sp_mean<-ddply(fish_sp, .(consumers, functionalgroup, genusspecies), summarize,
      meanbiomass=mean(biomass))

#herbivores
herbivore_sp<-ggplot()+
  geom_bar(data=subset(fish_sp_mean, functionalgroup!="Corallivore"), aes(x=functionalgroup, y=meanbiomass, fill=genusspecies), stat="identity")+
  #geom_errorbar(data=herbivore_means_no_nutrients, aes(x=consumers, y=meanbiomass, ymin=meanbiomass-sebiomass, ymax=meanbiomass+sebiomass), width = .25)+
  #geom_text(data=herbivore_means, aes(label = group, x=consumers, y= meanbiomass+sebiomass), vjust =-2) +
  #scale_fill_manual(values=c("#ffd166","#06d6a0","#118ab2","#073b4c"), name = "Functional group")+
  facet_wrap(~consumers, nrow=1)+
  labs(fill="Species")+
  guides(fill=guide_legend(ncol = 2))+
  xlab("Functional Group")+
  ylab("Mean biomass (g/exclosure)")+
  theme_minimal()+
  theme(legend.position = "right")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  theme(strip.text.x = element_text(size=12, color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
herbivore_sp

#corallivores
corallivore_sp<-ggplot()+
  geom_bar(data=subset(fish_sp_mean, functionalgroup=="Corallivore"), aes(x=functionalgroup, y=meanbiomass, fill=genusspecies), stat="identity")+
  #geom_errorbar(data=herbivore_means_no_nutrients, aes(x=consumers, y=meanbiomass, ymin=meanbiomass-sebiomass, ymax=meanbiomass+sebiomass), width = .25)+
  #geom_text(data=herbivore_means, aes(label = group, x=consumers, y= meanbiomass+sebiomass), vjust =-2) +
  #scale_fill_manual(values=c("#ffd166","#06d6a0","#118ab2","#073b4c"), name = "Functional group")+
  facet_wrap(~consumers, nrow=1)+
  labs(fill="Species")+
  #guides(fill=guide_legend(ncol = 2))+
  xlab("Functional Group")+
  ylab("Mean biomass (g/exclosure)")+
  theme_minimal()+
  theme(legend.position = "right")+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  theme(strip.text.x = element_text(size=12, color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
corallivore_sp
ggsave("figures/corallivore_species_plot.pdf", width=12, height=6, units="in")


# plot Figure S4 ----------------------------------------------------------------- 

## cowplot 
herbivore_corallivore_species_plot<-cowplot::plot_grid(herbivore_sp, corallivore_sp, nrow=2, ncol=1, align="vh",  labels=c("(a)","(b)"))
herbivore_corallivore_species_plot
ggsave("figures/herbivore_corallivore_species_plot.pdf", width=12, height=8.5, units="in")

# Bites on corals - a proxy for corallivory -------------------------------------- 

# read data
bites<-read.csv("data/data_clean/bites_on_corals.csv", header = TRUE, stringsAsFactors = TRUE)

bites_l<-pivot_longer(bites, cols=Bites_T10:Bites_T12, names_to = "Timepoints", values_to = "bites")


bites.mod<-glmmTMB(bites ~ Herb_Trt*Nutrient_Trt+ (1|Block_Plot/Coral), data=bites_l, family=binomial(link = "logit"))
car::Anova(bites.mod)
plot(simulateResiduals(fittedModel=bites.mod, quantileFunction = qnorm))

bites_predicted<-data.frame(emmeans(bites.mod, type="response", ~Herb_Trt))

bites_predicted$Herb_Trt<-factor(bites_predicted$Herb_Trt, levels=c("open", "3X3", "2X2", "1X1"))

bites_predicted$consumers <- as.factor(ifelse(bites_predicted$Herb_Trt=="open", "High", 
                                              ifelse(bites_predicted$Herb_Trt=="3X3" , "Medium",
                                                     ifelse(bites_predicted$Herb_Trt=="2X2", "Low", 
                                                            ifelse(bites_predicted$Herb_Trt== "1X1", "Very Low", "NA")))))
bites_predicted$consumers<-factor(bites_predicted$consumers, levels=c("High", "Medium", "Low", "Very Low"))

# plot Figure 1 ----------------------------------------------------------------- 

#combine into one dataframe
herbs_bites<-left_join(herbivore_means_no_nutrients, bites_predicted, by="consumers")

herbs_bites$prob_max<-herbs_bites$prob + herbs_bites$SE
herbs_bites$prob_min<-herbs_bites$prob - herbs_bites$SE
herbs_bites$meanbiomass_max<-herbs_bites$meanbiomass + herbs_bites$sebiomass
herbs_bites$meanbiomass_min<-herbs_bites$meanbiomass - herbs_bites$sebiomass

herbs_bites$herbs<-as.factor("herbivores")
herbs_bites$corallivores<-as.factor("corallivores")

# Calculate the scaling factor
# This ensures both variables occupy roughly the same vertical space
scale_factor <- max(herbs_bites$meanbiomass) / max(herbs_bites$prob)
scale_factor

#Create the dual-axis plot
herb_bites_plot<-ggplot(herbs_bites, aes(x = consumers)) +
  # Plot energy consumption as columns (using the scaled data)
  geom_point(aes(y = prob*scale_factor), color = "#005AB5", size=3, stroke=1, pch=1) +
  geom_errorbar(aes(ymax=prob_max*scale_factor, ymin=prob_min*scale_factor), color = "#005AB5", width=0.2) +
  geom_line(aes(y=prob*scale_factor, group=corallivores), color="#005AB5")+
  # Plot temperature as a line
  geom_point(aes(y = meanbiomass), color = "#DC3220", size = 3, stroke=1,pch=2) +
  geom_errorbar(aes(ymax=meanbiomass_max, ymin=meanbiomass_min), color = "#DC3220", width=0.2) +
  geom_line(aes(y=meanbiomass, group=herbs), color="#DC3220")+
  # Define the scales for the y-axes
  scale_y_continuous(
    # Name for the primary (left) axis
    name = "Herbivore biomass",
    # Define the secondary (right) axis
    sec.axis = sec_axis(
      # The inverse transformation to correctly label the secondary axis
      trans = ~ . / scale_factor,
      # Name for the secondary (right) axis
      name = "Corallivory\n (proportion of corals bitten)")) +
  # Add labels and a theme for better readability
  #labs(title = "")+
  xlab("Consumer Pressure")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "black", fill = NA))+
  theme(strip.text.x = element_text(size=12, color="black"))+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)
herb_bites_plot
ggsave("figures/herb_bites_plot.pdf", width=4, height=4, units="in")

