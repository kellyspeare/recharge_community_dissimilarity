
## packages
# load data

# Nitrate + Nitrite - analyzed through time ----------------------------------------------------


nut.mod.1<-lmer(Nitrite_plus_Nitrate ~ Nutrient_Trt*Days_since_nutrient_swap +(1|Sample_period), data=water_nutrients)
car::Anova(nut.mod.1)
plot(allEffects(nut.mod.1))
plot(simulateResiduals(fittedModel=nut.mod, quantileFunction = qnorm))

nut.mod.2<-lm(Nitrite_plus_Nitrate ~ Nutrient_Trt*Days_since_nutrient_swap, data=water_nutrients)
car::Anova(nut.mod)
plot(allEffects(nut.mod))
plot(simulateResiduals(fittedModel=nut.mod, quantileFunction = qnorm))

# including sampling timepoint as a random effect
# Response: Nitrite_plus_Nitrate
# Chisq Df Pr(>Chisq)    
# Nutrient_Trt                          20.9608  1  4.688e-06 ***
# Days_since_nutrient_swap               1.9252  1   0.165287    
# Nutrient_Trt:Days_since_nutrient_swap  7.4436  1   0.006366 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# no random effects
# Sum Sq Df F value    Pr(>F)    
# Nutrient_Trt                          0.35656  1 12.5841 0.0007159 ***
# Days_since_nutrient_swap              0.24934  1  8.8000 0.0041721 ** 
# Nutrient_Trt:Days_since_nutrient_swap 0.11512  1  4.0628 0.0478510 *  
# Residuals                             1.89839 67                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

nut_ref<-ref_grid(nut.mod.1, at=list(Days_since_nutrient_swap=c(0,25,50,75,100,125,150,175)))
nut_pred<-data.frame(emmip(nut_ref, Nutrient_Trt ~ Days_since_nutrient_swap, cov.reduce=range, plotit=FALSE)) #look at them visually

nit_plot<-ggplot()+
  geom_point(data=water_nutrients, aes(x=Days_since_nutrient_swap, y=Nitrite_plus_Nitrate, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 5),  size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=nut_pred, aes(x=Days_since_nutrient_swap, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=nut_pred, aes(x=Days_since_nutrient_swap, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Days since diffuser swap")+
  ylab(expression(paste("Nitrate + Nitrite (", mu, "M)", sep="")))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)

nut.mod.2<-lm(Nitrite_plus_Nitrate ~ Nutrient_Trt*Days_since_nutrient_swap, data=water_nutrients)
car::Anova(nut.mod)

# Nitrate + Nitrite - end of deployment data ----------------------------------------------------

nut.mod.3<-lm(Nitrite_plus_Nitrate ~ Nutrient_Trt, data=subset(water_nutrients, Days_since_nutrient_swap >80))
car::Anova(nut.mod.3)

# Response: Nitrite_plus_Nitrate
#                  Sum Sq Df F value Pr(>F)
#Nutrient_Trt 0.035623  1  3.7943 0.06117 .
#Residuals    0.272273 29   

# Phosphate --------------------------------------------------------------------

phos.mod.1<-lmer(Phosphate ~ Nutrient_Trt*Days_since_nutrient_swap +(1|Sample_period), data=water_nutrients)
car::Anova(phos.mod.1)
summary(phos.mod.1)
plot(allEffects(phos.mod.1))
plot(simulateResiduals(fittedModel=phos.mod.1, quantileFunction = qnorm))

# Chisq Df Pr(>Chisq)    
# Nutrient_Trt                          12.2211  1  0.0004725 ***
#   Days_since_nutrient_swap               0.3707  1  0.5426012    
# Nutrient_Trt:Days_since_nutrient_swap  2.2701  1  0.1318942    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


phos.mod.2<-lm(Phosphate ~ Nutrient_Trt*Days_since_nutrient_swap, data=water_nutrients)
car::Anova(phos.mod.2)
summary(phos.mod.2)
plot(allEffects(phos.mod.2))
plot(simulateResiduals(fittedModel=phos.mod.2, quantileFunction = qnorm))



phos_ref<-ref_grid(phos.mod.1, at=list(Days_since_nutrient_swap=c(0,25,50,75,100,125,150,175)))
phos_pred<-data.frame(emmip(phos_ref, Nutrient_Trt ~ Days_since_nutrient_swap, cov.reduce=range, plotit=FALSE)) #look at them visually

phos_plot<-ggplot()+
  geom_point(data=water_nutrients, aes(x=Days_since_nutrient_swap, y=Phosphate, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 5),  size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  #ylim(c(-0.01,.25))+
  geom_ribbon(data=phos_pred, aes(x=Days_since_nutrient_swap, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=phos_pred, aes(x=Days_since_nutrient_swap, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Days since diffuser swap")+
  ylab(expression(paste("Phosphorus (", mu, "M)", sep="")))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)

plot_legend<-ggplot()+
  geom_point(data=water_nutrients, aes(x=Days_since_nutrient_swap, y=Phosphate, color=Nutrient_Trt, group=Nutrient_Trt, shape=Nutrient_Trt), 
             position=position_jitter(w = 5),  size=2, stroke = 1, alpha=0.5)+
  scale_shape_manual(values=c(1,2), name="Nutrients")+
  geom_ribbon(data=phos_pred, aes(x=Days_since_nutrient_swap, ymin=yvar-SE, ymax=yvar+SE, fill=Nutrient_Trt),alpha= 0.3) +
  geom_line(data=phos_pred, aes(x=Days_since_nutrient_swap, y = yvar, color=Nutrient_Trt), size = 0.9)+
  geom_smooth(method="lm")+
  scale_color_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  scale_fill_manual(values=c("#A9C837", "#37C89E"), name = "Nutrients", labels = c("Ambient","Enriched"))+
  xlab("Days since diffuser swap")+
  ylab(expression(paste("Phosphorus (", mu, "M)", sep="")))+
  theme_bw()+
  theme(legend.position = "right")+
  theme(axis.text=element_text(colour="black", size=12), axis.title = element_text(size=14))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.border =  element_rect(colour="black",size=0.75), axis.ticks=element_line(color="black"))+
  theme(aspect.ratio = 6/6)

legend<-get_legend(plot_legend)

# Phosphate - end of deployment data ----------------------------------------------------

phos.mod.3<-lm(Phosphate ~ Nutrient_Trt, data=subset(water_nutrients, Days_since_nutrient_swap >80))
car::Anova(phos.mod.3)
# 
# Response: Phosphate
#                   Sum Sq Df F value Pr(>F)
# Nutrient_Trt 0.0019923  1  2.8564 0.1017
# Residuals    0.0202271 29      


# combine plots ----------------------------------------------------------------

nit_phos_plot<-cowplot::plot_grid(nit_plot, phos_plot, 
                   align="hv", ncol=2, 
                   axis = "bt",
                   labels=c("(a)","(b)"))
nit_phos_fig<-cowplot::plot_grid(nit_phos_plot, legend, align="hv", ncol=2, 
                   axis = "bt",
                   rel_widths = c(3,.5) )
ggsave("figures/nit_phos_fig.pdf", units="in", width=9, height=4)
