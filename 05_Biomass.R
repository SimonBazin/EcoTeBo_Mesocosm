################### Biomass  ###################
################################################
setwd("C:/Users/sibazin/Desktop/Simon/Doctorat/chapitre 5/article/Data")
#------------------------------------------------------- 
# This script is used to analyze community biomass data  

# 1) Full model
# 2) Macroinvertebrate
# 3) Zooplankton
# 4) Phytoplankton

#------------------------------------------------------- 
# Packages and data

# Packages
library(car)
library(lmerTest)
library(ggplot2)
library(plyr)

# Data
df_taxa<-read.table(file="df_taxa.csv",header=T,sep=";",dec=".")

# Rename
df_taxa$compartment <- factor(df_taxa$compartment, levels = c("inv", "zoo", "phyto"), 
                              labels = c("Invertebrate", "Zooplankton", "Phytoplankton"))

# Mesocosm_ID as character
df_taxa$mesocosm_ID<-as.character(df_taxa$mesocosm_ID)
# Mesocosm metadata
m<-aggregate(data=df_taxa, cbind(fish,warming) ~ mesocosm_ID + date + tmean15, unique)

# Biomass
b.full<-ddply(df_taxa, .(date,mesocosm_ID,fish,warming,compartment), summarise,
              C_mug_l = sum(C_mug_l),
              log.C_mug_l= log(sum(C_mug_l)),
              tmean15 = unique(tmean15))

b.full[b.full$compartment=="Invertebrate",]$log.C_mug_l<-b.full[b.full$compartment=="Invertebrate",]$log.C_mug_l-mean(b.full[b.full$compartment=="Invertebrate",]$log.C_mug_l)
b.full[b.full$compartment=="Zooplankton",]$log.C_mug_l<-b.full[b.full$compartment=="Zooplankton",]$log.C_mug_l-mean(b.full[b.full$compartment=="Zooplankton",]$log.C_mug_l)
b.full[b.full$compartment=="Phytoplankton",]$log.C_mug_l<-b.full[b.full$compartment=="Phytoplankton",]$log.C_mug_l-mean(b.full[b.full$compartment=="Phytoplankton",]$log.C_mug_l)


#------------------------------------------------------- 
# Full model

# Quadratic model
mod<-lmer(log(C_mug_l)~(tmean15+I(tmean15^2))*fish*warming*compartment+(1|date)+(1|mesocosm_ID),data=b.full,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(mod)
# Stepwise backward regression
step_res <- step(mod,alpha.random = 1)
final <- get_model(step_res)
Anova(final)
summary(final)
res_lme=residuals(final)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)

#------------------------------------------------------- 
# MacroInvertebrate compartment

Invertebrate<-b.full[b.full$compartment=="Invertebrate",]

# Quadratic model
mod<-lmer(log(C_mug_l)~(tmean15+I(tmean15^2))*warming+(1|date)+(1|mesocosm_ID),data=Invertebrate,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(mod)
# Stepwise backward regression
step_res <- step(mod,alpha.random = 1)
final.inv <- get_model(step_res)
Anova(final.inv)
summary(final.inv)
res_lme=residuals(final.inv)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)


#------------------------------------------------------- 
# Zooplankton compartment

Zooplankton<-b.full[b.full$compartment=="Zooplankton",]

# Quadratic model
mod<-lmer(log(C_mug_l)~(tmean15+I(tmean15^2))*warming+(1|date)+(1|mesocosm_ID),data=Zooplankton,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(mod)
# Stepwise backward regression
step_res <- step(mod,alpha.random = 1)
final.zoo <- get_model(step_res)
Anova(final.zoo)
summary(final.zoo)
res_lme=residuals(final.zoo)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)


#------------------------------------------------------- 
# Phytoplankton compartment

Phytoplankton<-b.full[b.full$compartment=="Phytoplankton",]

# Quadratic model
mod<-lmer(log(C_mug_l)~(tmean15+I(tmean15^2))*warming+(1|date)+(1|mesocosm_ID),data=Phytoplankton,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(mod)
summary(mod)
# Stepwise
step_res <- step(mod,alpha.random = 1)
final.phyto <- get_model(step_res)
Anova(final.phyto)
summary(final.phyto)
res_lme=residuals(final.phyto)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)



# Mean and CI for plot
df_plot<-ddply(b.full, .(compartment,date,warming), summarise,
               m.C_mug_l            = mean(log(C_mug_l),na.rm=T),
               CI.C_mug_l          = 1.96*sd(log(C_mug_l), na.rm=TRUE)/sqrt(length(na.omit(log(C_mug_l)))),
               tmean15            = mean(tmean15))

# Predict from final model
b.pred<-b.full
b.pred$pred<-NA
b.pred[b.pred$compartment=="Invertebrate",]$pred<-predict(final.inv,newdata=b.pred[b.pred$compartment=="Invertebrate",])
b.pred[b.pred$compartment=="Zooplankton",]$pred<-predict(final.zoo,newdata=b.pred[b.pred$compartment=="Zooplankton",])
b.pred[b.pred$compartment=="Phytoplankton",]$pred<-predict(final.phyto,newdata=b.pred[b.pred$compartment=="Phytoplankton",])


# Plot predict quadratic model
ggplot(b.full,aes(x=tmean15,y=log(C_mug_l),color=warming))+
  geom_jitter(width = 0.15,shape=1)+
  geom_pointrange(data=df_plot,aes(y=m.C_mug_l,ymin=m.C_mug_l-CI.C_mug_l,ymax=m.C_mug_l+CI.C_mug_l))+
  theme_classic()+
  stat_smooth(data=b.pred[-c(which(b.pred$compartment=="Invertebrate"),which(b.pred$compartment=="Phytoplankton" & b.pred$warming=="NW")),]
              ,aes(x=tmean15,y=pred),method = lm, formula = y ~ poly(x, 2, raw = TRUE),se=F)+
  facet_wrap(~factor(compartment,levels=c("Invertebrate","Zooplankton","Phytoplankton")),scales="free")+
  labs(y=expression(Log~(biomass~(µg~C.L^-1))),
       x=bquote(paste(italic(T[15])," (\u00B0C)")),
       fill="Temperature")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black","red"))+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"),
        strip.text.x = element_text(size=10, color="black",
                                    face="bold"),
        legend.position="none")

