################### Density  ###################
################################################

#------------------------------------------------------- 
# This script is used to analyze community density data  

# 1) Full model
# 2) Macroinvertebrate
# 3) Zooplankton
# 4) Phytoplankton

#------------------------------------------------------- 
# Packages and data

# Packages
library(plyr)
library(car)
library(lmerTest)
library(ggplot2)

# Data
df_taxa<-read.table(file="df_taxa.csv",header=T,sep=";",dec=".")

# Rename
df_taxa$compartment <- factor(df_taxa$compartment, levels = c("inv", "zoo", "phyto"), 
                             labels = c("Invertebrate", "Zooplankton", "Phytoplankton"))

# Mesocosm_ID as character
df_taxa$mesocosm_ID<-as.character(df_taxa$mesocosm_ID)
# Mesocosm metadata
m<-aggregate(data=df_taxa, cbind(fish,warming) ~ mesocosm_ID + date + tmean15, unique)


# Density
d.full<-ddply(df_taxa, .(date,mesocosm_ID,fish,warming,compartment), summarise,
              N_l = sum(N_l),
              log.N_l = log(sum(N_l)),
              tmean15 = unique(tmean15))

d.full[d.full$compartment=="Invertebrate",]$log.N_l<-d.full[d.full$compartment=="Invertebrate",]$log.N_l-mean(d.full[d.full$compartment=="Invertebrate",]$log.N_l)
d.full[d.full$compartment=="Zooplankton",]$log.N_l<-d.full[d.full$compartment=="Zooplankton",]$log.N_l-mean(d.full[d.full$compartment=="Zooplankton",]$log.N_l)
d.full[d.full$compartment=="Phytoplankton",]$log.N_l<-d.full[d.full$compartment=="Phytoplankton",]$log.N_l-mean(d.full[d.full$compartment=="Phytoplankton",]$log.N_l)


#------------------------------------------------------- 
# Full model

# Quadratic model 
mod<-lmer(log(N_l)~(tmean15+I(tmean15^2))*fish*warming*compartment+(1|date)+(1|mesocosm_ID),data=d.full,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
summary(mod)
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
# Macroinvertebrate compartment

Invertebrate<-d.full[d.full$compartment=="Invertebrate",]

# Quadratic model 
mod<-lmer(log(N_l)~(tmean15+I(tmean15^2))*warming+(1|date)+(1|mesocosm_ID),data=Invertebrate,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
summary(mod)
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

Zooplankton<-d.full[d.full$compartment=="Zooplankton",]

# Quadratic model
mod<-lmer(log(N_l)~(tmean15+I(tmean15^2))*warming+(1|date)+(1|mesocosm_ID),data=Zooplankton,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(mod)
# Stepwise backward regression
step_res <- step(mod,alpha.random = 1)
final.zoo <- get_model(step_res)
Anova(final.zoo)
final.zoo<-lmer(log(N_l)~tmean15+I(tmean15^2)+warming+tmean15:warming+(1|date)+(1|mesocosm_ID),data=Zooplankton,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(final.zoo)
summary(final.zoo)
res_lme=residuals(final.zoo)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)



#------------------------------------------------------- 
# Phytoplankton compartment

Phytoplankton<-d.full[d.full$compartment=="Phytoplankton",]

# Quadratic model
mod<-lmer(log(N_l)~(tmean15+I(tmean15^2))*warming+(1|date)+(1|mesocosm_ID),data=Phytoplankton,REML=F,
          control = lmerControl(optimizer ="Nelder_Mead"))
Anova(mod)
summary(mod)
# Stepwise backward regression
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
df_plot<-ddply(d.full, .(compartment,date,warming), summarise,
               m.N_l            = mean(log(N_l),na.rm=T),
               CI.N_l           = 1.96*sd(log(N_l), na.rm=TRUE)/sqrt(length(na.omit(log(N_l)))),
               tmean15            = mean(tmean15))

# Predict from final models
d.pred<-d.full
d.pred$pred<-NA
d.pred[d.pred$compartment=="Invertebrate",]$pred<-predict(final.inv,newdata=d.pred[d.pred$compartment=="Invertebrate",])
d.pred[d.pred$compartment=="Zooplankton",]$pred<-predict(final.zoo,newdata=d.pred[d.pred$compartment=="Zooplankton",])
d.pred[d.pred$compartment=="Phytoplankton",]$pred<-predict(final.phyto,newdata=d.pred[d.pred$compartment=="Phytoplankton",])


# Plot predict quadratic model
ggplot(d.pred,aes(x=tmean15,y=log(N_l),color=warming))+
  geom_jitter(width = 0.15,shape=1)+
  geom_pointrange(data=df_plot,aes(y=m.N_l,ymin=m.N_l-CI.N_l,ymax=m.N_l+CI.N_l))+
  theme_classic()+
  stat_smooth(data=d.pred[-which(d.pred$compartment=="Phytoplankton" & d.pred$warming=="NW"),],aes(x=tmean15,y=pred),method = lm, formula = y ~ poly(x, 2, raw = TRUE),se=F)+
  facet_wrap(~compartment,scales="free")+
  labs(y=expression(Log~(density~(N.L^-1))),x=bquote(paste(italic(T[15])," (\u00B0C)")),fill="Temperature")+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black","red"))+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"),
        strip.text.x = element_text(size=10, color="black",
                                    face="bold"),
        legend.position="none")
