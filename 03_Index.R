################### Ecological index  ###################
#########################################################

#------------------------------------------------------- 
# This script is used to analyze community ecological index S and H'

# For each part, S and H' 

# For each index :
# 1) Full model
# 2) Macroinvertebrate
# 3) Phytoplankton


#------------------------------------------------------- 
# Packages and data
library(ggplot2)
library(car)
library(lmerTest)
library(plyr)
library(gridExtra)

# Data
df_index<-read.table(file="df_index.csv",header=T,sep=";",dec=".")

# Mesocosm as character
df_index$mesocosm_ID<-as.character(df_index$mesocosm_ID)

# Scaled ecological index
df_index$Hs<-NA
df_index$Ss<-NA

df_index[df_index$compartment=="inv",]$Hs<-df_index[df_index$compartment=="inv",]$H-mean(df_index[df_index$compartment=="inv",]$H)
df_index[df_index$compartment=="zoo",]$Hs<-df_index[df_index$compartment=="zoo",]$H-mean(df_index[df_index$compartment=="zoo",]$H)
df_index[df_index$compartment=="phyto",]$Hs<-df_index[df_index$compartment=="phyto",]$H-mean(df_index[df_index$compartment=="phyto",]$H)

df_index[df_index$compartment=="inv",]$Ss<-df_index[df_index$compartment=="inv",]$S-mean(df_index[df_index$compartment=="inv",]$S)
df_index[df_index$compartment=="zoo",]$Ss<-df_index[df_index$compartment=="zoo",]$S-mean(df_index[df_index$compartment=="zoo",]$S)
df_index[df_index$compartment=="phyto",]$Ss<-df_index[df_index$compartment=="phyto",]$S-mean(df_index[df_index$compartment=="phyto",]$S)


###############################
# Species richness S

# Full lmer model
m1<-lmer(S~(tmean15+I(tmean15^2))*fish*warming*compartment+(1|date)+(1|mesocosm_ID),data=df_index,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(m1)
# Stepwise backward regression
step_res <- step(m1,alpha.random = 1)
final.full <- get_model(step_res)
Anova(final.full)
# I(tmean15^2):fish but no tmean15:fish 
# Manually add tmean15:fish
final.full<-lmer(S~tmean15+I(tmean15^2)+fish+compartment+tmean15:fish+tmean15:fish+I(tmean15^2):fish+tmean15:compartment+I(tmean15^2):compartment+(1|date)+(1|mesocosm_ID),data=df_index,REML=F,
                 control=lmerControl(optimizer="Nelder_Mead"))
Anova(final.full)
# Manual backward regression
df_index$fish<-factor(df_index$fish,levels = c("NF","LB","SB"))
final.full<-lmer(S~tmean15+I(tmean15^2)+fish+compartment+tmean15:fish+tmean15:fish+tmean15:compartment+I(tmean15^2):compartment+(1|date)+(1|mesocosm_ID),data=df_index,REML=F,
                 control=lmerControl(optimizer="Nelder_Mead"))
Anova(final.full)
summary(final.full)
res_lme=residuals(final.full)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)


#---------------------------------
# Macroinvertebrate compartment

inv<-df_index[df_index$compartment=="inv",]

# Full lmer model
m1<-lmer(S~(tmean15+I(tmean15^2))+(1|date)+(1|mesocosm_ID),data=inv,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(m1)
# Manual backward regression
final.inv<-lmer(S~tmean15+(1|date)+(1|mesocosm_ID),data=inv,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(final.inv)
summary(final.inv)
res_lme=residuals(final.inv)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)


#---------------------------------
# Phytoplankton compartment

phyto<-df_index[df_index$compartment=="phyto",]

# Full lmer model
m1<-lmer(S~(tmean15+I(tmean15^2))+(1|date)+(1|mesocosm_ID),data=phyto,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(m1)
# Stepwise backward regression
step_res <- step(m1,alpha.random = 1)
final.phyto <- get_model(step_res)
Anova(final.phyto)
summary(final.phyto)
res_lme=residuals(final.phyto)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)



#---------------------------------
# Plot

# Mean and CI for plot
df_plot<-ddply(df_index, .(date,fish,warming), summarise,
               m.S          = mean(S,na.rm=T),
               CI.S         = 1.96*sd(S, na.rm=TRUE)/sqrt(length(na.omit(S))),
               tmean15      = mean(tmean15))
# Predict from final model
df.pred<-df_index
df.pred$pred<-predict(final.full,newdata=df.pred)

# Plot 
df.pred$fish<-factor(df.pred$fish,levels=c("NF","LB","SB"))
df_plot$fish<-factor(df_plot$fish,levels=c("NF","LB","SB"))
fish.labs <- c("No Fish", "Large-bodied Fish","Small-bodied Fish")
names(fish.labs) <- c("NF", "LB","SB")
df.pred[df.pred$fish=="NF",]$pred<-NA
ggplot(df.pred,aes(x=tmean15,y=S))+
  geom_jitter(width = 0.15,aes(shape=compartment))+
  geom_pointrange(data=df_plot,aes(y=m.S,ymin=m.S-CI.S,ymax=m.S+CI.S))+
  theme_classic()+
  facet_wrap(~fish,labeller = labeller(fish=fish.labs))+
  stat_smooth(data=df.pred,aes(x=tmean15,y=pred),method = lm,se=F,color="darkgrey")+
  theme(legend.position = "none")+
  scale_shape_manual(values = c(1,2))+
  labs(y=expression(italic(S)),x=bquote(paste(italic(T[15])," (\u00B0C)")))+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"),
        # strip.background = element_rect(colour="white", fill="white", 
        #                                 size=1.5, linetype="solid"),
        strip.text.x = element_text(size=10, color="black",
                                    face="bold"),
        strip.text.y = element_text(size=10, color="black",
                                    face="bold"),
        legend.position="none")





# Mean and CI for plot
df_plot<-ddply(df_index, .(compartment,date,warming), summarise,
               m.S          = mean(S,na.rm=T),
               CI.S         = 1.96*sd(S, na.rm=TRUE)/sqrt(length(na.omit(S))),
               tmean15      = mean(tmean15))


# Predict from final model
df.pred<-df_index
df.pred$pred<-NA
df.pred[df.pred$compartment=="inv",]$pred<-predict(final.inv,newdata=df.pred[df.pred$compartment=="inv",])
df.pred[df.pred$compartment=="phyto",]$pred<-predict(final.phyto,newdata=df.pred[df.pred$compartment=="phyto",])

# Plot 
comp.labs <- c("Invertebrate", "Phytoplankton")
names(comp.labs) <- c("inv", "phyto")
ggplot(df.pred,aes(x=tmean15,y=S))+
  geom_jitter(width = 0.15,aes(shape=compartment))+
  geom_pointrange(data=df_plot,aes(y=m.S,ymin=m.S-CI.S,ymax=m.S+CI.S))+
  theme_classic()+
  facet_wrap(~compartment,labeller = labeller(compartment=comp.labs))+
  stat_smooth(data=df.pred[df.pred$compartment=="inv",],aes(x=tmean15,y=pred),method = lm,se=F,color="darkgrey")+
  stat_smooth(data=df.pred[df.pred$compartment=="phyto",],aes(x=tmean15,y=pred),method = lm, formula = y ~ poly(x, 2, raw = TRUE),se=F,color="darkgrey")+
  theme(legend.position = "none")+
  scale_shape_manual(values = c(1,2))+
  labs(y=expression(italic(S)),x=bquote(paste(italic(T[15])," (\u00B0C)")))+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"),
        # strip.background = element_rect(colour="white", fill="white", 
        #                                 size=1.5, linetype="solid"),
        strip.text.x = element_text(size=10, color="black",
                                    face="bold"),
        strip.text.y = element_text(size=10, color="black",
                                    face="bold"),
        legend.position="none")


###############################
# Shannon index

# Full lmer model
m1<-lmer(H~(tmean15+I(tmean15^2))*fish*warming*compartment+(1|date)+(1|mesocosm_ID),data=df_index,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
# Stepwise backward regression
step_res <- step(m1,alpha.random = 1)
final.full <- get_model(step_res)
Anova(final.full)
summary(final.full)
res_lme=residuals(final.full)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)

#---------------------------------
# Macroinvertebrate compartment

# lmer model
m1<-lmer(H~(tmean15+I(tmean15^2))*warming*fish+(1|date)+(1|mesocosm_ID),data=inv,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(m1)
# Stepwise backward regression
step_res <- step(m1,alpha.random = 1)
final.inv <- get_model(step_res)
Anova(final.inv)
inv$fish<-factor(inv$fish,levels = c("NF","LB","SB"))
final.inv<-lmer(H~tmean15+I(tmean15^2)+warming+fish+tmean15:warming+I(tmean15^2):warming+tmean15:fish+I(tmean15^2):fish+(1|date)+(1|mesocosm_ID),data=inv,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(final.inv)
summary(final.inv)
res_lme=residuals(final.inv)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)


#-------------------------------
# Phytoplankton compartment

# lmer model
m1<-lmer(H~(tmean15+I(tmean15^2))*warming*fish+(1|date)+(1|mesocosm_ID),data=phyto,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(m1)
summary(m1)
# Stepwise backward regression
step_res <- step(m1,alpha.random = 1)
final.phyto <- get_model(step_res)
Anova(final.phyto)
phyto$fish<-factor(phyto$fish,levels = c("NF","LB","SB"))
final.phyto<-lmer(H~(tmean15+I(tmean15^2))*fish-I(tmean15^2):fish+(1|date)+(1|mesocosm_ID),data=phyto,REML=F,
         control=lmerControl(optimizer="Nelder_Mead"))
Anova(final.phyto)
summary(final.phyto)
res_lme=residuals(final.phyto)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)


#-------------------------------
# Plot

# Mean and CI for plot
df_plot<-ddply(inv, .(compartment,date,warming), summarise,
               m.H          = mean(H,na.rm=T),
               CI.H         = 1.96*sd(H, na.rm=TRUE)/sqrt(length(na.omit(H))),
               tmean15      = mean(tmean15))

# Predict from final model
df.pred<-inv
df.pred$pred<-predict(final.inv,newdata=df.pred)

ggplot(df.pred,aes(x=tmean15,y=H, color=warming))+
  geom_jitter(width = 0.15,shape=1)+
  geom_pointrange(data=df_plot,aes(y=m.H,ymin=m.H-CI.H,ymax=m.H+CI.H))+
  theme_classic()+
  #facet_grid(,labeller = labeller(compartment=comp.labs,fish=fish.labs))+
  stat_smooth(data=df.pred,aes(x=tmean15,y=pred),method = lm, formula = y ~ poly(x, 2, raw = TRUE),se=F)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("black","red"))+
  labs(y=expression(italic("H'")),x=bquote(paste(italic(T[15])," (\u00B0C)")))+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"),
        # strip.background = element_rect(colour="white", fill="white", 
        #                                 size=1.5, linetype="solid"),
        strip.text.x = element_text(size=10, color="black",
                                    face="bold"),
        strip.text.y = element_text(size=10, color="black",
                                    face="bold"),
        legend.position="none")


# Mean and CI for plot
df_plot<-ddply(df_index, .(fish,date,warming), summarise,
               m.H          = mean(H,na.rm=T),
               CI.H         = 1.96*sd(H, na.rm=TRUE)/sqrt(length(na.omit(H))),
               tmean15      = mean(tmean15))

# Predict from final model
df.pred<-df_index
df.pred$pred<-NA
df.pred[df.pred$compartment=="inv",]$pred<-predict(final.inv,newdata=df.pred[df.pred$compartment=="inv",])
df.pred[df.pred$compartment=="phyto",]$pred<-predict(final.phyto,newdata=df.pred[df.pred$compartment=="phyto",])

df.pred$fish<-factor(df.pred$fish,levels=c("NF","LB","SB"))
df_plot$fish<-factor(df_plot$fish,levels=c("NF","LB","SB"))
df_index$fish<-factor(df_index$fish,levels=c("NF","LB","SB"))
fish.labs <- c("No Fish", "Large-bodied Fish","Small-bodied Fish")
names(fish.labs) <- c("NF", "LB","SB")
comp.labs <- c("Invertebrate", "Phytoplankton")
names(comp.labs) <- c("inv", "phyto")

ggplot(df_index,aes(x=tmean15,y=H))+
  geom_jitter(width = 0.15,shape=1)+
  geom_pointrange(data=df_plot,aes(y=m.H,ymin=m.H-CI.H,ymax=m.H+CI.H))+
  theme_classic()+
  facet_grid(compartment~fish,labeller = labeller(compartment=comp.labs,fish=fish.labs))+
  stat_smooth(data=df.pred,aes(x=tmean15,y=pred),method = lm, formula = y ~ poly(x, 2, raw = TRUE),se=F,color="darkgrey")+
  theme(legend.position = "none")+
  labs(y=expression(italic("H'")),x=bquote(paste(italic(T[15])," (\u00B0C)")))+
  theme(axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"),
        # strip.background = element_rect(colour="white", fill="white", 
        #                                 size=1.5, linetype="solid"),
        strip.text.x = element_text(size=10, color="black",
                                    face="bold"),
        strip.text.y = element_text(size=10, color="black",
                                    face="bold"),
        legend.position="none")
