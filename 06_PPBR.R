################### PPBR  ######################
################################################

#------------------------------------------------------- 
# This script is used to analyze Predator Prey Biomass Ratio  

#-------------------------------------------------------
# Packages and Data 

# Packages
library(car)
library(lmerTest)
library(plyr)
library(ggplot2)

# Data
ppbr<-read.table(file="df_ppbr.csv",header=T,sep=";",dec=".")

# Mesocosm_ID as character
ppbr$mesocosm_ID<-as.character(ppbr$mesocosm_ID)

# Mesocosm metadata
m<-aggregate(data=ppbr, cbind(fish,warming) ~ mesocosm_ID + date + tmean15, unique)


# Full model
m1<-lmer(log(ZP,10)~(tmean15+I(tmean15^2))*fish*warming+(1|date)+(1|mesocosm_ID),data=ppbr,REML=F,
         control = lmerControl(optimizer ="Nelder_Mead"))
Anova(m1)
# Stepwise backward regression
step_res <- step(m1,alpha.random = 1)
final <- get_model(step_res)
Anova(final)
summary(final)
res_lme=residuals(final)
plot(res_lme)
qqnorm(res_lme,label=TRUE)
qqline(res_lme)
hist(res_lme)

# Mean and CI for plot
df1<-ddply(ppbr, .(warming,date), summarise,
           m.ZP          = mean(log(ZP),na.rm=T),
           CI.ZP         = 1.96*sd(log(ZP), na.rm=TRUE)/sqrt(length(na.omit(log(ZP)))),
           tmean15       = mean(tmean15))

# Predict from final model
ppbr.pred<-ppbr
ppbr.pred$pred<-predict(final,newdata=ppbr.pred)

# Plot predict PPBR
ggplot(ppbr,aes(x=tmean15,y=log(ZP),color=warming))+
  geom_jitter(alpha=0.5,width = 0.1,shape=1)+
  scale_color_manual(values=c("black","red"))+
  geom_pointrange(data=df1,aes(x=tmean15,y=m.ZP,ymin=m.ZP-CI.ZP,ymax=m.ZP+CI.ZP))+
  stat_smooth(data=ppbr.pred,aes(x=tmean15,y=pred),method = lm, formula = y ~ poly(x, 2, raw = TRUE),se=F)+
  theme_classic()+
  labs(x=bquote(paste(italic(T[15])," (\u00B0C)")),y="log(PPBR)",color="Temperature",fill="Compartment")+
  theme(legend.position="none",
        axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=12,color="black"),
        legend.title = element_text(size=12,color="black"),
        legend.text = element_text(size=11,color="black"))

