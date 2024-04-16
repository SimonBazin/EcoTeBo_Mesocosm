################ Biomass Size Spectra  ################
#######################################################

#------------------------------------------------------- 
# This script is used to analyze Biomass Size Spectra 

# Packages
library(plyr)
library(car)
library(lmerTest)
library(gridExtra)
library(merTools)
library(ggplot2)

# Functions
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Data
df_bss<-read.table(file="df_bss.csv",header=T,sep=";",dec=".")

# Format 
df_bss$mesocosm_ID<-as.character(df_bss$mesocosm_ID)


# Lmer model
df_bss$warming<-factor(df_bss$warming,levels = c("NW","W"))
hLT.lm = lmer(log10BiomNorm ~ log10binMid.scaled*tmean15*fish*warming+(1|date)+(1|mesocosm_ID), na.action = na.omit,data=df_bss,REML=F,
              control = lmerControl(optimizer ="Nelder_Mead"))
Anova(hLT.lm)
summary(hLT.lm)
# Stepwise backward regression
step_res <- step(hLT.lm,alpha.random = 1)
final <- get_model(step_res)
Anova(final)
summary(final)
res_lme=residuals(final)
plot(res_lme)
abline(h=0,col="red")


# Newdata
newdata<-expand.grid(log10binMid.scaled=unique(df_bss$log10binMid.scaled),
                     date=unique(df_bss$date),
                     mesocosm_ID=unique(df_bss$mesocosm_ID),
                     warming=c("NW","W"),
                     tmean15=c(10,20,30))
newdata <- cbind(newdata, predictInterval(final, newdata,level=0.95))
unique(newdata$mesocosm_ID)

# Plot tmean15 effect on BSS slope
h1<-ggplot(df_bss,aes(x=log10binMid.scaled,y=log10BiomNorm,color=tmean15))+
  geom_point()+
  theme_classic()+
  scale_color_gradient(low="blue",high="red")+
  labs(color=bquote(paste(italic(T[15])," (\u00B0C)")),x=expression(Log[10](binMid.scaled)),y=expression(log[10](BiomNorm)))+
  theme(legend.position = "none")+
  guides(fill=F)+
  geom_smooth(data=newdata,method="lm",aes(y=fit,fill=as.character(tmean15)),size=1.5,se=F)+
  geom_smooth(data=newdata,method="lm",aes(y=upr,fill=as.character(tmean15)),size=0.3,lty=2,se=F)+
  geom_smooth(data=newdata,method="lm",aes(y=lwr,fill=as.character(tmean15)),size=0.3,lty=2,se=F)


# Plot BSS intercept

df_coef<-data.frame(warming=c("W",rep("NW",3)),
                    tmean15=c(20,10,20,30),
                    coef=rep("Intercept",4),
                    fit=NA,lwr=NA,upr=NA)


# 95 % confidence intervals of slope and intercept are predicted scaling tmean15
df_bss$t15.scaled10<-df_bss$tmean15-10
df_bss$t15.scaled20<-df_bss$tmean15-20
df_bss$t15.scaled30<-df_bss$tmean15-30

mod.s<-list(mod.s10 = lmer(log10BiomNorm ~ log10binMid.scaled+t15.scaled10+warming+log10binMid.scaled:t15.scaled10+t15.scaled10:warming+(1|date)+(1|mesocosm_ID), na.action = na.omit,data=df_bss,REML=F,
                           control = lmerControl(optimizer ="Nelder_Mead")),
            mod.s20 = lmer(log10BiomNorm ~ log10binMid.scaled+t15.scaled20+warming+log10binMid.scaled:t15.scaled20+t15.scaled20:warming+(1|date)+(1|mesocosm_ID), na.action = na.omit,data=df_bss,REML=F,
                           control = lmerControl(optimizer ="Nelder_Mead")),
            mod.s30 = lmer(log10BiomNorm ~ log10binMid.scaled+t15.scaled30+warming+log10binMid.scaled:t15.scaled30+t15.scaled30:warming+(1|date)+(1|mesocosm_ID), na.action = na.omit,data=df_bss,REML=F,
                           control = lmerControl(optimizer ="Nelder_Mead")))


t<-c(10,20,30)
for(i in 1:length(t)){
  df_coef[df_coef$warming=="NW" & df_coef$coef=="Intercept" & df_coef$tmean15==t[i],]$fit<-summary(mod.s[[i]])$coefficients[1,1]
  df_coef[df_coef$warming=="NW" & df_coef$coef=="Intercept" & df_coef$tmean15==t[i],]$upr<-summary(mod.s[[i]])$coefficients[1,1]+1.96*summary(mod.s[[i]])$coefficients[1,2]
  df_coef[df_coef$warming=="NW" & df_coef$coef=="Intercept" & df_coef$tmean15==t[i],]$lwr<-summary(mod.s[[i]])$coefficients[1,1]-1.96*summary(mod.s[[i]])$coefficients[1,2]
}

# Warming as intercept in lmer model
df_bss$warming<-factor(df_bss$warming,levels = c("W","NW"))
mod.s20 = lmer(log10BiomNorm ~ log10binMid.scaled+t15.scaled20+warming+log10binMid.scaled:t15.scaled20+t15.scaled20:warming+(1|date)+(1|mesocosm_ID), na.action = na.omit,data=df_bss,REML=F,
              control = lmerControl(optimizer ="Nelder_Mead"))
# Fill df_coef
df_coef[df_coef$warming=="W" & df_coef$coef=="Intercept" & df_coef$tmean15==20,]$fit<-summary(mod.s20)$coefficients[1,1]
df_coef[df_coef$warming=="W" & df_coef$coef=="Intercept" & df_coef$tmean15==20,]$upr<-summary(mod.s20)$coefficients[1,1]+1.96*(summary(final)$coefficients[1,2])
df_coef[df_coef$warming=="W" & df_coef$coef=="Intercept" & df_coef$tmean15==20,]$lwr<-summary(mod.s20)$coefficients[1,1]-1.96*(summary(final)$coefficients[1,2])



# Plot Intercept and Slope
h2<-ggplot(df_coef,aes(x=warming,y=fit,color=as.character(tmean15)))+
  geom_point(position = position_dodge(0.5),size=2)+
  geom_pointrange(aes(ymin=lwr,ymax=upr),position = position_dodge(0.5))+
  scale_color_manual(values=c("blue","purple","red"))+
  theme_classic()+
  labs(x="",y="BSS Intercept",color=expression(italic(T[15])))+
  theme(legend.position="none",
        axis.text = element_text(size=11,color="black"),
        axis.title = element_text(size=11,color="black"))

# Get legend from h2 plot
lgd<-get_legend(h1)

# Combine BSS and coefficients plot
grid.arrange(h1,h2,lgd,ncol=3,nrow=1,widths = c(1.5,1.5,0.5))

