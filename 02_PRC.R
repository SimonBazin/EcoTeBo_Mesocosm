########### Principal Response Curve ########### 
################################################

#------------------------------------------------------- 
# This script is used to analyze community composition 

# 1) Macroinvertebrate
# 2) Phytoplankton


#------------------------------------------------------- 
# Packages and data

# Packages
library(ggplot2)
library(vegan)

# Data
df_taxa<-read.table(file="df_taxa.csv",header=T,sep=";",dec=".")


# Mesocosm_ID as character
df_taxa$mesocosm_ID<-as.character(df_taxa$mesocosm_ID)
# Mesocosm metadata
m<-aggregate(data=df_taxa, cbind(fish,warming) ~ mesocosm_ID + date + date_c + tmean15, unique)
m$treatment<-paste(m$fish,m$warming,sep="_")

# Treatment and week as factors for PRC
treatment<-factor(m$treatment,levels = c("NF_NW","NF_W","LB_NW","LB_W","SB_NW","SB_W"))
warming<-factor(m$warming)
fish<-factor(m$fish)
week<-factor(rep(c(4,8,12),each=24))


#--------------------------------------------------------------------
# Invertebrate


# DF inv
t_inv<-unique(df_taxa[df_taxa$compartment=="inv",]$taxa)
inv<-as.data.frame(matrix(NA,nrow=nrow(m),ncol=length(t_inv)))
colnames(inv)<-t_inv

# Density for each taxa in each mesocosm at each date
for(i in 1:nrow(m)){
  for(j in 1:length(t_inv)){
    if(nrow(df_taxa[df_taxa$mesocosm_ID==m[i,]$mesocosm_ID & df_taxa$date==m[i,]$date & df_taxa$taxa==t_inv[j],])>0){
    inv[i,j]<-df_taxa[df_taxa$mesocosm_ID==m[i,]$mesocosm_ID & df_taxa$date==m[i,]$date & df_taxa$taxa==t_inv[j],]$N_l
    }
    else{
      inv[i,j]<-0
    }
  }
}

# Hellinger transformation
inv<-decostand(inv, method = "hellinger")

# PRC
mod <- prc(inv, treatment, week)
mod            # RDA
summary(mod)   # PRC
sum_mod<-summary(mod)[["sp"]]  # Taxa with scores 

# Permutation anova
ctrl <- how(plots = Plots(strata = treatment,type = "free"),
            within = Within(type = "series"), nperm = 999)
anova(mod, permutations = ctrl, first=TRUE)

# Plot PRC
col<-c("red","black","red","black","red")
lty<-c(1,2,2,3,3)
pch<-c(1,2,2,3,3)
plot(mod, select = abs(sum_mod) > 0.1,species = TRUE, scaling = "symmetric",
     axis = 1, correlation = FALSE, type = "b",
     lty = lty, col = col, pch=pch, legpos="none",xlab="Week") 
legend(x=0.415,y=0.2,lty = c(1,lty), pch=c(1,pch),col = c("black",col),legend=levels(treatment),
       horiz=TRUE,cex=0.57,xjust=0.5)



#--------------------------------------------------------------------
# Phytoplankton


# DF phyto
t_phyto<-unique(df_taxa[df_taxa$compartment=="phyto",]$taxa)
phyto<-as.data.frame(matrix(NA,nrow=nrow(m),ncol=length(t_phyto)))
colnames(phyto)<-t_phyto

# Density for each taxa in each mesocosm at each date
for(i in 1:nrow(m)){
  for(j in 1:length(t_phyto)){
    if(nrow(df_taxa[df_taxa$mesocosm_ID==m[i,]$mesocosm_ID & df_taxa$date==m[i,]$date & df_taxa$taxa==t_phyto[j],])>0){
      phyto[i,j]<-df_taxa[df_taxa$mesocosm_ID==m[i,]$mesocosm_ID & df_taxa$date==m[i,]$date & df_taxa$taxa==t_phyto[j],]$N_l
    }
    else{
      phyto[i,j]<-0
    }
  }
}

# Hellinger transformation
phyto<-decostand(phyto, method = "hellinger")

# PRC
mod <- prc(phyto, treatment, week)
mod            # RDA
summary(mod)   # PRC
sum_mod<-summary(mod)[["sp"]]  # Taxa with scores 
# Permutation anova
ctrl <- how(plots = Plots(strata = treatment,type = "free"),
            within = Within(type = "series"), nperm = 999)
anova(mod, permutations = ctrl, first=TRUE)

# Plot PRC
col<-c("red","black","red","black","red")
lty<-c(1,2,2,3,3)
pch<-c(1,2,2,3,3)
plot(mod, select = abs(sum_mod) > 0.1,species = TRUE, scaling = "symmetric",
     axis = 1, correlation = FALSE, type = "b",
     lty = lty, col = col, pch=pch, legpos="",xlab="Week")
legend(x=0.405,y=0.2,lty = c(1,lty), pch=c(1,pch),col = c("black",col),legend=levels(treatment),
       horiz=TRUE,cex=0.57,xjust=0.5)

