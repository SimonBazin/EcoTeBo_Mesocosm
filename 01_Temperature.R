############## Temperature ###############
##########################################

#-----------------------------------------
# Water temperature evolution over time during mesocosm experiment
# Water temperature difference between heated and non-heated mesocosms

# Package
library(ggplot2)
library(scales)
library(car)

# Data
load(file = "df_temp.RData")

# English date system
Sys.setlocale("LC_TIME", "English_United States")

# As date
t_hobo$date_time<-as.POSIXct(t_hobo$date_time)

# Plot temperature over time from HOBO data
lims <- as.POSIXct(strptime(c("2021-03-24 06:00","2021-06-15 18:00"), format = "%Y-%m-%d %H:%M"))
ggplot(t_hobo, aes(x=date_time, y=Temp,color=warming, group = mesocosm_ID)) + 
  geom_line()+
  theme_bw() +
  geom_smooth(method = "auto",aes(fill=warming),linetype="blank")+
  scale_fill_manual(values = c("black","red"))+
  scale_color_manual(values = c("black","red"))+
  labs(x="",y="Temperature (°C)",fill="Temperature")+
  theme(legend.position = "none", axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12))+
  scale_x_datetime(limits =lims,breaks = date_breaks("2 weeks"),
                   labels = date_format("%B %d"))

# Linear model
mod<-lm(Temp~warming,data=t_hobo)
summary(mod)
Anova(mod)
# Water temperature of heated mesocosms was +3.35 °C higher
