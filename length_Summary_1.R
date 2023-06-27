library(readr)
length_time <- read_csv("length_time.csv")
View(length_time)
length_time$temp<-as.factor(length_time$temp)
library(dplyr)
library(ggplot2)
data.summary= length_time%>%
  group_by(temp,time)%>%
  summarise_if(is.numeric, funs(mean,sd))

plot.1<- ggplot(data=data.summary, aes(x =time, y =mean, color=temp))+
  geom_point()+
  geom_errorbar(aes(x=time, ymin=mean-sd, ymax=mean+sd), width=0.2, alpha=0.5, size=0.5)+
    geom_smooth()+
  labs(x=("time (days)"), y = "length as squared root area (cm)")
plot.1

write.table(data.summary, file = "length_mean_summary.txt", sep = ',', row.names = T)
