##Celebrity.R##
library(tidyverse)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
data = read.csv(args[2], header=TRUE, sep=",")
outfolder = args[1]

outfolder='/home/benjamin/Desktop/ViTime'

#make a barplot of counts in each celebrity
data$Celeb_score<- as.factor(data$Celeb_score)
plot<- ggplot(data,aes(x=Celeb_score,fill=Celeb_score))+
geom_bar()+ ggtitle("Celebrity Counts")

# make a facet plot using a single representative from each celebrity
data_subset<- data %>% group_by(Celeb_score)%>% slice(1)

data_plot<- data_subset %>% gather(key = "Day", value = "Abundance",-Celeb_score)
data_plot<- data_plot %>% group_by(Day) %>% mutate(Day=group_indices())

x<- ggplot(data_plot,aes(x=Day,y=Abundance,color=Celeb_score))+
geom_point()+geom_line()+
facet_wrap(~Celeb_score,scales='free_y')+
labs(x='Time',y='Abundance (RPKM)')+ #make x axis angled
theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(gridExtra)

ggsave(paste(outfolder,"/Celebrity_membership.png",sep=""),grid.arrange(plot,x,ncol=2),width=20,height=10)
