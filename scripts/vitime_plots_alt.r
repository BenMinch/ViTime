##Create an Rscript to make plots of timecourse data

##Load libraries
library(ggplot2)
library(tidyverse)

##Load data
args = commandArgs(trailingOnly=TRUE)
data = read.csv(args[1], header=TRUE, sep=",")
outfolder = args[2]

###Stacked barplot of order level abundance over time

#group data by order and sum all the other rows
datas<- data[,-1]
data2<-datas %>% group_by(Order) %>% summarise_all(sum)
data2<-data2[,-2]
df_long<- gather(data2, key = "Day", value = "Abundance", -Order)


df_long2<- df_long %>% group_by(Day) %>% mutate(Proportion = Abundance/sum(Abundance))

plot<- ggplot(df_long2,aes(x=Day,y=Proportion,fill=Order))+
geom_col()+ ggtitle("Order level abundance over time")+
#make x axis angled
theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(outfolder,"/Order_abundance_over_time.png",sep=""),plot,width=10,height=10)

###Individual genome heatmap

library(pheatmap)
heatmap_long<- data %>% gather(key = "Day", value = "Abundance", -Order,-Genome)
heatmap_long_z<- heatmap_long %>% group_by(Day) %>% mutate(zscore= (Abundance-mean(Abundance))/sd(Abundance))
head(heatmap_long_z)

#convert back to wide format
library(tidyr)
heatmap_wide<- tidyr::pivot_wider(heatmap_long_z, id_cols = Genome, names_from = Day,values_from = zscore)

heatmap_wide$Order<- data$Order
heatmap_matrix<- heatmap_wide %>% select(-Genome,-Order) %>% as.matrix()
rownames(heatmap_matrix)<- heatmap_wide$Genome

png(paste(outfolder,"/Genome_heatmap.png",sep=""),width=10,height=10,units="in",res=300)
pheatmap(heatmap_matrix,scale='none',cluster_rows = TRUE,main='Genome abundance over time',cluster_cols = FALSE,show_rownames = TRUE) 
dev.off()

hc<-hclust(dist(heatmap_matrix))
clusters<- cutree(hc,k=5)
data$cluster<- clusters

#visualize the clusters as a sideways barplot with counts of orders on the y axis

plot<- ggplot(data,aes(x=cluster,fill=Order))+
geom_bar()+ ggtitle("Cluster Membership")

cluster_subset<- data %>% group_by(cluster)%>% slice(1)
cluster_subset$cluster<- as.factor(cluster_subset$cluster)


cluster_plot<- cluster_subset %>% gather(key = "Day", value = "Abundance", -Order,-Genome,-cluster)

cluster_plot<- cluster_plot %>% group_by(Day) %>% mutate(Day=group_indices())

x<- ggplot(cluster_plot,aes(x=Day,y=Abundance,color=cluster))+
geom_point()+geom_line()+
facet_wrap(~cluster,scales='free_y')+
labs(x='Time',y='Abundance (RPKM)')+ #make x axis angled
theme(axis.text.x = element_text(angle = 90, hjust = 1))


#arrange these two plots next to each other
library(gridExtra)

ggsave(paste(outfolder,"/Genome_cluster_membership.png",sep=""),grid.arrange(plot,x,ncol=2),width=20,height=10)
