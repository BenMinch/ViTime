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

###Environmental plots

environment= read.csv(args[3], header=TRUE, sep=",")
environment= read.csv('/home/benjamin/Desktop/ViTime/environment_data.csv')


#Change Date to day and use group_indices
environment_cols<- colnames(environment)
environment_cols<- environment_cols[-1]
environment<- environment %>% group_by(Date)%>%mutate(Day=group_indices())

#map this data onto cluster_plot by day
environment$Day<- as.factor(environment$Day)

cluster_plot$Day<- as.factor(cluster_plot$Day)

merged_df<- merge(cluster_plot,environment,by='Day')

#do overall spearman correlation for each environmental variable
#create a dataframe to store the results with p values correlation coefficient and variable

cor_df<- data.frame('variable'=character(),'correlation'=numeric(),'p.value'=numeric(),stringsAsFactors=FALSE)

for (i in environment_cols){
    cor<- cor.test(merged_df$Abundance,merged_df[,i],method='spearman')
    cor_df<- rbind(cor_df,data.frame('variable'=i,'correlation'=cor$estimate,'p.value'=cor$p.value,stringsAsFactors=FALSE))
}
#reset index
rownames(cor_df)<- NULL
write.csv(cor_df,paste(outfolder,"/Environmental_correlations_total.csv",sep=""))

cor<- ggplot(cor_df,aes(x=variable,y=correlation,fill=variable))+
geom_bar(stat='identity',position='dodge')+
geom_text(aes(label=ifelse(p.value<0.05,'*',NA)),position=position_dodge(width=0.9),vjust=-0.5)+
labs(x='Environmental Variable',y='Correlation')+
ggtitle('Correlation of Environmental Variables with Abundance')
ggsave(paste(outfolder,"/Environmental_correlations_total_barplot.png",sep=""),cor,width=10,height=10)
#Make a grid of plots for each environmental variable and abundance
plot_list<- list()
for (i in environment_cols){
    #make a list to store the plots
    plot<-ggplot(merged_df,aes(x=merged_df[,i],y=Abundance))+
    geom_point()+geom_smooth(method='lm')+
    labs(x=i,y='Abundance (RPKM)')+
    ggtitle(paste('Abundance vs',i))+
    #scale x axis to range of data
    scale_x_continuous(limits=c(min(merged_df[,i]),max(merged_df[,i])))
    #add text for correlation and p value
    plot_list[[i]]<- plot
}

grid<- grid.arrange(grobs=plot_list,ncol=2)
ggsave(paste(outfolder,"/Environmental_correlations_total.png",sep=""),grid,width=10,height=10)

#Do the same for each cluster

cluster1<- merged_df %>% filter(cluster==1)
cluster2<- merged_df %>% filter(cluster==2)
cluster3<- merged_df %>% filter(cluster==3)
cluster4<- merged_df %>% filter(cluster==4)
cluster5<- merged_df %>% filter(cluster==5)

corcluster_df<- data.frame('variable'=character(),'correlation'=numeric(),'p.value'=numeric(),'cluster'= character(),stringsAsFactors=FALSE)

for (i in environment_cols){
    cor<- cor.test(cluster1$Abundance,cluster1[,i],method='spearman')
    corcluster_df<- rbind(corcluster_df,data.frame('variable'=i,'correlation'=cor$estimate,'p.value'=cor$p.value, 'cluster'='one',stringsAsFactors=FALSE))
}

for (i in environment_cols){
    cor<- cor.test(cluster2$Abundance,cluster2[,i],method='spearman')
    corcluster_df<- rbind(corcluster_df,data.frame('variable'=i,'correlation'=cor$estimate,'p.value'=cor$p.value, 'cluster'='two',stringsAsFactors=FALSE))
}

for (i in environment_cols){
    cor<- cor.test(cluster3$Abundance,cluster3[,i],method='spearman')
    corcluster_df<- rbind(corcluster_df,data.frame('variable'=i,'correlation'=cor$estimate,'p.value'=cor$p.value, 'cluster'='three',stringsAsFactors=FALSE))
}

for (i in environment_cols){
    cor<- cor.test(cluster4$Abundance,cluster4[,i],method='spearman')
    corcluster_df<- rbind(corcluster_df,data.frame('variable'=i,'correlation'=cor$estimate,'p.value'=cor$p.value, 'cluster'= 'four',stringsAsFactors=FALSE))
}

for (i in environment_cols){
    cor<- cor.test(cluster5$Abundance,cluster5[,i],method='spearman')
    corcluster_df<- rbind(corcluster_df,data.frame('variable'=i,'correlation'=cor$estimate,'p.value'=cor$p.value, 'cluster'='five',stringsAsFactors=FALSE))
}
rownames(corcluster_df)<- NULL
write.csv(corcluster_df,paste(outfolder,"/Environmental_correlations_cluster.csv",sep=""))
#Make a barplot for each cluster showing variable on x axis and correlation on y axis with a star for sigificant p values

clusters<- ggplot(corcluster_df,aes(x=variable,y=correlation,fill=cluster))+
geom_bar(stat='identity',position='dodge')+
geom_text(aes(label=ifelse(p.value<0.05,'*',NA)),position=position_dodge(width=0.9),vjust=-0.5)+
labs(x='Environmental Variable',y='Correlation')+
ggtitle('Correlation of Environmental Variables with Abundance')

ggsave(paste(outfolder,"/Environmental_correlations_cluster.png",sep=""),clusters,width=10,height=10)

#check R version

