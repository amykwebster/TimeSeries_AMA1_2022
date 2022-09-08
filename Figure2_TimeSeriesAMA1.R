#Figure 2
#!!!All data used for these scripts should be found on sheets in Supplementary Files 1 and 2!!!
#write a script that loops through each of the clusters and overlaps with gene lists for transcriptional regulators of interest
#also includes script for heatmap of genes included in clustering

#SET working directory to where files are on your computer
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries/Cluster_EnrichmentAnalysis/")
library(gplots)
library(Hmisc)
library(ggstance)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(gridExtra)

#Pull in the file that is in the Clusters_SeqNames tab of Supp File 1
#Clusters<-read.csv("Clustering_6027genes_0.8_129clusters.csv",header = F)
head(Clusters)
rownames(Clusters)<-Clusters$V1
Clusters<-Clusters[,-1]

#Pull in the file that is in the Zscore_DiffExpressedGenes tab of Supp File 1
#Zscores<-read.csv("clustering_z_means.csv",header = T)
head(Zscores)
rownames(Zscores)<-Zscores$gene
Zscores<-Zscores[,-1]

#Pull in the file that is in the TimeSeries_DiffExpResults tab of Supp File 1
ExpressedGenes<-read.csv("StarvationTimeSeries_v2_resultspage.csv",header = T)
ExpressedGenes<-subset(ExpressedGenes,FDR!="NA")
ExpressedGenes<-ExpressedGenes[!duplicated(ExpressedGenes$primaryIdentifier),]
ExpressedGenes<-ExpressedGenes[complete.cases(ExpressedGenes),]
rownames(ExpressedGenes)<-ExpressedGenes$primaryIdentifier
head(ExpressedGenes)
head(Clusters)

#merge clusters and expressed genes to get clusters with WBIDs
ClusterExpressedMerge<-merge(Clusters,ExpressedGenes,by.x=0,by.y="sequence")
head(ClusterExpressedMerge)
rownames(ClusterExpressedMerge)<-ClusterExpressedMerge$primaryIdentifier
Cluster_WBID<-ClusterExpressedMerge[,-1]
head(Cluster_WBID)
Cluster_WBID<-Cluster_WBID[,-130:-195]

ZscoreExpressedMerge<-merge(Zscores,ExpressedGenes,by.x = 0,by.y="sequence")
head(ZscoreExpressedMerge)
rownames(ZscoreExpressedMerge)<-ZscoreExpressedMerge$primaryIdentifier
Zscore_WBID<-ZscoreExpressedMerge[,-1]
Zscore_WBID<-Zscore_WBID[,-13:-78]
head(Zscore_WBID)

Zscore_Clusters_WBID<-Zscore_WBID[rownames(Cluster_WBID),]
Zscore_Clusters_Seq<-Zscores[rownames(Clusters),]
Zscore_Clusters_Seq2<-as.matrix(Zscore_Clusters_Seq)
#heatmap.2(as.matrix(Zscore_Clusters_WBID[ClusterGeneList,]), Colv = NA,labRow = FALSE)

Zscore_Clusters_Seq[4830:4840,] #row 4838 has an NA
Zscore_Clusters_Seq<-Zscore_Clusters_Seq[-4838,]

#this chunk calculates the mean z-score for each cluster
Clusters_means<-c()
head(Zscore_Clusters_Seq)
for (k in 1:129){
  Cluster1<-as.data.frame(Clusters[,k])
  rownames(Cluster1)<-rownames(Clusters)
  colnames(Cluster1)<-c("Cluster")
  ClusterGeneList<-rownames(subset(Cluster1,Cluster=="TRUE"))
  ClusterMean<-colMeans(Zscore_Clusters_Seq[ClusterGeneList,],na.rm = TRUE)
  ClusterMean<-as.data.frame(ClusterMean)
  Clusters_means<-rbind(Clusters_means,t(ClusterMean))
}

Clusters_means<-as.data.frame(Clusters_means)
head(Clusters_means)
rownames(Clusters_means)<-c(1:129)
Cluster_means1<-Clusters_means

Cluster_means1$maxval<-do.call(pmax,Cluster_means1)
Cluster_means1$maxtime<-colnames(Cluster_means1)[max.col(Cluster_means1,ties.method="first")]
Cluster_means1$clustenumber<-rownames(Cluster_means1)
Cluster_means2<-Cluster_means1
Cluster_means2$minval<-apply(Cluster_means2[1:129,1:12],1,FUN = min)
Cluster_means2[,1:12]<-Cluster_means2[,1:12] - Cluster_means2$minval
Cluster_means2$sum<-apply(Cluster_means2[1:129,1:12],1,FUN = sum)
Cluster_means2$sumInHalf<-Cluster_means2$sum / 2
head(Cluster_means2)

for(i in 1:129){
  if (Cluster_means2[i,18]<Cluster_means2[i,1]){
    Cluster_means2[i,19]<-colnames(Cluster_means2)[1]
  }
  else{
    if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]){
      Cluster_means2[i,19]<-colnames(Cluster_means2)[2]
    }
    else{
      if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]){
        Cluster_means2[i,19]<-colnames(Cluster_means2)[3]
      }
      else{
        if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]){
          Cluster_means2[i,19]<-colnames(Cluster_means2)[4]
        }
        else{
          if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]){
            Cluster_means2[i,19]<-colnames(Cluster_means2)[5]
          }
          else{
            if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]){
              Cluster_means2[i,19]<-colnames(Cluster_means2)[6]
            }
            else{
              if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]+Cluster_means2[i,7]){
                Cluster_means2[i,19]<-colnames(Cluster_means2)[7]
              }
              else{
                if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]+Cluster_means2[i,7]+Cluster_means2[i,8]){
                  Cluster_means2[i,19]<-colnames(Cluster_means2)[8]
                }
                else{
                  if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]+Cluster_means2[i,7]+Cluster_means2[i,8]+Cluster_means2[i,9]){
                    Cluster_means2[i,19]<-colnames(Cluster_means2)[9]
                  }
                  else{
                    if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]+Cluster_means2[i,7]+Cluster_means2[i,8]+Cluster_means2[i,9]+Cluster_means2[i,10]){
                      Cluster_means2[i,19]<-colnames(Cluster_means2)[10]
                    }
                    else{
                      if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]+Cluster_means2[i,7]+Cluster_means2[i,8]+Cluster_means2[i,9]+Cluster_means2[i,10]+Cluster_means2[i,11]){
                        Cluster_means2[i,19]<-colnames(Cluster_means2)[11]
                      }
                      else{
                        if (Cluster_means2[i,18]<Cluster_means2[i,1]+Cluster_means2[i,2]+Cluster_means2[i,3]+Cluster_means2[i,4]+Cluster_means2[i,5]+Cluster_means2[i,6]+Cluster_means2[i,7]+Cluster_means2[i,8]+Cluster_means2[i,9]+Cluster_means2[i,10]+Cluster_means2[i,11]+Cluster_means2[i,12]){
                          Cluster_means2[i,19]<-colnames(Cluster_means2)[12]
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

head(Cluster_means2,20)
subset(Cluster_means2,maxtime=="g_12hr")


time_points<-c("a_.2hr","b_0hr","c_2hr","d_4hr","e_6hr","f_9hr","g_12hr","h_1d","i_2d","j_4d","k_8d","l_12d")
Clusters_ordered<-c()
for (n in 1:12){
  ClusterTime_subset<-subset(Cluster_means2,maxtime==time_points[n])
  for (m in 1:length(levels(factor(ClusterTime_subset$V19)))) {
    ClusterTime_subset2<-subset(ClusterTime_subset,V19==levels(factor(ClusterTime_subset$V19))[m])
    if (n < 7){
      ClusterTime_subset2<-ClusterTime_subset2[order(-ClusterTime_subset2$maxval),]
    }
    else{
      ClusterTime_subset2<-ClusterTime_subset2[order(ClusterTime_subset2$maxval),]
    }
    Clusters_ordered<-rbind(Clusters_ordered,ClusterTime_subset2)
  }
}

tail(Clusters_ordered,20)
dim(subset(Clusters_ordered,maxtime=="a_.2hr")) #19 clusters have peak at -2hr
dim(subset(Clusters_ordered,maxtime=="l_12d")) #36 clusters have peak at 12 d 


rearrange_order<-Clusters_ordered$clustenumber

paletteLength <- 100
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(paletteLength)

myBreaks <- c(seq(min(as.matrix(Clusters_means)), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(as.matrix(Clusters_means))/paletteLength, max(as.matrix(Clusters_means)), length.out=floor(paletteLength/2)))

#put all the clustered genes in order by cluster determined in the above loop
#clusterorder < 11 determines how many clusters I want to have their own color

rearrange_order

rowannotation2<-c()
Clusters_ordered<-c()
for (k in 1:129){
  clusterorder<-as.numeric(rearrange_order[k])
  if (clusterorder < 11){ #change this to change which clusters get their own color later
    colnameforcluster<-clusterorder
    Cluster1<-as.data.frame(Clusters[,clusterorder])
    rownames(Cluster1)<-rownames(Clusters)
    colnames(Cluster1)<-c("Cluster")
    ClusterGeneList<-rownames(subset(Cluster1,Cluster=="TRUE"))
    length(ClusterGeneList)
    Clusters_ordered<-rbind(Clusters_ordered,Zscore_Clusters_Seq[ClusterGeneList,])
    rowannotation<-data.frame(c(rep(paste("cluster", colnameforcluster), length(ClusterGeneList))))
    rowannotation2<-rbind(rowannotation2,rowannotation)
  }
  else {
    colnameforcluster<-c("other")
    Cluster1<-as.data.frame(Clusters[,clusterorder])
    rownames(Cluster1)<-rownames(Clusters)
    colnames(Cluster1)<-c("Cluster")
    ClusterGeneList<-rownames(subset(Cluster1,Cluster=="TRUE"))
    length(ClusterGeneList)
    Clusters_ordered<-rbind(Clusters_ordered,Zscore_Clusters_Seq[ClusterGeneList,])
    rowannotation<-data.frame(c(rep(paste("cluster", colnameforcluster), length(ClusterGeneList))))
    rowannotation2<-rbind(rowannotation2,rowannotation)
  }
}

colnames(rowannotation2)<-c("Clusters")
rownames(rowannotation2)<-rownames(Clusters_ordered)

colors_annotate<-colorRampPalette(brewer.pal(11,"Set3"))

mycolors <- colors_annotate(length(unique(rowannotation2$Clusters)))
mycolors<-c("#8DD3C7",  "#BEBADA" ,"#FFFFB3", "#80B1D3", "#FB8072","#FDB462" ,"#B3DE69", "#FCCDE5", "#BC80BD","#79AB94" ,"#CCEBC5")
names(mycolors) <- unique(rowannotation2$Clusters)
mycolors
mycolors <- list(Clusters = mycolors)
mycolors

pheatmap(as.matrix(Clusters_ordered),color = myColor,breaks = myBreaks, treeheight_row = 0,treeheight_col = 0,show_rownames = FALSE, cluster_cols = FALSE,cluster_rows = FALSE,annotation_row = rowannotation2,annotation_colors = mycolors,labels_col = c("-2 hr","0 hr","2 hr","4 hr","6 hr","9 hr","12 hr","1 d","2 d","4 d","8 d","12 d"),cellwidth = 15,main="Differentially expressed genes (FDR< 10^-30)") 


#this is the heatmap if I didn't do any of that ordering...
pheatmap(as.matrix(Zscore_Clusters_Seq),treeheight_row = 0,treeheight_col = 0,show_rownames = FALSE, cluster_cols = FALSE,cluster_rows = TRUE)


Zscore_Clusters_Seq2<-Zscore_Clusters_Seq
Zscore_Clusters_Seq2$maxval<-do.call(pmax,Zscore_Clusters_Seq2)
head(Zscore_Clusters_Seq2)
Zscore_Clusters_Seq2$maxtime<-colnames(Zscore_Clusters_Seq2)[max.col(Zscore_Clusters_Seq2,ties.method="first")]



#Pull in the gene lists which are found on the GeneLists_WBID tab of Supp File 1
#GeneLists<-read.csv("GeneLists_TranscriptionalRegulators_WBID.csv",header =T )


head(GeneLists)
new_df3<-c()
plotlist<-list()
for (i in c(40,41,10,11,3,4,7,8)){ #how many columns are in the file
  GeneList_sub<-GeneLists[,i]
  GeneList_sub<-as.data.frame(GeneList_sub)
  #remove NAs and empty values from the dataframe
  GeneList_sub <- GeneList_sub[!apply(is.na(GeneList_sub) | GeneList_sub == "", 1, all),]
  GeneList_sub<-as.data.frame(GeneList_sub)
  GeneList_sub<-distinct(GeneList_sub,GeneList_sub)
  rownames(GeneList_sub)<-GeneList_sub$GeneList_sub
  venn1<-venn(list(ClusteredGenes=rownames(Cluster_WBID),GeneList=rownames(GeneList_sub)),show.plot = FALSE)
  venn1_overlap<-attr(venn1,"intersections")$`ClusteredGenes:GeneList`
  venn2<-venn(list(ExpressedGenes=rownames(ExpressedGenes),GeneList=rownames(GeneList_sub)),show.plot = FALSE)
  venn2_overlap<-attr(venn2,"intersections")$`ExpressedGenes:GeneList`
  genes_of_interest_dataframe<-as.data.frame(Zscore_WBID[venn1_overlap,])
  row.has.na <- apply(genes_of_interest_dataframe, 1, function(x){any(is.na(x))})
  genes_of_interest_dataframe<-genes_of_interest_dataframe[!row.has.na,]
  head(genes_of_interest_dataframe)
  genes_of_interest_dataframe$genes<-rownames(genes_of_interest_dataframe)
  genes_of_interest_dataframe<-melt(genes_of_interest_dataframe)
  head(genes_of_interest_dataframe)
  genes_of_interest_dataframe<-genes_of_interest_dataframe[order(genes_of_interest_dataframe$genes),]
  gene_averages<-genes_of_interest_dataframe
  gene_averages<-unique(gene_averages)
  head(gene_averages,20)
  gene_averages$hours<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
  gene_averages$samplenumber<-c(1,2,3,4,5,6,7,8,9,10,11,12)
  gene_averages<-gene_averages[,1:4]
  colnames(gene_averages)<-c("genes","variable","Zscore","hours")
  gene_averages_cast<-dcast(gene_averages,genes ~ hours, value.var="Zscore")
  rownames(gene_averages_cast)<-gene_averages_cast$genes
  gene_averages_cast<-gene_averages_cast[,-1]
  ClusterPlot<-ggplot(gene_averages,aes(x=hours,y=Zscore))+
    #geom_line(aes(group=genes),alpha=0.1)+
    theme_classic(base_size = 6)+theme(aspect.ratio = 1)+
    stat_summary(fun.data=mean_se,geom="line")+
    stat_summary(fun.data=mean_cl_boot,fun.args=(conf.int=0.99),size=0.075)+
    labs(y="Z-score",x="time point")+ggtitle(colnames(GeneLists)[i])
  plotlist[[i]] = ClusterPlot
  
  for (j in 1:129){
    Cluster1<-as.data.frame(Cluster_WBID[,j])
    rownames(Cluster1)<-rownames(Cluster_WBID)
    colnames(Cluster1)<-c("Cluster")
    ClusterGeneList<-rownames(subset(Cluster1,Cluster=="TRUE"))
    venn3<-venn(list(ClusterGeneList,venn1_overlap),show.plot = FALSE)
    overlap_length<-length(attr(venn3,"intersections")$`A:B`)
    phyper_venn<-phyper(overlap_length-1,length(ClusterGeneList),16699-length(ClusterGeneList),length(venn2_overlap),lower.tail = F,log.p = F)
    new_df3<-rbind(new_df3,c(colnames(GeneLists)[i],length(venn2_overlap),length(venn1_overlap),length(ClusterGeneList),j[1],overlap_length,phyper_venn))
    
  }
  
}


grid.arrange(grobs=plotlist[c(40,41,10,11,3,4,7,8)],ncol=2)

#grid.arrange(
#  grobs = list(plot1,plot2,plot3,plot4,plot9),
#  widths = c(1,1,1,1,1),
#  layout_matrix = rbind(c(NA,NA,1,2),
#                       c(NA,3,4, 5))
#)

new_df3

new_df3<-as.data.frame(new_df3,stringsAsFactors = FALSE)

new_df4<-new_df3
new_df3[] <- lapply(new_df3, function(x) if(is.na(as.numeric(x[1]))) as.factor(x) else  as.numeric(x))
colnames(new_df3)<-c("GeneList","NumInBackground","NumInClusters","ClusterSize","ClusterNumber","NumInList","Phyper")
tail(new_df3)
new_df3_significant<-subset(new_df3,Phyper<0.00038) #bonferroni p-value of 0.05/129
new_df3_significant<-subset(new_df3,Phyper<0.01) 
tail(new_df3_significant,40)

tail(new_df3_significant,100)
write.table(new_df3_significant,"ClusterEnrichments_121020.txt",quote = F,sep = "\t")
write.table(new_df3,"ClusterEnrichments_all_122920.txt",quote = F,sep = "\t")
#ggplot(new_df3_significant,aes(x=reorder(as.factor(ClusterNumber),log10(Phyper)),y= -log10(Phyper),fill=GeneList))+
#  geom_bar(stat = "identity")+facet_grid(.~GeneList,scales="free_x")+theme_classic(base_size = 12)+
#  labs(x="Cluster Number")+coord_flip()


ggplot(new_df3_significant,aes(x= -log10(Phyper),y=reorder(as.factor(ClusterNumber),-ClusterNumber)))+
  geom_barh(stat="identity")+theme_classic(base_size = 12)+
  labs(x="-log10(p-value)",y="Cluster Number")+facet_wrap(vars(GeneList),nrow=2,scales = "free_y")


#color-coded plot for cluster enrichments
new_df3_10clust<-subset(new_df3,ClusterNumber<12)
colors_annotate<-c("#8DD3C7",  "#BEBADA" ,"#FFFFB3", "#80B1D3", "#FB8072","#FDB462" ,"#B3DE69", "#FCCDE5", "#BC80BD","#79AB94" ,"#CCEBC5")
#reorder the clusters based on heatmap ordering
new_df3_10clust$ClusterNumber<-factor(new_df3_10clust$ClusterNumber,levels = c("11","7","3","8","1","5","10","9","2","4","6"))

ggplot(new_df3_10clust,aes(x= -log10(Phyper),y= ClusterNumber,fill=ClusterNumber))+
  geom_barh(stat="identity")+theme_classic(base_size = 12)+scale_fill_manual(values=colors_annotate)+
  labs(x="-log10(p-value)",y="Cluster Number")+facet_wrap(vars(GeneList),nrow=3,scales = "free_y")+
  scale_y_discrete(limits = rev(levels(new_df3_10clust$ClusterNumber)))+xlim(0,45)

new_df3_10clust

#now to plot clusters 1-10 to go next to the heatmap

plotlist2<-list()
for (m in c(7,3,8,1,5,10,9,2,4,6)){ #how many clusters
  Cluster1<-as.data.frame(Clusters[,m])
  rownames(Cluster1)<-rownames(Clusters)
  colnames(Cluster1)<-c("Cluster")
  ClusterGeneList<-rownames(subset(Cluster1,Cluster=="TRUE"))
  genes_of_interest_dataframe<-as.data.frame(Zscores[ClusterGeneList,])
  row.has.na <- apply(genes_of_interest_dataframe, 1, function(x){any(is.na(x))})
  genes_of_interest_dataframe<-genes_of_interest_dataframe[!row.has.na,]
  head(genes_of_interest_dataframe)
  genes_of_interest_dataframe$genes<-rownames(genes_of_interest_dataframe)
  genes_of_interest_dataframe<-melt(genes_of_interest_dataframe)
  head(genes_of_interest_dataframe)
  genes_of_interest_dataframe<-genes_of_interest_dataframe[order(genes_of_interest_dataframe$genes),]
  gene_averages<-genes_of_interest_dataframe
  gene_averages<-unique(gene_averages)
  head(gene_averages,20)
  gene_averages$hours<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
  gene_averages$samplenumber<-c(1,2,3,4,5,6,7,8,9,10,11,12)
  gene_averages<-gene_averages[,1:4]
  colnames(gene_averages)<-c("genes","variable","Zscore","hours")
  gene_averages_cast<-dcast(gene_averages,genes ~ hours, value.var="Zscore")
  rownames(gene_averages_cast)<-gene_averages_cast$genes
  gene_averages_cast<-gene_averages_cast[,-1]
  ClusterPlot<-ggplot(gene_averages,aes(x=hours,y=Zscore))+
    geom_line(aes(group=genes),alpha=0.005)+
    theme_classic(base_size = 6)+theme(aspect.ratio = 1)+
    # stat_summary(fun.data=mean_se,geom="line")+
    #  stat_summary(fun.data=mean_cl_boot,fun.args=(conf.int=0.99),size=0.075)+
    labs(y="Z-score",x="time point")
  plotlist2[[m]] = ClusterPlot
}


grid.arrange(grobs=list(plotlist2[[7]],plotlist2[[3]],plotlist2[[8]],plotlist2[[1]],
                        plotlist2[[5]],plotlist2[[10]],plotlist2[[9]],plotlist2[[2]],plotlist2[[4]],plotlist2[[6]]),ncol=1)


