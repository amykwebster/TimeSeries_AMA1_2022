#Code to generate parts of Figures 3 and 6 (Fig 3B-E, Fig 6A-C) which are all color-coded plots (by cluster) of Z-scores for DE genes 
#!!!All data used for these scripts should be found on sheets in Supplementary Files 1 and 2!!!

#SET working directory to where files are on your computer
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries/Cluster_EnrichmentAnalysis/")
library(eulerr)
library(gplots)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(RColorBrewer)


#Pull in the Clusters, Z-scores, time series gene expression results, and gene lists for enrichment analyses as done in the Figure 2 script.
#These files are tabs in Supplementary File 1 
Clusters<-read.csv("Clustering_6027genes_0.8_129clusters.csv",header = F)
head(Clusters)
rownames(Clusters)<-Clusters$V1
Clusters<-Clusters[,-1]

Zscores<-read.csv("clustering_z_means.csv",header = T)
head(Zscores)
rownames(Zscores)<-Zscores$gene
Zscores<-Zscores[,-1]

ExpressedGenes<-read.csv("StarvationTimeSeries_v2_resultspage.csv",header = T)
ExpressedGenes<-subset(ExpressedGenes,FDR!="NA")
ExpressedGenes<-ExpressedGenes[!duplicated(ExpressedGenes$primaryIdentifier),]
ExpressedGenes<-ExpressedGenes[complete.cases(ExpressedGenes),]
rownames(ExpressedGenes)<-ExpressedGenes$primaryIdentifier


GeneLists<-read.csv("GeneLists_TranscriptionalRegulators_WBID.csv",header =T )
head(GeneLists)


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



head(GeneLists[,42:43])
new_df3<-c()
new_df5<-c()
for (i in c(12,13,22,30,31,32,33,34,42,43)){ #how many columns are in the file
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
#cluster plot would go here
  
  for (j in 1:129){
    Cluster1<-as.data.frame(Cluster_WBID[,j])
    rownames(Cluster1)<-rownames(Cluster_WBID)
    colnames(Cluster1)<-c("Cluster")
    ClusterGeneList<-rownames(subset(Cluster1,Cluster=="TRUE"))
    venn3<-venn(list(ClusterGeneList,venn1_overlap),show.plot = FALSE)
    attr(venn3,"intersections")$`A:B`
    overlap_length<-length(attr(venn3,"intersections")$`A:B`)
    phyper_venn<-phyper(overlap_length-1,length(ClusterGeneList),16699-length(ClusterGeneList),length(venn2_overlap),lower.tail = F,log.p = F)
    new_df3<-rbind(new_df3,c(colnames(GeneLists)[i],length(venn2_overlap),length(venn1_overlap),length(ClusterGeneList),j[1],overlap_length,phyper_venn))
 if (overlap_length > 0){
     if (phyper_venn < 0.01){
       if (j<=10){
         overlap_genes<-attr(venn3,"intersections")$`A:B`
         overlap_zscore<-Zscore_WBID[overlap_genes,]
         overlap_zscore$geneset<- colnames(GeneLists)[i]
         overlap_zscore$DEcluster<- j
         new_df5<-rbind(new_df5,overlap_zscore)
       }
     else{
       overlap_genes<-attr(venn3,"intersections")$`A:B`
       overlap_zscore<-Zscore_WBID[overlap_genes,]
       overlap_zscore$geneset<- colnames(GeneLists)[i]
       overlap_zscore$DEcluster<-"DE cluster >10"
       new_df5<-rbind(new_df5,overlap_zscore)
     }
     }
     else{
       overlap_genes<-attr(venn3,"intersections")$`A:B`
       overlap_zscore<-Zscore_WBID[overlap_genes,]
       overlap_zscore$geneset<- colnames(GeneLists)[i]
       overlap_zscore$DEcluster<-"not DE"
       new_df5<-rbind(new_df5,overlap_zscore)
     }
   
 }   

    
  }
  
}

new_df5$WBID<-rownames(new_df5)
melted_df5<-melt(new_df5)
head(melted_df5)
melted_df5<-melted_df5[order(melted_df5$WBID),]
head(melted_df5)

melted_df5$hours<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
melted_df5$samplenumber<-c(1,2,3,4,5,6,7,8,9,10,11,12)


colors_annotate<-colorRampPalette(brewer.pal(11,"Set3")[c(4,7,10,2,9,6,11,3,5,8,1)]) #same colors as heatmap clusters
#colors_annotate<-colorRampPalette(brewer.pal(11,"Set3")[c(4,6,2,7,11,1,5,3,10,9,8)]) #same colors as heatmap clusters
#colors_annotate<-colorRampPalette(brewer.pal(11,"Set3")[c(2,1,10,4,11,3,9,6,7,5,8)]) 


#NOTE: I include extra gene lists because the color coding oonly works if each cluster is enriched at least once in the set
#Also by choosing ncol=7, I get a consistent size across various figures...

#Fig 3B and 3E
#to plot germlineENRICH, somaENRICH, PolII active, docked, and no signal
# chose i in c(12,13,18,19,20,21,22) and set p<0.01
ggplot(subset(melted_df5,DEcluster!="not DE" & DEcluster!="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster))+
  scale_color_manual(values=colors_annotate(11))+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster),alpha=0.2)+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="not DE"),aes(hours,value,factor=WBID),alpha=0.01)+
  geom_line(alpha=0.5)+theme_classic(base_size = 14)+theme(aspect.ratio = 1)+
  facet_wrap(vars(geneset),ncol=7)

#Fig 3C and 3D
#c(12,13,28,29,35,36,37) to plot PGC Seydoux dataset, set p<0.05
ggplot(subset(melted_df5,DEcluster!="not DE" & DEcluster!="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster))+
  scale_color_manual(values=colors_annotate(11))+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster),alpha=0.2)+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="not DE"),aes(hours,value,factor=WBID),alpha=0.1)+
  geom_line(alpha=0.5)+theme_classic(base_size = 14)+theme(aspect.ratio = 1)+
  facet_wrap(vars(geneset),ncol=7)

#Fig 6A-C
#c(12,13,22,30,31,32,33,34,42,43) to plot ama-1 dependent and independent genes, set p<0.01
ggplot(subset(melted_df5,DEcluster!="not DE" & DEcluster!="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster))+
  scale_color_manual(values=colors_annotate(11))+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster),alpha=0.2)+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="not DE"),aes(hours,value,factor=WBID),alpha=0.01)+
  geom_line(alpha=0.5)+theme_classic(base_size = 14)+theme(aspect.ratio = 1)+
  facet_wrap(vars(geneset),ncol=7)


grid.arrange(grobs=plotlist[c(28,29,36,37)],ncol=2)




new_df3<-as.data.frame(new_df3,stringsAsFactors = FALSE)
new_df3[]
head(new_df3)
levels(factor(new_df3$GeneList))
 
new_df3[] <- lapply(new_df3, function(x) if(is.na(as.numeric(x[1]))) as.factor(x) else  as.numeric(x))
colnames(new_df3)<-c("GeneList","NumInBackground","NumInClusters","ClusterSize","ClusterNumber","NumInList","Phyper")
head(new_df3)
ama1_2679_enrichedclusters<-subset(new_df3,GeneList=="ama1_dep_downinaux" & Phyper<0.05)
highconfup_1091clusters<-subset(new_df3,GeneList=="HighConfUp_1091" & Phyper<0.05)
highconfup_1091clusters
commonclusters<-venn(list(ama1_2679_enrichedclusters$ClusterNumber,highconfup_1091clusters$ClusterNumber))
attr(commonclusters,"intersections")$`A:B`

#determine whether soma and germline genes overlap with active or docked pol II

head(GeneLists)
colnames(GeneLists)[42]
GeneList_sub2<-GeneLists[,21] #soma enriched
GeneList_sub2<-GeneLists[,22] #germline enriched
GeneList_sub2<-GeneLists[,43]
GeneList_sub2<-as.data.frame(GeneList_sub2)
#remove NAs and empty values from the dataframe
GeneList_sub2 <- GeneList_sub2[!apply(is.na(GeneList_sub2) | GeneList_sub2 == "", 1, all),]
GeneList_sub2<-as.data.frame(GeneList_sub2)

GeneList_sub2<-distinct(GeneList_sub2,GeneList_sub2)
rownames(GeneList_sub2)<-GeneList_sub2$GeneList_sub2


new_df8<-c()
for (i in 19:22){ #how many columns are in the file
  GeneList_sub<-GeneLists[,i]
  GeneList_sub<-as.data.frame(GeneList_sub)
  #remove NAs and empty values from the dataframe
  GeneList_sub <- GeneList_sub[!apply(is.na(GeneList_sub) | GeneList_sub == "", 1, all),]
  GeneList_sub<-as.data.frame(GeneList_sub)
  GeneList_sub<-distinct(GeneList_sub,GeneList_sub)
  rownames(GeneList_sub)<-GeneList_sub$GeneList_sub

  venn1<-venn(list(ClusteredGenes=rownames(Cluster_WBID),GeneList=rownames(GeneList_sub)),show.plot = FALSE)
  venn1_overlap<-attr(venn1,"intersections")$`ClusteredGenes:GeneList`
  venn2<-venn(list(ClusteredGenes=rownames(Cluster_WBID),GeneList2=rownames(GeneList_sub2)),show.plot = FALSE)
  venn2_overlap<-attr(venn2,"intersections")$`ClusteredGenes:GeneList2`
  venn3<-venn(list(venn1_overlap,venn2_overlap), show.plot = FALSE)
  venn3_overlap<-attr(venn3,"intersections")$`A:B`
  overlap_length<-length(venn3_overlap)
  phyper_venn<-phyper(overlap_length-1,length(venn1_overlap),dim(Cluster_WBID)[1]-length(venn1_overlap),length(venn2_overlap),lower.tail = F,log.p = F)
  new_df8<-rbind(new_df8,c(colnames(GeneLists)[i],dim(GeneList_sub2)[1],dim(GeneList_sub)[1],length(venn1_overlap),length(venn2_overlap),overlap_length,phyper_venn))
  
}


new_df8

#make four way venn for germline, soma, docked, and active
venn_numbers<-c(Germline=821,PolII_Docked=218,PolII_Active=382,"Germline&PolII_Active"=46,"Germline&PolII_Docked"=81)
venn_numbers2<-c(Soma=3923,PolII_Docked=83,PolII_Active=63,"Soma&PolII_Docked"=216,"Soma&PolII_Active"=365)

set.seed(5)
fit1<-euler(venn_numbers,fills="blue")
Overlap_Germline_DockedActive<-plot(fit1,quantities=TRUE)
Overlap_Germline_DockedActive

set.seed(2)
fit2<-euler(venn_numbers2,fills="blue")
Overlap_Soma_DockedActive<-plot(fit2,quantities=TRUE)
Overlap_Soma_DockedActive

grid.arrange(Overlap_Germline_DockedActive,Overlap_Soma_DockedActive,ncol=2)





#Define gene sets of high confidence up genes and ama-1 dependent genes that are up-regulated in the original time series
head(GeneLists[,30:31])
GeneList_sub<-GeneLists[,30] #ama-1 dependent up, 304 genes, change to 31 for ama-1 independent,296 genes
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

ama1DepUp<-subset(genes_of_interest_dataframe,j_4d > i_2d | k_8d > i_2d | l_12d > i_2d | i_2d > h_1d | j_4d > h_1d | k_8d > h_1d | l_12d > h_1d | k_8d > j_4d | l_12d > j_4d | l_12d > k_8d) #477 / 539 included in Z-score. 252/304

ama1IndUp<-subset(genes_of_interest_dataframe,j_4d > i_2d | k_8d > i_2d | l_12d > i_2d | i_2d > h_1d | j_4d > h_1d | k_8d > h_1d | l_12d > h_1d | k_8d > j_4d | l_12d > j_4d | l_12d > k_8d) #876 / 910 included in Z-score.  276/296 clustered

write.table(ama1IndUp,"ama1IndUp_Clustered_876.txt",quote = F,sep = "\t")
write.table(ama1DepUp,"ama1DepUp_Clustered_477.txt",quote = F,sep = "\t")

