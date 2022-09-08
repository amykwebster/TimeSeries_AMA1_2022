
#Supplementary Figures!
library(ggplot2)
library(gridExtra)
library(reshape2)
library(dplyr)
library(gplots)

#Supp Fig 1 - plot clusters
#Bring in the Cluster and Z-score files from Supp File 1, as in the Figure 2 script
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries")
clustering_data_1<-read.csv("clustering_z_means.csv",header = T,row.names = "gene")
head(clustering_data_1)
Clustering<-read.csv("Clustering_6027genes_0.8_129clusters.csv",header = F)
#Clustering<-Clustering[1:6027,]

newdataframe<-c()
plotlist2<-list()
#by doing 1-85, include all clusters with more than 4 genes
#to plot all 129 clusters, change base_size, ncol, and plotlist2 range...
for (i in 1:85) {
  Clustering_subset<-subset(Clustering,Clustering[,i]=="TRUE")
  newdataframe<-rbind(newdataframe,nrow(Clustering_subset))
  sequence_names<-rownames(Clustering_subset)
  clustering_data<-clustering_data_1[sequence_names,]
  row.has.na <- apply(clustering_data, 1, function(x){any(is.na(x))})
  clustering_data<-clustering_data[!row.has.na,]
  sig_genes_goterms<-rownames(clustering_data)
  clustering_mean<-colMeans(clustering_data[,1:12])
  clustering_mean<-data.frame(clustering_mean)
  clustering_mean$samplenumber<-c(1,2,3,4,5,6,7,8,9,10,11,12)
  clustering_mean$hours<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
  clustering_data$genes<-rownames(clustering_data)
  clustering_data1<-melt(clustering_data)
  clustering_data1<-clustering_data1[order(clustering_data1$genes),]
  clustering_data1$hours<-c(-2,0,2,4,6,9,12,24,48,96,192,288)
  clustering_data1$samplenumber<-c(1,2,3,4,5,6,7,8,9,10,11,12)
  clustering_data3<-subset(clustering_data1,samplenumber=="1"|samplenumber=="2"|samplenumber=="3"|samplenumber=="4"|samplenumber=="5"|samplenumber=="6"|samplenumber=="7"|samplenumber=="8")
  clustering_mean3<-subset(clustering_mean,samplenumber=="1"|samplenumber=="2"|samplenumber=="3"|samplenumber=="4"|samplenumber=="5"|samplenumber=="6"|samplenumber=="7"|samplenumber=="8")
  head(clustering_data3,12)
  title1<-paste("Cl.",i,":",nrow(Clustering_subset),"Genes",sep = " ")
  plot2<-(ggplot(data=clustering_data1,aes(x=hours,y=value,factor=genes))+
    geom_line(alpha=0.05,size=0.15)+
    labs(y="Z-score norm.",x="Hours of L1 arrest")+
    ggtitle(title1)+theme_classic(base_size = 6)+theme(aspect.ratio = 1))

plotlist2[[i]] = plot2
}

grid.arrange(grobs=plotlist2[1:85],ncol=8)


#Supp Fig 2 - gene regulatory lists not in Fig 2
#Again, bring in the Clusters, Z-scores, and Starvation Time Series RNA-seq results from Supp File 1, as in the Figure 2 script
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries/Cluster_EnrichmentAnalysis/")
Clusters<-read.csv("Clustering_6027genes_0.8_129clusters.csv",header = F)
head(Clusters)
rownames(Clusters)<-Clusters$V1
Clusters<-Clusters[,-1]

Zscores<-read.csv("clustering_z_means.csv",header = T)
head(Zscores)
rownames(Zscores)<-Zscores$gene
Zscores<-Zscores[,-1]

#bring in the expressed genes
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

GeneLists<-read.csv("GeneLists_TranscriptionalRegulators_WBID.csv",header =T )
head(GeneLists)


head(GeneLists)
new_df3<-c()
plotlist<-list()
for (i in c(1,2,38,39,6,9,12,13,14,15)){ #how many columns are in the file
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
    theme_classic(base_size = 8)+theme(aspect.ratio = 1)+
    stat_summary(fun.data=mean_se,geom="line")+
    stat_summary(fun.data=mean_cl_boot,fun.args=(conf.int=0.99),size=0.075)+
    labs(y="Z-score",x="Hours of L1 arrest")+ggtitle(colnames(GeneLists)[i])
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

tail(new_df3)


grid.arrange(grobs=plotlist[c(1,2,38,39,6,9,12,13,14,15)],ncol=2)



new_df3

new_df3<-as.data.frame(new_df3,stringsAsFactors = FALSE)
new_df3[]
subset(new_df3,GeneList=="SKN1_upRNAi") 
subset(new_df3,GeneList=="Germline_4fold") #only cluster 4 marginally significant
subset(new_df3,GeneList=="Germline_8fold") #none of the clusters are significant
new_df4<-new_df3
new_df3[] <- lapply(new_df3, function(x) if(is.na(as.numeric(x[1]))) as.factor(x) else  as.numeric(x))
colnames(new_df3)<-c("GeneList","NumInBackground","NumInClusters","ClusterSize","ClusterNumber","NumInList","Phyper")
head(new_df3[904:906,])



#Supp Fig 3 - higher fold germline enrichment
head(GeneLists)
new_df3<-c()
new_df5<-c()
for (i in c(1,2,21,22,24,25,26)){ #how many columns are in the file
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
      if (phyper_venn < 0.05){
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

#to plot germlineENRICH, germline 2x, 4x, 8x
# chose i in c(1,2,21,22,24,25,26) and set p<0.05
ggplot(subset(melted_df5,DEcluster!="not DE" & DEcluster!="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster))+
  scale_color_manual(values=colors_annotate(11))+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="DE cluster >10"),aes(hours,value,factor=WBID,color=DEcluster),alpha=0.2)+
  geom_line(inherit.aes = FALSE,data=subset(melted_df5,DEcluster=="not DE"),aes(hours,value,factor=WBID),alpha=0.01)+
  geom_line(alpha=0.5)+theme_classic(base_size = 14)+theme(aspect.ratio = 1)+
  labs(x="Hours of L1 arrest",y="Z-score")+
  facet_wrap(vars(geneset),ncol=7)



#This is now part of regular Fig 4.. 
#hatching efficiency of soma strain, part of Supp File 2 (Figure 4B)
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/FollowUp/")
ama1_dose<-read.csv("ama1_CA1200_doseresponse.csv",header = T)
head(ama1_dose,50)
ama1_dose2<-ggplot(ama1_dose,aes(x=log10(dose+0.0001),y=prop_L1s,color=strain,shape=strain,factor=replicate))+
  geom_point(size=3)+geom_line()+theme_classic(base_size=20)+
  scale_color_manual(values=c("#56B4E9", "#999999","#999999","#999999"))+
  labs(y="hatching efficiency 1 d post bleach",x="log10(auxin dose + 1e-4)")+
  theme(aspect.ratio = 1)+scale_x_continuous(breaks=c(-4,-3,-2,-1,0),labels = c("0.0001","0.001","0.01","0.1","1"))+
  labs(x="auxin dose (mM)")

ama1_dose2

ama1_1mMdose<-subset(ama1_dose,dose=="1")
ama1_1mMdose
pairwise.t.test(ama1_1mMdose$prop_L1s,ama1_1mMdose$strain,p.adjust.method = "none")

ama1_0.1mMdose<-subset(ama1_dose,dose=="0.1")
pairwise.t.test(ama1_0.1mMdose$prop_L1s,ama1_0.1mMdose$strain,p.adjust.method = "none")

ama1_0.01mMdose<-subset(ama1_dose,dose=="0.01")
pairwise.t.test(ama1_0.01mMdose$prop_L1s,ama1_0.01mMdose$strain,p.adjust.method = "none")

ama1_0.001mMdose<-subset(ama1_dose,dose=="0.001")
pairwise.t.test(ama1_0.001mMdose$prop_L1s,ama1_0.001mMdose$strain,p.adjust.method = "none")

ama1_0mMdose<-subset(ama1_dose,dose=="0")
pairwise.t.test(ama1_0mMdose$prop_L1s,ama1_0mMdose$strain,p.adjust.method = "none")

#Supplementary Figure 4C, data part of Supplementary File 2
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries/")
Amanitin_SS<-read.csv("AlphaAmanitin_HighDose_SS.csv",header = T)
Amanitin_SS
Protein_ama1$hours2<-factor(Protein_ama1$hours2,levels = c("36 hr", "132 hr +EtOH", "132 hr +aux"))
Amanitin_SS$amanitin<-factor(Amanitin_SS$amanitin,levels = c("none","alpha-amanitin"))
ggplot(Amanitin_SS,aes(x=amanitin,y=prop_alive))+
  geom_point(size=3)+facet_grid(.~day)+theme_classic(base_size = 14)+
  labs(x="Treatment",y="Proportion alive",tag="C")+theme(aspect.ratio = 1)+
  ggtitle("25 ug/mL alpha-amanitin added 12 hours after hatch to N2")



#Supplementary Fig 5

setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq/")

RNAyield<-read.csv("RNAyield_RNAseq.csv",header = T)
head(RNAyield,12)
RNAyield$Day_Treatment<-factor(RNAyield$Day_Treatment,levels = c("d2", "d6_EtOH", "d6_aux"))
TotalRNA1<-ggplot(RNAyield,aes(x=Day_Treatment,y=Mass.L1.pg.))+
  geom_boxplot()+labs(x="Hours of L1 arrest",y="Total RNA per L1 (pg)")+labs(tag="B")+
  geom_point(size=2,alpha=0.5)+
  theme_classic(base_size=14)+theme(aspect.ratio = 1)+ylim(0,45)+theme(aspect.ratio = 1)+
  facet_grid(.~Strain)+
  theme(aspect.ratio=0.75,axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")+
  scale_x_discrete(labels=c("d2" = "36", "d6_EtOH" = "132 hr +EtOH", "d6_aux" = "132 hr +aux"))

TotalRNA1

RNA_model7<-lme(Mass.L1.pg.~Day,random=~1|Rep,data = subset(RNAyield,Total.Q.!="NA"))
RNA_tval<-summary(RNA_model7)
RNA_tval$tTable

pairwise.t.test(RNAyield$Total.Q.,RNAyield$Day_Treatment)

#protein
Protein_ama1<-read.csv("Protein_Conc_ama1.csv",header = T)
head(Protein_ama1,9)
Protein_ama1$hours<-c("36 h","132 h","132 h","36 h","132 h","132 h","36 h","132 h","132 h")
Protein_ama1$hours2<-c("36 hr","132 hr +EtOH","132 hr +aux","36 hr","132 hr +EtOH","132 hr +aux","36 hr","132 hr +EtOH","132 hr +aux")


Protein_ama1$hours<-factor(Protein_ama1$hours,levels = c("36 h", "132 h"))
Protein_ama1$hours2<-factor(Protein_ama1$hours2,levels = c("36 hr", "132 hr +EtOH", "132 hr +aux"))

Protein_ama1$strain<-"Peft-3::TIR1;ama-1::AID"
ProteinYield1<-ggplot(Protein_ama1,aes(x=hours2,y=Protein_perL1_ng))+
  geom_boxplot()+labs(x="Hours of L1 arrest",y="Protein per L1 (ng)")+labs(tag="C")+
  geom_point(size=2,alpha=0.5)+facet_grid(.~strain)+
  theme_classic(base_size=14)+ylim(0,5)+
  theme(aspect.ratio=0.75,axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
  
ProteinYield1

protein_model<-lme(Protein_perL1_ng~ID,random=~1|Rep,data = Protein_ama1)
summary(protein_model)

head(Protein_ama1)
Protein_ama1_2<-subset(Protein_ama1,hours!="36 h")
pairwise.t.test(Protein_ama1_2$Protein_perL1_ng,Protein_ama1_2$ID,paired = TRUE)
protein_model2<-lme(Protein_perL1_ng~ID,random=~1|Rep,data = subset(Protein_ama1,hours!="36 h"))
summary(protein_model2)

setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq")
ama1_d2d6survival<-read.csv("ama1_survivald2d6_RNAseq.csv",header = T)
head(ama1_d2d6survival)
ama1_d2d6survival$hours<-factor(ama1_d2d6survival$hours,levels = c("36 hr", "132 hr +EtOH", "132 hr +aux"))
d2d6Survival<-ggplot(ama1_d2d6survival,aes(hours,prop_alive))+
  geom_point(size=2,alpha=0.5)+facet_grid(.~strain)+
  theme_classic(base_size=14)+theme(aspect.ratio = 0.75,axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(x="Hours of L1 arrest",y="Proportion alive",tag="A")+ylim(0,1)
d2d6Survival

grid.arrange(d2d6Survival,TotalRNA1,ProteinYield1,ncol=2)

