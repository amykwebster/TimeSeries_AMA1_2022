#Spike-in analysis with absolute normalization, Figure 6H-I and Supplementary Figure 6
library(dplyr)
library(reshape2)
library(gridExtra)
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq/")

#Spike-in concentrations available from the manufacturer
SpikeIn_conc<-read.csv("cms_095046.csv",header = T)

head(SpikeIn_conc)

setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq/WS273_counts_namesort/")
#import all count files, then put them in a list. Assign original names to the objects in the list.
ama1_counts= list.files(pattern= ".counts$",full.names = TRUE)
ama1_counts_files=lapply(ama1_counts,read.table, header=F)
names(ama1_counts_files)<-ama1_counts

for (i in 1:48){
  colnames(ama1_counts_files[[i]])<-c("gene",names(ama1_counts_files[i]))
}
colnames(ama1_counts_files[[1]])

head(ama1_counts_files$`./ama1_CA1200_d2_14_S14_merge_name_sorted.counts`)
ama1_counts_df<-data.frame(Reduce(cbind,ama1_counts_files))
head(ama1_counts_df)

ama1_counts_df<-ama1_counts_df[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96)]
head(ama1_counts_df)
rownames(ama1_counts_df)<-ama1_counts_df$gene
ama1_counts_df<-ama1_counts_df[,-1]
ERCC<-read.csv("ERCC.csv",header = F)
rownames(ERCC)<-ERCC$V1
ama1_counts_ERCC<-ama1_counts_df[rownames(ERCC),]
ama1_counts_ERCC
ama1_counts_withoutERCC = ama1_counts_df[which(rownames(ama1_counts_df) %in% rownames(ERCC)), ]
SumCounts_ERCC<-data.frame(colSums(ama1_counts_ERCC),colSums(ama1_counts_withoutERCC))
SumCounts_ERCC$Ratio<-SumCounts_ERCC$colSums.ama1_counts_ERCC. / (SumCounts_ERCC$colSums.ama1_counts_withoutERCC. + SumCounts_ERCC$colSums.ama1_counts_ERCC.)

counts_strains<-c("ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200","ama1_CA1200",
                  "ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352","ama1_CA1352",
                  "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200", "CA1200",
                  "CA1352",  "CA1352", "CA1352", "CA1352", "CA1352", "CA1352", "CA1352", "CA1352", "CA1352", "CA1352", "CA1352", "CA1352")

counts_timepoints<-c("d2","d2","d2","d2","d6_aux","d6_aux","d6_aux","d6_aux","d6_EtOH","d6_EtOH","d6_EtOH","d6_EtOH",
                     "d2","d2","d2","d2","d6_aux","d6_aux","d6_aux","d6_aux","d6_EtOH","d6_EtOH","d6_EtOH","d6_EtOH",
                     "d2","d2","d2","d2","d6_aux","d6_aux","d6_aux","d6_aux","d6_EtOH","d6_EtOH","d6_EtOH","d6_EtOH",
                     "d2","d2","d2","d2","d6_aux","d6_aux","d6_aux","d6_aux","d6_EtOH","d6_EtOH","d6_EtOH","d6_EtOH")


counts_timepoints2<-c("d2","d2","d2","d2","d6","d6","d6","d6","d6","d6","d6","d6",
                      "d2","d2","d2","d2","d6","d6","d6","d6","d6","d6","d6","d6",
                      "d2","d2","d2","d2","d6","d6","d6","d6","d6","d6","d6","d6",
                      "d2","d2","d2","d2","d6","d6","d6","d6","d6","d6","d6","d6")
replicates<-c("rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4")
counts_groups<-factor(paste(counts_strains,counts_timepoints,replicates,sep = "."))
counts_groups<-as.character(counts_groups)
counts_groups
SumCounts_ERCC$Condition<-counts_groups
SumCounts_ERCC$Replicate<-c("rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4")
SumCounts_ERCC$Strains<-counts_strains
SumCounts_ERCC$Timepoint<-counts_timepoints
SumCounts_ERCC$Timepoint2<-counts_timepoints2


colnames(ama1_counts_ERCC)<-counts_groups
head(ama1_counts_ERCC)
head(SpikeIn_conc)

ama1_counts_df_MergeAbsSpike<-merge(SpikeIn_conc,ama1_counts_ERCC,by.x = "ERCC.ID",by.y = 0)
head(ama1_counts_df_MergeAbsSpike)
ama1_df_absSpike_melt<-melt(ama1_counts_df_MergeAbsSpike,id.vars = c("ERCC.ID","Re.sort.ID","subgroup","concentration.in.Mix.1..attomoles.ul.","concentration.in.Mix.2..attomoles.ul.","expected.fold.change.ratio","log2.Mix.1.Mix.2."))
head(ama1_df_absSpike_melt)
ama1_df_absSpike_melt_rmvZero<-subset(ama1_df_absSpike_melt,value>0)
head(ama1_df_absSpike_melt_rmvZero)
ggplot(ama1_counts_df_MergeAbsSpike,aes(x=log10(ama1_CA1200.d6_aux),y=log10(concentration.in.Mix.1..attomoles.ul.)))+
  geom_point()+theme(aspect.ratio = 1)+xlim(-5,12)+ylim(-5,12)

ggplot(ama1_counts_df_MergeAbsSpike,aes(x=log10(ama1_CA1200.d6_EtOH),y=log10(concentration.in.Mix.1..attomoles.ul.)))+
  geom_point()+theme(aspect.ratio = 1)+xlim(-5,12)+ylim(-5,12)

ggplot(ama1_df_absSpike_melt,aes(x=log10(value),y=log10(concentration.in.Mix.1..attomoles.ul.)))+
  geom_point()+theme(aspect.ratio = 1)+xlim(-5,12)+ylim(-5,12)+facet_wrap(vars(variable),ncol=8)

head(ama1_df_absSpike_melt)

#Supplementary Figure 6
ggplot(ama1_df_absSpike_melt_rmvZero,aes(x=log10(value),y=log10(concentration.in.Mix.1..attomoles.ul.)))+
  geom_point(shape=1)+xlim(-2,6)+ylim(-2,6)+facet_wrap(vars(variable),ncol=8)+
  theme_classic(base_size = 5)+theme(aspect.ratio = 1)+
  labs(x="log10(counts)",y="log10(conc Mix 1 attomoles/uL)")

#To obtain the graphs in Figure 6H-I:
#Fit a line to the graphs in supplementary figure 6
summary1<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d2.rep1")))
summary2<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d2.rep2")))
summary3<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d2.rep3")))
summary4<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d2.rep4")))

summary5<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_aux.rep1")))
summary6<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_aux.rep2")))
summary7<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_aux.rep3")))
summary8<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_aux.rep4")))

summary9<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_EtOH.rep1")))
summary10<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_EtOH.rep2")))
summary11<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_EtOH.rep3")))
summary12<-summary(lm(log10(concentration.in.Mix.1..attomoles.ul.)~log10(value),data = subset(ama1_df_absSpike_melt_rmvZero,variable=="ama1_CA1200.d6_EtOH.rep4")))


colnames(ama1_counts_df)<-counts_groups
head(ama1_counts_df)
head(ama1_counts_df_filter)
ama1_counts_df<-ama1_counts_df[rownames(ama1_counts_df_filter),]

summary1
summary1$coefficients[1,1]
summary1$coefficients[2,1]

#Pull the coefficients for slope and intercept from the lines above
#To normalize count data: y=counts*10^(slope+intercept) 
ama1_counts_df$ama1_CA1200.d2.rep1_NORM<-(ama1_counts_df$ama1_CA1200.d2.rep1*10^(summary1$coefficients[1,1]+summary1$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d2.rep2_NORM<-(ama1_counts_df$ama1_CA1200.d2.rep2*10^(summary2$coefficients[1,1]+summary2$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d2.rep3_NORM<-(ama1_counts_df$ama1_CA1200.d2.rep3*10^(summary3$coefficients[1,1]+summary3$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d2.rep4_NORM<-(ama1_counts_df$ama1_CA1200.d2.rep4*10^(summary4$coefficients[1,1]+summary4$coefficients[2,1]))

ama1_counts_df$ama1_CA1200.d6_aux.rep1_NORM<-(ama1_counts_df$ama1_CA1200.d6_aux.rep1*10^(summary5$coefficients[1,1]+summary5$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d6_aux.rep2_NORM<-(ama1_counts_df$ama1_CA1200.d6_aux.rep2*10^(summary6$coefficients[1,1]+summary6$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d6_aux.rep3_NORM<-(ama1_counts_df$ama1_CA1200.d6_aux.rep3*10^(summary7$coefficients[1,1]+summary7$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d6_aux.rep4_NORM<-(ama1_counts_df$ama1_CA1200.d6_aux.rep4*10^(summary8$coefficients[1,1]+summary8$coefficients[2,1]))

ama1_counts_df$ama1_CA1200.d6_EtOH.rep1_NORM<-(ama1_counts_df$ama1_CA1200.d6_EtOH.rep1*10^(summary9$coefficients[1,1]+summary9$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d6_EtOH.rep2_NORM<-(ama1_counts_df$ama1_CA1200.d6_EtOH.rep2*10^(summary10$coefficients[1,1]+summary10$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d6_EtOH.rep3_NORM<-(ama1_counts_df$ama1_CA1200.d6_EtOH.rep3*10^(summary11$coefficients[1,1]+summary11$coefficients[2,1]))
ama1_counts_df$ama1_CA1200.d6_EtOH.rep4_NORM<-(ama1_counts_df$ama1_CA1200.d6_EtOH.rep4*10^(summary12$coefficients[1,1]+summary12$coefficients[2,1]))
head(ama1_counts_df)

setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq/")
WS273_genenames<-read.csv("WS273_geneNames.csv",header = T)
ama1_counts_merge_biotype<-merge(ama1_counts_df,WS273_genenames,by.x = 0,by.y = "WB_id")
rownames(ama1_counts_merge_biotype)<-ama1_counts_merge_biotype$Row.names
GC_counts_byBiotype<-ama1_counts_merge_biotype[,-1]
head(ama1_counts_merge_biotype)
GC_counts_byBiotype<-GC_counts_byBiotype[,-61:-63]
head(GC_counts_byBiotype)
ama1_proteincoding<-subset(GC_counts_byBiotype,type=="protein_coding_gene") #use this! 
GC_counts_byBiotype_melt<-melt(GC_counts_byBiotype)
GC_byBiotype<-aggregate(GC_counts_byBiotype_melt$value,list(variable=GC_counts_byBiotype_melt$variable,type=GC_counts_byBiotype_melt$type),sum)


setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/RNAseq_timeseries/Cluster_EnrichmentAnalysis/")
GeneLists<-read.csv("GeneLists_TranscriptionalRegulators_WBID.csv",header = T)
head(GeneLists[,21:22])
SomaGenes<-as.data.frame(GeneLists[,21])
rownames(SomaGenes)<-SomaGenes$`GeneLists[, 21]`
GermlineGenes<-as.data.frame(GeneLists[,22])
rownames(GermlineGenes)<-GermlineGenes$`GeneLists[, 22]`


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

library(gplots)
GeneList_sub<-GeneLists[,21]
GeneList_sub<-as.data.frame(GeneList_sub)
#remove NAs and empty values from the dataframe
GeneList_sub <- GeneList_sub[!apply(is.na(GeneList_sub) | GeneList_sub == "", 1, all),]
GeneList_sub<-as.data.frame(GeneList_sub)
GeneList_sub<-distinct(GeneList_sub,GeneList_sub)
rownames(GeneList_sub)<-GeneList_sub$GeneList_sub
venn1<-venn(list(ProteinCodingGenes=rownames(ama1_proteincoding),GeneList=rownames(GeneList_sub)))
venn1_overlap<-attr(venn1,"intersections")$`ProteinCodingGenes:GeneList`
venn3<-venn(list(ClusteredGenes=rownames(Cluster_WBID),GeneList_background=venn1_overlap))
venn3_overlap<-attr(venn3,"intersections")$`ClusteredGenes:GeneList_background` #soma 4491

head(GeneLists)
GeneList_sub<-GeneLists[,22]
GeneList_sub<-as.data.frame(GeneList_sub)
#remove NAs and empty values from the dataframe
GeneList_sub <- GeneList_sub[!apply(is.na(GeneList_sub) | GeneList_sub == "", 1, all),]
GeneList_sub<-as.data.frame(GeneList_sub)
GeneList_sub<-distinct(GeneList_sub,GeneList_sub)
rownames(GeneList_sub)<-GeneList_sub$GeneList_sub
venn2<-venn(list(ProteinCodingGenes=rownames(ama1_proteincoding),GeneList=rownames(GeneList_sub)))
venn2_overlap<-attr(venn2,"intersections")$`ProteinCodingGenes:GeneList`
venn4<-venn(list(ClusteredGenes=rownames(Cluster_WBID),GeneList_background=venn2_overlap))
venn4_overlap<-attr(venn4,"intersections")$`ClusteredGenes:GeneList_background` #germline 948



ama1_proteincoding$ama1_CA1200.d2_NORM_Mean<-rowMeans(ama1_proteincoding[,49:52])
ama1_proteincoding$ama1_CA1200.d6_EtOH_NORM_Mean<-rowMeans(ama1_proteincoding[,57:60])
ama1_proteincoding$ama1_CA1200.d6_aux_NORM_Mean<-rowMeans(ama1_proteincoding[,53:56])

ama1_proteincoding2<-subset(ama1_proteincoding,ama1_CA1200.d2_NORM_Mean>10 | ama1_CA1200.d6_EtOH_NORM_Mean>10 | ama1_CA1200.d6_aux_NORM_Mean>10)
head(ama1_proteincoding2)
head(ama1_proteincoding2[venn3_overlap,63:65])

#Figure 6H,I
#Multiply by the volume of spike-in and divide by average number of worms needed to acquire the amount of total RNA used for each condition
AllReps_d2d6EtOH2<-ggplot()+
  geom_point(data=ama1_proteincoding2[venn3_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_EtOH_NORM_Mean*5/1873.6),color='soma enriched'),alpha=0.1)+
  geom_point(data=ama1_proteincoding2[venn4_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_EtOH_NORM_Mean*5/1873.6),color='germline enriched'),alpha=0.1)+
  #geom_abline(slope=1, intercept = 0,color="black")+
  geom_smooth(data=ama1_proteincoding2[venn3_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_EtOH_NORM_Mean*5/1873.6),color='soma enriched'),method='lm', formula= y~x)+
  geom_smooth(data=ama1_proteincoding2[venn4_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_EtOH_NORM_Mean*5/1873.6),color='germline enriched'),method='lm', formula= y~x)+
  theme_classic(base_size = 8)+theme(aspect.ratio = 1)+
  labs(x="log10(attomoles of transcript at 36 hr)",y='log10(attomoles of transcript at 132 hr)')

AllReps_d2d6EtOH2

AllReps_d2d6aux2<-ggplot()+
  geom_point(data=ama1_proteincoding2[venn3_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_aux_NORM_Mean*5/3901.4),color='soma enriched'),alpha=0.1)+
  geom_point(data=ama1_proteincoding2[venn4_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_aux_NORM_Mean*5/3901.4),color='germline enriched'),alpha=0.1)+
 # geom_abline(slope=1, intercept = 0,color="black")+
  geom_smooth(data=ama1_proteincoding2[venn3_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_aux_NORM_Mean*5/3901.4),color='soma enriched'),method='lm', formula= y~x)+
  geom_smooth(data=ama1_proteincoding2[venn4_overlap,],aes(x=log10(ama1_CA1200.d2_NORM_Mean*5/854.5),y=log10(ama1_CA1200.d6_aux_NORM_Mean*5/3901.4),color='germline enriched'),method='lm', formula= y~x)+
  theme_classic(base_size = 8)+theme(aspect.ratio = 1)+
  labs(x="log10(attomoles of transcript at 36 hr)",y='log10(attomoles of transcript at 132 hr with AMA-1/Pol II Depletion)')
AllReps_d2d6aux2



grid.arrange(AllReps_d2d6EtOH2,AllReps_d2d6aux2,ncol=2)


