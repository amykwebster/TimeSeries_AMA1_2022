#Code to generate parts of Figures 5 and 6 (Fig 5B,D,F, Fig 6E-I) which include analysis of the ama-1 RNA-seq data
#Includes the analysis to generate AMA-1 RNA-seq results in Supp File 1

library(gplots)
library(reshape2)
library(ggplot2)
library(RUVSeq)
library(Hmisc)
library(eulerr)
library(gridExtra)

#Bring in the Counts_ama-1 data available in Supp File 1. I have already merged that data into a single sheet, whereas in this script I imported each count file separately then merged
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
ama1_counts_withoutERCC = ama1_counts_df[which(rownames(ama1_counts_df) %nin% rownames(ERCC)), ]
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

counts_groups<-factor(paste(counts_strains,counts_timepoints,sep = "."))
counts_groups<-as.character(counts_groups)
SumCounts_ERCC$Condition<-counts_groups
SumCounts_ERCC$Replicate<-c("rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4","rep1","rep2","rep3","rep4")
SumCounts_ERCC$Strains<-counts_strains
SumCounts_ERCC$Timepoint<-counts_timepoints
SumCounts_ERCC$Timepoint2<-counts_timepoints2



#total RNA, Figure 6F (Supp File 2)
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq/")

RNAyield<-read.csv("RNAyield_RNAseq.csv",header = T)
head(RNAyield)
TotalRNA<-ggplot(RNAyield,aes(x=Day,y=Mass.L1.pg.))+
  geom_boxplot()+labs(x="Hours of L1 arrest",y="Total RNA per L1 (picograms)")+labs(tag="E")+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=0.5,binwidth=1.5,stackdir = "center")+
  theme_classic(base_size=18)+theme(aspect.ratio = 1)+ylim(0,45)+theme(aspect.ratio = 1)+
  scale_x_discrete(labels=c("d2" = "36", "d6" = "132"))

TotalRNA

RNA_model7<-lme(Mass.L1.pg.~Day,random=~1|Rep,data = subset(RNAyield,Total.Q.!="NA"))
RNA_tval<-summary(RNA_model7)
RNA_tval$tTable

#protein, Figure 6G (Supp FIle 2)
Protein_ama1<-read.csv("Protein_Conc_ama1.csv",header = T)
head(Protein_ama1,9)
Protein_ama1$hours<-c("36 h","132 h","132 h","36 h","132 h","132 h","36 h","132 h","132 h")

Protein_ama1$hours<-factor(Protein_ama1$hours,levels = c("36 h", "132 h"))


ProteinYield<-ggplot(Protein_ama1,aes(x=hours,y=Protein_perL1_ng,color=Protein_ama1$Rep))+
  geom_boxplot()+labs(x="Hours of L1 arrest",y="Protein per L1 (ng)")+labs(tag="E")+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=0.5,binwidth=0.25,stackdir = "center")+
  theme_classic(base_size=18)+theme(aspect.ratio = 1)+ylim(0,5)+theme(aspect.ratio = 1)

ProteinYield

protein_model<-lme(Protein_perL1_ng~ID,random=~1|Rep,data = Protein_ama1)
summary(protein_model)


#mRNA (spike-in proportion), Figure 6E continued from analysis above
SumCounts_ERCC$Condition<-factor(SumCounts_ERCC$Condition,levels = c("ama1_CA1200.d2","ama1_CA1200.d6_EtOH","ama1_CA1200.d6_aux","CA1200.d2","CA1200.d6_EtOH","CA1200.d6_aux","ama1_CA1352.d2","ama1_CA1352.d6_EtOH","ama1_CA1352.d6_aux","CA1352.d2","CA1352.d6_EtOH","CA1352.d6_aux"))
Prop_ERCC<-ggplot(SumCounts_ERCC,aes(x=Condition,y=Ratio,color=Condition))+
  geom_boxplot()+
  geom_point(size=1,alpha=0.5)+
  theme_classic(base_size = 10)+ylim(0,0.05)+
  labs(y="Prop ERCC reads out of mapped reads")+labs(tag = "B")+
  scale_colour_manual(values=c("#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999","#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999"))+
  theme(aspect.ratio=0.5,axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")+
  scale_x_discrete(labels=c("36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux"))

Prop_ERCC

pairwise.t.test(SumCounts_ERCC$Ratio,SumCounts_ERCC$Condition,pool.sd = FALSE,p.adjust="none")

head(SumCounts_ERCC)
Prop_ERCC2<-ggplot(SumCounts_ERCC,aes(x=Timepoint2,y=Ratio))+
  geom_boxplot(aes(color=Condition))+
  geom_dotplot(inherit.aes = FALSE,data=SumCounts_ERCC,aes(x=Timepoint2,y=Ratio,fill=Condition),binaxis = "y",alpha=0.5,binwidth=0.001,stackdir = "center",stackgroups = TRUE,binpositions = "all")+
  theme_classic(base_size = 18)+ylim(0,0.05)+
  labs(y="Prop ERCC reads out of mapped reads")+labs(tag = "B")+
  #scale_colour_manual(values=c("#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999","#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999"))+
  theme(aspect.ratio=1,axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")+
  scale_x_discrete(labels=c("36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux"))

Prop_ERCC2

SpikeInModel<-lme(Ratio~Timepoint2,random=~1|Replicate,data = SumCounts_ERCC)
summary(SpikeInModel)

head(SumCounts_ERCC)


#PCA, Figure 5B
setwd("/Users/amywebster/Documents/BaughLab/StarvationTimeCourse/ama1_RNAseq/")
WS273_genenames<-read.csv("WS273_geneNames.csv",header = T) #Version WS273 of the womm genome
ama1_counts_merge_biotype<-merge(ama1_counts_df,WS273_genenames,by.x = 0,by.y = "WB_id")
rownames(ama1_counts_merge_biotype)<-ama1_counts_merge_biotype$Row.names
GC_counts_byBiotype<-ama1_counts_merge_biotype[,-1]
GC_counts_byBiotype<-GC_counts_byBiotype[,-49:-51]
ama1_proteincoding<-subset(GC_counts_byBiotype,type=="protein_coding_gene") #use this! 
GC_counts_byBiotype_melt<-melt(GC_counts_byBiotype)
GC_byBiotype<-aggregate(GC_counts_byBiotype_melt$value,list(variable=GC_counts_byBiotype_melt$variable,type=GC_counts_byBiotype_melt$type),sum)

head(ama1_proteincoding)
ama1_proteincoding2<-ama1_proteincoding[,-49]
filter <- apply(ama1_proteincoding2, 1, function(x) length(x[x>10])>=4) #approximately 1 cpm since there are roughly 10 million reads...
ama1_filtered <- ama1_proteincoding2[filter,]

genes <- rownames(ama1_filtered)[grep("^WBGene", rownames(ama1_filtered))]
spikes <- rownames(ama1_counts_df)[grep("^ERCC", rownames(ama1_counts_df))]

ama1_counts_df_filter<-ama1_counts_df[rownames(ama1_filtered),]
ama1_counts_df_spikes<-ama1_counts_df[spikes,]
ama1_counts_df_filterspikes<-rbind(ama1_counts_df_filter,ama1_counts_df_spikes) #15742 including genes and spikes

cpm_d2<-cpm(ama1_counts_df_filter)
cpm_d2_df<-data.frame(cpm_d2)
cpm_d2_df$mean<-rowMeans(cpm_d2_df) #calculate mean CPM across all libraries
#change the 1:10 depending on how many libraries there are
cpm_d2_df2<-cpm_d2_df[,1:48]/cpm_d2_df$mean #mean normalize
cpm_d2_df2<-log2(cpm_d2_df2+1) #log2 transform 
pca = prcomp(t(cpm_d2_df2)) #principal component analysis (PCA) on the log2 mean normalized CPM values
summary(pca)
pca$x
pca_genes<-pca$x
pca_genes_dataframe<-as.data.frame(pca_genes)
#change the conditions
conditions<-counts_groups
pca_genes_dataframe<-data.frame(conditions,pca_genes_dataframe)
pca_genes_dataframe2<-data.frame(pca_genes_dataframe,counts_strains)
pca_genes_dataframe2<-data.frame(pca_genes_dataframe2,counts_timepoints)
counts_fullcondition<-paste(counts_strains,counts_timepoints)
replicates<-c("rep1","rep2","rep3","rep4","rep5","rep1","rep2","rep3","rep4","rep5")
pca_genes_dataframe2$fullcondition<-counts_fullcondition

colorBlindGrey8   <- c("#56B4E9", "#E69F00", "#999999", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
scales::show_col(colorBlindGrey8)

PCA<-(ggplot(pca_genes_dataframe2,aes(x=PC1,y=PC2,colour=fullcondition,shape=counts_strains))+
         geom_point(size=3,alpha=0.8)+
         ggtitle("PCA of ama-1 RNAseq, 15,650 genes")+
         scale_color_manual(values=c("#009E73","#56B4E9", "#E69F00","#009E73","#56B4E9","#E69F00","#009E73","#999999","#CC79A7","#009E73","#999999","#CC79A7"))+
         labs(x="PC1 (25.2% of variance)",y="PC2 (17.8% of variance)")+
         theme_classic(base_size = 10)+
         theme(aspect.ratio = 1,legend.position = "none")+labs(tag="C"))
PCA


#table with number of DE genes



#boxplots of transcription-dependent genes over all 12 conditions 

counts_reps<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19",
               "20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36",
               "37","38","39","40","41","42","43","44","45","46","47","48")

counts_groups2<-factor(paste(counts_groups,counts_reps,sep = "."))
counts_groups2

colnames(ama1_counts_df_filterspikes)<-counts_groups2
x <- as.factor(rep(c("ama1_CA1200_d2","ama1_CA1200_aux_d6", "ama1_CA1200_EtOH_d6","ama1_CA1352_d2","ama1_CA1352_aux_d6","ama1_CA1352_EtOH_d6","CA1200_d2","CA1200_aux_d6","CA1200_EtOH_d6","CA1352_d2","CA1352_aux_d6","CA1352_EtOH_d6"), each=4))

head(ama1_counts_df_filterspikes)
set <- newSeqExpressionSet(as.matrix(ama1_counts_df_filterspikes),phenoData = data.frame(x, row.names=counts_groups2))
set1 <- betweenLaneNormalization(set, which="upper")
set3 <- RUVg(set1, spikes, k=1) #using the spike-ins to normalize the data

design <- model.matrix(~0 + x + W_1, data=pData(set3))
design
y <- DGEList(counts=counts(set3), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
design


#day 6 aux vs etoh in same strain
#start with ama-1/CA1200
design
lrt <- glmLRT(fit, contrast = c(1,0,-1,0,0,0,0,0,0,0,0,0,0))
lrt_sort<-topTags(lrt,n=nrow(lrt$table))$table
significant<-subset(lrt_sort,FDR<0.05) #5081 genes
sig_down<-subset(significant,logFC<0) #2679
cpm_sigdown<-as.data.frame(cpm_d2)[rownames(sig_down),]
head(cpm_sigdown)
cpm_sigdown$WBID<-rownames(cpm_sigdown)
head(cpm_sigdown)
cpm_sigdown_melted<-melt(cpm_sigdown)
tail(cpm_sigdown_melted,30)
cpm_sigdown_melted<-cpm_sigdown_melted[order(cpm_sigdown_melted$WBID),]
cpm_sigdown_melted$condition<-counts_groups
cpm_sigdown_melted_means<-aggregate(cpm_sigdown_melted$value,list(condition=cpm_sigdown_melted$condition,WBID=cpm_sigdown_melted$WBID),mean)
cpm_sigdown_melted_means$condition <- factor(cpm_sigdown_melted_means$condition,levels = c("ama1_CA1200.d2","ama1_CA1200.d6_EtOH","ama1_CA1200.d6_aux","CA1200.d2","CA1200.d6_EtOH","CA1200.d6_aux","ama1_CA1352.d2","ama1_CA1352.d6_EtOH","ama1_CA1352.d6_aux","CA1352.d2","CA1352.d6_EtOH","CA1352.d6_aux"))
head(cpm_sigdown_melted_means,100)

cpm_sigdown_melted_means <- cpm_sigdown_melted_means[grep("^WBGene", cpm_sigdown_melted_means$WBID),]
head(cpm_sigdown_melted_means) #Only include genes with WBGene start in WBID column! 2675



#Figure 5D
TranscriptionDependent_2679<-ggplot(cpm_sigdown_melted_means,aes(x=condition,y=log2(x+0.1),color=factor(condition)))+
  geom_boxplot(alpha=0.1)+
  theme_classic(base_size = 10)+labs(x="Condition",y="log2(CPM)")+
  ggtitle("Transcription-dependent (2,679 genes)")+
  scale_colour_manual(values=c("#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999","#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(aspect.ratio = 0.5,legend.position = "none")+labs(tag="E")
TranscriptionDependent_2679

pairwise.t.test(cpm_sigdown_melted_means$x,cpm_sigdown_melted_means$condition,p.adjust.method = "none",paired = TRUE)
cpm_sigdown_melted_means_subset1<-cpm_sigdown_melted_means_subset1[complete.cases(cpm_sigdown_melted_means_subset1),]
cpm_sigdown_melted_means_subset1<-cpm_sigdown_melted_means[complete.cases(cpm_sigdown_melted_means),]
cpm_sigdown_melted_means_subset2<-subset(cpm_sigdown_melted_means_subset1,condition=="CA1200.d6_EtOH")
cpm_sigdown_melted_means_subset3<-subset(cpm_sigdown_melted_means_subset1,condition=="CA1200.d6_aux")
cpm_sigdown_melted_means_subset4<-subset(cpm_sigdown_melted_means_subset1,condition=="ama1_CA1200.d6_EtOH")
cpm_sigdown_melted_means_subset5<-subset(cpm_sigdown_melted_means_subset1,condition=="ama1_CA1200.d6_aux")
cpm_sigdown_melted_means_subset6<-subset(cpm_sigdown_melted_means_subset1,condition=="ama1_CA1352.d6_EtOH")
cpm_sigdown_melted_means_subset7<-subset(cpm_sigdown_melted_means_subset1,condition=="ama1_CA1352.d6_aux")
cpm_sigdown_melted_means_subset8<-subset(cpm_sigdown_melted_means_subset1,condition=="CA1352.d6_EtOH")
cpm_sigdown_melted_means_subset9<-subset(cpm_sigdown_melted_means_subset1,condition=="CA1352.d6_aux")
cpm_sigdown_melted_means_subset10<-subset(cpm_sigdown_melted_means_subset1,condition=="CA1200.d2")
cpm_sigdown_melted_means_subset11<-subset(cpm_sigdown_melted_means_subset1,condition=="ama1_CA1200.d2")
cpm_sigdown_melted_means_subset12<-subset(cpm_sigdown_melted_means_subset1,condition=="CA1352.d2")
cpm_sigdown_melted_means_subset13<-subset(cpm_sigdown_melted_means_subset1,condition=="ama1_CA1352.d2")



ks.test(cpm_sigdown_melted_means_subset3$x,cpm_sigdown_melted_means_subset2$x,alternative = "two.sided")#ns control soma
ks.test(cpm_sigdown_melted_means_subset4$x,cpm_sigdown_melted_means_subset5$x,alternative = "two.sided") #p<2.2e-16, soma d6
ks.test(cpm_sigdown_melted_means_subset6$x,cpm_sigdown_melted_means_subset7$x,alternative = "two.sided") #ns, germline d6
ks.test(cpm_sigdown_melted_means_subset8$x,cpm_sigdown_melted_means_subset9$x,alternative = "two.sided") #ns, control gline
ks.test(cpm_sigdown_melted_means_subset10$x,cpm_sigdown_melted_means_subset2$x,alternative = "two.sided") #p<1.3e-13, CA1200 d2 d6
ks.test(cpm_sigdown_melted_means_subset11$x,cpm_sigdown_melted_means_subset4$x,alternative = "two.sided") #p<2.2e-16, soma d2 d6
ks.test(cpm_sigdown_melted_means_subset12$x,cpm_sigdown_melted_means_subset8$x,alternative = "two.sided") #p<2.2e-16, CA1352 d2 d6
ks.test(cpm_sigdown_melted_means_subset13$x,cpm_sigdown_melted_means_subset6$x,alternative = "two.sided") #p<2.2e-16 gline d2 d6

#venn of high confidence overlap and transcription-dependent genes. Sig overlap

#d2 vs d6 etoh comparisons within strain background
lrt5 <- glmLRT(fit, contrast = c(0,-1,1,0,0,0,0,0,0,0,0,0,0))
lrt_sort5<-topTags(lrt5,n=nrow(lrt5$table))$table
significant5<-subset(lrt_sort5,FDR<0.05) #5080 genes
sig_down5<-subset(significant5,logFC<0) #2560
sig_up5<-subset(significant5,logFC>0) #2520

lrt6 <- glmLRT(fit, contrast = c(0,0,0,0,-1,1,0,0,0,0,0,0,0))
lrt_sort6<-topTags(lrt6,n=nrow(lrt6$table))$table
significant6<-subset(lrt_sort6,FDR<0.05) #5839 genes
sig_down6<-subset(significant6,logFC<0) #2896
sig_up6<-subset(significant6,logFC>0) #2943

lrt7 <- glmLRT(fit, contrast = c(0,0,0,0,0,0,0,-1,1,0,0,0,0))
lrt_sort7<-topTags(lrt7,n=nrow(lrt7$table))$table
significant7<-subset(lrt_sort7,FDR<0.05) #7457 genes
sig_down7<-subset(significant7,logFC<0) #4025
sig_up7<-subset(significant7,logFC>0) #3432

lrt8 <- glmLRT(fit, contrast = c(0,0,0,0,0,0,0,0,0,0,-1,1,0))
lrt_sort8<-topTags(lrt8,n=nrow(lrt8$table))$table
significant8<-subset(lrt_sort8,FDR<0.05) #6458 genes
sig_down8<-subset(significant8,logFC<0) #3290
sig_up8<-subset(significant8,logFC>0) #3168


#define high confidence "up genes" by looking at the overlap of the four gene sets going up d2 to d6
#(did this before for the typical analysis)
venn_UP<-venn(list(ama1_soma=rownames(sig_up5),ama1_germline=rownames(sig_up6),somaDriver=rownames(sig_up7),germlineDriver=rownames(sig_up8))) #1684 genes
plot(venn_UP)
print(venn_UP)
HighConfUp<-attr(venn_UP,"intersections")$`ama1_soma:ama1_germline:somaDriver:germlineDriver`

ama1Dep_Up<-venn(list(HighConfUp=HighConfUp,ama1_Down=rownames(sig_down)))#593 genes - enrichment

venn_numbers<-c(HighConfUp=1091,ama1Down=2086,"HighConfUp&ama1Down"=593)
set.seed(4)
fit1<-euler(venn_numbers,fills="blue")
Overlap_593<-plot(fit1,quantities=TRUE)
Overlap_593


ama1Dep_UpGenes<-attr(ama1Dep_Up,"intersections")$`HighConfUp:ama1_Down`
ama1Dep_UpGenes<-as.data.frame(ama1Dep_UpGenes)
rownames(ama1Dep_UpGenes)<-ama1Dep_UpGenes$ama1Dep_UpGenes
head(ama1Dep_UpGenes)
head(cpm_d2)
cpm_ama1_593genes<-cpm_d2[rownames(ama1Dep_UpGenes),]
counts_groups
head(cpm_ama1_593genes)
cpm_ama1_593genes<-as.data.frame(cpm_ama1_593genes)
cpm_ama1_593genes$WBID<-rownames(cpm_ama1_593genes)
cpm_ama1_593genes_melted<-melt(cpm_ama1_593genes)
head(cpm_ama1_593genes_melted,50)
cpm_ama1_593genes_melted<-cpm_ama1_593genes_melted[order(cpm_ama1_593genes_melted$WBID),]
cpm_ama1_593genes_melted$condition<-counts_groups
head(cpm_ama1_593genes_melted)
cpm_ama1_593genes_melted_means<-aggregate(cpm_ama1_593genes_melted$value,list(condition=cpm_ama1_593genes_melted$condition,gene=cpm_ama1_593genes_melted$WBID),mean)
head(cpm_ama1_593genes_melted_means,20)
cpm_ama1_593genes_melted_means_subsetSoma<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1200.d2"|condition=="ama1_CA1200.d6_aux"|condition=="ama1_CA1200.d6_EtOH"|condition=="CA1200.d2"|condition=="CA1200.d6_aux"|condition=="CA1200.d6_EtOH")
cpm_ama1_593genes_melted_means_subsetSoma$condition <- factor(cpm_ama1_593genes_melted_means_subsetSoma$condition,levels = c("ama1_CA1200.d2","ama1_CA1200.d6_EtOH","ama1_CA1200.d6_aux","CA1200.d2","CA1200.d6_EtOH","CA1200.d6_aux"))
cpm_ama1_593genes_melted_means$condition <- factor(cpm_ama1_593genes_melted_means$condition,levels = c("ama1_CA1200.d2","ama1_CA1200.d6_EtOH","ama1_CA1200.d6_aux","CA1200.d2","CA1200.d6_EtOH","CA1200.d6_aux","ama1_CA1352.d2","ama1_CA1352.d6_EtOH","ama1_CA1352.d6_aux","CA1352.d2","CA1352.d6_EtOH","CA1352.d6_aux"))



#ama1 INDEPENDENT up genes
ama1Ind_UpGenes<-attr(ama1Dep_Up,"intersections")$`HighConfUp`
ama1Ind_UpGenes<-as.data.frame(ama1Ind_UpGenes)
rownames(ama1Ind_UpGenes)<-ama1Ind_UpGenes$ama1Ind_UpGenes
head(ama1Ind_UpGenes)
#how many of these have small fold changes
subset_up<-sig_up5[rownames(ama1Ind_UpGenes),]#mean 1.62, median 1.19
subset_up<-sig_up6[rownames(ama1Ind_UpGenes),] #mean 1.75, median 1.303
subset_up<-sig_up7[rownames(ama1Ind_UpGenes),]#mean 2.18, median 1.76
subset_up<-sig_up8[rownames(ama1Ind_UpGenes),] #1.78, 1.36
colMeans(subset_up)
colMedians(as.matrix(subset_up))
subset_up_smallFC<-subset(subset_up,logFC<1 & logFC>-1 ) #399,333,149,316of 1091 have |FC|<1



cpm_ama1_1091genes<-cpm_d2[rownames(ama1Ind_UpGenes),]
head(cpm_ama1_1091genes)
cpm_ama1_1091genes<-as.data.frame(cpm_ama1_1091genes)
cpm_ama1_1091genes$WBID<-rownames(cpm_ama1_1091genes)
cpm_ama1_1091genes_melted<-melt(cpm_ama1_1091genes)
head(cpm_ama1_1091genes_melted)
cpm_ama1_1091genes_melted<-cpm_ama1_1091genes_melted[order(cpm_ama1_1091genes_melted$WBID),]
cpm_ama1_1091genes_melted$condition<-counts_groups
cpm_ama1_1091genes_melted_means<-aggregate(cpm_ama1_1091genes_melted$value,list(condition=cpm_ama1_1091genes_melted$condition,gene=cpm_ama1_1091genes_melted$WBID),mean)
head(cpm_ama1_1091genes_melted_means)
cpm_ama1_1091genes_melted_means_subsetSoma<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1200.d2"|condition=="ama1_CA1200.d6_aux"|condition=="ama1_CA1200.d6_EtOH"|condition=="CA1200.d2"|condition=="CA1200.d6_aux"|condition=="CA1200.d6_EtOH")
cpm_ama1_1091genes_melted_means_subsetSoma$condition <- factor(cpm_ama1_1091genes_melted_means_subsetSoma$condition,levels = c("ama1_CA1200.d2","ama1_CA1200.d6_EtOH","ama1_CA1200.d6_aux","CA1200.d2","CA1200.d6_EtOH","CA1200.d6_aux"))
cpm_ama1_1091genes_melted_means$condition <- factor(cpm_ama1_1091genes_melted_means$condition,levels = c("ama1_CA1200.d2","ama1_CA1200.d6_EtOH","ama1_CA1200.d6_aux","CA1200.d2","CA1200.d6_EtOH","CA1200.d6_aux","ama1_CA1352.d2","ama1_CA1352.d6_EtOH","ama1_CA1352.d6_aux","CA1352.d2","CA1352.d6_EtOH","CA1352.d6_aux"))



cpm_ama1_593genes_melted_means$genegroup<-"TranscriptionDepUp"
cpm_ama1_1091genes_melted_means$genegroup<-"TranscriptionIndUp"
cpm_1684genes_merge<-rbind(cpm_ama1_593genes_melted_means,cpm_ama1_1091genes_melted_means)

#Figure 5F
Up_1684<-ggplot(cpm_1684genes_merge,aes(condition,log2(x+0.1),color=factor(condition)))+
  geom_boxplot()+theme_classic(base_size = 10)+labs(x="Condition",y="log2(CPM)")+
  ggtitle("High confidence up genes")+
  facet_grid(genegroup ~ .)+
  scale_colour_manual(values=c("#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999","#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(aspect.ratio = 0.5,legend.position = "none")+labs(tag = "G")
Up_1684
#plot of the ~500 genes in the highly significant overlap across 12 conditions

head(cpm_ama1_593genes_melted_means)
cpm_ama1_593genes_melted_means_subset1<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1200.d2")
cpm_ama1_593genes_melted_means_subset2<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1200.d6_EtOH")
cpm_ama1_593genes_melted_means_subset3<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1200.d6_aux")
cpm_ama1_593genes_melted_means_subset4<-subset(cpm_ama1_593genes_melted_means,condition=="CA1200.d2")
cpm_ama1_593genes_melted_means_subset5<-subset(cpm_ama1_593genes_melted_means,condition=="CA1200.d6_EtOH")
cpm_ama1_593genes_melted_means_subset6<-subset(cpm_ama1_593genes_melted_means,condition=="CA1200.d6_aux")
cpm_ama1_593genes_melted_means_subset7<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1352.d2")
cpm_ama1_593genes_melted_means_subset8<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1352.d6_EtOH")
cpm_ama1_593genes_melted_means_subset9<-subset(cpm_ama1_593genes_melted_means,condition=="ama1_CA1352.d6_aux")
cpm_ama1_593genes_melted_means_subset10<-subset(cpm_ama1_593genes_melted_means,condition=="CA1352.d2")
cpm_ama1_593genes_melted_means_subset11<-subset(cpm_ama1_593genes_melted_means,condition=="CA1352.d6_EtOH")
cpm_ama1_593genes_melted_means_subset12<-subset(cpm_ama1_593genes_melted_means,condition=="CA1352.d6_aux")

ks.test(cpm_ama1_593genes_melted_means_subset1$x,cpm_ama1_593genes_melted_means_subset2$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset2$x,cpm_ama1_593genes_melted_means_subset3$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset4$x,cpm_ama1_593genes_melted_means_subset5$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset5$x,cpm_ama1_593genes_melted_means_subset6$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset7$x,cpm_ama1_593genes_melted_means_subset8$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset8$x,cpm_ama1_593genes_melted_means_subset9$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset10$x,cpm_ama1_593genes_melted_means_subset11$x,alternative = "two.sided")
ks.test(cpm_ama1_593genes_melted_means_subset11$x,cpm_ama1_593genes_melted_means_subset12$x,alternative = "two.sided")


cpm_ama1_1091genes_melted_means_subset1<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1200.d2")
cpm_ama1_1091genes_melted_means_subset2<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1200.d6_EtOH")
cpm_ama1_1091genes_melted_means_subset3<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1200.d6_aux")
cpm_ama1_1091genes_melted_means_subset4<-subset(cpm_ama1_1091genes_melted_means,condition=="CA1200.d2")
cpm_ama1_1091genes_melted_means_subset5<-subset(cpm_ama1_1091genes_melted_means,condition=="CA1200.d6_EtOH")
cpm_ama1_1091genes_melted_means_subset6<-subset(cpm_ama1_1091genes_melted_means,condition=="CA1200.d6_aux")
cpm_ama1_1091genes_melted_means_subset7<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1352.d2")
cpm_ama1_1091genes_melted_means_subset8<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1352.d6_EtOH")
cpm_ama1_1091genes_melted_means_subset9<-subset(cpm_ama1_1091genes_melted_means,condition=="ama1_CA1352.d6_aux")
cpm_ama1_1091genes_melted_means_subset10<-subset(cpm_ama1_1091genes_melted_means,condition=="CA1352.d2")
cpm_ama1_1091genes_melted_means_subset11<-subset(cpm_ama1_1091genes_melted_means,condition=="CA1352.d6_EtOH")
cpm_ama1_1091genes_melted_means_subset12<-subset(cpm_ama1_1091genes_melted_means,condition=="CA1352.d6_aux")

ks.test(cpm_ama1_1091genes_melted_means_subset1$x,cpm_ama1_1091genes_melted_means_subset2$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset2$x,cpm_ama1_1091genes_melted_means_subset3$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset4$x,cpm_ama1_1091genes_melted_means_subset5$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset5$x,cpm_ama1_1091genes_melted_means_subset6$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset7$x,cpm_ama1_1091genes_melted_means_subset8$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset8$x,cpm_ama1_1091genes_melted_means_subset9$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset10$x,cpm_ama1_1091genes_melted_means_subset11$x,alternative = "two.sided")
ks.test(cpm_ama1_1091genes_melted_means_subset11$x,cpm_ama1_1091genes_melted_means_subset12$x,alternative = "two.sided")



head(cpm_1684genes_merge)
cpm_1684genes_merge_d6aux<-subset(cpm_1684genes_merge,condition=="ama1_CA1200.d6_aux"|condition=="ama1_CA1352.d6_aux"|condition=="CA1200.d6_aux"|condition=="CA1352.d6_aux")
cpm_1684genes_merge_d6EtOH<-subset(cpm_1684genes_merge,condition=="ama1_CA1200.d6_EtOH"|condition=="ama1_CA1352.d6_EtOH"|condition=="CA1200.d6_EtOH"|condition=="CA1352.d6_EtOH")
cpm_1684genes_merge_d2<-subset(cpm_1684genes_merge,condition=="ama1_CA1200.d2"|condition=="ama1_CA1352.d2"|condition=="CA1200.d2"|condition=="CA1352.d2")
cpm_1684genes_merge_d6EtOH$norm_x<-cpm_1684genes_merge_d6EtOH$x/TotalRNANormFactor
cpm_1684genes_merge_d6EtOH$norm_x<-cpm_1684genes_merge_d6EtOH$x/TotalRNANormFactor
cpm_1684genes_merge_d2$norm_x<-cpm_1684genes_merge_d2$x
cpm_1684genes_merge_normx<-rbind(cpm_1684genes_merge_d2,cpm_1684genes_merge_d6EtOH)

head(RNAyield_means)
RNAyield_means<-aggregate(RNAyield1$Mass.L1.pg.,list(RNAyield1$Strain,RNAyield1$Day),mean)
RNAyield_means$norm<-c(1,1,1,1,RNAyield_means[5,3]/RNAyield_means[1,3],RNAyield_means[6,3]/RNAyield_means[2,3],RNAyield_means[7,3]/RNAyield_means[3,3],RNAyield_means[8,3]/RNAyield_means[4,3])
RNAyield_means
RNAyield_means$norminverse<-1/RNAyield_means$norm
RNAyield_means<-rbind(RNAyield_means,RNAyield_means[5:8,])
RNAyield_means$conditions<-c("ama1_CA1200.d2","ama1_CA1352.d2","CA1200.d2","CA1352.d2","ama1_CA1200.d6_aux","ama1_CA1352.d6_aux","CA1200.d6_aux","CA1352.d6_aux","ama1_CA1200.d6_EtOH","ama1_CA1352.d6_EtOH","CA1200.d6_EtOH","CA1352.d6_EtOH")

RNAyield_means
mergeRNAyield<-merge(cpm_1684genes_merge,RNAyield_means,by.x="condition",by.y="conditions")
head(mergeRNAyield)
mergeRNAyield$cpmNorm<-mergeRNAyield$x.x/mergeRNAyield$norminverse
head(mergeRNAyield)
Up_1684_normbyTotalRNA<-ggplot(mergeRNAyield,aes(condition,log2(cpmNorm+0.1),color=factor(condition)))+
  geom_boxplot()+theme_classic(base_size = 10)+labs(x="Condition",y="log2(CPM)")+
  ggtitle("High confidence up genes")+
  facet_grid(genegroup~.)+
  scale_colour_manual(values=c("#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999","#009E73", "#E69F00","#56B4E9","#009E73","#CC79A7","#999999"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(aspect.ratio = 0.5,legend.position = "none")+labs(tag = "F")+
  scale_x_discrete(labels=c("36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux","36 hr","132 hr + EtOH","132 hr + aux"))

Up_1684_normbyTotalRNA

head(mergeRNAyield)
mergeRNAyield_dcast <- dcast(mergeRNAyield, gene + genegroup ~ condition, value.var = "cpmNorm")
head(mergeRNAyield_dcast)
mergeRNAyield_dcast$ama1_CA1200_FC<-log2(mergeRNAyield_dcast$ama1_CA1200.d6_EtOH/mergeRNAyield_dcast$ama1_CA1200.d2)
mergeRNAyield_dcast$CA1200_FC<-log2(mergeRNAyield_dcast$CA1200.d6_EtOH/mergeRNAyield_dcast$CA1200.d2)
mergeRNAyield_dcast$ama1_CA1352_FC<-log2(mergeRNAyield_dcast$ama1_CA1352.d6_EtOH/mergeRNAyield_dcast$ama1_CA1352.d2)
mergeRNAyield_dcast$CA1352_FC<-log2(mergeRNAyield_dcast$CA1352.d6_EtOH/mergeRNAyield_dcast$CA1352.d2)
mergeRNAyield_dcast_2<-mergeRNAyield_dcast[,c(1,2,15,16,17,18)]
head(mergeRNAyield_dcast_2)
mergeRNAyield_dcast_2<-melt(mergeRNAyield_dcast_2)
head(mergeRNAyield_dcast_2)


mergeRNAyield_dcast3 <- dcast(mergeRNAyield, gene + genegroup ~ condition, value.var = "x.x")
head(mergeRNAyield_dcast3)
mergeRNAyield_dcast3$ama1_CA1200_FC<-log2(mergeRNAyield_dcast3$ama1_CA1200.d6_EtOH/mergeRNAyield_dcast3$ama1_CA1200.d2)
mergeRNAyield_dcast3$CA1200_FC<-log2(mergeRNAyield_dcast3$CA1200.d6_EtOH/mergeRNAyield_dcast3$CA1200.d2)
mergeRNAyield_dcast3$ama1_CA1352_FC<-log2(mergeRNAyield_dcast3$ama1_CA1352.d6_EtOH/mergeRNAyield_dcast3$ama1_CA1352.d2)
mergeRNAyield_dcast3$CA1352_FC<-log2(mergeRNAyield_dcast3$CA1352.d6_EtOH/mergeRNAyield_dcast3$CA1352.d2)
mergeRNAyield_dcast_4<-mergeRNAyield_dcast3[,c(1,2,15,16,17,18)]
head(mergeRNAyield_dcast_4)
mergeRNAyield_dcast_4<-melt(mergeRNAyield_dcast_4)
head(mergeRNAyield_dcast_4)
head(mergeRNAyield_dcast_2)
mergeRNAyield_dcast_2$normalized<-"norm by total RNA"
mergeRNAyield_dcast_4$normalized<-"edgeR output"
mergeRNAyield_dcast_all<-rbind(mergeRNAyield_dcast_2,mergeRNAyield_dcast_4)
head(mergeRNAyield_dcast_all)
mergeRNAyield_dcast_avg<-aggregate(mergeRNAyield_dcast_all$value,list(genegroup=mergeRNAyield_dcast_all$genegroup,gene=mergeRNAyield_dcast_all$gene,normalized=mergeRNAyield_dcast_all$normalized),mean)
tail(mergeRNAyield_dcast_avg)

highFC<-subset(mergeRNAyield_dcast_avg,genegroup=="TranscriptionIndUp" & x>1 & normalized=="norm by total RNA")
highFC2<-subset(mergeRNAyield_dcast_avg,genegroup=="TranscriptionDepUp" & x>1 & normalized=="norm by total RNA")

mergeRNAyield_avgs
mergeRNAyield_avgs<-aggregate(mergeRNAyield_dcast_avg$x,list(genegroup=mergeRNAyield_dcast_avg$genegroup,normalized=mergeRNAyield_dcast_avg$normalized),median)
(mergeRNAyield_avgs$x[1]-mergeRNAyield_avgs$x[3])/mergeRNAyield_avgs$x[1] #reduction of 53%

(mergeRNAyield_avgs$x[2]-mergeRNAyield_avgs$x[4])/mergeRNAyield_avgs$x[2] #reduction of 67%

NormTotalRNA_plot<-ggplot(mergeRNAyield_dcast_avg,aes(x=genegroup,y=x,color=normalized))+
  geom_boxplot(alpha=0.1)+theme_classic(base_size = 14)+
  theme(aspect.ratio = 0.75,legend.position = "none")+
  labs(x="Gene group",y="log2(132 hr EtOH / 36 hr)",tag = "H")

NormTotalRNA_plot                                                        


cpm_1684_TransIndUp<-subset(mergeRNAyield,genegroup=="TranscriptionIndUp")
pairwise.t.test(cpm_1684_TransIndUp$cpmNorm,cpm_1684_TransIndUp$condition,pool.sd = FALSE,p.adjust="none") #normalizing by total RNA abrogates the effect. No significant difference

cpm_1684_TransDepUp<-subset(mergeRNAyield,genegroup=="TranscriptionDepUp")
pairwise.t.test(cpm_1684_TransDepUp$cpmNorm,cpm_1684_TransDepUp$condition,pool.sd = FALSE,p.adjust="none")


phyper(593-1,2679,15650-2679,1684,lower.tail = F,log.p = F) #8.57e-82
#This is the calculation for Figure 5E
phyper(592,1684,15650-1684,2675,lower.tail = FALSE,log.p = FALSE) #4.5e-82 after removing 4 ERCC from the 2679
#two ways of calculating the hypergeometric p-value.
#first is the number of genes in both sets minus 1
#next is the total number of DE genes in one list
#subtract that number (above) from the total number of genes
#the last number is the number of DE genes in the other list
#lower tail is false means the p-value is the prob of getting that many genes genes or greater overlap




grid.arrange(Prop_ERCC,PCA,TranscriptionDependent_2679,Up_1684,NormTotalRNA_plot,TotalRNA,ProteinYield,nrow=2)


