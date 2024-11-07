#----Libraries ----
(setwd("C:/checkmate238"))
library(tidyverse)
library(dplyr)
library(viridis)
library(viridisLite)
library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(reshape2)
library(scDblFinder) #new doublet rate since using the HT kit...assumed to be 0.4% per thousand cells 
library(pheatmap)
library(SingleCellExperiment)
library(tidyr)
library(scRepertoire)
library(readxl)
library(DESeq2)
library(SummarizedExperiment)
library(GSEABase)
library(GSVA)
library(ggbiplot)
library(glmGamPoi)
library(scuttle)
library(limma)
library(SingleR)
library(celldex)
library(scRNAseq)
library(rstatix)
library(magrittr)
library(cetcolor)
library(gridExtra)
library(scuttle)
library(ComplexHeatmap)
library(ggplot2)

sessionInfo()

#----Subset to CD8 T cells ----
mergedObject.integrated.sct <- readRDS(file ="checkmateSNP_line141_azimuth_mergedObject.integrated.sct.rds")

Idents(mergedObject.integrated.sct) <- 'predicted.celltype.l1'
cd8t <- subset(x = mergedObject.integrated.sct, idents = c("CD8 T"))

#saveRDS(cd8t, "cd8t.rds")
#cd8t <- readRDS(file="cd8t.rds")

#add batch data
meta.data <- cd8t@meta.data
cd8t@meta.data <- meta.data %>%
  mutate(runDate = case_when(
    endsWith(subject2, "C")  ~ "2023_07",endsWith(subject2, "D")  ~ "2023_07", endsWith(subject2, "E")  ~ "2023_07",
    endsWith(subject2, "G")  ~ "2023_07", endsWith(subject2, "B")  ~ "2022_07",endsWith(subject2, "H")  ~ "2023_07",
    endsWith(subject2, "A")  ~ "2023_07",endsWith(subject2, "F")  ~ "2022_07"))
rm(meta.data)

#NOTE: Seurat cluster 13 and subsequent clusters have been renamed in Illustrator due to the absence of cluster 12 in CD8 T cells
#table(cd8t$seurat_clusters)

#NOTE: For Figure 5 D and Extended Figure 5 A, subjects were renamed to letters A-H in Illustrator

#----Figure 5 C ----
cd8t.w3 <- subset(x=cd8t, subset = visit=="W3D1")
Idents(cd8t.w3) <- "seurat_clusters"
cd8t$trt_visit<-paste(cd8t$description_of_actual_arm, cd8t$visit, sep = "_")

p1 <- DimPlot(cd8t, pt.size=.3,label=F, label.color="blue",label.size=4, repel=T,group.by = "seurat_clusters", raster=F,ncol=1) + ylim(-10,20) + xlim(-15,15)

p1 & ylim(-9,16.5) & xlim(-12.5,11) & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "after_stat(level)"), 
                  linewidth = 0.5, geom = "density_2d_filled",
                  colour = "black", alpha = 0) & scale_fill_gradient(low="white", high="white", labels=NULL, guide=NULL) & theme(panel.border = element_rect(colour = "black", fill="NA",size=1))

ggsave("c238_umap_ipi&nivo_fig.pdf", width=12, height=8, units=c("in"))

##----Extended Figure 5 B ----
p2 <- DimPlot(cd8t, pt.size=.3,label=F, label.color="blue",label.size=4, repel=T, split.by = "trt_visit",group.by = "seurat_clusters", raster=F,ncol=2) + ylim(-10,20) + xlim(-15,15)

p2 & ylim(-9,16.5) & xlim(-12.5,11) & 
  stat_density_2d(aes_string(x = "UMAP_1", y = "UMAP_2", fill = "after_stat(level)"), 
                  linewidth = 0.5, geom = "density_2d_filled",
                  colour = "black", alpha = 0) & scale_fill_gradient(low="white", high="white", labels=NULL, guide=NULL) & theme(panel.border = element_rect(colour = "black", fill="NA",size=1))

rm(p1, p2)

#----Figure 5 D ----
all <- as.data.frame(table(cd8t$subject3,cd8t$seurat_clusters))
all <- subset(all, Var2 != "12") #cluster 12 removed because it contains no cells in the cd8t population
trt <- as.data.frame(table(cd8t$description_of_actual_arm,cd8t$subject3))
trt <- subset(trt, Freq > 0)
comb <- merge(trt,all, by.x="Var2", by.y="Var1")
comb[c('subject1','subject2','visit')] <- str_split_fixed(comb$Var2, '_',3)
comb$subject <- paste(comb$subject1, comb$subject2, sep="_")
colnames(comb) <- c("Var2", "treatment","total","cluster","cluster_total","subject1","subject2","visit","subject")
comb[c('Week2','day')] <- str_split_fixed(comb$visit, 'D',2)
comb[c('Week3','Week')] <- str_split_fixed(comb$Week2, 'W', 2)
comb$Participant <- comb$subject
comb2 <- subset(comb, select = -c(Var2,subject1,subject2,day,Week3))
comb2$frequency <- (((comb2$cluster_total)/(comb2$total))*100)

f.i <- subset(comb2, subset= treatment=="IPILIMUMAB")
f.n <- subset(comb2, subset= treatment=="NIVOLUMAB")

f.i$cluster_Participant<-paste(f.i$cluster, f.i$Participant, sep = "_")
f.n$cluster_Participant<-paste(f.n$cluster, f.n$Participant, sep = "_")

ci <- ggplot(f.i) &
  geom_line(aes(x=Week,y=frequency,group=cluster_Participant,color=Participant), linewidth=1) & geom_point(aes(x=Week,y=frequency,color=Participant), size=2)&
  facet_grid(.~cluster,scales="free") & theme_bw() & theme(panel.border = element_rect(colour = "darkgrey", fill="NA",size=1)) & ylab("% CD8 T cells") & ggtitle("Ipilimumab")
cn <- ggplot(f.n) &
  geom_line(aes(x=Week,y=frequency,group=cluster_Participant,color=Participant), linewidth=1) & geom_point(aes(x=Week,y=frequency,color=Participant), size=2)&
  facet_grid(.~cluster,scales="free") & theme_bw() & theme(panel.border = element_rect(colour = "darkgrey", fill="NA",size=1)) & ylab("% CD8 T cells") & ggtitle("Nivolumab")

library(gridExtra)
cin <- grid.arrange(ci, cn, ncol=1); cin

ggsave("cluster-frequency-plot.pdf", height=8, width=12, units=c("in"))

rm(all, ci, cin, cn, comb, f.i, f.n, trt)

##----Figure 5 D: statistics----
comb2 <- subset(comb2, select = -c(Week2,Week))

clusteripi <- subset(comb2, treatment=="IPILIMUMAB")
clusternivo <- subset(comb2, treatment=="NIVOLUMAB")

comb2$cluster_visit_treatment<-paste(comb2$cluster,comb2$visit,comb2$treatment,sep = "_")

#anova; cluster:visit is significant in ipi only 
anovaipi <- anova_test(data=clusteripi, dv= frequency, wid= subject, within=c(cluster,visit))
anovanivo <- anova_test(data=clusternivo, dv= frequency, wid= subject, within=c(cluster,visit))

#post-hoc Tukey's HSD
tukey <- tukey_hsd(comb2, frequency ~ cluster_visit_treatment)
#write.csv(tukey,"tukeyipinivo.csv")

rm(anovaipi, anovanivo, clusteripi, clusternivo, comb2, tukey)

#----Extended Figure 5 A ----
#only timepoint W3D1
cd8t.sub.3 <- subset(x=cd8t, subset = visit=="W3D1")

#pseudobulk aggregate

library(scuttle)
DefaultAssay(cd8t.sub.3) <- "RNA"
cd8t.sub.3.sce <- cd8t.sub.3 %>% as.SingleCellExperiment()

#determine the number of cells per sample
table(cd8t.sub.3.sce$subject3)
groups <- colData(cd8t.sub.3.sce)[, c("subject3")]
cd8t.sub.3.sce
cd8t.sub.3.sce <- removeAltExps(cd8t.sub.3.sce) 
#aggregate across cluster-sample groups
pseudo_bulk_CD8_3 <- scuttle::aggregateAcrossCells(cd8t.sub.3.sce, ids = colData(cd8t.sub.3.sce)[, c("subject3")])

class(pseudo_bulk_CD8_3)
dim(pseudo_bulk_CD8_3)

#create metaData file 
metaData <- colnames(pseudo_bulk_CD8_3) %>% as.data.frame()
colnames(metaData) <- "subject3"

metaData <- metaData %>%
  mutate(trt = case_when(
    endsWith(subject3, "C_W3D1")  ~ "ipi",endsWith(subject3, "D_W3D1")  ~ "ipi", endsWith(subject3, "E_W3D1")  ~ "nivo",
    endsWith(subject3, "G_W3D1")  ~ "nivo", endsWith(subject3, "B_W3D1")  ~ "ipi",endsWith(subject3, "H_W3D1")  ~ "nivo",
    endsWith(subject3, "A_W3D1")  ~ "ipi",endsWith(subject3, "F_W3D1")  ~ "nivo"))

metaData <- metaData %>%
  mutate(sex = case_when(
    endsWith(subject3, "C_W3D1")  ~ "male",endsWith(subject3, "D_W3D1")  ~ "female", endsWith(subject3, "E_W3D1")  ~ "female",
    endsWith(subject3, "G_W3D1")  ~ "female", endsWith(subject3, "B_W3D1")  ~ "female",endsWith(subject3, "H_W3D1")  ~ "male",
    endsWith(subject3, "A_W3D1")  ~ "male",endsWith(subject3, "F_W3D1")  ~ "female"))

metaData <- metaData %>%
  mutate(age2 = case_when(
    endsWith(subject3, "C_W3D1")  ~ "48",endsWith(subject3, "D_W3D1")  ~ "45", endsWith(subject3, "E_W3D1")  ~ "63",
    endsWith(subject3, "G_W3D1")  ~ "60", endsWith(subject3, "B_W3D1")  ~ "68",endsWith(subject3, "H_W3D1")  ~ "56",
    endsWith(subject3, "A_W3D1")  ~ "72",endsWith(subject3, "F_W3D1")  ~ "41"))

metaData <- metaData %>%
  mutate(stage = case_when(
    endsWith(subject3, "C_W3D1")  ~ "3B",endsWith(subject3, "D_W3D1")  ~ "3C", endsWith(subject3, "E_W3D1")  ~ "3C",
    endsWith(subject3, "G_W3D1")  ~ "3C", endsWith(subject3, "B_W3D1")  ~ "3C",endsWith(subject3, "H_W3D1")  ~ "3C",
    endsWith(subject3, "A_W3D1")  ~ "3C",endsWith(subject3, "F_W3D1")  ~ "3C"))

metaData <- metaData %>%
  mutate(runDate = case_when(
    endsWith(subject3, "C_W3D1")  ~ "2023_07",endsWith(subject3, "D_W3D1")  ~ "2023_07", endsWith(subject3, "E_W3D1")  ~ "2023_07",
    endsWith(subject3, "G_W3D1")  ~ "2023_07", endsWith(subject3, "B_W3D1")  ~ "2022_07",endsWith(subject3, "H_W3D1")  ~ "2023_07",
    endsWith(subject3, "A_W3D1")  ~ "2023_07",endsWith(subject3, "F_W3D1")  ~ "2022_07"))

metaData <- metaData %>%
  unite('subgroup', trt, sex, runDate, remove=FALSE)
metaData

#create DESeq2 object        
pseudo_bulk_CD8_3 <- assay(pseudo_bulk_CD8_3) 
cluster_counts <- as.data.frame(as.matrix(pseudo_bulk_CD8_3[, which(colnames(pseudo_bulk_CD8_3) %in% metaData$subject3)])) 
#View(cluster_counts) #just making sure metaData and counts have same info
#check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
all(metaData$subject3 == colnames(cluster_counts))         

dds <- DESeqDataSetFromMatrix(cluster_counts, 
                              colData = metaData, 
                              design = ~ subgroup) 

#variance stabilizing transformation 
vsd <- vst(dds, blind=FALSE)
#vsd$subgroup <- factor(vsd$subgroup, levels = c("group1","group2"))
z <- plotPCA(vsd, intgroup = "trt")
#z <- plotPCA(vsd)
z + geom_point( size = 6)+ theme_classic()  

z <- plotPCA(vsd, intgroup = "subgroup")
z + geom_point( size = 6)+ theme_classic()  

#batch variation removed using removeBatchEffect 
#removed any shifts in the log2-scale expression data that can be explained by batch (runDate) 
mat <- assay(vsd)
mm <- model.matrix(~ trt, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$runDate, design=mm)
assay(vsd) <- mat
plotPCA(vsd, intgroup = "trt")

p5 <- plotPCA(vsd, intgroup = "trt")
p5 + geom_point( size = 7)+ theme_classic()+
  geom_text(
    label=vsd$subject3,
    nudge_x = 0.5, nudge_y = 0.5, 
    check_overlap = F
  )

rm(cd8t.sub.3, cd8t.sub.3.sce, cluster_counts, dds, mat, mm, p5, pseudo_bulk_CD8_3, vsd, z, groups, metaData)

#----Extended Figure 5 C ----
cd8t.3 <- subset(cd8t, subset =  seurat_clusters=="3")
cd8t.3.w3 <- subset(cd8t.3, subset = visit=="W3D1")

DefaultAssay(cd8t.3.w3) <- "integrated"
VlnPlot(cd8t.3.w3, features= c('GZMB', 'GZMH', 'NKG7', 'PRF1', 'GNLY', 'IL21'), group.by='description_of_actual_arm', ncol=2, pt.size=0)

ggsave("vlnplot_w3_cluster3_ipi_v_nivo_fig.pdf", width=6, height=8, units="in")

#----Extended Figure 5 D ----
DefaultAssay(cd8t.3.w3) <- "SCT"
VlnPlot(cd8t.3.w3, features= c('GZMB', 'GZMH', 'NKG7', 'PRF1', 'GNLY'), group.by='seurat_clusters',ncol=5, pt.size=0)

#----Figure 5 E ----
cd8t.week3 <- subset(cd8t, subset = visit=="W3D1")

DefaultAssay(cd8t.week3) <- "RNA"
a <- VlnPlot(cd8t.week3, features= c('IL21R'), group.by='seurat_clusters',pt.size=.7,split.plot=F, log=F ) + ggtitle("IL21R") + stat_ydensity(drop=FALSE) + theme(legend.position="none"); a

b <- VlnPlot(cd8t.week3, features= c('GZMB'), group.by='seurat_clusters',pt.size=.7,split.plot=F, log=F ) + ggtitle("GZMB") + stat_ydensity(drop=FALSE) + theme(legend.position="none"); b

c <- VlnPlot(cd8t.week3, features= c('PRF1'), group.by='seurat_clusters',pt.size=.7,split.plot=F, log=F ) + ggtitle("PRF1") + stat_ydensity(drop=FALSE) + theme(legend.position="none"); c

abc <- grid.arrange(a, b, c, ncol=1); abc & stat_ydensity(drop=FALSE) & theme(legend.position="none")

ggsave("IL21R_GZMB_PRF1_vlnplot_fig.pdf", width=12, height=8, units=c("in"))

rm(a, b, c, abc, cd8t.3.w3)

#----Figure 5 F G H ----
#packages <- c("nbpMatching", "singscore")
#BiocManager::install(packages)

#devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
#devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
#devtools::install_github("jackbibby1/SCPA", force=TRUE)

library(SCPA)
library(msigdbr)
library(Seurat)
library(dplyr)
library(ggplot2)

#cd8t- week 3 only
cd8t.w3 <- subset(x=cd8t, subset = visit=="W3D1")

#import genesets
pathways <- "C:/checkmate238/IL21gmt.csv"

Idents(cd8t.w3) <- "description_of_actual_arm"
DimPlot(cd8t.w3) + theme(aspect.ratio=1)

ipi <- seurat_extract(cd8t.w3,
                      meta1 = "description_of_actual_arm", value_meta1 = "IPILIMUMAB")
nivo <- seurat_extract(cd8t.w3,
                       meta1 = "description_of_actual_arm", value_meta1 = "NIVOLUMAB")

scpa_out <- compare_pathways(samples = list(ipi, nivo), 
                             pathways = pathways,
                             downsample=5000)

##----Figure 5 G ----
p2 <- plot_rank(scpa_out = scpa_out, 
          pathway = "IL21", highlight_point_color = "#60c5f7",base_point_size=3,
          highlight_point_size = 4, label_pathway=T); p2 & xlim(-5,10) & ggtitle("Week 3 Pathway Enrichment");p2

ggsave("pathwayrank_fig2.pdf",width=12, height=8, units=c("in"))

##----Figure 5 H ----
cd8t.w1 <- subset(x=cd8t, subset = visit=="W1D1")

cd8t.w1 <- seurat_extract(cd8t, meta1 = "visit", value_meta1 = "W1D1")

scpa_out.w1 <- compare_pathways(samples = list(cd8t.w1,cd8t.w1), pathways = pathways, downsample=7000)

scpa_out.ipi <- compare_pathways(samples = list(cd8t.w1,ipi), pathways = pathways, downsample=7000)

scpa_out.nivo <- compare_pathways(samples = list(cd8t.w1, nivo), pathways = pathways, downsample=7000)

scpa_out.w1 <- scpa_out.w1 %>% mutate(w1_qval= qval)
scpa_out.ipi <- scpa_out.ipi %>% mutate(ipiw3_qval= qval)
scpa_out.nivo <- scpa_out.nivo %>% mutate(nivow3_qval= qval)
scpa_out.6 <- merge(x = scpa_out.w1, y = c(scpa_out.ipi,scpa_out.nivo), by="Pathway")
scpa_out.7 <- scpa_out.6[,c("Pathway","w1_qval","nivow3_qval","ipiw3_qval")]

plot_heatmap(scpa_out.7,row_fontsize=10, column_fontsize=10,
             column_names = c("Week 1 ", "Nivolumab Week 3", "Ipilimumab Week 3"),
             cluster_columns = F,
             show_row_names = T, scale_breaks=3) + ggtitle("all cd8t's")

ggsave(a, "plotheatmap-w1-w3n-w3i.pdf", width=8, height=12, units=c("in"))

#save heatmap
#pdf command
#pdf(file = "C:/checkmate238/plotheatmap-w1-w3n-w3i.pdf",   #directory to save the file in
#    width = 8, # The width of the plot in inches
#    height = 12) # The height of the plot in inches
#plot_heatmap(scpa_out.7,row_fontsize=10, column_fontsize=10,
#             column_names = c("Week 1 ", "Nivolumab Week 3", "Ipilimumab Week 3"),
#             cluster_columns = F,
#             show_row_names = T, scale_breaks=3) + ggtitle("all cd8t's")
#create the file
#dev.off()

rm(scpa_out.w1, scpa_out.nivo, scpa_out.ipi, scpa_out.6, scpa_out.7, nivo, ipi, cd8t.w1, pathways)

#----Figure 5 F ----
#color based on value
color2 <- ifelse(scpa_out$FC < 0, "skyblue", "tomato")

ggplot(scpa_out, aes(y=reorder(Pathway,FC), x= FC, fill=color2)) +
  geom_col() +
  guides(fill=guide_legend(title="Enriched For")) +
  scale_fill_hue(labels=c("nivo","ipi")) +
  xlab("Fold Change") +
  ylab("Pathway") +
  ggtitle("Week 3 Pathway Enrichment") + theme_classic()

ggsave("pathwayenrichment-classic_fig.pdf",width=12, height=8, units=c("in"))

rm(scpa_out)

#----Extended Figure 5 E ----
cd8t.week3 <- subset(cd8t, visit == "W3D1")
cd8t.week3.s <- ScaleData(cd8t.week3, assay = "RNA")

Idents(cd8t.week3.s) <- 'seurat_clusters'
pbmc.markers <- FindAllMarkers(cd8t.week3.s, test.use="t", logfc.threshold=1)

Idents(cd8t.week3) <- 'seurat_clusters'
pbmc.markers <- FindAllMarkers(cd8t.week3, test.use="t", logfc.threshold=1)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DefaultAssay(cd8t.week3) <- "RNA"
cd8t.week3 <- ScaleData(cd8t.week3, assay = "RNA")
DoHeatmap(cd8t.week3, features = top10$gene) + NoLegend()

all(top10$gene %in% rownames(cd8t.week3))

features <- intersect(top10$gene, rownames(top10$gene[['RNA']]))

cd8t.week3.s <- ScaleData(cd8t.week3, assay = "RNA")

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n=10) %>%
  ungroup() -> top10
DoHeatmap(cd8t.week3.s, features = top10$gene, slot='scale.data',size=2, angle=0, draw.lines=T) + theme(axis.text.y = element_text(size = 5)) +  scale_fill_gradientn(colors = c("blue", "white","red"))

rm(list= ls())