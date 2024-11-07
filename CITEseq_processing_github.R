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

sessionInfo()

#----Metadata ----
#set up the metadata; only needs to be run once- will be used for all sets in that batch

#compile meta.data files; alternative to using preexisting meta.data file (meta.csv)
#metaC <- read.csv(file = "metaC.csv")
#metaB <- read.csv(file = "metaB.csv")
#meta <- merge(x=metaC, y=metaB, by.x=c("usubjid"), by.y=c("Unique.Subject.Identifier"))

#meta.data compiled with extra info, titles adjusted
meta <- read.csv(file = "meta.csv")
sampleSet <- read.csv(file = "c238_sample_set.csv"); sampleSet$X <- NULL

meltedManifest <- melt(sampleSet, id.vars = "Biolegend.link", measure.vars = c("Samples.L12.L34", "Samples.ER1.W3"))
meltedManifest <- meltedManifest[which(meltedManifest$value != ""),]
rm(meltedManifest)

#PreIntegration 2022-08-22 ----
##----Merge ----
readData <- function(object.out, object)
{
  #new assay to store ADT information
  adt_object <- object.out$`Antibody Capture`
  adt_object@Dimnames[[1]]
  
  #remove HTOs from ADT panel 
  adt_assay_object <- adt_object[setdiff(rownames(x = adt_object),grep("Hashtag",adt_object@Dimnames[[1]], value=T )), ] #list HTOs so that they are removed from ADT assay 
  adt_assay_object <- CreateAssayObject(counts = adt_assay_object)
  
  #add assay to previously created Seurat object
  object[["ADT"]] <- adt_assay_object
  
  #add HTO as separate assay
  HTO_assay_object <- adt_object[intersect(rownames(x = adt_object),grep("Hashtag",adt_object@Dimnames[[1]], value=T )), ] #list HTOs 
  HTO_assay_object <- CreateAssayObject(counts = HTO_assay_object)
  
  #add assay to previously created Seurat object
  object[["HTO"]] <- HTO_assay_object
  return(object)
}

object.out <- Read10X(data.dir = "C:/checkmate238/2022-08-22/multi-L12-GEX/sample_filtered_feature_bc_matrix/") #folder containing filtered barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
object <- CreateSeuratObject(counts = object.out$`Gene Expression`, Project = "xxx")
object <- readData(object.out, object)

object@meta.data$orig.ident <- "L12"

object.out <- Read10X(data.dir = "C:/checkmate238/2022-08-22/multi-L34-GEX/sample_filtered_feature_bc_matrix/") #folder containing filtered barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
object.1 <- CreateSeuratObject(counts = object.out$`Gene Expression`, Project = "xxx")
object.1 <- readData(object.out, object.1)

object.1@meta.data$orig.ident <- "L34"

object <- merge(object, y = c(object.1), add.cell.ids = c("L12", "L34"), project = "xxx")
object

#check you see the full list of IDs
head(colnames(object))
tail(colnames(object))

#check ADT and HTO names
rownames(object[["ADT"]]);  rownames(object[["HTO"]])

#quantity of each identity
table(object$orig.ident)

setwd("C:/checkmate238")
saveRDS(object, file="object_filtered_multi-L12-L34.rds")
saveRDS(meta, file="meta_total_checkmate.rds")
#saveRDS(sampleSet, file="sampleSet_checkmate.rds")

rm(object.1, object.out, readData)

##----Demultiplex ----
object_filtered <- object

#if importing
#object_filtered <- readRDS(file="object_filtered_multi-L12-L34.rds")

#normalize HTO data, here we use centered log-ratio (CLR) transformation
object <- NormalizeData(object_filtered, assay = "HTO", normalization.method = "CLR")

#demultiplex
object <- HTODemux(object, assay = "HTO", positive.quantile = 0.99)

table(object$HTO_classification.global) 
table(object$HTO_maxID) 
Idents(object) <- "HTO_maxID"

#add in the correct values for HTOs based on sample
metaDat <- object@meta.data
metaDat$subject <- "blank"
metaDat$subject[grep("Hashtag-1", metaDat$HTO_maxID)] <- "F_W1D1"
metaDat$subject[grep("Hashtag-2", metaDat$HTO_maxID)] <- "F_W1D1"
metaDat$subject[grep("Hashtag-3", metaDat$HTO_maxID)] <- "F_W3D1"
metaDat$subject[grep("Hashtag-4", metaDat$HTO_maxID)] <- "F_W3D1"
metaDat$subject[grep("Hashtag-5", metaDat$HTO_maxID)] <- "B_W1D1"
metaDat$subject[grep("Hashtag-6", metaDat$HTO_maxID)] <- "B_W1D1"
metaDat$subject[grep("Hashtag-7", metaDat$HTO_maxID)] <- "B_W3D1"
metaDat$subject[grep("Hashtag-8", metaDat$HTO_maxID)] <- "B_W3D1"
object <- AddMetaData(object, metaDat)

VlnPlot(object, features = "nCount_HTO", pt.size = 0.1, log = TRUE)  + labs(title = "object") +  NoLegend() + 
  geom_hline(yintercept = 6)+ scale_y_continuous(limits = c(1,100000), trans = "log10")

##----Quality Control ----
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(object@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(object@assays[["RNA"]][mito.genes, ])/Matrix::colSums(object@assays[["RNA"]])

object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident" ) 

#filter out rare subset of cells with an outlier level of high mitochondrial percentage and low UMI content
FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.15)
FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") +  geom_hline(yintercept = c(110, 5400))

#do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(object) <- "HTO_classification.global"
RidgePlot(object, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

#definition for gating
object <- subset(x = object, subset = nFeature_RNA > 110 & nFeature_RNA < 5400 & percent.mito >  -Inf & percent.mito < 0.15) 

#visualization and clustering steps 
#normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins to calculate z-score for dispersion within each bin
object <- FindVariableFeatures(object, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
object <- ScaleData(object, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
object <- RunPCA(object, verbose = FALSE)

#examine and visualize PCA results a few different ways
DimPlot(object = object, reduction = "pca")
DimHeatmap(object = object, reduction = "pca", dims = 1:6, cells = 500, balanced = TRUE)
ElbowPlot(object)

object <- FindNeighbors(object, dims = 1:20) #adjust based on elbow plot
object <- FindClusters(object, resolution = 0.8, verbose = FALSE)
object <- RunUMAP(object, dims = 1:20) #adjust based on elbow plot

DimPlot(object, label = TRUE, group.by = "seurat_clusters")

#normalize ADT data
object <- NormalizeData(object, normalization.method = "CLR", margin = 2, assay = "ADT") 
object <- ScaleData(object, assay = "ADT")

##----Find doublets ----
DefaultAssay(object) <- "RNA"

#remove negatives from object
Idents(object) <- "HTO_classification.global"
object <- subset(object, idents = "Negative", invert = TRUE)
table(object$HTO_classification.global) 

saveRDS(object, "object-L12-L34.rds")

##----SNPs ----

#if importing
#object <- readRDS(file="object-L12-L34.rds")

orig.meta <- object@meta.data

library(readr)

#load in clusters.tsv files from souporcell
batch1a <-readr::read_tsv("batchL12_clusters.tsv")
batch1b <-readr::read_tsv("batchL34_clusters.tsv")

#for each tsv, add sample prefix, combine with barcode, and overwrite barcode
batch1a$ident <- "L12"
batch1a <- batch1a %>%
  unite('barcode', ident, barcode, remove=FALSE)

batch1b$ident <- "L34"
batch1b <- batch1b %>%
  unite('barcode', ident, barcode, remove=FALSE)

#combine tsv files, subset on singlets, and remove duplicate barcodes (though there shouldn't be)
batch1 <- rbind(batch1a,batch1b)
batch1 <- subset(batch1, subset =  status=="singlet")
batch1 <- batch1[!duplicated(batch1$barcode),]

batch1.2 <- data.frame(batch1$barcode, batch1$status, batch1$assignment, batch1$ident)
colnames(batch1.2) <- c("barcode","status","HTO_assignment","ident")

#add information to object metadata dataframe
orig.meta$barcode <- rownames(orig.meta)
orig.meta$toname <- rownames(orig.meta)
comb3 <- merge(x = orig.meta, y = batch1.2, by = "barcode")
rownames(comb3) <- comb3$toname

comb3 <- subset(comb3, select = -c(toname))

#use merged cells to subset object
b <- comb3$barcode
object<- subset(object, cells=b)

#reassign metadata
object@meta.data <- comb3

saveRDS(object, file="L12-L34_SNP.rds")

#scDblFinder
sce.object <- as.SingleCellExperiment(object)
sce.object <- scDblFinder(sce.object, dbr = (0.004 * ncol(sce.object)/1000), samples="HTO_maxID") #doublet rate of 0.4% per 1000 cells
table(sce.object$scDblFinder.class, sce.object$HTO_maxID)

identical(colnames(object),colnames(sce.object))

object$scDblFinder.class <- sce.object$scDblFinder.class
object$scDblFinder.score <- sce.object$scDblFinder.score
object$scDblFinder.weighted <- sce.object$scDblFinder.weighted
rm(sce.object)

object$scDblFinder.class <- factor(object$scDblFinder.class, levels=c("doublet", "singlet"))

#look at numbers of singlets and doublets using scDblFinder
object@meta.data %>%
  ggplot(aes(x= scDblFinder.class)) +
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = after_stat(count)), stat = "count", vjust = -1.5)

##----Filter out doublets----
table(object$scDblFinder.class)

#CD19 vs CD3 gene expression before filtering out doublets
FeatureScatter(object, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class") + NoLegend() 

#subset to singlets 
object_singlet <- subset(object, subset =  scDblFinder.class=="singlet")

#subset singlets round 2 by excluding cells that are double positive for CD19 and CD3E genes 
FeatureScatter(object_singlet, feature1 = "CD3E", feature2 = "CD19") + NoLegend() + labs(title="stim1") + geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1)

object_singlet.2 <- subset(object_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)

#look at number of singlets
object_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

saveRDS(object_singlet.2, file = "multi-L12-L34_postDblt_no_negative.rds")
saveRDS(meta, file="meta_total_checkmate.rds")
rm(object, object_filtered, object_singlet, object_singlet.2, metaDat, mito.genes, percent.mito)
rm(batch1, batch1.2, batch1a, batch1b, comb3, orig.meta, b)

##----TCR/BCR ----
setwd("C:/checkmate238/contig_annotations")
suppressMessages(library(scRepertoire))

###----TCR ----
object1_TCR <- read.csv("2022-08-22/TCR/multi_L12_TCR_filtered_contig_annotations.csv")
object2_TCR <- read.csv("2022-08-22/TCR/multi_L34_TCR_filtered_contig_annotations.csv")

contig_list_TCR <- list(object1_TCR, object2_TCR)

combined_TCR <- combineTCR(contig_list_TCR, 
                           samples =  c("object1","object2"),
                           ID = c("L12", "L34"))

#adjust connector and num_connects based on barcode format
for (i in seq_along(combined_TCR)) {
  combined_TCR[[i]] <- stripBarcode(combined_TCR[[i]], 
                                    column = 1, connector = "_", num_connects = 3)
}

#see barcode format
View(combined_TCR$object1_L12$barcode)

(setwd("C:/checkmate238"))
mergedObject.integrated.sct.TCR <- readRDS("multi-L12-L34_postDblt_no_negative.rds")
view(mergedObject.integrated.sct.TCR)

setwd("C:/checkmate238/contig_annotations")

#adding the barcode_prefix to be the same as in mergedObject.integrated.sct
combined_TCR$object1_L12$barcode_prefix <- "L12"
combined_TCR$object2_L34$barcode_prefix <- "L34"

View(combined_TCR)

for (i in seq_along(combined_TCR)) {
  combined_TCR[[i]] <- combined_TCR[[i]]  %>% unite('barcode',barcode_prefix, barcode, remove=T)
}
View(combined_TCR)

#adding more meta data to combined
data_to_add <- rownames_to_column(mergedObject.integrated.sct.TCR@meta.data, "barcodes")
data_to_add <- data_to_add[,c("barcodes","HTO_maxID")] %>% as.data.frame()
colnames(data_to_add) <- c("barcode","HTO_maxID")
for (i in seq_along(combined_TCR)) {
  combined_TCR[[i]] <- left_join(combined_TCR[[i]], data_to_add, by = "barcode")}

saveRDS(combined_TCR, file="combined_TCR_L12-L34.rds")

rm(contig_list_TCR, object1_TCR, object2_TCR, combined_TCR, data_to_add, i, mergedObject.integrated.sct.TCR)

###----BCR ----
object1_BCR <- read.csv("2022-08-22/BCR/multi_L12_BCR_filtered_contig_annotations.csv")
object2_BCR <- read.csv("2022-08-22/BCR/multi_L34_BCR_filtered_contig_annotations.csv")

contig_list_BCR <- list(object1_BCR, object2_BCR)

combined_BCR <- combineBCR(contig_list_BCR, 
                           samples =  c("object1","object2"),
                           ID = c("L12", "L34"))

#adjust connector and num_connects based on barcode format
for (i in seq_along(combined_BCR)) {
  combined_BCR[[i]] <- stripBarcode(combined_BCR[[i]],
                                    column = 1, connector = "_", num_connects = 3)
}

#see barcode format
View(combined_BCR$object1_L12$barcode)

(setwd("C:/checkmate238"))
mergedObject.integrated.sct.BCR <- readRDS("multi-L12-L34_postDblt_no_negative.rds")
view(mergedObject.integrated.sct.BCR)

setwd("C:/checkmate238/contig_annotations")

#adding the barcode_prefix to be the same as in mergedObject.integrated.sct
combined_BCR$object1_L12$barcode_prefix <- "L12"
combined_BCR$object2_L34$barcode_prefix <- "L34"

View(combined_BCR)

for (i in seq_along(combined_BCR)) {
  combined_BCR[[i]] <- combined_BCR[[i]]  %>% unite('barcode',barcode_prefix, barcode, remove=T)
}
View(combined_BCR)

#adding more meta data to combined
data_to_add <- rownames_to_column(mergedObject.integrated.sct.BCR@meta.data, "barcodes")
data_to_add <- data_to_add[,c("barcodes","HTO_maxID")] %>% as.data.frame()
colnames(data_to_add) <- c("barcode","HTO_maxID")
for (i in seq_along(combined_BCR)) {
  combined_BCR[[i]] <- left_join(combined_BCR[[i]], data_to_add, by = "barcode")}

saveRDS(combined_BCR, file="combined_BCR_L12-L34.rds")

rm(contig_list_BCR, object1_BCR, object2_BCR, combined_BCR, data_to_add, i, mergedObject.integrated.sct.BCR)

#PreIntegration 2023-07-27 ----
(setwd("C:/checkmate238"))
##----Merge ----
readData <- function(object.out, object)
{
  #new assay to store ADT information
  adt_object <- object.out$`Antibody Capture`
  adt_object@Dimnames[[1]]
  
  #remove HTOs from ADT panel 
  adt_object@Dimnames[[1]][grep("CD8_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0080_anti-human_CD8a_Antibody"
  adt_object@Dimnames[[1]][grep("CD19_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0050_anti-human_CD19_Antibody"  
  adt_object@Dimnames[[1]][grep("CD2_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0367_anti-human_CD2_Antibody"  
  adt_object@Dimnames[[1]][grep("CD5_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0138_anti-human_CD5_Antibody"  
  adt_object@Dimnames[[1]][grep("IsoIgG1_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0090_Mouse_IgG1___isotype_Ctrl_Antibody"  
  adt_object@Dimnames[[1]][grep("IsoIgG2A_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0091_Mouse_IgG2a___isotype_Ctrl_Antibody"  
  adt_object@Dimnames[[1]][grep("CD4_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0072_anti-human_CD4_Antibody"  
  adt_object@Dimnames[[1]][grep("CD27_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0154_anti-human_CD27_Antibody"  
  adt_object@Dimnames[[1]][grep("CD45RA_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0063_anti-human_CD45RA_Antibody"  
  adt_object@Dimnames[[1]][grep("CD152_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0151_anti-human_CD152__CTLA-4_Antibody"  
  adt_object@Dimnames[[1]][grep("CD279_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0088_anti-human_CD279__PD-1_Antibody"  
  adt_object@Dimnames[[1]][grep("CD38_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0410_anti-human_CD38_Antibody"  
  adt_object@Dimnames[[1]][grep("IgG1_ADT", adt_object@Dimnames[[1]])] <- "Mouse_Anti-Human_IgG1_Fc-UNLB__HP6001"   
  adt_object@Dimnames[[1]][grep("IgG4_ADT", adt_object@Dimnames[[1]])] <- "Mouse_Anti-Human_IgG4_pFc-UNLB__HP6023"  
  adt_object@Dimnames[[1]][grep("CXCR5_ADT", adt_object@Dimnames[[1]])] <- "TotalSeq-C0144_anti-human_CD185__CXCR5_Antibody"    
  
  adt_object@Dimnames[[1]][grep("E_W1D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0251_anti-human_Hashtag_1_Antibody"  
  adt_object@Dimnames[[1]][grep("E_W3D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0252_anti-human_Hashtag_2_Antibody"  
  adt_object@Dimnames[[1]][grep("D_W1D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0253_anti-human_Hashtag_3_Antibody"  
  adt_object@Dimnames[[1]][grep("D_W3D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0254_anti-human_Hashtag_4_Antibody"  
  adt_object@Dimnames[[1]][grep("G_W1D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0255_anti-human_Hashtag_5_Antibody"  
  adt_object@Dimnames[[1]][grep("G_W3D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0256_anti-human_Hashtag_6_Antibody"  
  adt_object@Dimnames[[1]][grep("H_W1D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0257_anti-human_Hashtag_7_Antibody"  
  adt_object@Dimnames[[1]][grep("H_W3D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0258_anti-human_Hashtag_8_Antibody"  
  adt_object@Dimnames[[1]][grep("A_W1D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0259_anti-human_Hashtag_9_Antibody"  
  adt_object@Dimnames[[1]][grep("A_W3D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0260_anti-human_Hashtag_10_Antibody"  
  adt_object@Dimnames[[1]][grep("C_W1D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0262_anti-human_Hashtag_12_Antibody" 
  adt_object@Dimnames[[1]][grep("C_W3D1", adt_object@Dimnames[[1]])] <- "TotalSeq-C0263_anti-human_Hashtag_13_Antibody" 
  
  adt_object@Dimnames[[1]]
  
  adt_assay_object <- adt_object[setdiff(rownames(x = adt_object),grep("Hashtag",adt_object@Dimnames[[1]], value=T )), ] #list HTOs so that they are removed from ADT assay 
  adt_assay_object <- CreateAssayObject(counts = adt_assay_object)
  
  #add assay to previously created Seurat object
  object[["ADT"]] <- adt_assay_object
  
  #add HTO as separate assay
  HTO_assay_object <- adt_object[intersect(rownames(x = adt_object),grep("Hashtag",adt_object@Dimnames[[1]], value=T )), ] #list HTOs 
  HTO_assay_object <- CreateAssayObject(counts = HTO_assay_object)
  
  #add assay to previously created Seurat object
  object[["HTO"]] <- HTO_assay_object
  return(object)
}
object.out <- Read10X(data.dir = "C:/checkmate238/2023-07-27/multi-ER1-GEX/sample_filtered_feature_bc_matrix/") #folder containing filtered barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
object <- CreateSeuratObject(counts = object.out$`Gene Expression`, Project = "xxx")
object <- readData(object.out, object)

object@meta.data$orig.ident <- "ER1"

object.out <- Read10X(data.dir = "C:/checkmate238/2023-07-27/multi-W1-GEX/sample_filtered_feature_bc_matrix/") #folder containing filtered barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
object.1 <- CreateSeuratObject(counts = object.out$`Gene Expression`, Project = "xxx")
object.1 <- readData(object.out, object.1)

object.1@meta.data$orig.ident <- "W1"

object.out <- Read10X(data.dir = "C:/checkmate238/2023-07-27/multi-W2-GEX/sample_filtered_feature_bc_matrix/") #folder containing filtered barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
object.2 <- CreateSeuratObject(counts = object.out$`Gene Expression`, Project = "xxx")
object.2 <- readData(object.out, object.2)

object.2@meta.data$orig.ident <- "W2"

object.out <- Read10X(data.dir = "C:/checkmate238/2023-07-27/multi-W3-GEX/sample_filtered_feature_bc_matrix/") #folder containing filtered barcodes.tsv.gz; features.tsv.gz; matrix.mtx.gz
object.3 <- CreateSeuratObject(counts = object.out$`Gene Expression`, Project = "xxx")
object.3 <- readData(object.out, object.3)

object.3@meta.data$orig.ident <- "W3"

object <- merge(object, y = c(object.1, object.2, object.3), add.cell.ids = c("ER1", "W1", "W2", "W3"), project = "xxx")
object

#check you see the full list of IDs
head(colnames(object))
tail(colnames(object))

#check ADT's and HTO's
rownames(object[["ADT"]]);  rownames(object[["HTO"]])

#see how many of each identity
table(object$orig.ident) 

saveRDS(object, file="object_filtered_ER1-W3.rds")
saveRDS(meta, file="meta_total_checkmate.rds")

rm(object.1, object.2, object.3, object.out)

##----Demultiplex ----
object_filtered <- object

#if importing
#object_filtered <- readRDS(file="object_filtered_ER1-W3.rds")

#normalize HTO data, here we use centered log-ratio (CLR) transformation
object <- NormalizeData(object_filtered, assay = "HTO", normalization.method = "CLR")

#demultiplex
object <- HTODemux(object, assay = "HTO", positive.quantile = 0.99)

table(object$HTO_classification.global) 
table(object$HTO_maxID) 
Idents(object) <- "HTO_maxID"

#add in the correct values for HTOs based on sample
metaDat <- object@meta.data
metaDat$subject <- "blank"
metaDat$subject[grep("Hashtag-1", metaDat$HTO_maxID)] <- "E_W1D1"
metaDat$subject[grep("Hashtag-2", metaDat$HTO_maxID)] <- "E_W3D1"
metaDat$subject[grep("Hashtag-3", metaDat$HTO_maxID)] <- "D_W1D1"
metaDat$subject[grep("Hashtag-4", metaDat$HTO_maxID)] <- "D_W3D1"
metaDat$subject[grep("Hashtag-5", metaDat$HTO_maxID)] <- "G_W1D1"
metaDat$subject[grep("Hashtag-6", metaDat$HTO_maxID)] <- "G_W3D1"
metaDat$subject[grep("Hashtag-7", metaDat$HTO_maxID)] <- "H_W1D1"
metaDat$subject[grep("Hashtag-8", metaDat$HTO_maxID)] <- "H_W3D1"
metaDat$subject[grep("Hashtag-9", metaDat$HTO_maxID)] <- "A_W1D1"
metaDat$subject[grep("Hashtag-10", metaDat$HTO_maxID)] <- "A_W3D1"
metaDat$subject[grep("Hashtag-12", metaDat$HTO_maxID)] <- "C_W1D1"
metaDat$subject[grep("Hashtag-13", metaDat$HTO_maxID)] <- "C_W3D1"
object <- AddMetaData(object, metaDat)

VlnPlot(object, features = "nCount_HTO", pt.size = 0.1, log = TRUE, raster=F)  + labs(title = "object") +  NoLegend() + 
  geom_hline(yintercept = 6)+ scale_y_continuous(limits = c(1,100000), trans = "log10")

##----Quality Control ----
#mtDNA% (fraction of mitochondrial transcript counts of total transcript counts) threshold is used to filter out dead, stressed, low-quality cells in data 
mito.genes <- grep(pattern = "^MT-", x = rownames(object@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(object@assays[["RNA"]][mito.genes, ])/Matrix::colSums(object@assays[["RNA"]])

object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito") 
VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident" ) 

#filter out rare subset of cells with an outlier level of high mitochondrial percentage and also low UMI content
FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "orig.ident" ) + geom_hline(yintercept = 0.17) 
FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") +  geom_hline(yintercept = c(250, 5500))

#do doublets have higher nCount and nFeature and negatives have lower nCount and nFeature?
Idents(object) <- "HTO_classification.global"
RidgePlot(object, features = c("nCount_RNA", "nFeature_RNA"), log = TRUE)

#definition for gating
object <- subset(x = object, subset = nFeature_RNA > 250 & nFeature_RNA < 5500 & percent.mito >  -Inf & percent.mito < 0.17 ) 

#perform visualization and clustering steps 
#normalize gene expression measurements for each cell by total expression and multiply this by scale factor 10000
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)

#calculate the average expression and dispersion for each gene and place these into bins to calculate z-score for dispersion within each bin
object <- FindVariableFeatures(object, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))

#removing unwanted sources of variation (technical noise, batch effects, biological sources of variation (cell cycle stages_))
object <- ScaleData(object, vars.to.regress = c("nCount_RNA", "percent.mito"))

#perform linear dimensional reduction
object <- RunPCA(object, verbose = FALSE)

#examine and visualize PCA results a few different ways
DimPlot(object = object, reduction = "pca",raster=F)
DimHeatmap(object = object, reduction = "pca", dims = 1:6, cells = 500, balanced = TRUE)
ElbowPlot(object)

object <- FindNeighbors(object, dims = 1:20) #adjust based on elbow plot
object <- FindClusters(object, resolution = 0.8, verbose = FALSE)
object <- RunUMAP(object, dims = 1:20) #adjust based on elbow plot

DimPlot(object, label = TRUE, group.by = "seurat_clusters",raster=F)

#normalize ADT data
object <- NormalizeData(object, normalization.method = "CLR", margin = 2, assay = "ADT") 
object <- ScaleData(object, assay = "ADT")

##----Find doublets----
DefaultAssay(object) <- "RNA"

#remove negatives from object
Idents(object) <- "HTO_classification.global"
object <- subset(object, idents = "Negative", invert = TRUE)
table(object$HTO_classification.global) 

saveRDS(object, file="object-ER1-W3.rds")

##----SNPs ----
#if importing
#object <- readRDS(file="object-ER1-W3.rds")

orig.meta <- object@meta.data

library(readr)

#load in clusters.tsv files from souporcell
batch1a <-readr::read_tsv("batchER1_clusters.tsv")
batch1b <-readr::read_tsv("batchW1_clusters.tsv")
batch1c <-readr::read_tsv("batchW2_clusters.tsv")
batch1d <-readr::read_tsv("batchW3_clusters.tsv")

#for each tsv, add sample prefix, combine with barcode, and overwrite barcode
batch1a$ident <- "ER1"
batch1a <- batch1a %>%
  unite('barcode', ident, barcode, remove=FALSE)

batch1b$ident <- "W1"
batch1b <- batch1b %>%
  unite('barcode', ident, barcode, remove=FALSE)

batch1c$ident <- "W2"
batch1c <- batch1c %>%
  unite('barcode', ident, barcode, remove=FALSE)

batch1d$ident <- "W3"
batch1d <- batch1d %>%
  unite('barcode', ident, barcode, remove=FALSE)

#combine tsv files, subset on singlets, and remove duplicate barcodes (though there shouldn't be)
batch1 <- rbind(batch1a,batch1b,batch1c,batch1d)
batch1 <- subset(batch1, subset =  status=="singlet")
batch1 <- batch1[!duplicated(batch1$barcode),]

batch1.2 <- data.frame(batch1$barcode, batch1$status, batch1$assignment, batch1$ident)
colnames(batch1.2) <- c("barcode","status","HTO_assignment","ident")

#add information to object metadata df
orig.meta$barcode <- rownames(orig.meta)
orig.meta$toname <- rownames(orig.meta)
comb3 <- merge(x = orig.meta, y = batch1.2, by = "barcode")
rownames(comb3) <- comb3$toname

comb3 <- subset(comb3, select = -c(toname))

#use merged cells to subset object
b <- comb3$barcode
object<- subset(object, cells=b)

#reassign metadata
object@meta.data <- comb3

saveRDS(object, file="ER1-W3_SNP.rds")

sce.object <- as.SingleCellExperiment(object)
sce.object <- scDblFinder(sce.object, dbr = (0.004 * ncol(sce.object)/1000), samples="HTO_maxID") #doublet rate of 0.4% per 1000 cells
table(sce.object$scDblFinder.class, sce.object$HTO_maxID)

identical(colnames(object),colnames(sce.object))

object$scDblFinder.class <- sce.object$scDblFinder.class
object$scDblFinder.score <- sce.object$scDblFinder.score
object$scDblFinder.weighted <- sce.object$scDblFinder.weighted
rm(sce.object)

object$scDblFinder.class <- factor(object$scDblFinder.class, levels=c("doublet", "singlet"))

#look at numbers of singlets and doublets using scDblFinder
object@meta.data %>%
  ggplot(aes(x= scDblFinder.class)) +
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -1.5)

##----Filter out doublets ----
table(object$scDblFinder.class)

#CD19 vs CD3 gene expression before filtering out doublets
FeatureScatter(object, feature1 = "CD3E", feature2 = "CD19", group.by = "scDblFinder.class",raster=F) + NoLegend() 

#subset to singlets 
object_singlet <- subset(object, subset =  scDblFinder.class=="singlet")

#subset singlets round 2 by excluding cells that are double positive for CD19 and CD3E genes 
FeatureScatter(object_singlet, feature1 = "CD3E", feature2 = "CD19") + NoLegend() + geom_vline(xintercept = 0.1) +  geom_hline(yintercept = 0.1)

object_singlet.2 <- subset(object_singlet, subset = CD3E < 0.1 & CD19 > 0.1 |CD3E > 0.1 & CD19 < 0.1 | CD3E < 0.1 & CD19 < 0.1)

#look at number of singlets
object_singlet.2@meta.data %>% 
  ggplot(aes(x=orig.ident)) + 
  geom_bar(aes(fill = HTO_maxID)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  geom_text(aes(label = ..count..), stat = "count", vjust = -0.5)

saveRDS(object_singlet.2, file = "multi-ER1-W3_postDblt_no_negative.rds")

rm(object, object_filtered, object_singlet.2, metaDat, s1, mito.genes, percent.mito, object_singlet)
rm(batch1, batch1.2, batch1a, batch1b, batch1c, batch1d, comb2, comb3, orig.meta, orig.meta.2, b, readData)

##----TCR/BCR----
setwd("C:/checkmate238/contig_annotations")
suppressMessages(library(scRepertoire))

###----TCR ----
object1_TCR <- read.csv("2023-07-27/TCR/multi_ER1_TCR_filtered_contig_annotations.csv")
object2_TCR <- read.csv("2023-07-27/TCR/multi_W1_TCR_filtered_contig_annotations.csv")
object3_TCR <- read.csv("2023-07-27/TCR/multi_W2_TCR_filtered_contig_annotations.csv")
object4_TCR <- read.csv("2023-07-27/TCR/multi_W3_TCR_filtered_contig_annotations.csv")

contig_list_TCR <- list(object1_TCR, object2_TCR, object3_TCR, object4_TCR)

combined_TCR <- combineTCR(contig_list_TCR, 
                           samples =  c("object1","object2","object3", "object4"),
                           ID = c("ER1", "W1", "W2", "W3"))

#adjust connector and num_connects based on barcode format
for (i in seq_along(combined_TCR)) {
  combined_TCR[[i]] <- stripBarcode(combined_TCR[[i]], 
                                    column = 1, connector = "_", num_connects = 3)
}

(setwd("C:/checkmate238"))
mergedObject.integrated.sct.TCR <- readRDS("multi-ER1-W3_postDblt_no_negative.rds")
view(mergedObject.integrated.sct.TCR)

setwd("C:/checkmate238/contig_annotations")

#adding the barcode_prefix to be the same as in mergedObject.integrated.sct
combined_TCR$object1_ER1$barcode_prefix <-"ER1"
combined_TCR$object2_W1$barcode_prefix <- "W1"
combined_TCR$object3_W2$barcode_prefix <- "W2"
combined_TCR$object4_W3$barcode_prefix <-"W3"

View(combined_TCR)

for (i in seq_along(combined_TCR)) {
  combined_TCR[[i]] <- combined_TCR[[i]]  %>% unite('barcode',barcode_prefix, barcode, remove=T)
}
View(combined_TCR)

#adding more meta data to combined
data_to_add <- rownames_to_column(mergedObject.integrated.sct.TCR@meta.data, "barcodes")
data_to_add <- data_to_add[,c("barcodes","HTO_maxID")] %>% as.data.frame()
colnames(data_to_add) <- c("barcode","HTO_maxID")
for (i in seq_along(combined_TCR)) {
  combined_TCR[[i]] <- left_join(combined_TCR[[i]], data_to_add, by = "barcode")}

saveRDS(combined_TCR, file="combined_TCR_ER1-W3.rds")

rm(i, contig_list_TCR, object1_TCR, object2_TCR, object3_TCR, object4_TCR, combined_TCR, data_to_add, mergedObject.integrated.sct.TCR)

###----BCR ----
object1_BCR <- read.csv("2023-07-27/BCR/multi_ER1_BCR_filtered_contig_annotations.csv")
object2_BCR <- read.csv("2023-07-27/BCR/multi_W1_BCR_filtered_contig_annotations.csv")
object3_BCR <- read.csv("2023-07-27/BCR/multi_W2_BCR_filtered_contig_annotations.csv")
object4_BCR <- read.csv("2023-07-27/BCR/multi_W3_BCR_filtered_contig_annotations.csv")

contig_list_BCR <- list(object1_BCR, object2_BCR, object3_BCR, object4_BCR)

combined_BCR <- combineBCR(contig_list_BCR, 
                           samples =  c("object1","object2","object3","object4"),
                           ID = c("ER1", "W1", "W2", "W3"))

#adjust connector and num_connects based on barcode format
for (i in seq_along(combined_BCR)) {
  combined_BCR[[i]] <- stripBarcode(combined_BCR[[i]],
                                    column = 1, connector = "_", num_connects = 3)
}

#see barcode format
View(combined_BCR)

(setwd("C:/checkmate238"))
mergedObject.integrated.sct.BCR <- readRDS("multi-ER1-W3_postDblt_no_negative.rds")
view(mergedObject.integrated.sct.BCR)

setwd("C:/checkmate238/contig_annotations")

#adding the barcode_prefix to be the same as in mergedObject.integrated.sct
combined_BCR$object1_ER1$barcode_prefix <- "ER1"
combined_BCR$object2_W1$barcode_prefix <- "W1"
combined_BCR$object3_W2$barcode_prefix <- "W2"
combined_BCR$object4_W3$barcode_prefix <- "W3"

View(combined_BCR)

for (i in seq_along(combined_BCR)) {
  combined_BCR[[i]] <- combined_BCR[[i]]  %>% unite('barcode',barcode_prefix, barcode, remove=T)
}
View(combined_BCR)

#adding more meta data to combined
data_to_add <- rownames_to_column(mergedObject.integrated.sct.BCR@meta.data, "barcodes")
data_to_add <- data_to_add[,c("barcodes","HTO_maxID")] %>% as.data.frame()
colnames(data_to_add) <- c("barcode","HTO_maxID")
for (i in seq_along(combined_BCR)) {
  combined_BCR[[i]] <- left_join(combined_BCR[[i]], data_to_add, by = "barcode")}

saveRDS(combined_BCR, file="combined_BCR_ER1-W3.rds")

rm(contig_list_BCR, object1_BCR, object2_BCR, object3_BCR, object4_BCR, combined_BCR, data_to_add, i, mergedObject.integrated.sct.BCR)

#----Integration ----
##----Reading in pre-processed objects and identifying each lane ----
(setwd("C:/checkmate238"))

#load in all pre-process objects per lane
object1 <- readRDS(file = "multi-L12-L34_postDblt_no_negative.rds")
object2 <- readRDS(file = "multi-ER1-W3_postDblt_no_negative.rds")

##--- Integrating lanes ----
#SOURCE: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
#using SCTransform and RPCA to integrate datasets quickly...default k.anchors = 5 for FindIntegrationFeatures
mergedObject <- merge(object1, y = c(object2)) #merge all lanes 
rm(object1, object2)

mergedObject$orig.ident <- as.factor(mergedObject$orig.ident)
levels(mergedObject$orig.ident)
DefaultAssay(mergedObject) <- "RNA"

#split the dataset into a list of the seurat objects
mergedObject <- SplitObject(mergedObject, split.by = "orig.ident")

#normalize using SCTransform and identify variable features for each dataset independently; this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
mergedObject <- lapply(X = mergedObject, FUN = SCTransform, method = "glmGamPoi")

#select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
#run the PrepSCTIntegration() function prior to identifying anchors
features <- SelectIntegrationFeatures(object.list = mergedObject, nfeatures = 5000)
mergedObject <- PrepSCTIntegration(object.list = mergedObject, anchor.features = features)
mergedObject <- lapply(X = mergedObject, FUN = RunPCA, features = features)

DefaultAssay(mergedObject$`L12`) #just to check; should be SCT

#identify anchors using FindIntegrationAnchors()(takes list of Seurat objects as input) and use these anchors to integrate the two datasets together with IntegrateData()
immune.anchors <- FindIntegrationAnchors(object.list = mergedObject, normalization.method = "SCT", 
                                         reference = c(1),  #reference datasets help speed up process, will reference first dataset in list 
                                         anchor.features = features, 
                                         dims = 1:20, 
                                         reduction = "rpca", 
                                         k.anchor = 20)

rm(mergedObject) #save space

mergedObject.integrated.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
mergedObject.integrated.sct <- RunPCA(mergedObject.integrated.sct)
mergedObject.integrated.sct <- RunUMAP(mergedObject.integrated.sct,reduction = "pca", dims = 1:20)
mergedObject.integrated.sct <- FindNeighbors(mergedObject.integrated.sct, reduction = "pca", dims = 1:20)
mergedObject.integrated.sct <- FindClusters(mergedObject.integrated.sct, resolution = 0.5)

#final merged dataset
#saveRDS(mergedObject.integrated.sct, file = "mergedObject.integrated.sct.rds") 

#mergedObject.integrated.sct <- readRDS(file ="C:/contain_cdrive/CONTAIN_ES/20230728_mergedObject.integrated.sct.rds")

p1 <- DimPlot(mergedObject.integrated.sct, label=T, raster = FALSE);p1

##----Azimuth ----
library(Azimuth)
library(SeuratData)
library(patchwork)

mergedObject.integrated.sct <- RunAzimuth(mergedObject.integrated.sct, reference = "pbmcref")

Idents(mergedObject.integrated.sct) <- "predicted.celltype.l2"

p2 <- DimPlot(mergedObject.integrated.sct, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, raster=F); p2
p1 + p2

metaDat <- mergedObject.integrated.sct@meta.data
metaDat$rownames <- rownames(metaDat)
mergeMeta <- merge(meta, metaDat, by.x = "subject3", by.y="subject")
rownames(mergeMeta) <- mergeMeta$rownames; mergeMeta$rownames <- NULL
mergedObject.integrated.sct <- AddMetaData(mergedObject.integrated.sct, mergeMeta)

#make sure all samples are listed
table(mergedObject.integrated.sct$orig.ident)

#final merged dataset
saveRDS(mergedObject.integrated.sct, file = "checkmateSNP_line141_azimuth_mergedObject.integrated.sct.rds") 

rm(list=ls())

#----PostIntegration ----
##----TCR ----
setwd("C:/checkmate238/contig_annotations")
combined_TCR_1 <- readRDS(file = "combined_TCR_L12-L34.rds")
combined_TCR_2 <- readRDS(file = "combined_TCR_ER1-W3.rds")

list.receptors.TCR <- c(combined_TCR_1, combined_TCR_2)

#use the object of interest (integrated object, Bcells, cd8t, etc.) to query into the TCR's; pulls out the TCR sequences for cells with matching barcodes in the object of interest

#change directory if the object of interest is stored separately from the contig_annotations
setwd("C:/checkmate238")
seurat_obj <- readRDS(file = "checkmate_line141_azimuth_mergedObject.integrated.sct.rds")
setwd("C:/checkmate238/contig_annotations")

combined.object.TCR <-combineExpression(list.receptors.TCR,
                                        seurat_obj,
                                        group.by="ID",
                                        cloneCall="strict",
                                        proportion=T)

#see TCR information in the right-hand columns
View(combined.object.TCR@meta.data[,34:40])

#final merged dataset with TCRs added 
saveRDS(combined.object.TCR, file = "checkmate.mergedObject.integrated.sct.TCR.rds") 

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

DimPlot(combined.object.TCR, group.by = "cloneType", raster=F) +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

ggplot(combined.object.TCR@meta.data, aes(fill=factor(cloneType), x=subject,), raster=F) + geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") + 
  scale_fill_manual(values = colorblind_vector(5), na.value="grey")

rm(list= ls())

##----BCR ----
setwd("C:/checkmate238/contig_annotations")

combined_BCR_1 <- readRDS(file = "combined_BCR_L12-L34.rds")
combined_BCR_2 <- readRDS(file = "combined_BCR_ER1-W3.rds")
 
list.receptors.BCR <- c(combined_BCR_1, combined_BCR_2)
 
#use the object of interest (integrated object, Bcells, cd8t, etc.) to query into the BCR's; pulls out the BCR sequences for cells with matching barcodes in the object of interest

#change directory if the object of interest is stored separately from the contig_annotations
setwd("C:/checkmate238")
seurat_obj <- readRDS(file="checkmate_line141_azimuth_mergedObject.integrated.sct.rds")
setwd("C:/checkmate238/contig_annotations")

combined.object.BCR <-combineExpression(list.receptors.BCR,
                                        seurat_obj,
                                        group.by="ID",
                                        cloneCall="strict",
                                         proportion=T)

#see BCR information in the right-hand columns
View(combined.object.BCR@meta.data[,34:40])

#final merged dataset with BCRs added 
saveRDS(combined.object.BCR, file = "checkmate.mergedObject.integrated.sct.BCR.rds") 
 
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                             "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                             "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
 
DimPlot(combined.object.BCR, group.by = "cloneType", raster=F) +
   scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
   theme(plot.title = element_blank())
 
ggplot(combined.object.BCR@meta.data, aes(fill=factor(cloneType), x=subject,), raster=F) + geom_bar(position="fill")+theme_classic()+labs(y = "Count", x = "Cluster") + 
   scale_fill_manual(values = colorblind_vector(5), na.value="grey")

rm(list= ls())