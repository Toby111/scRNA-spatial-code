rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)

#####load data####
Sham<-Read10X("./Ctrl/")
AKI<-Read10X("./IRD3/")

###labeling###
colnames(Sham) <- paste(colnames(Sham),"Sham",sep = "_")
colnames(AKI) <- paste(colnames(AKI),"AKI",sep = "_")

#####CreateSeuratObject########
Sham <- CreateSeuratObject(
  Sham,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
Sham[["percent.mt"]] <- PercentageFeatureSet(Sham,pattern = "^mt-")
Sham<-NormalizeData(Sham,verbose = FALSE)
Sham<- FindVariableFeatures(Sham, selection.method = "vst",nfeatures = 2000)

AKI <- CreateSeuratObject(
  AKI,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
AKI[["percent.mt"]] <- PercentageFeatureSet(AKI,pattern = "^mt-")
AKI<-NormalizeData(AKI,verbose = FALSE)
AKI<- FindVariableFeatures(AKI, selection.method = "vst", nfeatures = 2000)


sce = merge(x = Sham, y = AKI)

##QC####
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
ggsave2(filename = "QC-VlnPlot.pdf",width = 8,height = 4.5)
dev.off()

#### filter####
cat("Before filter :",nrow(sce@meta.data),"cells\n")
sce <- subset(sce, subset = 
                nFeature_RNA > 500 & 
                nCount_RNA > 1000 & 
                nCount_RNA < 20000 &
                percent.mt < 40)
cat("After filter :",nrow(sce@meta.data),"cells\n")

sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
DimPlot(sce, reduction = "pca")

####batch effects correction####
library(harmony)
seuratObj <- RunHarmony(sce, "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "harmony",label=T) 
DimPlot(seuratObj,reduction = "umap",label=T) 

sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = 1:15) 
sce <- FindClusters(sce, resolution = 0.8)
DimPlot(sce,reduction = "umap",label=T) 

save(sce,file="sce_AKI_2S.RData")


dim.use <- 1:20


library("scales")
mypal = pal_d3("category20", alpha = 1)(20)
mypal
show_col(mypal)

sce <- RunUMAP(sce, dims = dim.use, do.fast = TRUE)
DimPlot(object = sce, pt.size=0.1,label = T,cols = mypal)
ggsave2(filename = "CellCluster-UMAP.pdf", width = 7,height = 6)
dev.off()

DimPlot(object = sce, 
        split.by="orig.ident", label = T,
        pt.size=0.1,reduction = "umap", cols = mypal)
ggsave2(filename = "CellCluster-UMAP-split.pdf", width = 18,height = 6)
dev.off()


###########mouse kidney annotation#########
####marker gene expression####
DotPlot(sce, features = c("Itgax","Ptprc","Cd68","Adgre1","Cd79a","Cd79b", "Cd3d","Cd3e","Lrp2","Pck1","Slc12a1","Umod",
                           "Atp6v0d2","Atp6v1g3","Aqp2","Fxyd4","S100a8","S100a9","Uncx", "Cited1",
                           "Pecam1","Egfl7","Nphs1", "Nphs2","Col1a1","Col3a1","Trpv5","Slc12a3","Pvalb","Slc34a1",
                          "Slc5a2","Kap","Epcam", "Havcr1", "Krt8", "Spp1", "Spp2", "Vcam1"), cols = c("white", "#66101F"))
ggsave2(filename = "DotPlot-markergene_zs.pdf", width = 25,height = 8)
dev.off()

###annotation######
sce <- RenameIdents(object = sce,
  "0" = "MAC","1" = "PT","2" = "PT","3" = "Neutrophil","4" = "PT",
  "5" = "MAC","6" = "CDPC","7" = "LOH","8" = "PT",
  "9" = "PT","10" = "PT","11" = "LOH",
  "12" = "T cell","13" = "PT","14" = "Podocyte","15" = "CDIC","16" = "PT",
  "17" = "EC","18" = "MAC",  "19" = "PT","20" = "B cell",
  "21" = "Fibroblast","22" = "Podocyte","23" = "PT","24" = "MAC",
  "25" = "T cell","26" = "Fibroblast")
DimPlot(object = sce, pt.size= 0.1,
        reduction = "umap",label = T,cols = mypal)
ggsave2(filename = "CellCluster-UMAPPlot_Rename_.pdf", width = 6,height = 4)
dev.off()
DimPlot(object = sce, pt.size=0.1,
        reduction = "umap",label = T, label.size = 4,cols = mypal,split.by = "orig.ident")
ggsave2(filename = "CellCluster-UMAPPlot_Rename_split.pdf", width = 12,height = 4)
dev.off()

######cell percentage#######
meta_data <- sce@meta.data 
plot_data <- data.frame(table(meta_data$orig.ident,meta_data$celltype))
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
pic_percentage <- ggplot(plot_data,aes(x = Var1,y = Percentage,fill = Var2)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + labs(fill = "Cluster")+
  theme_bw() + scale_fill_d3()
ggsave(pic_percentage, file='pic_percentage_other.pdf', width=12, height=10)


####################enrichment################
sce_PT = subset(sce, celltype == "PT")
table(Idents(sce_PT))
Idents(sce_PT) = sce_PT$orig.ident

deg<- FindMarkers(sce_PT,
                  ident.1 = "AKI",
                  ident.2 = "Sham")

sig_deg<- subset(deg, p_val_adj<0.05 & abs(avg_log2FC)>1)

write.csv(deg,'all_degM.csv')
write.csv(sig_deg,'asig_degM.csv')

library(data.table)
deg=fread('all_degM.csv')
sig_deg=fread('asig_degM.csv')

## GSEA
###############################################################################
library(Seurat)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(fgsea)
library(enrichplot)


alldiff <- deg[order(deg$avg_log2FC,decreasing = T),]
id <- alldiff$avg_log2FC
names(id) <- rownames(alldiff)
id
###gmtfile <- "h.all.v7.4.symbols.gmt"
gmtfile <- "mh.all.v2023.1.Mm.symbols.gmt"

hallmark <- read.gmt(gmtfile)

y <- GSEA(id,TERM2GENE = hallmark)

save(y,file = "gsea_mac.Rdata")

dotplot(y,showCategory=20,split=".sign")+facet_grid(~.sign)

yd <- data.frame(y)

write.table(yd,
            file="AKI_mac_GSEA_hallmark.txt",
            sep="\t",quote = F,row.names = T)




######spatial transcriptomics#####


rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)
library(dplyr)
library(hdf5r)
library(hdf5r)
packageVersion("Seurat")
################################################################################
Sham <- Load10X_Spatial(data.dir ="./data/Sham/SpaceRanger/", 
                            filename = "filtered_feature_bc_matrix.h5",
                               slice ="Sham")
Sham@project.name <-"Sham"
Idents(Sham) <-"Sham"
Sham$orig.ident <-"Sham"
Sham


AKI <- Load10X_Spatial(data.dir ="./data/AKI/SpaceRanger/", 
                            filename = "filtered_feature_bc_matrix.h5",
                               slice ="AKI")
AKI@project.name <-"AKI"
Idents(AKI) <-"AKI"
AKI$orig.ident <-"AKI"
AKI

##################################################
Sham <- SCTransform(Sham, assay = "Spatial", verbose = FALSE)
Sham <- RunPCA(Sham, assay = "SCT", verbose = FALSE) 
plot1 <- DimPlot(Sham, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(Sham, ndims=20, reduction="pca") 
plot1+plot2
ggsave("Sham_PCA.pdf",width = 10,height=5)

pc.num=1:15
Sham <- FindNeighbors(Sham, reduction = "pca", dims = pc.num)
Sham <- FindClusters(Sham, verbose = FALSE,resolution = 0.3)
# UMAP####
Sham <- RunUMAP(Sham, reduction = "pca", dims =  pc.num)
p1 <- DimPlot(Sham, reduction = "umap", label = TRUE)
# SpatialDimPlot####
p2 <- SpatialDimPlot(Sham, label = TRUE, label.size = 3)
p1 + p2
ggsave("Sham_Spatial.pdf",width = 10,height=5)


###############################################################################
position <-data.table::fread("./data/Sham/SpaceRanger/spatial/tissue_positions_list.csv")
position <- as.data.frame(position)

colnames(position) <-c("barcode",
                       "in_tissue",
                       "row",
                       "col",
                       "pxl_row",
                       "pxl_col")
position <- position[which(position$in_tissue=="1"),]

str(position)

position$row <-as.numeric(position$row)
position$col <-as.numeric(position$col)

ggplot(position, aes(x = row, y = col)) + 
  geom_point(data = position,color="blue")+
  theme_bw()


meta <-Sham@meta.data
rownames(position) <-position$barcode
position <-position[rownames(meta),]
meta$x=position$row
meta$y=position$col

p3 <-ggplot(meta, aes(x = x, y = y)) + 
     geom_point(data = meta, aes(color=seurat_clusters))+
     theme_bw()

save(position,file = "Sham_position.Rds")
save(IRD3,file = "Sham_stRNA.Rds")
###################################################################################

AKI <- Load10X_Spatial(data.dir ="./data/AKI/SpaceRanger/", 
                      filename = "filtered_feature_bc_matrix.h5",
                      slice ="AKI")
AKI@project.name <-"AKI"
Idents(AKI) <-"AKI"
AKI$orig.ident <-"AKI"
AKI

##################################################
AKI <- SCTransform(AKI, assay = "Spatial", verbose = FALSE)
AKI <- RunPCA(AKI, assay = "SCT", verbose = FALSE) 

pc.num=1:15

AKI <- FindNeighbors(AKI, reduction = "pca", dims = pc.num)
AKI <- FindClusters(AKI, verbose = FALSE,resolution = 0.3)

AKI <- RunUMAP(AKI, reduction = "pca", dims =  pc.num)
p1 <- DimPlot(AKI, reduction = "umap", label = TRUE)

################################################################################
position <-data.table::fread("./data/AKI/SpaceRanger/spatial/tissue_positions_list.csv")
position <- as.data.frame(position)

colnames(position) <-c("barcode",
                       "in_tissue",
                       "row",
                       "col",
                       "pxl_row",
                       "pxl_col")
position <- position[which(position$in_tissue=="1"),]

str(position)

position$row <-as.numeric(position$row)
position$col <-as.numeric(position$col)

ggplot(position, aes(x = row, y = col)) + 
  geom_point(data = position,color="blue")+
  theme_bw()


meta <-AKI@meta.data
rownames(position) <-position$barcode
position <-position[rownames(meta),]
meta$x=position$row
meta$y=position$col

p3 <-ggplot(meta, aes(x = x, y = y)) + 
  geom_point(data = meta, aes(color=seurat_clusters))+
  theme_bw()

save(position,file = "AKI_position.Rds")
save(AKI,file = "AKI_stRNA.Rds")
###################################################################################

load("./Sham_stRNA.Rds")
load("./AKI_stRNA.Rds")

p1 <- SpatialFeaturePlot(Sham, features ="Lcn2",min.cutoff = 0,
                         max.cutoff = 2.5)

p2 <- SpatialFeaturePlot(AKI, features ="Lcn2",min.cutoff = 0,
                         max.cutoff = 2.5)

p1 | p2 

p3 <- VlnPlot(Sham, features ="Lcn2",
              pt.size = 0, group.by = "orig.ident", y.max = 3
              )
p4 <- VlnPlot(AKI, features ="Lcn2",
              pt.size = 0, group.by = "orig.ident", y.max = 3
              )

p1 | p2 | p3 | p4 

ggsave("Compare_Lcn2.pdf",width = 12,height=5)



p1 <- SpatialDimPlot(Sham, label = T, label.size = 3)
p2 <- SpatialDimPlot(AKI, label = T, label.size = 3)
p1 | p2 
ggsave("Compare-total.pdf",width = 12,height=5)




FeaturePlot(IRD3,features = c("Otud5","Gpx4"),blend=T, blend.threshold = 0.1, cols = 
              c("lightgrey","#A20643","#5E4FA2"))
ggsave2("Merge_IRD3_Otud5-Gpx4.pdf",width = 15,height=5)













