rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(patchwork)
library(harmony)

S<-Read10X("./sham/")
A<-Read10X("./aki/")

colnames(A) <- paste(colnames(A),"A",sep = "_")
colnames(S) <- paste(colnames(S),"S",sep = "_")

A <- CreateSeuratObject(
  A,
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
A[["percent.MT"]] <- PercentageFeatureSet(A,pattern = "^mt-")
A<-NormalizeData(A,verbose = FALSE)
A<- FindVariableFeatures(A, selection.method = "vst",nfeatures = 2000)

S <- CreateSeuratObject(
  S,
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
S[["percent.MT"]] <- PercentageFeatureSet(S,pattern = "^MT-")
S<-NormalizeData(S,verbose = FALSE)
S<- FindVariableFeatures(S, selection.method = "vst",nfeatures = 2000)
sce = merge(x = A, y = S)
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
sce <- RunUMAP(sce, dims = 1:15, reduction = "pca")
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:15) 
save(sce,file="sce.RData")
cat("Before filter :",nrow(sce@meta.data),"cells\n")
sce <- subset(sce, subset = 
                nFeature_RNA > 500 & 
                nCount_RNA > 1000 & 
                nCount_RNA < 20000 &
                percent.MT < 25)
cat("After filter :",nrow(sce@meta.data),"cells\n")
all.markers <- FindAllMarkers(sce, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)



