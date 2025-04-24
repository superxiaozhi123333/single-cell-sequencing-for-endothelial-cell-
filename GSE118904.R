#GSE118904
library(Seurat)
library(dplyr)
library(patchwork)
library(BiocManager)
library(cowplot) 
library(ggplot2) 
#Tumor
count1 <- read.table("GSE118904_rawcounts tumor.csv.gz",  header=T,stringsAsFactors = F,row.names = 1,sep=",")
Tumor<- CreateSeuratObject(count1,min.cells = 5,project = "Tumor",min.features = 200,names.field = 1,names.delim = "_")
Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^MT-")
Tumor <- subset(Tumor, subset = nCount_RNA > 3000 & nCount_RNA < 40000 & percent.mt < 10)
Tumor <- NormalizeData(Tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor <- FindVariableFeatures(Tumor, selection.method = "vst", nfeatures = 2000)
Tumor
#Normal
count2 <- read.table("GSE118904_rawcounts normal.csv.gz",  header=T,stringsAsFactors = F,row.names = 1,sep=",")
Normal<- CreateSeuratObject(count2,min.cells = 5,project = "Normal",min.features = 200,names.field = 1,names.delim = "_")
Normal[["percent.mt"]] <- PercentageFeatureSet(Normal, pattern = "^MT-")
Normal <- subset(Normal, subset = nCount_RNA > 3000 & nCount_RNA < 40000 & percent.mt < 10)
Normal <- NormalizeData(Normal, normalization.method = "LogNormalize", scale.factor = 10000)
Normal <- FindVariableFeatures(Normal, selection.method = "vst", nfeatures = 2000)
Normal
#merge
ALL.anchors <- FindIntegrationAnchors(object.list = list(Normal, Tumor), dims = 1:20)
total <- IntegrateData(anchorset = ALL.anchors, dims = 1:20)

#Perform an integrated analysis
DefaultAssay(total) <- "integrated"

total <- ScaleData(total)
total <- RunPCA(total, features = VariableFeatures(object = total))
total <- FindNeighbors(total, dims = 1:10)
total <- FindClusters(total, resolution = 0.5)
total <- RunUMAP(total, dims = 1:10)
total <- RunTSNE(total, dims = 1:10)
DimPlot(total, reduction = "umap", label = T, label.size = 5)
DimPlot(total, reduction = "tsne", label = T, label.size = 5)

DimPlot(total, reduction = "umap", group.by = "orig.ident")
DimPlot(total, reduction = "umap", split.by = "orig.ident")

EC.Marker <- FindMarkers(total, ident.1 = "Tumor", ident.2 = "Normal", verbose = FALSE,logfc.threshold = 0,min.pct = 0)
EC.Marker["Pdlim5",]



#GSE110501
library(Seurat)
library(dplyr)
library(patchwork)
library(BiocManager)
#Tumor
count1 <- read.table("GSM2994872_898096.colo205.xenograft.hFc.24hr.mm10.UMI.txt.gz",header = T, stringsAsFactors = F,row.names = 1)
Tumor<- CreateSeuratObject(count1,min.cells = 10,project = "Tumor",min.features = 200,names.field = 1,names.delim = "_")
Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^MT-")
Tumor <- subset(Tumor, subset = nCount_RNA > 3000 & nCount_RNA < 40000 & percent.mt < 10)
Tumor <- NormalizeData(Tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Tumor <- FindVariableFeatures(Tumor, selection.method = "vst", nfeatures = 2000)
#Heart
count2 <- read.table("GSM2994876_898100.heart.hFc.24hr.UMI.txt.gz",header = T, stringsAsFactors = F,row.names = 1)
Heart<- CreateSeuratObject(count2,min.cells = 10,project = "Heart",min.features = 200,names.field = 1,names.delim = "_")
Heart[["percent.mt"]] <- PercentageFeatureSet(Heart, pattern = "^MT-")
Heart <- subset(Heart, subset = nCount_RNA > 3000 & nCount_RNA < 40000 & percent.mt < 10)
Heart <- NormalizeData(Heart, normalization.method = "LogNormalize", scale.factor = 10000)
Heart <- FindVariableFeatures(Heart, selection.method = "vst", nfeatures = 2000)
#merge EC=extracted cells 
EC.anchors <- FindIntegrationAnchors(object.list = list(Heart, Tumor), dims = 1:20)
EC <- IntegrateData(anchorset = EC.anchors, dims = 1:20)
#Perform an integrated analysis
DefaultAssay(EC) <- "integrated"
EC <- ScaleData(EC)
EC <- RunPCA(EC, features = VariableFeatures(object = EC))
EC <- FindNeighbors(EC, dims = 1:10)
EC <- FindClusters(EC, resolution = 0.5)
EC <- RunUMAP(EC, dims = 1:10)
EC <- RunTSNE(EC, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = T, label.size = 5)
DimPlot(EC, reduction = "umap", label = T, label.size = 5)
DimPlot(EC, reduction = "tsne", label = T, label.size = 5)
View(EC)
DimPlot(EC, reduction = "umap", label = T, label.size = 5,split.by = orig.ident)
DimPlot(EC, reduction = "umap", label = T, label.size = 5,split.by = "orig.ident")
#Endothelial cells
VlnPlot(EC, features = c("Pecam1"))
VlnPlot(EC, features =   c("Pecam1"), slot = "counts", log = TRUE)
FeaturePlot(EC, features = "Pecam1", pt.size = 0.5, cols = c("grey80", "red", "blue"))
#Fibroblasts
VlnPlot(EC, features = c("Thy1"))
FeaturePlot(EC, features = "Thy1", pt.size = 0.5, cols = c("grey80", "red", "blue"))
#Immune cells
VlnPlot(EC, features = c("Ptprc"))
FeaturePlot(EC, features = "Ptprc", pt.size = 0.5, cols = c("grey80", "red", "blue"))
DefaultAssay(EC) <- "RNA"
#annotation
new.cluster.ids <- c( "Fibroblasts",   # Thy1
                      "Fibroblasts",    # Thy1
                      "Fibroblasts",  # Thy1
                      "Fibroblasts",   # Thy1
                      "Endothelial cells",         # Pecam1
                      "Endothelial cells",         # Pecam1
                      "Fibroblasts",             # Thy1
                      "Immune cells",           # Ptprc Cd68
                      "Fibroblasts",          # Thy1
                      "Endothelial cells",         # Pecam1
                      "Immune cells"    )      # Ptprc Cd68     # PPBP
names(new.cluster.ids) <- levels(EC)
new.cluster.ids <- c( "Fibroblasts",   # Thy1
                      "Fibroblasts",    # Thy1
                      "Fibroblasts",  # Thy1
                      "Fibroblasts",   # Thy1
                      "Endothelial cells",         # Pecam1
                      "Endothelial cells",         # Pecam1
                      "Fibroblasts",             # Thy1
                      "Immune cells",           # Ptprc Cd68
                      "Fibroblasts",          # Thy1
                      "Endothelial cells",         # Pecam1
                      "Immune cells" ,
                      "Endothelial cells")      # Ptprc Cd68     # PPBP
names(new.cluster.ids) <- levels(EC)
new.cluster.ids
EC <- RenameIdents(EC, new.cluster.ids)
Idents(EC)
head(EC@meta.data)
table(EC@meta.data$orig.ident)
table(EC@meta.data$celltype)
Idents(EC)
EC@metadata$celltype=Idents(EC)
View(EC)
EC@meta.data$celltype=Idents(EC)
table(EC@meta.data$celltype)
EC@meta.data$group=EC@meta.data$orig.ident
Total=EC

#EC=Endothelial cells
EC= Total[,Total@meta.data$celltype%in% c("Endothelial cells")]
table(EC@meta.data$celltype)
table(EC@meta.data$group)
cluster1.markers <- FindMarkers(EC, ident.1 = "Tumor", ident.2 = "Heart", verbose = FALSE,assay = 'RNA',logfc.threshold = 0,min.pct = 0)
Ident(EC)=EC@meta.data$group
Idents(EC)=EC@meta.data$group
cluster1.markers <- FindMarkers(EC, ident.1 = "Tumor", ident.2 = "Heart", verbose = FALSE,assay = 'RNA',logfc.threshold = 0,min.pct = 0)
TECvsNEC_deg=cluster1.markers
#export
write.csv(TECvsNEC_deg,file = "TECvsNEC_deg.csv")


#GSE131907
library(Seurat)
library(dplyr)
library(patchwork)
library(BiocManager)
library(cowplot) 
library(ggplot2) 
library(stringr)
library(paletteer)
library(MySeuratWrappers)
library(SCpubr)
library(viridis)
library(ComplexHeatmap)

counts <- Read10X_h5("NSCLC_GSE131907_expression.h5")
Total<- CreateSeuratObject(counts,min.cells = 10,min.features = 200,names.field = 1,names.delim = "_")
Total
Total[["percent.mt"]] <- PercentageFeatureSet(Total, pattern = "^MT-")
Total <- subset(Total, subset =  nCount_RNA > 200 & nCount_RNA < 10000 & percent.mt < 20)
Total <- NormalizeData(Total, normalization.method = "LogNormalize", scale.factor = 10000)
Total <- FindVariableFeatures(Total, selection.method = "vst", nfeatures = 2000)

Total <- ScaleData(Total)
Total<- RunPCA(Total, features = VariableFeatures(object =Total))

Total <- FindNeighbors(Total, dims = 1:10)
Total <- FindClusters(Total, resolution = 0.5)
table(Total@meta.data$seurat_clusters)

Total <- RunUMAP(Total, dims = 1:10)
Total <- RunTSNE(Total, dims = 1:10)
DimPlot(Total, reduction = "umap", label = T, label.size = 3)
DimPlot(Total, reduction = "tsne", label = T, label.size = 3)
table(Total@meta.data$sample)

#grepl("LUNG_N06",table$cust_id)
Total$group <- "NA"
Total$probe_sample<- "NA"
Total$group <- "NA"
Total$probe_sample<- "NA"
for(i in 165854:203298){
  t=colnames(Total[,i])
  if (grepl("LUNG_N",t)) {Total$sample[i] = "LUNG_N"} else { 
    if (grepl("LUNG_T",t)) {Total$sample[i] = "LUNG_T"} else{
      if (grepl("EBUS",t)) {Total$sample[i] = "EBUS"} else{  
        if (grepl("BRONCHO",t)) {Total$sample[i] = "BRONCHO"} else{ 
          if (grepl("LN",t)) {Total$sample[i] = "LN"} else{ 
            if (grepl("EFFUSION",t)) {Total$sample[i] = "EFFUSION"} else{ 
              if (grepl("NS",t)) {Total$sample[i] = "NS"} }}}}}}
  print(i)
}

Total$sample <- "NA"
for(i in 1:203298){
  if (Total$orig.ident[i]%in% c("Tumor1","Tumor2","Tumor3","Tumor4")) {Total$group[i] = "Tumor"} 
  if (Total$orig.ident[i]%in% c("Normal1","Normal2","Normal3","Normal4")) {Total$group[i] = "Normal"} 
  
}

table(Total@meta.data$active.ident)

Lung_Normal = Lung_total[,Lung_total@meta.data$sample%in% c("LUNG_N")]


Lung_total = Total[,Total@meta.data$sample%in% c("LUNG_N","LUNG_T")]
table(Lung_total@meta.data$sample)
save(Lung_total,file = 'Lung_total.Rdata') 
load(file = 'Lung_total.Rdata')

Lung_total <- NormalizeData(Lung_total, normalization.method = "LogNormalize", scale.factor = 10000)
Lung_total <- FindVariableFeatures(Lung_total, selection.method = "vst", nfeatures = 2000)

Lung_total<- ScaleData(Lung_total)
Lung_total<- RunPCA(Lung_total, features = VariableFeatures(object =Lung_total))

Lung_total <- FindNeighbors(Lung_total, dims = 1:10)
Lung_total<- FindClusters(Lung_total, resolution = 0.5)
table(Lung_total@meta.data$seurat_clusters)

Lung_total<- RunUMAP(Lung_total, dims = 1:10)
Lung_total <- RunTSNE(Lung_total, dims = 1:10)
DimPlot(Lung_total, reduction = "umap", label = T, label.size = 3)
DimPlot(Lung_total, reduction = "tsne", label = T, label.size = 3)
table(Lung_total@meta.data$sample)
table(Lung_total@meta.data$celltype)

DimPlot(Lung_total, reduction = "umap", group.by = "sample")
DimPlot(Lung_total, reduction = "tsne", group.by = "sample")
DimPlot(Lung_total, reduction = "umap", split.by = "sample")
DimPlot(Lung_total, reduction = "tsne", split.by = "sample")

Lung_total.markers <- FindAllMarkers(Lung_total, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(Lung_total.markers)
abc= Lung_total.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(abc,file = "Lung_total.markers10.csv")

new.cluster.ids <- c("CD4+T cells",   #
                     "Macrophages",    # 
                     "CD8+T cells",  # 
                     "NK cells",   # 
                     "Macrophages",    # 
                     "Monocyte-DCs",  #  
                     "B cells",    # 
                     "CD8+T cells",     # 
                     "Alveolar",             #
                     "Fibroblasts",     # 
                     "Mast cells",     # 
                     "Cancer",     # 
                     "Alveolar",     # 
                     "Endothelial",     # 
                     "Macrophages",     # 
                     "Epithelial cells",     # 
                     "Alveolar",     # 
                     "pDCs")      # 
names(new.cluster.ids) <- levels(Lung_total)
new.cluster.ids
Lung_total<- RenameIdents(Lung_total, new.cluster.ids)

DimPlot(Lung_total, reduction = "umap", label = T, label.size = 3)
DimPlot(Lung_total, reduction = "tsne", label = T, label.size = 3)

DimPlot(Lung_total, reduction = "umap", group.by = "sample")
DimPlot(Lung_total, reduction = "tsne", group.by = "sample")
DimPlot(Lung_total, reduction = "umap", split.by = "sample")
DimPlot(Lung_total, reduction = "tsne", split.by = "sample")

Lung_total$celltype <- Idents(Lung_total)
Lung_total$celltype.sample <- paste(Idents(Lung_total),Lung_total$sample, sep = "_")
Idents(Lung_total) <- "celltype.sample"
Idents(Lung_total) <- "celltype"
table(Lung_total@meta.data$celltype)
My_levels <- c("CD4+T cells","Macrophages","CD8+T cells","NK cells","Monocyte-DCs","B cells","Alveolar","Fibroblasts","Mast cells","Endothelial","Epithelial cells","pDCs","Cancer")
Idents(Lung_total) <- factor(Idents(Lung_total), levels= My_levels)
#Except for tumor cells
Lung_total2= Lung_total[,Lung_total@meta.data$celltype!= c("Cancer")]
Idents(Lung_total2) <- "celltype"
VlnPlot(Lung_total2, split.by = "sample",features = c("PDLIM5"))

table(Lung_total@meta.data$celltype)
Lung_EC= Lung_total[,Lung_total@meta.data$celltype%in% c("Endothelial")]
table(Lung_EC@meta.data$sample)
table(Lung_EC@meta.data$celltype)
save(Lung_EC,file = 'Lung_EC.Rdata') 
load(file = 'Lung_EC.Rdata')
head(Lung_EC@meta.data)
Idents(Lung_EC)<- "sample"
cluster1.markers <- FindMarkers(Lung_EC, ident.1 = "LUNG_T", ident.2 = "LUNG_N", verbose = FALSE,assay = 'RNA',logfc.threshold = 0,min.pct = 0)
TECvsNEC_deg=cluster1.markers

write.csv(TECvsNEC_deg,file = "TECvsNEC_deg.csv")


deg <- read.table("TECvsNEC_deg.csv",  header=T,stringsAsFactors = F,row.names = 1,sep=",")
genelist <- read.table("genelist.csv",  header=T,stringsAsFactors = F,sep=",", fileEncoding = 'utf-8')
exp<-deg

genelist$ID=toupper(genelist$ID)

exp$ID <- rownames(exp)
colnames(genelist) <- c('ID')
genelist$ID=toupper(genelist$ID)
express <- merge(x = exp, y =genelist , by = "ID")
rownames(express)=express$ID
express$ID =NULL

write.csv(express,file = "selectedgene_deg.csv")

deg=express
head(deg)

deg <-  read.table("selectedgene_deg.csv",  header=T,stringsAsFactors = F,sep=",", fileEncoding = 'utf-8')

library(limma)

logFC=0.1
P.Value = 0.05
k1 = (deg$p_val < P.Value)&(deg$avg_log2FC < -logFC)
k2 = (deg$p_val < P.Value)&(deg$avg_log2FC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)

k1 = (express$p_val < P.Value)&(express$avg_log2FC < -logFC)
k2 = (express$p_val < P.Value)&(express$avg_log2FC > logFC)
express$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(express$change)


logfc.cut=0.1
p.cut=0.05
library(ggplot2)
p<-ggplot(data = deg, aes(x = avg_log2FC,y = -log10(p_val)))+
  geom_point(alpha=0.8, size=2, aes(color=change)) +
  xlab(expression('log'[2]*'(FoldChange)'))+ 
  ylab(expression('-log'[10]*'(pvalue)'))+
  geom_vline(xintercept=c(-logfc.cut,logfc.cut),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(p.cut),lty=2,col="black",lwd=0.5) +
  theme_test() 
p

p+scale_color_manual(values=c("#6697ea", "grey","#b02428"))



library(ggrepel)
options(ggrepel.max.overlaps = Inf)

interested_genes<-c("PDLIM1","PDLIM2","PDLIM3","PDLIM4","PDLIM5","PDLIM7","CXCR4","MYO1B","SLC6A4","S100A4","SYNPO")
rownames(deg)=deg$X
deg$Label<-ifelse(rownames(deg)%in%interested_genes,rownames(deg),'')


ggplot(deg,aes(avg_log2FC,-log10(p_val),fill=change))+
  geom_point(shape=21,alpha=1,size=2,aes(color=change))+
  scale_fill_manual(values=c("#6697ea", "grey60","#b02428"))+ 
  scale_color_manual(values=c("#6697ea", "grey60","#b02428"))+ 
  guides(size=F)+  
  geom_label_repel(aes(x = avg_log2FC,y = -log10(p_val),size=2,label =Label),fill = "white")+
  xlab(expression('log'[2]*'(FoldChange)'))+
  ylab(expression('-log'[10]*'(pvalue)'))+
  geom_vline(xintercept=c(-logfc.cut,logfc.cut),lty=2,col="black",lwd=0.5) +
  geom_hline(yintercept = -log10(p.cut),lty=2,col="black",lwd=0.5) +
  theme_test()+
  xlim(-1.5,1.5)+ylim(0,28)+
  annotate('text',x=-4.5,y=27,fontface="bold.italic",size=7,color='red',
           label='Genes \n Down Regulated')+
  annotate('text',x=4.5,y=27,fontface="bold.italic",size=7,color='red',
           label='Genes \n  Up Regulated')
