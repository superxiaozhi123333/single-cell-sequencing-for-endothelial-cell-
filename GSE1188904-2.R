.libPaths("/home/GZY/R/x86_64-conda-linux-gnu-library/4.2")
library(Seurat)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(tidyr)
library(viridis)
library(MySeuratWrappers)
#############build folders for results############
dir.create('1.QC')
dir.create('2.Cluster')
my36colors <-c('#E95C59','#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
############ load data#####################
data_sample <- read.csv(file="GSE118904_rawcounts_negs.csv", header= T,row.names = 1) #导入数据
# 数据类型为行是细胞，列是基因
dat <- as.matrix(data_sample) #dataframe转换为矩阵类型
dat[1:3,1:3]
allcell_raw <- CreateSeuratObject(dat,project = "GSE118904")
Idents(allcell_raw)<-allcell_raw$orig.ident
allcell_raw<-subset(allcell_raw,idents='tumor')
############ QC#############
#mitochondrial QC metrics
allcell_raw[["percent.mt"]] <- PercentageFeatureSet(allcell_raw, pattern = "^mt-")
# Visualize QC metrics as a violin plot
Idents(allcell_raw)<-allcell_raw$orig.ident
pdf('1.QC/QC_metrics.pdf',width=10,height = 4)
dev.off()
VlnPlot(allcell_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3,pt.size = 0.00)
# feature-feature relationships
pdf('1.QC/feature-feature relationships.pdf',width=10,height = 5)
plot1 <- FeatureScatter(allcell_raw, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(allcell_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#ggplot2_nUMI_percent.mito
p1 <- ggplot(allcell_raw@meta.data,aes(x=nCount_RNA,y=percent.mt,fill=nCount_RNA)) +
  geom_point(alpha=0.5,shape=21,colour="black") +
  scale_fill_gradientn(colours=c("yellow","green","purple1","purple2","purple3","purple4")) +
  geom_hline(yintercept=5,linetype="dashed",color="red") +
  geom_vline(xintercept=1000,linetype="dashed",color="blue") +
  xlim(0,15000) +
  theme_bw()
ggsave(p1,filename="1.QC/ggplot2_nUMI_percent.mito.pdf",width=5,height=4)
saveRDS(allcell_raw,'allcell_raw.rds')

##filter
allcell <- subset(allcell_raw, subset = 
                    nFeature_RNA > 400 &
                    nCount_RNA >1000 &
                    nCount_RNA <25000 & 
                    percent.mt < 10 
)
table(allcell$orig.ident)
#control   tumor 
#1802     849 
table(allcell_raw$orig.ident)
##control   tumor 
##1802     849 
saveRDS(allcell,'allcell.rds')

######标准化######
DefaultAssay(allcell)<-"RNA"
allcell <- NormalizeData(allcell)
allcell <- FindVariableFeatures(allcell, selection.method = "vst", nfeatures = 2000)
allcell <- ScaleData(allcell)
allcell <- RunPCA(allcell)
ElbowPlot(allcell)
pc.num=1:30
allcell <- RunUMAP(allcell, dims=pc.num)
allcell <- RunTSNE(allcell, dims=pc.num)
allcell <- FindNeighbors(allcell, dims = pc.num) %>% FindClusters(resolution = 0.5)
Idents(allcell) <- allcell@meta.data$seurat_cluster
new.cluster.ids <- c(1:length(levels(allcell)))
names(x = new.cluster.ids) <- levels(x = allcell)
allcell <- RenameIdents(object = allcell, new.cluster.ids)
allcell[["Clusters"]]<-Idents(allcell)
dev.off()
############ find markers###################
DefaultAssay(allcell)<-"RNA"
allcell.markers_res0.5 <- FindAllMarkers(allcell, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allcell.markers_res0.5 <- allcell.markers_res0.5
allcell.markers_res0.5$nLog_padj <- -log10(allcell.markers_res0.5$p_val_adj) 
allcell.markers_res0.5$ratio <- allcell.markers_res0.5$pct.1/allcell.markers_res0.5$pct.2
out <- cbind(gene=allcell.markers_res0.5$gene,allcell.markers_res0.5)
write.table(out,'2.Cluster/allcell_marker_res0.5.xls',sep = '\t',quote = F,row.names = F)
write.csv(allcell.markers_res0.5,"2.Cluster/allcell_markers_res0.5.csv")
top10markers <- allcell.markers_res0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DotPlot(allcell,features =rev(unique(top10markers$gene)),dot.scale = 4)+RotatedAxis()+
  theme(axis.text.x = element_text(face="italic",angle = 90,vjust = 0.5),axis.title = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  scale_color_viridis(option = "D")#20X5.1

top10ratio <- allcell.markers_res0.5 %>% group_by(cluster) %>% top_n(n = 10, wt = ratio)
DotPlot(allcell,features =rev(unique(top10ratio$gene)),dot.scale = 4)+RotatedAxis()+
  theme(axis.text.x = element_text(face="italic",angle = 90,vjust = 0.5),axis.title = element_blank())+
  theme(panel.border = element_rect(colour = "black", fill=NA))+
  scale_color_viridis(option = "D")#20X5.1
saveRDS(allcell,'allcell.rds')
###########
Idents(allcell)<-allcell$Clusters
Idents(allcell)<-allcell$orig.ident
VlnPlot(allcell, features = c('Pdlim1','Pdlim2','Pdlim3', 'Pdlim4','Pdlim5','Pdlim6','Pdlim7',
                              'Cdh5','Pecam1','Vwf'),cols = my36colors,pt.size = 0.01,stack=T)

allcell$PD1 <- allcell@assays$RNA@data['Pdlim1',]
allcell$PD2 <- allcell@assays$RNA@data['Pdlim2',]
allcell$PD4 <- allcell@assays$RNA@data['Pdlim4',]
allcell$PD5 <- allcell@assays$RNA@data['Pdlim5',]
allcell$PD7 <- allcell@assays$RNA@data['Pdlim7',]

############作图###########################
#1
pdf("umap.pdf", width = 5, height = 4)
DimPlot(object = allcell , reduction = "umap", label = T,label.size = 6,repel = T,pt.size = 1)  + NoLegend()+NoAxes()#分组umap
dev.off()
#2
pdf("tsne.pdf", width = 5, height = 4)
DimPlot(object = allcell , reduction = "tsne", label = T,label.size = 6,repel = T,pt.size = 1)  + NoLegend()+NoAxes()#分组umap
dev.off()
#3
library(ggplot2)
pdlim_data <- allcell@meta.data[, c("PD1", "PD2", "PD4", "PD5", "PD7")]
pdlim_data_melted <- reshape2::melt(pdlim_data)
pdlim_data_melted$variable <- gsub("PD", "Pdlim", pdlim_data_melted$variable)
# Define custom colors for each Pdlim marker
custom_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
# Create a box plot with unique colors for each Pdlim marker and color-coded points
pdf("Box_plot1.pdf", width = 5, height = 4)
ggplot(pdlim_data_melted, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(color = "#2b2b2b", alpha = 0.7, width = 0.7) +
  geom_jitter(aes(color = variable), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Expression of Pdlim Markers",
       x = " ",
       y = "Expression Level"
       ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        plot.caption = element_text(size = 10, color = "gray"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()
pdf("Box_plot2.pdf", width = 5, height = 4)
ggplot(pdlim_data_melted, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(color = "#2b2b2b", alpha = 0.7, width = 0.7) +
  geom_point(position = position_dodge(width = 0.7), aes(color = variable), alpha = 0.5) +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Expression of Pdlim Markers",
       x = "",
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        plot.caption = element_text(size = 10, color = "gray"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
dev.off()
#4
pdf("VlnPlot.pdf", width = 6, height = 3)
Idents(allcell)=allcell$Clusters
VlnPlot(allcell, features = c('Pdlim5'),cols = my36colors,pt.size = 0)
dev.off()
#5 pie plot看比例####
dim(allcell)#849
dim(subset(allcell,subset = PD1>0))#155
dim(subset(allcell,subset = PD2>0))#173
dim(subset(allcell,subset = PD4>0))#172
dim(subset(allcell,subset = PD5>0))#219
dim(subset(allcell,subset = PD7>0))#177
library(ggplot2)
pdlim_data <- allcell@meta.data[, c("PD1", "PD2", "PD4", "PD5", "PD7")]
custom_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
create_pie_chart <- function(pd_data, color, title) {
  positive_cells <- sum(pd_data > 0)
  negative_cells <- ncol(allcell) - positive_cells
  positive_percentage <- round((positive_cells / ncol(allcell)) * 100, 1)
  
  data <- data.frame(label = c("Positive", "Negative"),
                     value = c(positive_cells, negative_cells))
  
  pie <- ggplot(data, aes(x = "", y = value, fill = label)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = c(color, "lightgray")) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(paste(title, "\n", positive_percentage, "% Positive"))
  
  return(pie)
}
pie_pd1 <- create_pie_chart(allcell$PD1, custom_colors[1], "Pdlim1")
pie_pd2 <- create_pie_chart(allcell$PD2, custom_colors[2], "Pdlim2")
pie_pd4 <- create_pie_chart(allcell$PD4, custom_colors[3], "Pdlim4")
pie_pd5 <- create_pie_chart(allcell$PD5, custom_colors[4], "Pdlim5")
pie_pd7 <- create_pie_chart(allcell$PD7, custom_colors[5], "Pdlim7")
pdf("pie_charts.pdf", width = 10, height = 8)
grid.arrange(pie_pd1, pie_pd2, pie_pd4, pie_pd5, pie_pd7, ncol = 2)
dev.off()



saveRDS(allcell,'allcell.rds')
