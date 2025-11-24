library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(EnhancedVolcano)

currentDirectory = "/Users/williamjarassier/TRAVAIL/Analyse-transcriptomique-NR/"
celpath = paste0(currentDirectory,"/DATA/RawData_E18_035_NRose/")
filescel = dir(celpath,pattern=".CEL",full.names = T)

library(oligo)
data = read.celfiles(filescel)
#eset = rma(data)

sns <- sampleNames(data)
sns <- gsub('\\.CEL$', '', sns)
sns <- gsub('\\d{2}_', '', sns)
sampleNames(data) = sns

exp_raw <- log2(Biobase::exprs(data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Sample = c("control","control","control","control","control","control",
                                "Wnt","Wnt","Wnt","Wnt","Wnt","Wnt"),
                     MW = c("40PA","40PA","40PA","12PA","12PA","12PA",
                            "40PA","40PA","40PA","12PA","12PA","12PA")) 

library(ggplot2)
ggplot(dataGG, aes(x=PC1, y=PC2,colour = Sample, shape = MW)) +
  geom_point() +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15))

oligo::boxplot(data, target = "core",main = "Boxplot of log2-intensitites for the raw data")
#library(reshape2)
#b = melt(data.frame(eset@assayData$exprs))
#ggplot(b,aes(x=variable,y=value,colour=variable)) + 
#  geom_boxplot() +
#  theme_classic() + 
#  theme(axis.text.x = element_text(angle=90)) +
#  scale_x_discrete(name="Sample") +
#  scale_y_discrete(name="Log2-intensities")

library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = data,
                    outdir = celpath,
                    force = TRUE, do.logtransform = TRUE)

RMA_calibration = function(eset){
  row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(eset)))
  RLE_data <- sweep(Biobase::exprs(eset), 1, row_medians_assayData)
  RLE_data <- as.data.frame(RLE_data)
  RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)
  print((ggplot2::ggplot(RLE_data_gathered, aes(patient_array,log2_expression_deviation)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(c(-2, 2)) + 
    theme(axis.text.x = element_text(colour = "aquamarine4", 
                                     angle = 60, size = 6.5, hjust = 1 ,
                                     face = "bold"))))
}

p1=RMA_calibration(oligo::rma(data, normalize = F))
p2=RMA_calibration(oligo::rma(data, normalize = T))
p1 + p2

############PROBLEME
#QC_PCA = function(data) {
eset = oligo::rma(data, normalize = T)
exprs <- Biobase::exprs(eset)
PCA <- prcomp(t(exprs), scale = FALSE)
ercentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Sample = c("control","control","control","control","control","control",
                                "Wnt","Wnt","Wnt","Wnt","Wnt","Wnt"),
                     MW = c("40PA","40PA","40PA","12PA","12PA","12PA",
                            "40PA","40PA","40PA","12PA","12PA","12PA")) 
ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Sample, colour = MW)) +
  ggtitle("PCA plot of the calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
}
#p1=QC_PCA(data)
#p2=QC_PCA(data)
#grid.arrange(  p1,  p2,  nrow = 2)

eset_norm = oligo::rma(data, normalize = T)
expr <- Biobase::exprs(eset_norm)
annotation_for_heatmap <- 
  data.frame(Sample = c("control","control","control","control","control","control",
                        "Wnt","Wnt","Wnt","Wnt","Wnt","Wnt"),
             MW = c("40PA","40PA","40PA","12PA","12PA","12PA",
                    "40PA","40PA","40PA","12PA","12PA","12PA"))
row.names(annotation_for_heatmap) <- row.names(pData(eset_norm))
dists <- as.matrix(dist(t(expr), method = "manhattan"))
rownames(dists) <- row.names(pData(eset_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))
colnames(dists) <- NULL
diag(dists) <- NA
ann_colors <- list(
  Sample = c(control = "chartreuse4", Wnt = "burlywood3"),
  MW = c('40PA' = "blue4", '12PA' = "cadetblue2")
)
library(pheatmap)
pheatmap(dists, col = (hmcol), 
         annotation_row = annotation_for_heatmap,
         annotation_colors = ann_colors,
         legend = TRUE, 
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm = TRUE), 
                           max(dists, na.rm = TRUE)), 
         legend_labels = (c("small distance", "large distance")),
         main = "Clustering heatmap for the calibrated samples")


data_medians <- rowMedians(Biobase::exprs(eset_norm))
hist_res <- hist(data_medians, 100, col = "cornsilk1", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "antiquewhite4",
                 xlab = "Median intensities")
abline(v = 2.5, col = "coral4", lwd = 2)

idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){sum(x > 2.4) >= samples_cutoff})
table(idx_man_threshold)


pData(eset_norm)$Sample = c("control","control","control","control","control","control",
                            "Wnt","Wnt","Wnt","Wnt","Wnt","Wnt")
pData(eset_norm)$MW = c("40PA","40PA","40PA","12PA","12PA","12PA",
                        "40PA","40PA","40PA","12PA","12PA","12PA")
no_of_samples <- table(paste0(pData(eset_norm)$Sample, "_", pData(eset_norm)$MW))
no_of_samples 
samples_cutoff <- min(no_of_samples)
idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){
                             sum(x > 2.5) >= samples_cutoff})
table(idx_man_threshold)
data_filtered <- subset(eset_norm, idx_man_threshold)

BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)
head(ls("package:hugene10sttranscriptcluster.db"))
head( keys(hugene10sttranscriptcluster.db) )
BiocManager::install("hgu95av2.db")
data_annotation <- AnnotationDbi::select(hgu95av2.db,
                                       keys = (featureNames(eset_norm)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")
anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))

