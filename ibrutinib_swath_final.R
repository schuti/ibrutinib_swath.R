################################################################################################################################
# ibrutinib_swath.R: A custom R script for automatic data preprocessing, analysis and visualization of ibrutinib_swath dataset #
# Author: Somchai Chutipongtanate                                                                                              #
# Last update: April 9, 2019                                                                                                   #
# Source: https://github.com/schuti/ibrutinib_swath.R                                                                          #
################################################################################################################################

# Load: R packages
library(readxl)
library(dplyr)
library(tidyr)
library(biomaRt)
library(preprocessCore)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(FactoMineR)
library(pheatmap)

# Load: ibrutinib-SWATH dataset (PXD013402)
setwd("~/Desktop")
data_path <- "~/Desktop/ibrutinib_SWATH.xlsx"

# Start: Data preprocess -----------------------------------------------------------------------------
## loading
group <- as.factor(c("W", "W", "W", "iW", "iW", "iW", "Q", "Q", "Q", "iQ", "iQ", "iQ"))
sample_label <- as.character(c("W1", "W2", "W3", "iW1", "iW2", "iW3", "Q1", "Q2", "Q3", "iQ1", "iQ2", "iQ3"))
areaPept <- read_excel(data_path, sheet = "Area - peptides")
areaProt <- read_excel(data_path, sheet = "Area - proteins")

## gene mapping using biomaRt package (ref#1)
df <- areaProt[ ,1] %>% 
  tidyr::separate(Protein, c("sp", "uniProtID", "entry_name"), sep = "\\|") %>%
  tidyr::separate(entry_name, c("entry_names", "species"), sep = "_") %>%
  dplyr::select(uniProtID, entry_names, species) 
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
tmp <- getBM(attributes = c('uniprotswissprot', 'external_gene_name'), 
             filters = 'uniprotswissprot', 
             values = df$uniProtID, 
             mart = ensembl) 
colnames(tmp) <- c('uniProtID', "gene.SYMBOL")
df <- left_join(df, tmp[!duplicated(tmp$uniProtID), ], by = "uniProtID") 
ind <- is.na(df$gene.SYMBOL)
df$gene.SYMBOL[ind] <- df$entry_names[ind]
id_all <- df

## Quantile normalization using preprocessCore package (ref#2)
expr_raw <- areaProt[ , 2:length(areaProt)]
colnames(expr_raw) <- sample_label
Quantile <- as.data.frame(normalize.quantiles(log2(as.matrix(expr_raw))))
colnames(Quantile) <- sample_label

## Missing values replaced by zero
ind <- which(is.na(Quantile), arr.ind = TRUE)
Quantile[ind] <- 0
expr_processed <- Quantile

## Collect datasets
raw_ds <- cbind(id_all, expr_raw)
process_ds <- cbind(id_all, expr_processed) 
df <- t(expr_processed)
colnames(df) <- id_all$gene.SYMBOL
log_ds <- data.frame(group, df)

# End: Data preprocess -----------------------------------------------------------------------------

# Start: Data analysis and visualization -----------------------------------------------------------------------------
## Group average
tmp <- data.frame(group = log_ds[ , 1], 2^log_ds[ , 2:length(log_ds)]) %>% 
  gather(gene.SYMBOL, expression, -group) %>%
  dplyr::group_by(group, gene.SYMBOL) %>% 
  dplyr::summarize(group_mean = mean(expression)) %>%
  spread(gene.SYMBOL, group_mean)
gr_avr <- as.data.frame(tmp[ , 2:length(tmp)])
rownames(gr_avr) <- tmp$group
gr_pair <- combn(unique(tmp$group), 2) 	
fc <- (gr_avr[gr_pair[1, ], ] / gr_avr[gr_pair[2, ], ]) %>% log2()  
rownames(fc) <- paste0('log2', '(', gr_pair[1, ], '/', gr_pair[2, ], ')')
log2fc_ds <- fc

## Group SD
tmp <-  data.frame(group = log_ds[ , 1], 2^log_ds[ , 2:length(log_ds)]) %>% 
  gather(gene.SYMBOL, expression, -group) %>%
  dplyr::group_by(group, gene.SYMBOL) %>% 
  dplyr::summarize(group_sd = sd(expression)) %>%
  spread(gene.SYMBOL, group_sd)
gr_sd <- as.data.frame(tmp[ , 2:length(tmp)])
rownames(gr_sd) <- tmp$group

## Coefficient of variation
qc <- 100 *gr_sd/gr_avr 
qc <- data.frame(group = tmp$group, qc)
QC <- qc %>% gather(gene, CV, -group)

# Calculate median-CV of each group
medianCV <- QC %>% dplyr::group_by(group) %>% summarise(CV = round(median(CV), 1))

# Violin plot of inter-group CV 
plot.qc <- ggplot(QC, aes(x=group, y=CV)) + 
              geom_violin(aes(fill = group), trim=FALSE, width = 0.8, 
                          na.rm = TRUE, position = "dodge")+
              geom_boxplot(width=0.1, fill = 'white', outlier.size = 0, 
                          na.rm = TRUE, position = "dodge")+
              geom_text(data = medianCV, aes(label = CV), position = position_dodge(width = 1), 
                          hjust = -0.5, vjust = -0.5, size = 5) +
              xlab("") + ylab("% Coefficient of Variation") +
              scale_y_continuous(breaks=c(0, 10, 20, 50, ceiling(max(QC$CV, na.rm=TRUE)))) +
              theme_light(base_size = 12)
plot.qc
print(paste0("Median-CV: iQ, ", medianCV[1,2], "%; iW, ", medianCV[2,2], "%; Q, ", medianCV[3,2], "%; W, ", medianCV[4,2], "%"))
pdf("QC_violinPlot.pdf", width = 6, height = 4)
print(plot.qc)
dev.off()

## Correlation heatmap
corr <- 2^expr_processed
corr <- round(cor(corr, method = "pearson"),3)
corr[lower.tri(corr)] <- NA
melted_corr <- melt(corr, na.rm = TRUE)
plot_corrHM <- ggplot(data = melted_corr, aes(x = Var2, y = Var1, fill = value))+  
                  geom_tile(color = "white")+
                  scale_fill_gradient2(low = "white", high = "red", mid = "yellow",   
                                       midpoint = 0.9, limit = c(0.8, 1), space = "Lab", 
                                       name= paste("Pearson", "\ncorrelation") ) +
                  labs(x = "", y = "") +
                  theme_minimal() +
                  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                                   size = 10, hjust = 1)) +
                  coord_fixed() + 
                  geom_text(aes(label = value), color = "black", size = 2) +  
                  theme(axis.text.y = element_text(color = "black", size=10),
                        panel.grid.major = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.ticks = element_blank(),
                        legend.justification = c(1, 0),
                        legend.position = c(0.5, 0.7),
                        legend.direction = "horizontal")+
                  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

plot_corrHM
pdf("plot_corrHM.pdf", width = 6, height = 4)
print(plot_corrHM)
dev.off()

## nPP plot
n_pept_prot <- areaPept %>% 
  dplyr::group_by(Protein) %>% 
  dplyr::summarize(n_pept = n()) %>% 
  arrange(desc(n_pept))
nPP <- data.frame(n_pept = c("1", "2-5", "6-10"), 
                    n_prot = rbind(n_pept_prot %>% filter(n_pept ==1) %>% nrow(),
                                   n_pept_prot %>% filter(n_pept >=2 & n_pept <= 5) %>% nrow(), 
                                   n_pept_prot %>% filter(n_pept >=6) %>% nrow()))
nPP_plot <- ggplot(nPP, aes(x = n_pept, y= n_prot)) + 
                    geom_bar(stat = "identity", fill = "steelblue") + 
                    ylim(0, max(nPP$n_prot)+50) +
                    geom_text(aes(label= n_prot), vjust=-0.3, color="black", size=4.5) +
                    geom_text(aes(label= paste0(round(100*n_prot/sum(n_prot), 1), "%")), vjust=1.6, color="white", size=4.5) +
                    xlab("Numbers of peptide") + ylab("Numbers of protein") +
                    theme_light(base_size = 12)
nPP_plot
pdf("number_pept_prot.pdf", width = 4, height = 3)
print(nPP_plot)
dev.off()

## PCA individual plot using FactorMineR package (ref#3)
fit_pca <- PCA(log_ds[ , 2:length(log_ds)], graph = FALSE, scale.unit = TRUE)
percentage <- fit_pca$eig[ , 2]
PCs <- data.frame(fit_pca$ind$coord)
PCs$group <- group  
plotPCA <- ggplot(data = PCs, aes(x = Dim.1, y = Dim.2)) +
              geom_point(aes(colour = group), size = 3) +
              xlab(paste0('PC1', ' ', '(', round(percentage[1], 2), '%)')) + 
              ylab(paste0('PC2', ' ', '(', round(percentage[2], 2), '%)')) +
              scale_fill_hue(l=40) + 
              coord_fixed(ratio=1, xlim=range(PCs$Dim.1), ylim=range(PCs$Dim.2)) +
              geom_text_repel(label = rownames(PCs)) +
              theme_light(base_size = 15)
plotPCA
pdf("plotPCA.pdf", width = 6, height = 4)
print(plotPCA)
dev.off()

## Protein abundance heatmap by pheatmap package (ref#4)
qc_hm <- 2^expr_processed
rownames(qc_hm) <- process_ds$gene.SYMBOL 
for(i in seq_along(qc_hm)){
  if(qc_hm[i] != 0){
    qc_hm[i] <- log10(qc_hm[i])
  } else {
    qc_hm[i] <- 0
  }}    
n_missing <- sum(qc_hm == 0)
n_total <- dim(qc_hm)[1] * dim(qc_hm)[2]

qc_hm_plot <- pheatmap(t(qc_hm),  
                         breaks = seq(0, max(qc_hm), length.out=101), 
                         legend_breaks = seq(0, round(max(qc_hm), 0), length.out=8), 
                         legend_labels = c("1e+00", "1e+01", "1e+02", "1e+03", "1e+04", "1e+05", "1e+06", "1e+07"),
                         color = colorRampPalette(c("black", "#8ea1ff", "#14ff57", "yellow", 	"orange", "#ea4444"))(100), 
                         border_color = "gray",
                         clustering_distance_cols = "maximum", 
                         clustering_method_columns = "complete",
                         cluster_rows = FALSE, 
                         fontsize_row = 8, fontsize_col = 1, 
                         scale = "none",
                         main = paste0("Total ", n_total, " data points; ", n_missing, " missing values (", round(100*n_missing/n_total, 2), "%) showed in black)") )
qc_hm_plot
pdf("qc_heatmap.pdf", width = 8, height = 4)
print(qc_hm_plot)
dev.off()

## ANOVA with Tukey's post-hoc
tmp <- as.matrix(log_ds[, 2:length(log_ds)])
fit.aov <- aov(tmp ~ group) 
output.aov <- summary.aov(fit.aov) 
anova.pVal <- numeric(length = ncol(tmp))
for (i in 1:length(output.aov)){
  anova.pVal[i] <- output.aov[[i]][1, 5]
}
adj.pVal <- matrix(nrow = ncol(tmp), ncol = nrow(log2fc_ds))
colnames(adj.pVal) <-  paste(gr_pair[1, ], " vs ", gr_pair[2, ])
rownames(adj.pVal) <- colnames(tmp)
for (i in 1:ncol(tmp)){ 
  adj.pVal[i, ] <- (TukeyHSD((aov(tmp[, i] ~ group))))[[1]][ ,4] 
}
anova_ds <- cbind(anova.pVal, adj.pVal)

## Pairwise-Volcano plot
tmp <- data.frame(gene = rownames(anova_ds), anova_ds)   
colnames(tmp) <- c("gene", "anova.pVal", paste0(gr_pair[1, ], "/", gr_pair[2, ]) )
long_ano <- gather(tmp, compare, adj_pVal, -gene, -anova.pVal)     
fc.vp <- t(log2fc_ds)
fc.vp <- data.frame(gene = colnames(log2fc_ds), fc.vp)
colnames(fc.vp) <- c("gene", paste0(gr_pair[1, ], "/", gr_pair[2, ]) )
long_fc <- gather(fc.vp, compare, log2FC, -gene)
long_ano.fc <- long_ano %>% 
  left_join(long_fc, by = c("gene", "compare"))
long_ano.fc$gene <- as.character(long_ano.fc$gene)
volcano_all <- ggplot(data = long_ano.fc, aes(x= log2FC, y=-log10(adj_pVal))) +
                    geom_point(aes(color = as.factor(abs(log2FC) >= log2(1.5) &	anova.pVal < 0.05 & adj_pVal < 0.05)), size = 3, alpha = 0.5, show.legend = FALSE) +
                    scale_color_manual(values = c("grey", "red")) +
                    xlab("log2 (fold change)") + ylab("-log10 (adjusted p-value)") +
                    ggtitle(label = paste0("Volcano plot at ", 1.5, 
                                           "x fold change and adjusted P-value < ", 0.05)) + 
                    theme_grey(base_size = 15) +
                    facet_wrap(~ compare)
volcano_all
pdf("valcano_all.pdf", width = 12, height = 6)
print(volcano_all)
dev.off()

## Significant protein heatmap by pheatmap (ref#4)
tmp <- as.matrix(log_ds[ , 2:length(log_ds)])
med <- apply(t(tmp), 1, mean)
medScale <- (t(tmp) - med)
tmp <- anova_ds[, 1]
medScale <- data.frame(medScale, 
                      anova_pVal = tmp, 
                      gene = rownames(medScale))
medScale_sig <- medScale %>% filter(anova_pVal < 0.05)
rownames(medScale_sig) <- medScale_sig$gene
medScale_sig <- medScale_sig[, 1: (length(medScale_sig) - 2)]
nprot_sig <- nrow(medScale_sig)
hm_sig <- pheatmap(medScale_sig, silent = FALSE,
                     breaks = seq(-(max(round(medScale_sig, 0))), max(round(medScale_sig, 0)), length.out=101), 
                     legend_breaks = seq(-(max(round(medScale_sig, 0))), max(round(medScale_sig, 0)), length.out=5),
                     color = colorRampPalette(c("darkblue", "blue", "white", "orangered", "red"))(100),  
                     border_color = NA,
                     annotation_col = data.frame(group = factor(group), 
                                                 row.names = sample_label),
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation", 
                     clustering_method = "average", 
                     fontsize_row = 2, fontsize_col = 10, 
                     scale = "none",
                     main = paste0(nprot_sig, " significant proteins (ANOVA p-value < ", 0.05, ")", "\nScale: Log2(fold change) with mean center", "\nClustering: Correlation distance and average linkage"))
hm_sig
pdf("heatmap_sig.pdf", width = 6, height = 9)
print(hm_sig)
dev.off()

# End: Data analysis and visualization -----------------------------------------------------------------------------

# References
## 1. Durinck S, Spellman P, Birney E, Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, 1184–1191.
## 2. Bolstad B (2018). preprocessCore: A collection of pre-processing functions. R package version 1.44.0, 
## 3. Lê, S., Josse, J. & Husson, F. (2008). FactoMineR: An R Package for Multivariate Analysis. Journal of Statistical Software. 25(1). pp. 1-18.
## 4. Raivo Kolde (2018). pheatmap: Pretty Heatmaps. R package version 1.0.10. 


