# WGCNA on gene expression data

# set up
setwd()

library(WGCNA)
library(DESeq2)
library(tidyverse)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)


# 1. load data (row header is gene name, column header is sample name)
data <- read.csv("GE_WGCNA_input.csv", header=TRUE)
rownames(data)<-data$Gene_ID
data <- data[,-1]
data <- round(data, digits=0)
head(data)

# 2. load metadata (Pre-treatment, treatment, sympiont density, PAM, Mortality, etc)
phenodata <- read.csv("Metadata.csv", header=TRUE)
head(phenodata)
phenodata$Treatment <- as.factor(phenodata$Treatment)
phenodata$Acute_heat <- as.factor(phenodata$Acute_heat)
phenodata$Heat_history <- as.factor(phenodata$Heat_history)
phenodata$History_type <- as.factor(phenodata$History_type)


# 3. Outlier detection
# this function requires the data to be transposed first
gsg <- goodSamplesGenes(t(data))
gsg$allOK # if it says FALSE, genes and/or samples have been detected as outliers. It says FALSE.
# to detect outliers:
table(gsg$goodGenes) # 1059 outliers
table(gsg$goodSamples) # 0 outliers
# to remove outliers
data <- data[gsg$goodGenes==TRUE,] 


# another way to check for outliers is with a PCA plot. Do this before and after removing outliers like above.
pca <-prcomp(t(data))
pca_data <- pca$x
pca_var <- pca$sdev^2
pca_var_percent <- round(pca_var/sum(pca_var)*100,digits=2)
# before
pdf("PCA_OutlierCheck.pdf")
ggplot(pca_data,aes(PC1, PC2))+geom_point()+geom_text(label=rownames(pca_data))+labs(x=paste0('PC1: ', pca_var_percent[1], ' %'), y=paste0('PC2: ',pca_var_percent[2], ' %'))
dev.off()
# after
pdf("PCA_OutliersRemoved.pdf")
ggplot(pca_data,aes(PC1, PC2))+geom_point()+geom_text(label=rownames(pca_data))+labs(x=paste0('PC1: ', pca_var_percent[1], ' %'), y=paste0('PC2: ',pca_var_percent[2], ' %'))
dev.off()
# if you want to make the graph cleaner:
pdf("PCA_nodots.pdf")
ggplot(pca_data,aes(PC1, PC2))+geom_text(label=rownames(pca_data))+labs(x=paste0('PC1: ', pca_var_percent[1], ' %'), y=paste0('PC2: ',pca_var_percent[2], ' %'))
dev.off()

# based on the PCA plots I decided to remove STP_H3. I will reload the data, remove STP_H3 and run
# everything up to here again. Lines used to remove STP_H3:
data<-data[,-25]
phenodata <- phenodata[-25,]
# Make a new PCA plot
pdf("PCA_STP_H3_Removed.pdf")
ggplot(pca_data,aes(PC1, PC2))+geom_text(label=rownames(pca_data))+labs(x=paste0('PC1: ', pca_var_percent[1], ' %'), y=paste0('PC2: ',pca_var_percent[2], ' %'))
dev.off()


# 4. Normalization
# exclude outliers if there are any.
# make a list of the sample names:
sample <- c('STP_C1', 'STP_C2', 'STP_C3', 'STP_C4', 'LTP_H1', 'LTP_H2', 'LTP_H3', 'N_H1', 'N_H2', 'N_H3', 'A_H1', 'A_H2', 'A_H3', 'LTP_C1', 'LTP_C2', 'LTP_C3', 'A_C1', 'A_C2', 'A_C3', 'N_C1', 'N_C2', 'N_C3', 'STP_H1', 'STP_H2', 'STP_H4')
ColData<-data.frame(sample)

# create dds, not specifying a model (design =~1)
dds<-DESeqDataSetFromMatrix(countData=data, colData=ColData,design=~1)

# To remove genes with overall low expression, but not exclude those highly expressed in specific treatments, I set a total row sum cut-off of 100
dds2 <- dds[rowSums(counts(dds)) >= 100, ]
nrow(dds2) # 13203

# perform variance stabilization
dds_norm<-vst(dds2)

# get normalized counts
norm.counts <- assay(dds_norm) %>% t()


# 5. Network Construction
# Choose the soft-thresholding power
powers = c(c(1:10), seq(from = 12, to=50, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                  powerVector = powers,
                  networkType = "signed", # this means we also want to consider the direction of expression, not just if it's sig different
                  verbose = 5)
sft.data <- sft$fitIndices

# visualizations to pick powers
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
 

pdf("WGCNA_Powers.pdf")
grid.arrange(a1, a2, nrow = 2)
dev.off()

# choose a soft threshold that is above the red line in the R^2 value, and retains the highest mean connectivity. I will choose 9. 
soft_power <- 9


# 5. Network construction
# for the blockwise function, data needs to be numeric:
norm.counts[]<-sapply(norm.counts, as.numeric)

# assign correlation function to WGCNA correlation function to avoid namespace error
temp_cor <- cor
cor <- WGCNA::cor

# blockwiseModules function provides a relatively fast and computationally inexpensive clustering
# these are mostly default settings, can try to change them if needed
bwnet2<-blockwiseModules(norm.counts,
				checkMissingData=TRUE,
				maxBlockSize= 60000, # this specifies how many genes we can max have per cluster (default 5000, depends on memory of your workstation. 32 GB can handle up to 30000 genes. I am choosing 60000 because it is recommended to go for as high as your machine can handle)
				blockSizePenaltyPower=Inf, # I really don't want only one block so I'm setting a high penalty
				loadTOM=FALSE, # I'm asking it to calculate TOM here
				corType="pearson",
				TOMType="signed", #determines how to deal with genes on the edge between blocks
				power=soft_power,
				mergeCutHeight=0.25,
				numericLabels=FALSE,
				randomSeed=1234,
				verbose=3) 
cor<-temp_cor


# save WGCNA results to a file
readr::write_rds(bwnet2, file=file.path("GE_WGCNA_results.RDS"))

# make a dendrogram:
moduleColors <- bwnet2$colors
geneTree <- bwnet2$dendrograms[[1]]

pdf("Module_Dendrogram.pdf")
plotDendroAndColors(geneTree, moduleColors,"Module",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

# 6. Module Eigengenes 
module_eigengenes <- bwnet2$MEs

# print out a preview
head(module_eigengenes)

# make a dendrogram of how close the modules are:
ME.dissimilarity = 1-cor(module_eigengenes, use="complete") # calculate eigengene dissimilarity
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 

pdf("MEs_Dendrogram.pdf")
plot(METree)
abline(h=.25, col = "red") # a height of .25 corresponds to correlation of .75, see mergeCutHeight in the blockwiseModules function
dev.off()


# get number of genes for each module
genes_per_module<-table(bwnet2$colors)
write.csv(genes_per_module, "GE_WGCNA_genes_per_module.csv", row.names=FALSE)

# What genes are a part of which module?
gene_module_key <- tibble::enframe(bwnet2$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))
# here I can ask which genes are part of this specific module
gene_module_key %>%
  dplyr::filter(module == "MEturquoise")

readr::write_tsv(gene_module_key,
  file = file.path("GE_WGCNA_gene_to_module.tsv")
)

# make lists for each module:
# Split by module
gene_lists <- split(gene_module_key, gene_module_key$module)

# Write one .csv per module
for (mod in names(gene_lists)) {
  filename <- paste0("GE_WGCNA_", mod, "_genes.csv")
  readr::write_csv(gene_lists[[mod]], file = filename)
}


# get an idea if our eigengenes relate to our metadata
# first check that samples are still in order:
all.equal(phenodata$Sample_ID, rownames(module_eigengenes)) # should say TRUE

#Create a design matrix from the different variables
des_mat_t<-model.matrix(~phenodata$Treatment) # 8 different treatments
des_mat_ah<-model.matrix(~phenodata$Acute_heat) # Yes or no
des_mat_hh<-model.matrix(~phenodata$Heat_history) # Yes or no
des_mat_ht<-model.matrix(~phenodata$History_type) # 4 pre-treatments


# Run linear model on each module. Limma wants our tests to be per row, so we need to transpose
fit<-limma::lmFit(t(module_eigengenes),design=des_mat_t)
# apply empirical Bayes to smooth standard errors:
fit<-limma::eBayes(fit)
# obtain multiple testing correction and obtain stats
stats_df_t<-limma::topTable(fit, number=ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
head(stats_df_t) #most significantly differential modules across the specified pheno observation are at the top
write.csv(stats_df_t, "GE_WGCNA_modules_by_treatment.csv", row.names=FALSE)

# repeat for all variables
fit<-limma::lmFit(t(module_eigengenes),design=des_mat_ah)
fit<-limma::eBayes(fit)
stats_df_ah<-limma::topTable(fit, number=ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
head(stats_df_ah)
write.csv(stats_df_ah, "GE_WGCNA_modules_by_acute_heat.csv", row.names=FALSE)

fit<-limma::lmFit(t(module_eigengenes),design=des_mat_hh)
fit<-limma::eBayes(fit)
stats_df_hh<-limma::topTable(fit, number=ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
head(stats_df_hh)
write.csv(stats_df_hh, "GE_WGCNA_modules_by_heat_history.csv", row.names=FALSE)

fit<-limma::lmFit(t(module_eigengenes),design=des_mat_ht)
fit<-limma::eBayes(fit)
stats_df_ht<-limma::topTable(fit, number=ncol(module_eigengenes)) %>% tibble::rownames_to_column("module")
head(stats_df_ht)
write.csv(stats_df_ht, "GE_WGCNA_modules_by_history_type.csv", row.names=FALSE)


# can do this on PC by loading csv files
# Make a heatmap by joining them all into one matrix:
df_t  <- stats_df_t  %>% select(module, adj.P.Val) %>% rename(Treatment = adj.P.Val)
df_ah <- stats_df_ah %>% select(module, adj.P.Val) %>% rename(AcuteHeat = adj.P.Val)
df_hh <- stats_df_hh %>% select(module, adj.P.Val) %>% rename(HeatHistory = adj.P.Val)
df_ht <- stats_df_ht %>% select(module, adj.P.Val) %>% rename(HistoryType = adj.P.Val)

combined_pvals <- df_t %>%
  full_join(df_ah, by = "module") %>%
  full_join(df_hh, by = "module") %>%
  full_join(df_ht, by = "module")

# remove MEgrey
combined_pvals <- combined_pvals[-14,]

# log transform p-values and turn non-sig ones into NAs
pvals_long <- pvals_long %>%
  mutate(neglog10_p = -log10(adj.P.Val),
         neglog10_p = ifelse(adj.P.Val > 0.05, NA, neglog10_p))

# order modules
combined_pvals <- combined_pvals %>%
  mutate(module = factor(module, levels = c("MEgreenyellow", "MEsalmon", "MEblue", "MEturquoise", "MEpink",
             "MEgreen", "MEblack", "MEpurple", "MEbrown", "MEmidnightblue", 
             "MEyellow", "MEtan", "MEmagenta", "MEred", "MEcyan", "MElightcyan"))) %>%
           arrange(module)

# convert to matrix:
pval_mat <- combined_pvals %>%
  column_to_rownames("module") %>%
  as.matrix()

# convert to -log10 and grey out any non-significant ones
neglog10_p <- -log10(pval_mat)
neglog10_p[pval_mat > 0.05] <- NA



textMatrix =  paste(signif(pval_mat, 2))
dim(textMatrix) = dim(pval_mat)

# Create a color palette: white to red
my_palette <- colorRampPalette(c("white", "red"))(n = 50)

pdf("ExpVariables_Correlations.pdf", width = 6, height = 8)
par(mar = c(5, 10, 2, 2))  # bottom, left, top, right
labeledHeatmap(Matrix = neglog10_p,
               xLabels = colnames(neglog10_p),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = my_palette,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               main = paste("Module-variable relationships"))
dev.off()




# make plots for all modules by all variables
# make heatmap functions for each variable (slightly different because I need to specify different colours/number of colours depending on how many levels each variable has), also see here for more explanation: for explanation:https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#413_Make_a_custom_heatmap_function

# treatment (8 levels):
make_module_heatmap_tr <-function(module_name,
						enrichment_mat=norm.counts,
						metadata_df=phenodata,
						gene_module_key_df=gene_module_key,
						module_eigengenes_df=module_eigengenes) {
	module_eigengene<-module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("Sample_ID")
	col_annot_df<-metadata_df %>% dplyr::select(Sample_ID, Treatment) %>% dplyr::inner_join(module_eigengene, by="Sample_ID") %>% dplyr::arrange(Treatment, Sample_ID)%>%tibble::column_to_rownames("Sample_ID")
	col_annot<-ComplexHeatmap::HeatmapAnnotation(Treatment=col_annot_df$Treatment,
		module_eigengene=ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df,module_name)),
		col=list(Treatment=c("STP_Control" = "#f1a340", "LTP_Heat" = "#998ec3", "N_Heat" = "#b8c38e", "A_Heat" = "#8eb8c3", "LTP_Control" = "#c3998e", "A_Control" = "#c38e9e", "N_Control" = "#8ec3b3", "STP_Heat" = "#c3b38e")))
	module_genes<-gene_module_key_df %>% dplyr::filter(module==module_name) %>% dplyr::pull(gene)
	mod_mat<-enrichment_mat %>% t() %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% module_genes) %>% dplyr::select(rownames(col_annot_df)) %>% as.matrix()
	mod_mat<-mod_mat %>% t() %>% scale() %>% t()
	color_func <- circlize::colorRamp2(
		c(-2, 0, 2),
		c("#67a9cf", "#f7f7f7", "#ef8a62"))
	heatmap<-ComplexHeatmap::Heatmap(mod_mat, 
		name=module_name, 
		col=color_func, 
		bottom_annotation=col_annot, 
		cluster_columns=FALSE,
		show_row_names=TRUE,
		show_column_names=TRUE)
	return(heatmap)
}

# Acute heat and heat history (2 levels):
make_module_heatmap_ah <-function(module_name,
						enrichment_mat=norm.counts,
						metadata_df=phenodata,
						gene_module_key_df=gene_module_key,
						module_eigengenes_df=module_eigengenes) {
	module_eigengene<-module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("Sample_ID")
	col_annot_df<-metadata_df %>% dplyr::select(Sample_ID, Acute_heat) %>% dplyr::inner_join(module_eigengene, by="Sample_ID") %>% dplyr::arrange(Acute_heat, Sample_ID)%>%tibble::column_to_rownames("Sample_ID")
	col_annot<-ComplexHeatmap::HeatmapAnnotation(Acute_heat=col_annot_df$Acute_heat,
		module_eigengene=ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df,module_name)),
		col=list(Acute_heat=c("Control" = "#f1a340", "Heat" = "#b8c38e")))
	module_genes<-gene_module_key_df %>% dplyr::filter(module==module_name) %>% dplyr::pull(gene)
	mod_mat<-enrichment_mat %>% t() %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% module_genes) %>% dplyr::select(rownames(col_annot_df)) %>% as.matrix()
	mod_mat<-mod_mat %>% t() %>% scale() %>% t()
	color_func <- circlize::colorRamp2(
		c(-2, 0, 2),
		c("#67a9cf", "#f7f7f7", "#ef8a62"))
	heatmap<-ComplexHeatmap::Heatmap(mod_mat, 
		name=module_name, 
		col=color_func, 
		bottom_annotation=col_annot, 
		cluster_columns=FALSE,
		show_row_names=TRUE,
		show_column_names=TRUE)
	return(heatmap)
}

make_module_heatmap_hh <-function(module_name,
						enrichment_mat=norm.counts,
						metadata_df=phenodata,
						gene_module_key_df=gene_module_key,
						module_eigengenes_df=module_eigengenes) {
	module_eigengene<-module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("Sample_ID")
	col_annot_df<-metadata_df %>% dplyr::select(Sample_ID, Heat_history) %>% dplyr::inner_join(module_eigengene, by="Sample_ID") %>% dplyr::arrange(Heat_history, Sample_ID)%>%tibble::column_to_rownames("Sample_ID")
	col_annot<-ComplexHeatmap::HeatmapAnnotation(Heat_history=col_annot_df$Heat_history,
		module_eigengene=ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df,module_name)),
		col=list(Heat_history=c("TRUE" = "#f1a340", "FALSE" = "#b8c38e")))
	module_genes<-gene_module_key_df %>% dplyr::filter(module==module_name) %>% dplyr::pull(gene)
	mod_mat<-enrichment_mat %>% t() %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% module_genes) %>% dplyr::select(rownames(col_annot_df)) %>% as.matrix()
	mod_mat<-mod_mat %>% t() %>% scale() %>% t()
	color_func <- circlize::colorRamp2(
		c(-2, 0, 2),
		c("#67a9cf", "#f7f7f7", "#ef8a62"))
	heatmap<-ComplexHeatmap::Heatmap(mod_mat, 
		name=module_name, 
		col=color_func, 
		bottom_annotation=col_annot, 
		cluster_columns=FALSE,
		show_row_names=TRUE,
		show_column_names=TRUE)
	return(heatmap)
}

# History type (4 colours):
make_module_heatmap_ht <-function(module_name,
						enrichment_mat=norm.counts,
						metadata_df=phenodata,
						gene_module_key_df=gene_module_key,
						module_eigengenes_df=module_eigengenes) {
	module_eigengene<-module_eigengenes_df %>% dplyr::select(all_of(module_name)) %>% tibble::rownames_to_column("Sample_ID")
	col_annot_df<-metadata_df %>% dplyr::select(Sample_ID, History_type) %>% dplyr::inner_join(module_eigengene, by="Sample_ID") %>% dplyr::arrange(History_type, Sample_ID)%>%tibble::column_to_rownames("Sample_ID")
	col_annot<-ComplexHeatmap::HeatmapAnnotation(History_type=col_annot_df$History_type,
		module_eigengene=ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df,module_name)),
		col=list(History_type=c("ST_Primed" = "#8eb8c3", "LT_Primed" = "#c3998e", "Naive" = "#c38e9e", "Acclimated" = "#8ec3b3")))
	module_genes<-gene_module_key_df %>% dplyr::filter(module==module_name) %>% dplyr::pull(gene)
	mod_mat<-enrichment_mat %>% t() %>% as.data.frame() %>% dplyr::filter(rownames(.) %in% module_genes) %>% dplyr::select(rownames(col_annot_df)) %>% as.matrix()
	mod_mat<-mod_mat %>% t() %>% scale() %>% t()
	color_func <- circlize::colorRamp2(
		c(-2, 0, 2),
		c("#67a9cf", "#f7f7f7", "#ef8a62"))
	heatmap<-ComplexHeatmap::Heatmap(mod_mat, 
		name=module_name, 
		col=color_func, 
		bottom_annotation=col_annot, 
		cluster_columns=FALSE,
		show_row_names=TRUE,
		show_column_names=TRUE)
	return(heatmap)
}



# make loops for each variable
# get a list of all modules:
MEs <- colnames(module_eigengenes)

# Treatment:
box_data <-module_eigengenes %>%
	tibble::rownames_to_column("Sample_ID") %>%
	dplyr::inner_join(phenodata %>% dplyr::select(Sample_ID, Treatment),
  by = c("Sample_ID" = "Sample_ID"))
# loop for all MEs:
for (ME in MEs) {
	# first the boxplot
  filename <- paste0(ME, "_Box_Treatment.pdf")
  pdf(filename)
  print(
  	ggplot(box_data, aes(x = Treatment, y = .data[[ME]], color = Treatment)) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic() +
    labs(y = ME, title = paste(ME, "by Treatment"))
    )
  dev.off()
  # then the heatmap:
  module_heatmap <- make_module_heatmap_tr(module_name=ME)
  filename <- paste0(ME, "_Heatmap_Treatment.pdf")
  pdf(filename)
  print(module_heatmap)
  dev.off()
}



# Acute heat:
box_data <-module_eigengenes %>%
	tibble::rownames_to_column("Sample_ID") %>%
	dplyr::inner_join(phenodata %>% dplyr::select(Sample_ID, Acute_heat),
  by = c("Sample_ID" = "Sample_ID"))
# loop for all MEs:
for (ME in MEs) {
	# first the boxplot
  filename <- paste0(ME, "_Box_AcuteHeat.pdf")
  pdf(filename)
  print(
  	ggplot(box_data, aes(x = Acute_heat, y = .data[[ME]], color = Acute_heat)) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic() +
    labs(y = ME, title = paste(ME, "by Acute Heat"))
    )
  dev.off()
  # then the heatmap:
  module_heatmap <- make_module_heatmap_ah(module_name=ME)
  filename <- paste0(ME, "_Heatmap_AcuteHeat.pdf")
  pdf(filename)
  print(module_heatmap)
  dev.off()
}



# Heat history:
box_data <-module_eigengenes %>%
	tibble::rownames_to_column("Sample_ID") %>%
	dplyr::inner_join(phenodata %>% dplyr::select(Sample_ID, Heat_history),
  by = c("Sample_ID" = "Sample_ID"))
# loop for all MEs:
for (ME in MEs) {
	# first the boxplot
  filename <- paste0(ME, "_Box_HeatHistory.pdf")
  pdf(filename)
  print(
  	ggplot(box_data, aes(x = Heat_history, y = .data[[ME]], color = Heat_history)) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic() +
    labs(y = ME, title = paste(ME, "by Heat History"))
    )
  dev.off()
  # then the heatmap:
  module_heatmap <- make_module_heatmap_hh(module_name=ME)
  filename <- paste0(ME, "_Heatmap_HeatHistory.pdf")
  pdf(filename)
  print(module_heatmap)
  dev.off()
}



# Heat type:
box_data <-module_eigengenes %>%
	tibble::rownames_to_column("Sample_ID") %>%
	dplyr::inner_join(phenodata %>% dplyr::select(Sample_ID, History_type),
  by = c("Sample_ID" = "Sample_ID"))
# loop for all MEs:
for (ME in MEs) {
	# first the boxplot
  filename <- paste0(ME, "_Box_HistoryType.pdf")
  pdf(filename)
  print(
  	ggplot(box_data, aes(x = History_type, y = .data[[ME]], color = History_type)) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    ggforce::geom_sina(maxwidth = 0.3) +
    theme_classic() +
    labs(y = ME, title = paste(ME, "by History Type"))
    )
  dev.off()
  # then the heatmap:
  module_heatmap <- make_module_heatmap_hh(module_name=ME)
  filename <- paste0(ME, "_Heatmap_HistoryType.pdf")
  pdf(filename)
  print(module_heatmap)
  dev.off()
}





# 7. Associate modules with phenodata (instructions from https://rpubs.com/natmurad/WGCNA)
#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

#Recalculating MEs with label colors
MEs0 = moduleEigengenes(norm.counts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs <- MEs[,-17] # remove MEgrey

# prepare the phenotypic data:
# remove non-numerical variables
pheno <- phenodata[,-c(1:6)]
rownames(pheno) <- sample
pheno$Mortality <- as.integer(pheno$Mortality > 0)


moduleTraitCor = cor(MEs, pheno, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
p_adj <- apply(moduleTraitPvalue, 2, p.adjust, method = "BH")


#sizeGrWindow(8,4)

#Displaying correlations and its p-values

textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(p_adj, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 20, 3, 3))

#Displaying the correlation values in a heatmap plot
pdf("TraitCorrelations_new.pdf", width = 6, height = 8)
par(mar = c(5, 10, 2, 2))  # bottom, left, top, right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(pheno),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



# I want to find out about intramodular connectivity of nodes
# https://github.com/kpatel427/YouTubeTutorials/blob/main/WGCNA.R

# module membership is calcualted as the correlation of the eigengene and the gene expression profile
module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')

nSamples <- nrow(norm.counts)
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]


# this gives a matrix of p-values for every gene in every module. 
# Next I want to extract the correlation and adjusted p-values of the genes sorted by modules. 
names(gene_lists)

path <- "./"

# List all module gene files (skip the overall mapping files)
module_files <- list.files(path, pattern = "GE_WGCNA_ME.*_genes\\.csv", full.names = TRUE)

# Loop through each module file
for (file in module_files) {
  
  # Extract module name, e.g. "MEgreen" from "GE_WGCNA_MEgreen_genes.csv"
  module_name <- sub("GE_WGCNA_(ME\\w+)_genes\\.csv", "\\1", basename(file))
  
  # Read gene list (assuming a single column of gene names)
  genes <- read.csv(file, header = TRUE)[,1]
  
  # Match these genes to your correlation and p-value matrices
  corr_values <- module.membership.measure[module_name, genes, drop = TRUE]
  pvals <- module.membership.measure.pvals[module_name, genes, drop = TRUE]
  
  # Adjust p-values for multiple testing (FDR)
  pvals_adj <- p.adjust(pvals, method = "fdr")
  
  # Create dataframe
  df <- data.frame(
    Gene = genes,
    Correlation = corr_values,
    P.Value = pvals,
    Adj.P.Value = pvals_adj
  )
  
  # Define output file name
  out_file <- file.path(path, paste0(module_name, "_hubgenes.csv"))
  
  # Save as CSV
  write.csv(df, out_file, row.names = FALSE)
  
  cat("Saved:", out_file, "\n")
}




# Finally I would like to calculate the gene significance and associated p-values in correlation with trait data
# e.g., which genes are significantly correlated with acclimated samples
# for this, I first need to binarize my trait information:

# for acute heat and heat history this is easy because they are already basically binary:
traits <- phenodata %>% mutate(acute_heat_bin = ifelse(grepl('Heat', Acute_heat), 1, 0)) %>% select(9)
traits$heat_history_bin <- phenodata %>% mutate(heat_history_bin = ifelse(grepl('TRUE', Heat_history), 1, 0)) %>% select(9)

# for treatment and history type I have to make separate columns where each category is 1 and all others are 0
levels(phenodata$Treatment)
phenodata$Treatment <- factor(phenodata$Treatment, levels = c("N_Control", "N_Heat", "STP_Control", "STP_Heat", "LTP_Control", "LTP_Heat", "A_Control", "A_Heat"))

treatment.out <- binarizeCategoricalColumns(phenodata$Treatment,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)


levels(phenodata$History_type)
phenodata$History_type <- factor(phenodata$History_type, levels = c("Naive", "ST_primed", "LT_primed", "Acclimated"))
history_type.out <- binarizeCategoricalColumns(phenodata$History_type,
                                            includePairwise = FALSE,
                                            includeLevelVsAll = TRUE,
                                            minCount = 1)

traits <- cbind(traits, treatment.out, history_type.out)


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.ST_primed.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


# make another heatmap with correlation values and without MEgrey:
module_eigengenes <- bwnet2$MEs
moduleColors <- bwnet2$colors
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)
MEs0 = moduleEigengenes(norm.counts, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEs <- MEs[,-17]
pheno <- phenodata[,-c(1:6)]
rownames(pheno) <- sample
traits <- cbind(traits, pheno) # only include binary categories (acute heat and heat history)
moduleTraitCor = cor(MEs, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
cor_text <- signif(moduleTraitCor, 2)
p_text <- ifelse(moduleTraitPvalue < 0.05,
                 paste0("\n(", signif(moduleTraitPvalue, 2), ")"),
                 "") # only put p-value if it's significant
textMatrix <- paste0(cor_text, p_text)
dim(textMatrix) <- dim(moduleTraitCor)
pdf("TraitCorrelations_new.pdf", width = 6, height = 8)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



# Intramodular connectivity
soft_power = 10 # as above
adjacency <- adjacency(norm.counts, power = soft_power)

# Compute TOM (topological overlap matrix) if not done
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Now calculate intramodular connectivity:
IMConn <- intramodularConnectivity(adjacency, moduleColors)

# add module membership to it:
IMConn$moduleColor <- moduleColors

# hub genes per module:
hubGenes <- by(IMConn, IMConn$moduleColor, function(df) {
  df[order(-df$kWithin), ]
})

# save as .csv files:
for (moduleColor in names(hubGenes)) {
  filename <- paste0("ME", moduleColor, "_intraCon.csv")
  write.csv(hubGenes[[moduleColor]], file = filename, row.names = TRUE)
}






