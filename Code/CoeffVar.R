# Coefficient of variance between replicates

# load packages
library(tidyverse)

# 1. Load TPM data
# load kallisto abundance estimates
setwd("/Users/dixmf/Documents/Projects/EnvironmentalMemory/DataAnalysis/RNAData/Old_genome")

# make sample IDs
sample_id<-dir(file.path("Kallisto"))
sample_id

# a list of paths to the kallisto results indexed by the sample IDs is collated with
kal_dirs<-file.path("Kallisto",sample_id, "abundance.tsv")
kal_dirs

# read in counts for all genes
g <- read.delim(kal_dirs[1], header = TRUE, row.names = 1)
genes <- g[,4]
genes <- as.data.frame(genes)
rownames(genes) <- rownames(g)
colnames(genes)[1] <- sample_id[1]

for(i in 2:length(kal_dirs)){
  g <- read.delim(kal_dirs[i], header = TRUE, row.names = 1)
  genes[,i] <- g[,4]
  colnames(genes)[i] <- sample_id[i]
}

head(genes)


# will delete same STP_H sample as in WGCNA:
genes <- genes[-25]


# rename the columns to my naming convention:
colnames(genes) <- c("STP_C_1", "STP_C_2", "STP_C_3", "STP_C_4",
                     "LTP_H_1", "LTP_H_2", "LTP_H_3",
                     "N_H_1", "N_H_2", "N_H_3",
                     "A_H_1", "A_H_2", "A_H_3",
                     "LTP_C_1", "LTP_C_2", "LTP_C_3",
                     "A_C_1", "A_C_2", "A_C_3",
                     "N_C_1", "N_C_2", "N_C_3",
                     "STP_H_1", "STP_H_2", "STP_H_3")


# I don't want low-expression genes to skew the analysis, so will remove them:
keep <- rowSums(genes) >= 100
genes_filtered <- genes[keep, ] # this still keeps 13179 genes


# 2. Calculate Coefficients of Variance:
# make a long data frame
genes_long <- genes_filtered %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "TPM"
  )


# add treatments
genes_long <- genes_long %>%
  mutate(
    Treatment = sub("_[0-9]+$", "", sample)
  )


# calculate coefficient of variance for each gene in each treatment as SD/mean:
cv_gene <- genes_long %>%
  group_by(Treatment, gene) %>%
  summarise(
    mean_TPM = mean(TPM, na.rm = TRUE),
    sd_TPM = sd(TPM, na.rm = TRUE),
    CV = sd_TPM / mean_TPM,
    .groups = "drop"
  )
head(cv_gene)


# 3. Stats
# Need to restructure the tibble so I can perform statistics:
cv_friedman <- cv_gene %>%
  select(gene, Treatment, CV) %>%
  mutate(
    Treatment = factor(
      Treatment,
      levels = c(
        "N_C", "N_H",
        "STP_C", "STP_H",
        "LTP_C", "LTP_H",
        "A_C", "A_H"
      )
    )
  ) 
cv_complete <- cv_friedman %>%
  group_by(gene) %>%
  filter(all(is.finite(CV))) %>%
  ungroup()
friedman.test(CV ~ Treatment | gene,
              data = cv_complete)
pairwise.wilcox.test(cv_complete$CV, cv_complete$Treatment, p.adjust.method="BH")


# 4. Make a bar plot
# get mean CV values per treatment:
mean_cv <- cv_gene %>%
  group_by(Treatment) %>%
  summarise(mean_CV = mean(CV, na.rm = TRUE),
            sd_CV = sd(CV, na.rm = TRUE)
  )
mean_cv <- as.data.frame(mean_cv)
mean_cv$Treatment <- factor(c("A_C", "A_H", "LTP_C", "LTP_H", "N_C", "N_H", "STP_C", "STP_H"), levels = c("N_C", "N_H", "STP_C", "STP_H", "LTP_C", "LTP_H", "A_C", "A_H"))

pdf("CVs_means.pdf", height = 4, width = 5)
ggplot(mean_cv, aes(x= Treatment, y = mean_CV, fill = Treatment)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean_CV - sd_CV, ymax = mean_CV + sd_CV),
                width = 0.2, color = "#404040") +
  scale_y_continuous(limits = c(0, 0.9)) +
  theme_bw(base_size = 14) +
  labs(x = "Treatment",
       y = "Coefficient of Variation") +
  scale_fill_manual(values = c(
    "A_C" = "#CC79A7",
    "A_H" = "#F2A7D3",
    "LTP_C" = "#D55E00",
    "LTP_H" = "#E69F00",
    "N_C" = "#0072B2",
    "N_H" = "#56B4E9",
    "STP_C" = "#009E73",
    "STP_H" = "#4DCFA0")) +
  theme(axis.text = element_text(size = 14, color = "#202020"),
        axis.title = element_text(size = 16, color = "#202020"))
dev.off()


# 5. Correlation tests
# get metadata
meta <- read.csv("Metadata.csv", header = TRUE) # either delete the outlier now or exclude it later

# make one df with Treatment, Module, mean change PAM and SD, and mean CV in that module
mean_CV_module <- map_dfr(
  CV_subsets,
  ~ data.frame(
    Treatment = colnames(.x),
    mean_CV = colMeans(.x, na.rm = TRUE)
  ),
  .id = "Module"
)

CV_df <- data.frame(
  Treatment = mean_CV_module$Treatment,
  Module = mean_CV_module$Module,
  CV = mean_CV_module$mean_CV,
  Change_PAM = rep(c(mean(meta$Change_PAM[17:19]), mean(meta$Change_PAM[11:13]), # Acclimated
                     mean(meta$Change_PAM[14:16]), mean(meta$Change_PAM[5:7]), # LTP
                     mean(meta$Change_PAM[20:22]), mean(meta$Change_PAM[8:10]), # Naive
                     mean(meta$Change_PAM[c(1,2,3,4)]), mean(meta$Change_PAM[c(23,24,26)])), 17), # STP
  Change_SD = rep(c(mean(meta$Change_SD[17:19]), mean(meta$Change_SD[11:13]), # Acclimated
                    mean(meta$Change_SD[14:16]), mean(meta$Change_SD[5:7]), # LTP
                    mean(meta$Change_SD[20:22]), mean(meta$Change_SD[8:10]), # Naive
                    mean(meta$Change_SD[c(1,2,3,4)]), mean(meta$Change_SD[c(23,24,26)])), 17) # STP
)

module_cor <- CV_df %>%
  group_by(Module) %>%
  summarise(
    cor_PAM = cor(CV, Change_PAM, method = "spearman"),
    p_PAM   = cor.test(CV, Change_PAM, method = "spearman")$p.value,
    cor_SD  = cor(CV, Change_SD, method = "spearman"),
    p_SD    = cor.test(CV, Change_SD, method = "spearman")$p.value
  )


# now also add correlation with plasticity, but for this I need to calculate change in CV from control to heat:
CV_change <- CV_df %>%
  mutate(
    History = sub("_[CH]$", "", Treatment),
    Condition = sub("^.*_", "", Treatment)
  ) %>%
  select(Module, History, Condition, CV) %>%
  pivot_wider(
    names_from = Condition,
    values_from = CV
  ) %>%
  mutate(
    CV_change = H - C
  ) %>%
  select(Module, History, CV_change)

Plastic <- Plastic_df[,-c(17:20)] # created this in Plasticity.R
CV_change <- CV_change[-c(25:28),] # remove MEgrey

# make a tibble of all the plasticity values per treatment in long format
Plastic_summary <- Plastic %>%
  rownames_to_column("Sample") %>%
  mutate(
    Treatment = sub("[0-9]+$", "", Sample)
  ) %>%
  pivot_longer(
    cols = -c(Sample, Treatment),
    names_to = "Module",
    values_to = "Plasticity"
  ) %>%
  group_by(Treatment, Module) %>%
  summarise(
    Plasticity = mean(Plasticity, na.rm = TRUE),
    .groups = "drop"
  )

# make the module names consistent 
CV_change <- CV_change %>%
  mutate(
    Module = sub("^GE_WGCNA_(.*)_genes$", "\\1", Module)
  )

# need to rename the treatments here to just N, STP, LTP and A
Plastic_history <- Plastic_summary %>%
  filter(grepl("_H$", Treatment)) %>%
  mutate(
    History = sub("_H$", "", Treatment)
  ) %>%
  select(History, Module, Plasticity)
# join the two data frames so I have one with CV_change and Plasticity
CV_change <- CV_change %>%
  left_join(
    Plastic_history,
    by = c(
      "History",
      "Module"
    )
  )


# make one data frame with all correlation tests
module_cor_CV_Plast <- CV_change %>%
  group_by(Module) %>%
  summarise(
    cor = cor(CV_change, Plasticity, method = "spearman"),
    p_value   = cor.test(CV_change, Plasticity, method = "spearman")$p.value
  )

module_cor <- module_cor[-7,] # remove MEgrey
module_cor$Module <- module_cor_CV_Plast$Module # shorter names
module_cor <- cbind(module_cor, module_cor_CV_Plast) # combine data frames
module_cor <- module_cor[,-6] # remove extra module name column

# adjust p-values:
module_cor <- module_cor %>%
  mutate(
    p_PAM_adj = p.adjust(p_PAM, method = "BH"),
    p_SD_adj  = p.adjust(p_SD,  method = "BH"),
    p_Plast_adj  = p.adjust(p_value,  method = "BH"))


# also make a plot showing the spearman's rhos and p-values for every module and every comparison
# take correlations and adjusted p-values only
module_cor_plot <- data.frame(module_cor)[,c(1,2,4,6,8,9, 10)]

# log transform p-values
module_cor_plot$p_PAM_adj <- -log10(module_cor_plot$p_PAM_adj)
module_cor_plot$p_SD_adj <- -log10(module_cor_plot$p_SD_adj)
module_cor_plot$p_Plast_adj <- -log10(module_cor_plot$p_Plast_adj)

# turn into a long data frame 
cor_long <- module_cor_plot %>%
  # rename into nicer terms for the graph
  transmute(
    Module,
    "Fv/Fm" = cor_PAM,
    "Symbiont Density" = cor_SD,
    "Transcriptomic Plasticity" = cor,
    p_PAM = p_PAM_adj,
    p_SymDens = p_SD_adj,
    p_Plast = p_Plast_adj
  ) %>%
  # make a long data frame out of it
  pivot_longer(
    cols = c("Fv/Fm", "Symbiont Density", "Transcriptomic Plasticity"),
    names_to = "Variable",
    values_to = "rho"
  ) %>%
  # add the p_adj to each row, assigned to the correct variable
  mutate(
    p_adj = case_when(
      Variable == "Fv/Fm" ~ p_PAM,
      Variable == "Symbiont Density" ~ p_SymDens,
      Variable == "Transcriptomic Plasticity" ~ p_Plast
    )
  )


pdf("CV_correlations_all.pdf", height = 4, width = 6)
ggplot(cor_long, aes(x = Module, y = rho, color = p_adj, shape = Variable)) +
  geom_point(size = 4) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 14) + 
  labs(x = "Module",
       y = "Spearman's rho",
       color = "-log10 adj. p-value") +
  scale_color_viridis_c(option = "plasma", limits = c(0, 0.4)) +
  scale_shape_manual(
    values = c(
      "Fv/Fm" = 16,
      "Symbiont Density" = 17,
      "Transcriptomic Plasticity" = 15
    )
  ) +
  theme(axis.text = element_text(size = 10, color = "#202020"), 
        axis.title = element_text(size = 12, color = "#202020"))
dev.off()
