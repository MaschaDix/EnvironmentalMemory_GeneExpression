# Sleuth
# R version 4.3.0

# BiocManager::install("pachterlab/sleuth")
suppressMessages({
  library(sleuth)
})


setwd()



# 1. Specify file paths,load metadata, make subsets (follow Pachterlab tutorial)
sample_id<-dir(file.path("Kallisto"))

# a list of paths to the kallisto results indexed by the sample IDs is collated with
kal_dirs<-file.path("Kallisto",sample_id)

# load an auxiliary table that describes the experimental design and relationship between the kallisto directories and samples:
s2c<-read.table(file.path("Metadata","metadata_sleuth.txt"),header=TRUE,stringsAsFactors=FALSE)
s2c<-dplyr::mutate(s2c, path=kal_dirs)
print(s2c)

# because I deleted STP_H3 in the WGCNA analysis, I will also remove it here
kal_dirs <- kal_dirs[-25]
s2c <- s2c[-25,]

# order the treatments
new_order<-c("25Control","32Control","Naive","LongPrimed","LongPrimedControl","ShortPrimed","ShortPrimedControl","Acclimatized")
s2c$treatment<-factor(s2c$treatment, levels=new_order)
s2c$tanks <- factor(s2c$tanks)


# 2. Pairwise comparisons using wald test to obtain log fold change

# make pairwise subsets, and add levels otherwise the loop later won't work
# everything vs naive control:
C25vsN <- subset(s2c, s2c$treatment %in% c("Naive", "25Control"))
C25vsN$treatment <- factor(C25vsN$treatment, levels=c("25Control", "Naive"))

C25vsLP <- subset(s2c, s2c$treatment %in% c("25Control", "LongPrimed"))
C25vsLP$treatment <- factor(C25vsLP$treatment, levels=c("25Control", "LongPrimed"))

C25vsLPC <- subset(s2c, s2c$treatment %in% c("25Control", "LongPrimedControl"))
C25vsLPC$treatment <- factor(C25vsLPC$treatment, levels=c("25Control", "LongPrimedControl"))

C25vsSP <- subset(s2c, s2c$treatment %in% c("25Control", "ShortPrimed"))
C25vsSP$treatment <- factor(C25vsSP$treatment, levels=c("25Control", "ShortPrimed"))

C25vsSPC <- subset(s2c, s2c$treatment %in% c("25Control", "ShortPrimedControl"))
C25vsSPC$treatment <- factor(C25vsSPC$treatment, levels=c("25Control", "ShortPrimedControl"))

C25vsC32 <- subset(s2c, s2c$treatment %in% c("25Control", "32Control"))
C25vsC32$treatment <- factor(C25vsC32$treatment, levels=c("25Control", "32Control"))

C25vsA <- subset(s2c, s2c$treatment %in% c("25Control", "Acclimatized"))
C25vsA$treatment <- factor(C25vsA$treatment, levels=c("25Control", "Acclimatized"))

pair_list <- list(
  C25vsN = C25vsN,
  C25vsLP = C25vsLP,
  C25vsLPC = C25vsLPC,
  C25vsSP = C25vsSP,
  C25vsSPC = C25vsSPC,
  C25vsC32 = C25vsC32,
  C25vsA = C25vsA)


# run the test
for (pair_name in names(pair_list)) {
  # get the data
  pair_data <- pair_list[[pair_name]]
  
  # print out which pairwise comparison we're working on
  cat("Running Sleuth for:", pair_name, "\n")
  
  # construct sleuth object
  so_w <- sleuth_prep(pair_data, extra_bootstrap_summary=TRUE)
  
  # fit models
  so_w <- sleuth_fit(so_w, ~treatment+tanks,'full')
  so_w <- sleuth_fit(so_w, ~1, 'reduced')
  
  # run Wald test using the second level specified earlier as the beta. b values indicate up and downregulation in the beta.
  levels_treatment <- levels(pair_data$treatment)
  target_level <- paste0("treatment", levels_treatment[2])
  so_w <- sleuth_wt(so_w,which_beta= target_level)
  
  # get results and filter for significance
  wald_table <- sleuth_results(so_w, target_level, test_type = "wt")
  
  # export
  out_file <- paste0("Sleuth_", pair_name, "_Wald_test_tanks.txt")
  write.table(wald_table, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
