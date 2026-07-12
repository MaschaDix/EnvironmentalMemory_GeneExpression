# Comparing log2 fold change between treatments within WGCNA modules

# 1. set-up
library(ggplot2)
library(tidyr)
library(dplyr)
library(emmeans)

setwd()

# 2. Get data
files <- list.files(pattern = ".csv")

for (f in files) {
  name <- tools::file_path_sans_ext(f)   # strip .csv
  assign(name, read.csv(f, stringsAsFactors = FALSE)) 
}


# 3. Do the thing
# Make a loop that includes creating a linear model of the data and plotting 
# estimated marginal means, and extracting significance tests from it. 
# get the module names
files <- list.files(pattern = ".csv")
module_names <- tools::file_path_sans_ext(files)
modules <- list(MEblack, MEblue, MEbrown, MEcyan, MEgreen, MEgreenyellow,
                MEgrey, MElightcyan, MEmagenta, MEmidnightblue, MEpink, MEpurple,
                MEred, MEsalmon, MEtan, MEturquoise, MEyellow)
names(modules) <- module_names


for (n in 1:length(module_names)) {
  df <- as.data.frame(modules[[n]])
  
  colnames(df) <- c("Gene_ID", "N_H", "STP_C", "STP_H", "LTP_C", "LTP_H", "A_C", "A_H")
  
  # I need to add one column for the naive control being 0
  df$N_C <- 0
  
  # Reshape long
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = -Gene_ID,
      names_to = "Treatment",
      values_to = "logFC"
    ) %>%
    mutate(
      Group = ifelse(grepl("_C$", Treatment), "Ambient", "Heat"),
      Condition = case_when(
        grepl("^A_", Treatment)  ~ "Acclimated",
        grepl("^LTP_", Treatment) ~ "Long-term primed",
        grepl("^STP_", Treatment) ~ "Short-term primed",
        grepl("^N_", Treatment)   ~ "Naive"
      )
    )
  
  # add direction within each treatment, based on direction
  df_long <- df_long %>%
    group_by(Gene_ID, Condition) %>%
    mutate(
      delta_HA = logFC[Group == "Heat"] - logFC[Group == "Ambient"],
      Direction = case_when(
        delta_HA > 0 ~ "UP",
        delta_HA < 0 ~ "DOWN",
        TRUE ~ NA_character_
      )
    ) %>%
    tidyr::fill(Direction, .direction = "downup") %>%
    ungroup()
  
  # order factors:
  df_long$Condition <- factor(
    df_long$Condition,
    levels = c("Naive", "Short-term primed",
               "Long-term primed", "Acclimated"))
  df_long$Group <- factor(df_long$Group, levels = c("Ambient", "Heat"))
  df_long$Direction <- factor(df_long$Direction, levels = c("UP", "DOWN"))
  
  # loop stats test for UP and DOWN groups:
  # I also need to save the emm output for the ggplot, so prepare a list:
  emm_list <- list()
  for (dire in levels(df_long$Direction)) {
    
    df_sub <- subset(df_long, Direction == dire)
    
    # run lm
    model <- lm(logFC ~ Condition * Group, data = df_sub)
    
    # Export model summary
    model_summary <- summary(model)
    sink(file = paste0("Stats_output/", module_names[n], "_", dire, "_model_summary.txt"))
    print(model_summary)
    sink()
    
    # Estimated marginal means
    emm <- emmeans(model, ~ Condition * Group)
    
    # Add to the list for ggplot
    emm_list[[dire]] <- as.data.frame(emm)
    
    # Test for differences between Naive and other treatments in each group
    tr <- contrast(emm, method = "pairwise", by = "Group")
    
    # Export test results
    sink(file = paste0("Stats_output/", module_names[n], "_", dire, "_sig_diff.txt"))
    print(tr)
    sink()
  }
  
  # add direction and level it
  emm_list$UP$Direction <- "Up"
  emm_list$DOWN$Direction <- "Down"
  emm_df <- do.call(rbind, emm_list)
  emm_df$Direction <- factor(emm_df$Direction, levels = c("Up", "Down"))
  
  # Plot
  p <- ggplot(emm_df, 
              aes(x = Group, y = emmean, color = Condition, group = interaction(Condition, Direction), 
                  linetype = Direction, shape = Direction)) +
    geom_point(size = 5, position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                  linewidth = 0.6, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_bw(base_size = 14) +
    scale_color_manual(values = c(
      "Naive" = "#2B93CE",
      "Short-term primed" = "#27B78A",
      "Long-term primed" = "#DE7F00",
      "Acclimated" = "#DF90BD")) +
    labs(title = module_names[n],
         y = "Estimated marginal mean log2FC",
         x = "",
         linetype = "Direction",
         shape = "Direction",
         color = "Condition") +
    theme(axis.text = element_text(size = 14, color = "#202020"), 
          axis.title = element_text(size = 16, color = "#202020"))
  
  pdf(file = paste0("Reaction_norm_graphs/", module_names[n], "_reaction_norm.pdf"), width = 6, height = 6)
  print(p)
  dev.off()
  
  # Message which module was done
  message("Finished ", module_names[n])
}


