# Plasticity correlation with physiological variables

library(tidyverse)
library(ggpubr) # to include stats in ggplots
library(car)
library(lme4)
library(lmerTest)
library(emmeans)

setwd()
bwnet2 <- readRDS("GE_WGCNA_results.RDS")

# Get Module Eigengenes 
module_eigengenes <- bwnet2$MEs
head(module_eigengenes)


# to calculate plasticity, I will average the control value of each pre-treatment, then calculate the distance to each heat sample.
Control_eigengenes <- matrix(NA, nrow = 4, ncol = 16)
rownames(Control_eigengenes) <- c("N", "STP", "LTP", "A")
colnames(Control_eigengenes) <- colnames(module_eigengenes)[-17]

for(i in 1:16) {
  # Naive
  Control_eigengenes[1,i] <- mean(module_eigengenes[20:22,i])
  
  # STP
  Control_eigengenes[2,i] <- mean(module_eigengenes[1:4,i])
  
  # LTP
  Control_eigengenes[3,i] <- mean(module_eigengenes[14:16,i])
  
  # A
  Control_eigengenes[4,i] <- mean(module_eigengenes[17:19,i])
}



Plasticities <- matrix(NA, nrow = 12, ncol = 16)
rownames(Plasticities) <- c("N_H1", "N_H2","N_H3",
                            "STP_H1", "STP_H2","STP_H3",
                            "LTP_H1", "LTP_H2", "LTP_H3", 
                            "A_H1", "A_H2", "A_H3")
colnames(Plasticities) <- colnames(module_eigengenes)[-17]


for(i in 1:16) {
  # Naive
  Plasticities[1:3,i] <- module_eigengenes[8:10, i] - Control_eigengenes[1,i]
  
  # STP
  Plasticities[4:6,i] <- module_eigengenes[23:25, i] - Control_eigengenes[2,i]
  
  # LTP
  Plasticities[7:9,i] <- module_eigengenes[5:7, i] - Control_eigengenes[3,i]
  
  # A
  Plasticities[10:12,i] <- module_eigengenes[11:13, i] - Control_eigengenes[4,i]
}



# Now I need to get the physiological variables. 
list.files()
meta <- read.csv("Metadata.csv", header = TRUE)

Plastic_df <- as.data.frame(Plasticities)
Plastic_df$Treatment <- c("N", "N", "N",
                          "STP", "STP", "STP",
                          "LTP", "LTP", "LTP",
                          "A", "A", "A")

Plastic_df$Change_PAM <- c(meta$Change_PAM[8:10], meta$Change_PAM[c(23,24,26)], 
                           meta$Change_PAM[5:7], meta$Change_PAM[11:13])
Plastic_df$Change_SD <- c(meta$Change_SD[8:10], meta$Change_SD[c(23,24,26)], 
                           meta$Change_SD[5:7], meta$Change_SD[11:13])
Plastic_df$Mortality <- c(meta$Mortality[8:10], meta$Mortality[c(23,24,26)], 
                          meta$Mortality[5:7], meta$Mortality[11:13])


# mean plasticities per sample
summary(t(Plastic_df[,-c(17:20)]))


# convert to long data frame
Plastic_long <- Plastic_df %>% 
  rownames_to_column("Sample") %>%
  pivot_longer(cols = starts_with("ME"),
               names_to = "Module",
               values_to = "Plasticity")


# make a barplot of absolute mean and sd of each treatment
mean_plast <- data.frame(Treatment = factor(c("N", "STP", "LTP", "A"), levels = c("N", "STP", "LTP", "A")),
                         Mean = c(mean(abs(Plast_N$Plasticity)), mean(abs(Plast_STP$Plasticity)), mean(abs(Plast_LTP$Plasticity)), mean(abs(Plast_A$Plasticity))),
                         SD = c(sd(abs(Plast_N$Plasticity)), sd(abs(Plast_STP$Plasticity)), sd(abs(Plast_LTP$Plasticity)), sd(abs(Plast_A$Plasticity))))

pdf("Plasticities.pdf", height = 5, width = 4.5)
ggplot(mean_plast, aes(x= Treatment, y = Mean, fill = Treatment)) +
  geom_col() +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, color = "#404040") +
  theme_bw(base_size = 14) +
  labs(x = "Treatment",
       y = "Transcriptomic Plasticity") +
  scale_fill_manual(values = c(
    "N" = "#2B93CE",
    "STP" = "#27B78A",
    "LTP" = "#DE7F00",
    "A" = "#DF90BD")) +
  theme(axis.text = element_text(size = 18, color = "#202020"),
        axis.title = element_text(size = 20, color = "#202020"))
dev.off()



# use a linear mixed model to test for significant differences in plasticity between treatments.
# I'm using a LMM because my data are hierarchical and not independent due to the module structure.
# Alternatively, I could create averages across all modules https://stats.oarc.ucla.edu/other/mult-pkg/introduction-to-linear-mixed-models/
# fixed effects are variables that do not vary, so treatment. 
# modules are assumed random variables
# https://medium.com/@TingyuZou/understanding-and-interpreting-linear-mixed-models-702bf1b38500

Plastic_long$Treatment <- factor(Plastic_long$Treatment, levels = c("N", "STP", "LTP", "A"))

Plastic_long_abs <- Plastic_long
Plastic_long_abs$Plasticity <- abs(Plastic_long_abs$Plasticity)

m_plast <- lmer(
  Plasticity ~ Treatment + (1 | Module),
  data = Plastic_long_abs
)

# check assumptions
# normally distributed residuals:
res <- resid(m_plast)
qqnorm(res)
qqline(res) # Q-Q plot looks ok, indicating some skewness but that's expected because I forced absolute values

# homoscedasticity
plot(m_plast, which = 1)

anova(m_plast)
summary(m_plast)
confint(m_plast)

pb <- bootMer(
  m_plast,
  FUN = function(x) fixef(x),
  nsim = 1000
)

apply(pb$t, 2, quantile, probs = c(0.025, 0.975))

# post-hoc using emm
emm <- emmeans(m_plast, ~Treatment)
pairs(emm, adjust = "tukey")


# confidence intervals should not overlap zero, t-value > 2 or < -2 indicates significance


# there is one LTP treatment that sometimes appears as an outlier, will run the analysis with and without outlier
Plast_long_Outlier_rem <- Plastic_df[-8,] %>% 
  rownames_to_column("Sample") %>%
  pivot_longer(cols = starts_with("ME"),
               names_to = "Module",
               values_to = "Plasticity")

# now test the correlation by module:
module_cor <- Plastic_long %>%
  group_by(Module) %>%
  summarise(
    cor_PAM = cor(Plasticity, Change_PAM, method = "spearman"),
    p_PAM   = cor.test(Plasticity, Change_PAM, method = "spearman")$p.value,
    cor_SD  = cor(Plasticity, Change_SD, method = "spearman"),
    p_SD    = cor.test(Plasticity, Change_SD, method = "spearman")$p.value
  )

# with outlier removed
module_cor_outlier_rem <- Plast_long_Outlier_rem %>%
  group_by(Module) %>%
  summarise(
    cor_PAM = cor(Plasticity, Change_PAM, method = "spearman"),
    p_PAM   = cor.test(Plasticity, Change_PAM, method = "spearman")$p.value,
    cor_SD  = cor(Plasticity, Change_SD, method = "spearman"),
    p_SD    = cor.test(Plasticity, Change_SD, method = "spearman")$p.value
  )

# multiple test correction:
module_cor <- module_cor %>%
  mutate(
    p_PAM_adj = p.adjust(p_PAM, method = "BH"),
    p_SD_adj  = p.adjust(p_SD,  method = "BH"))

# outlier removed
module_cor_outlier_rem <- module_cor_outlier_rem %>%
  mutate(
    p_PAM_adj = p.adjust(p_PAM, method = "BH"),
    p_SD_adj  = p.adjust(p_SD,  method = "BH"))



# Only the correlation of MEblack with SD is signficant using the more robust spearman correlation and p-value adjustment.
# make a nice plot:
Plastic_long$Treatment <- factor(
  Plastic_long$Treatment,
  levels = c("N", "STP",
             "LTP", "A"))

Plastic_black <- Plastic_long %>%
  filter(Module == "MEblack")

pdf("MEblack_corr.pdf", height = 5, width = 4.8)
ggplot(Plastic_black, aes(x = Plasticity, y = Change_SD, fill = Treatment)) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "#404040", linewidth = 0.5) + # no dashed line because actually significant
  geom_point(shape = 21, size = 5, color = "#404040", alpha = 0.9) +
  theme_bw(base_size = 14) + 
  labs(y = "Change in symbiont density (%)",
       x = "Transcriptomic plasticity") +
  scale_fill_manual(values = c(
    "N" = "#2B93CE",
    "STP" = "#27B78A",
    "LTP" = "#DE7F00",
    "A" = "#DF90BD")) +
  theme(axis.text = element_text(size = 14, color = "#202020"), 
        axis.title = element_text(size = 16, color = "#202020"))
dev.off()


# also make a plot showing the spearman's rhos and p-values for every module
module_cor_plot <- data.frame(module_cor_outlier_rem)[,c(1,2,4,6,7)]
# log transform p-values
module_cor_plot$p_PAM_adj <- -log10(module_cor_plot$p_PAM_adj)
module_cor_plot$p_SD_adj <- -log10(module_cor_plot$p_SD_adj)

# PAM
pdf("PAM_correlations_all.pdf", height = 5, width = 6)
ggplot(module_cor_plot, aes(x = Module, y = cor_PAM, color = p_PAM_adj)) +
  geom_point(size = 5) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 14) + 
  labs(x = "Module",
       y = "Spearman's rho") +
  scale_color_viridis_c(option = "plasma", limits = c(0, 1.32)) +
  theme(axis.text = element_text(size = 14, color = "#202020"), 
        axis.title = element_text(size = 16, color = "#202020"))
dev.off()

# SD
pdf("SD_correlations_all.pdf", height = 5, width = 6)
ggplot(module_cor_plot, aes(x = Module, y = cor_SD, color = p_SD_adj)) +
  geom_point(size = 5) +
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(-1, 1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_bw(base_size = 14) + 
  labs(x = "Module",
       y = "Spearman's rho") +
  scale_color_viridis_c(option = "plasma", limits = c(0, 1.32)) +
  theme(axis.text = element_text(size = 14, color = "#202020"), 
        axis.title = element_text(size = 16, color = "#202020"))
dev.off()


# Option 2, p-values on y-axis and modules as fill (preferred the other ones)
# PAM
ggplot(module_cor_plot, aes(x = cor_PAM, y = p_PAM_adj, color = Module)) +
  geom_point(size = 5) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(
    "MEblack" = "black",
    "MEblue" = "blue",
    "MEbrown" = "brown",
    "MEcyan" = "cyan",
    "MEgreen" = "green",
    "MEgreenyellow" = "greenyellow",
    "MElightcyan" = "lightcyan",
    "MEmidnightblue" = "midnightblue",
    "MEpink" = "pink",
    "MEpurple" = "purple",
    "MEred" = "red",
    "MEsalmon" = "salmon",
    "MEtan" = "tan",
    "MEturquoise" = "turquoise",
    "MEyellow" = "yellow")) +
  theme_bw(base_size = 14) + 
  labs(x = "log10 of adj. p-value",
       x = "Spearman's rho") +
  theme(axis.text = element_text(size = 14, color = "#202020"), 
        axis.title = element_text(size = 16, color = "#202020"))

# SD
ggplot(module_cor_plot, aes(x = cor_SD, y = p_SD_adj, color = Module)) +
  geom_point(size = 5) +
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(
    "MEblack" = "black",
    "MEblue" = "blue",
    "MEbrown" = "brown",
    "MEcyan" = "cyan",
    "MEgreen" = "green",
    "MEgreenyellow" = "greenyellow",
    "MElightcyan" = "lightcyan",
    "MEmidnightblue" = "midnightblue",
    "MEpink" = "pink",
    "MEpurple" = "purple",
    "MEred" = "red",
    "MEsalmon" = "salmon",
    "MEtan" = "tan",
    "MEturquoise" = "turquoise",
    "MEyellow" = "yellow")) +
  theme_bw(base_size = 14) + 
  labs(y = "log10 of adj. p-value",
       x = "Spearman's rho") +
  theme(axis.text = element_text(size = 14, color = "#202020"), 
        axis.title = element_text(size = 16, color = "#202020"))

# My mortality data are not appropriate for correlation tests like pearson or spearman. 
# Are changes in plasticity associated with occurrence of mortality?
# make mortality binary:
Plastic_long$Mortality <- as.integer(Plastic_long$Mortality > 0)

glm_by_module <- Plastic_long %>%
  group_by(Module) %>%
  do({
    m <- glm(Mortality ~ Plasticity,
             family = binomial,
             data = .)
    
    tidy(m)
  }) %>%
  ungroup()

