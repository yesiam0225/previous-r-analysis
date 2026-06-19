# ============================================================
# ANCOVA (Analysis of Covariance) — Age-Adjusted Comparisons
#
# This script tests group and time effects on gait outcomes
# after adjusting for participant age. Age control is essential
# in my work because I am comparing children to adults, and
# absolute gait parameters scale strongly with body size and
# developmental stage.
#
# Two datasets and two ANCOVA designs are used:
#
#   1. Two-way ANCOVA on the joint kinematic / spatiotemporal
#      dataset (`spt`), with Age as covariate and
#      Group * Time as the factorial design
#
#   2. One-way ANCOVA on a separate growth dataset (`grw`)
#      with Age as covariate and Group as the only factor
#
# Post-hoc analysis uses emmeans_test with Bonferroni correction
# for pairwise comparisons within Group or within Time.
#
# Author : Yeon-Joo Kang
# Inst.  : Georgia State University
# Period : 2021-2024
# Status : Archived dissertation analysis
# ============================================================

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

# Load joint kinematic / spatiotemporal dataset
joint <- read.csv("jointGC.csv", header = TRUE)

# Set categorical factors. The working dataframe in the original
# file is `spt`, which holds the joint kinematic columns side by
# side with spatiotemporal columns.
spt$Group <- as.factor(spt$Group)
spt$Time  <- as.factor(spt$Time)
spt$Group <- factor(spt$Group, levels = c("Training", "Control"))


# ============================================================
# Part 1: Two-way ANCOVAs
# Group * Time, adjusting for Age
# ============================================================

# ------------------------------------------------------------
# Knee flexion percentage of gait cycle
# ------------------------------------------------------------
res.aov1 <- joint %>%
  anova_test(KneeFlex. ~ Age + Group * Time)
get_anova_table(res.aov1)

# ------------------------------------------------------------
# Stride length
# ------------------------------------------------------------
res.aov2 <- spt %>%
  anova_test(Stride_Length ~ Age + Group * Time)
get_anova_table(res.aov2)

# ------------------------------------------------------------
# Step time
# ------------------------------------------------------------
res.aov3 <- spt %>%
  anova_test(Step_time ~ Age + Group * Time)
get_anova_table(res.aov3)

# ------------------------------------------------------------
# Step length
# ------------------------------------------------------------
res.aov4 <- spt %>%
  anova_test(Step_length ~ Age + Group * Time)
get_anova_table(res.aov4)

# ------------------------------------------------------------
# Step width
# ------------------------------------------------------------
res.aov5 <- spt %>%
  anova_test(Step_width ~ Age + Group * Time)
get_anova_table(res.aov5)

# ------------------------------------------------------------
# Swing time
# ------------------------------------------------------------
res.aov6 <- spt %>%
  anova_test(Swing ~ Age + Group * Time)
get_anova_table(res.aov6)

# ------------------------------------------------------------
# Single support
# ------------------------------------------------------------
res.aov7 <- spt %>%
  anova_test(Single ~ Age + Group * Time)
get_anova_table(res.aov7)

# ------------------------------------------------------------
# Post-hoc: time effect within each group, age-adjusted
# ------------------------------------------------------------

# Stride length: effect of Time within each Group
spt %>%
  group_by(Group) %>%
  anova_test(Stride_Length ~ Age + Time)

# Stride length: effect of Group within each Time
spt %>%
  group_by(Time) %>%
  anova_test(Stride_Length ~ Age + Group)

# ------------------------------------------------------------
# Post-hoc: pairwise comparisons with Bonferroni correction
# ------------------------------------------------------------

# Single support: Time differences within Training group
pwc <- spt %>%
  group_by(Group) %>%
  emmeans_test(Single ~ Time, covariate = Age,
               p.adjust.method = "bonferroni")
pwc %>% filter(Group == "Training")

# Single support: Time differences within Control group
pwc <- spt %>%
  group_by(Group) %>%
  emmeans_test(Single ~ Time, covariate = Age,
               p.adjust.method = "bonferroni")
pwc %>% filter(Group == "Control")

# Single support: Group differences at Time 1
pwc <- spt %>%
  group_by(Time) %>%
  emmeans_test(Single ~ Group, covariate = Age,
               p.adjust.method = "bonferroni")
pwc %>% filter(Time == "1")

# Single support: Group differences at Time 2
pwc <- spt %>%
  group_by(Time) %>%
  emmeans_test(Single ~ Group, covariate = Age,
               p.adjust.method = "bonferroni")
pwc %>% filter(Time == "2")


# ============================================================
# Part 2: One-way ANCOVAs on growth dataset
# Group only, adjusting for Age
# ============================================================

grw <- read.csv("growth.csv", header = TRUE)
grw$Group <- as.factor(grw$Group)

# ------------------------------------------------------------
# Model diagnostics for one representative outcome (step width)
# ------------------------------------------------------------

# Fit the linear model with the covariate first
model <- lm(Step_width ~ Age + Group, data = grw)

# Diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted, -.se.fit)
head(model.metrics, 3)

# Identify potential outliers
model.metrics %>%
  filter(abs(.std.resid) > 3) %>%
  as.data.frame()

# ------------------------------------------------------------
# One-way ANCOVAs across spatiotemporal outcomes
# ------------------------------------------------------------

# Stride time
res.aov8 <- grw %>% anova_test(Stride_time ~ Age + Group)
get_anova_table(res.aov8)

# Stride length
res.aov9 <- grw %>% anova_test(Stride_Length ~ Age + Group)
get_anova_table(res.aov9)

# Step time
res.aov10 <- grw %>% anova_test(Step_time ~ Age + Group)
get_anova_table(res.aov10)

# Step length
res.aov11 <- grw %>% anova_test(Step_length ~ Age + Group)
get_anova_table(res.aov11)

# Step width
res.aov12 <- grw %>% anova_test(Step_width ~ Age + Group)
get_anova_table(res.aov12)

# Swing time
res.aov13 <- grw %>% anova_test(Swing ~ Age + Group)
get_anova_table(res.aov13)

# Single support
res.aov14 <- grw %>% anova_test(Single ~ Age + Group)
get_anova_table(res.aov14)

# ------------------------------------------------------------
# Linearity assumption check
#
# ANCOVA assumes a linear relationship between the covariate
# (Age) and the dependent variable within each combination of
# Group and Time. The scatter plot with a smoothed loess curve
# is the standard visual diagnostic.
# ------------------------------------------------------------
ggscatter(spt, x = "Age", y = "Step_width",
          facet.by = c("Group", "Time"),
          short.panel.labs = FALSE) +
  stat_smooth(method = "loess", span = 0.9)
