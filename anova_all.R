# ============================================================
# Mixed ANOVAs and Repeated-Measures ANOVAs
#
# Three analyses are bundled in this single file (the same way
# they lived in my original working file), covering:
#
#   1. COP balance metrics from a Rocker Board (RB) and Wobble
#      Board (WB) intervention study with adults and children
#   2. Joint kinematic peak values and percentages across the
#      gait cycle (ankle, hip, knee) from the unilateral ankle
#      loading study
#   3. Spatiotemporal parameters from the unilateral ankle
#      loading study
#
# Statistical methods used:
#   - 3-way mixed ANOVA: aov_car with group * intervention * time
#     and Error(ID/(intervention*time))
#   - 2-way mixed ANOVA: anova_test with Group (between) and
#     Time (within)
#   - Post-hoc interaction contrasts via emmeans + pairs
#
# Inputs:
#   - cop.csv   : COP metrics produced by my Python cop_analysis.ipynb
#   - jointGC1.csv (or similar): joint kinematics
#                 from my Python joint kinematics pipeline
#   - spt data (variable name "spt"): spatiotemporal parameters
#
# Author : Yeon-Joo Kang
# Inst.  : Georgia State University
# Period : 2021-2024
# Status : Archived dissertation analysis
# ============================================================

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
library(Matrix)
library(foreign)
library(lme4)
library(sjPlot)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(afex)
library(tidyr)
library(emmeans)

# ============================================================
# Analysis 1: COP intervention study
#
# 3-way mixed ANOVA on Center-of-Pressure (COP) outcomes from
# an RB/WB balance training intervention. The between-subjects
# factor is group (adults vs children). The within-subjects
# factors are intervention (RB, WB) and time (pre, post).
#
# Outcomes tested:
#   - Ellipse area, height, width, angle (95% confidence ellipse)
#   - Sway range AP, sway range ML
#   - Sway velocity AP, sway velocity ML
# ============================================================

cop <- read.csv("cop.csv", header = TRUE)

# Set factor levels explicitly so the plots and contrasts read
# in a deliberate order
cop$time <- factor(cop$time, levels = c("pre", "pos"))
cop$intervention <- factor(cop$intervention, levels = c("RB", "WB"))
cop$group <- factor(cop$group, levels = c("adults", "children"))

# ------------------------------------------------------------
# Quick exploratory boxplot
# ------------------------------------------------------------
bxp <- ggboxplot(cop, x = "session", y = "ellipse_area_95.",
                 color = "group", palette = "jco")
bxp

# ------------------------------------------------------------
# Ellipse area
# ------------------------------------------------------------
Mixed.EA <- aov_car(ellipse_area_95. ~ group * intervention * time
                    + Error(ID / (intervention * time)),
                    data = cop)
knitr::kable(nice(Mixed.EA))

# Intervention effect within each group
Mixed_EA_Interaction <- emmeans(Mixed.EA, ~ intervention | group)
pairs(Mixed_EA_Interaction)

bxp <- ggboxplot(cop, x = "session", y = "ellipse_area_95.",
                 color = "group", palette = "jco")
bxp

# ------------------------------------------------------------
# Ellipse height (major axis)
# ------------------------------------------------------------
Mixed.h <- aov_car(ellipse_height ~ group * intervention * time
                   + Error(ID / (intervention * time)),
                   data = cop)
knitr::kable(nice(Mixed.h))

# Time effect within each group
Mixed_h_Interaction <- emmeans(Mixed.h, ~ time | group)
pairs(Mixed_h_Interaction)

bxp <- ggboxplot(cop, x = "session", y = "ellipse_height",
                 color = "group", palette = "jco")
bxp

# ------------------------------------------------------------
# Ellipse width (minor axis)
# ------------------------------------------------------------
Mixed.w <- aov_car(ellipse_width ~ group * intervention * time
                   + Error(ID / (intervention * time)),
                   data = cop)
knitr::kable(nice(Mixed.w))

# Intervention effect within each group
Mixed_w_Interaction <- emmeans(Mixed.w, ~ intervention | group)
pairs(Mixed_w_Interaction)

bxp <- ggboxplot(cop, x = "session", y = "ellipse_width",
                 color = "group", palette = "jco")
bxp

# ------------------------------------------------------------
# Ellipse rotation angle
# ------------------------------------------------------------
Mixed.ag <- aov_car(ellipse_angle ~ group * intervention * time
                    + Error(ID / (intervention * time)),
                    data = cop)
knitr::kable(nice(Mixed.ag))

bxp <- ggboxplot(cop, x = "session", y = "ellipse_angle",
                 color = "group", palette = "jco")
bxp

# ------------------------------------------------------------
# Sway range AP (anterior-posterior)
# ------------------------------------------------------------
bxp <- ggboxplot(cop, x = "session", y = "sway_range_AP",
                 color = "group", palette = "jco")
bxp

Mixed.ap <- aov_car(sway_range_AP ~ group * intervention * time
                    + Error(ID / (intervention * time)),
                    data = cop)
knitr::kable(nice(Mixed.ap))

Mixed_ap_Interaction <- emmeans(Mixed.ap, ~ time | group)
pairs(Mixed_ap_Interaction)

# ------------------------------------------------------------
# Sway range ML (medial-lateral)
# ------------------------------------------------------------
bxp <- ggboxplot(cop, x = "session", y = "sway_range_ML",
                 color = "group", palette = "jco")
bxp

Mixed.ml <- aov_car(sway_range_ML ~ group * intervention * time
                    + Error(ID / (intervention * time)),
                    data = cop)
knitr::kable(nice(Mixed.ml))

# Time effect within each group
Mixed_ml_Interaction <- emmeans(Mixed.ml, ~ time | group)
pairs(Mixed_ml_Interaction)

# Intervention effect within each group
Mixed_ml_Interaction <- emmeans(Mixed.ml, ~ intervention | group)
pairs(Mixed_ml_Interaction)

# ------------------------------------------------------------
# Sway velocity AP
# ------------------------------------------------------------
bxp <- ggboxplot(cop, x = "session", y = "sway_velocity_AP",
                 color = "group", palette = "jco")
bxp

Mixed.vap <- aov_car(sway_velocity_AP ~ group * intervention * time
                     + Error(ID / (intervention * time)),
                     data = cop)
knitr::kable(nice(Mixed.vap))

Mixed_vap_Interaction <- emmeans(Mixed.vap, ~ time | group)
pairs(Mixed_vap_Interaction)

# ------------------------------------------------------------
# Sway velocity ML
# ------------------------------------------------------------
bxp <- ggboxplot(cop, x = "session", y = "sway_velocity_ML",
                 color = "group", palette = "jco")
bxp

Mixed.vml <- aov_car(sway_velocity_ML ~ group * intervention * time
                     + Error(ID / (intervention * time)),
                     data = cop)
knitr::kable(nice(Mixed.vml))


# ============================================================
# Analysis 2: Joint kinematic 2-way mixed ANOVAs
#
# Joint angle peaks and percentages of the gait cycle compared
# between groups (between-subjects) across two time points
# (within-subjects). Uses anova_test from rstatix.
#
# Note: the working dataframe in the original analysis is `spt`,
# which holds the joint kinematic and spatiotemporal columns
# side by side. Set the factors before running the ANOVAs.
# ============================================================

spt$Group <- as.factor(spt$Group)
spt$Time  <- as.factor(spt$Time)

# ------------------------------------------------------------
# Knee ROM (sanity check using "joint" dataframe, if separate)
# ------------------------------------------------------------
res.aov <- anova_test(data = joint, dv = KneeROM, wid = Subject,
                      between = Group, within = Time)
res.aov

# ------------------------------------------------------------
# Ankle
# ------------------------------------------------------------

# Plantarflexion peak
res.aov1 <- anova_test(data = spt, dv = AnklePLantar, wid = Subject,
                       between = Group, within = Time)
res.aov1

# Dorsiflexion peak
res.aov2 <- anova_test(data = spt, dv = AnkleDorsi, wid = Subject,
                       between = Group, within = Time)
res.aov2

# Ankle ROM
res.aov3 <- anova_test(data = spt, dv = AnkleROM, wid = Subject,
                       between = Group, within = Time)
res.aov3

# Plantarflexion timing as percent of gait cycle
res.aov4 <- anova_test(data = spt, dv = AnklePLantarTime, wid = Subject,
                       between = Group, within = Time)
res.aov4

# Dorsiflexion timing as percent of gait cycle
res.aov5 <- anova_test(data = spt, dv = AnkleDorsiTime, wid = Subject,
                       between = Group, within = Time)
res.aov5

# ------------------------------------------------------------
# Hip
# ------------------------------------------------------------

bxp <- ggboxplot(spt, x = "Time", y = "HipExt.",
                 color = "Group", palette = "jco")
bxp

# Hip extension peak
res.aov6 <- anova_test(data = spt, dv = HipExt, wid = Subject,
                       between = Group, within = Time)
res.aov6

# Hip flexion peak
res.aov7 <- anova_test(data = spt, dv = HipFlex, wid = Subject,
                       between = Group, within = Time)
res.aov7

# Hip ROM
res.aov8 <- anova_test(data = spt, dv = HipROM, wid = Subject,
                       between = Group, within = Time)
res.aov8

# Hip extension timing as percent of gait cycle
res.aov9 <- anova_test(data = spt, dv = HipExt., wid = Subject,
                       between = Group, within = Time)
res.aov9

# Hip flexion timing as percent of gait cycle
res.aov10 <- anova_test(data = spt, dv = HipFlex., wid = Subject,
                        between = Group, within = Time)
res.aov10

# ------------------------------------------------------------
# Knee
# ------------------------------------------------------------

bxp <- ggboxplot(spt, x = "Time", y = "KneeFlex.",
                 color = "Group", palette = "jco")
bxp

# Knee extension peak
res.aov11 <- anova_test(data = spt, dv = KneeExt, wid = Subject,
                        between = Group, within = Time)
res.aov11

# Knee flexion peak
res.aov12 <- anova_test(data = spt, dv = KneeFlex, wid = Subject,
                        between = Group, within = Time)
res.aov12

# Knee ROM
res.aov13 <- anova_test(data = spt, dv = KneeROM, wid = Subject,
                        between = Group, within = Time)
res.aov13

# Knee extension timing as percent of gait cycle
res.aov14 <- anova_test(data = spt, dv = KneeExt., wid = Subject,
                        between = Group, within = Time)
res.aov14

# Knee flexion timing as percent of gait cycle
res.aov15 <- anova_test(data = spt, dv = KneeFlex., wid = Subject,
                        between = Group, within = Time)
res.aov15

# ------------------------------------------------------------
# Post hoc analysis (representative example: ankle plantarflexion)
# ------------------------------------------------------------

Mixed.SG <- aov_car(AnklePLantar ~ Group * Time + Error(Subject / Time),
                    data = spt)
knitr::kable(nice(Mixed.SG))

# Group differences at each Time
Single_Interaction <- emmeans(Mixed.SG, ~ Group | Time)
pairs(Single_Interaction)

# Time differences within each Group
Single_Interaction <- emmeans(Mixed.SG, ~ Time | Group)
pairs(Single_Interaction)


# ============================================================
# Analysis 3: Spatiotemporal between-group ANOVAs
#
# One-way ANOVAs comparing groups on spatiotemporal parameters.
# Uses anova_test from rstatix with Group as the between-subjects
# factor only (no within-subjects factor here).
# ============================================================

# Stride time
res.aov1 <- anova_test(data = spt, dv = Stride_time, wid = Subject,
                       between = Group)
res.aov1

# Stride length
res.aov2 <- anova_test(data = spt, dv = Stride_Length, wid = Subject,
                       between = Group)
res.aov2

# Step time
res.aov3 <- anova_test(data = spt, dv = Step_time, wid = Subject,
                       between = Group)
res.aov3

# Step length
res.aov4 <- anova_test(data = spt, dv = Step_length, wid = Subject,
                       between = Group)
res.aov4

# Step width
res.aov5 <- anova_test(data = spt, dv = Step_width, wid = Subject,
                       between = Group)
res.aov5

# Swing time
res.aov6 <- anova_test(data = spt, dv = Swing, wid = Subject,
                       between = Group)
res.aov6

# Single support time
res.aov7 <- anova_test(data = spt, dv = Single, wid = Subject,
                       between = Group)
res.aov7
