# ============================================================
# Independent-Samples t-Tests on Joint Kinematic Peaks and ROM
#
# Pairwise group comparison of joint kinematic outcomes
# (ankle, hip, knee peaks and ROM) from the joint kinematics
# pipeline output. Each test compares the two levels of `Group`
# assuming equal variances.
#
# Input:
#   - jointGC1.csv : joint kinematic descriptive measures
#                    produced by my Python joint kinematics
#                    pipeline (Stage 3 output)
#
# Author : Yeon-Joo Kang
# Inst.  : Georgia State University
# Period : 2021-2024
# Status : Archived dissertation analysis
# ============================================================

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------
library(ggpubr)

# Load joint kinematic descriptive dataset
joint <- read.csv("jointGC1.csv", header = TRUE)

# Set Group as factor
joint$Group <- as.factor(joint$Group)

# Inspect structure
str(joint)

# ------------------------------------------------------------
# Ankle
# ------------------------------------------------------------

# Plantarflexion percentage of gait cycle
t.test(AnklePLantar. ~ Group, var.equal = TRUE, data = joint)

# Dorsiflexion percentage of gait cycle
t.test(AnkleDorsix. ~ Group, var.equal = TRUE, data = joint)

# Ankle range of motion
t.test(AnkleROM ~ Group, var.equal = TRUE, data = joint)

# ------------------------------------------------------------
# Hip
# ------------------------------------------------------------

# Hip flexion peak
t.test(HipFlex ~ Group, var.equal = TRUE, data = joint)

# Hip extension percentage of gait cycle
t.test(HipExt. ~ Group, var.equal = TRUE, data = joint)

# Hip range of motion
t.test(HipROM ~ Group, var.equal = TRUE, data = joint)

# ------------------------------------------------------------
# Knee
# ------------------------------------------------------------

# Knee flexion peak
t.test(KneeFlex ~ Group, var.equal = TRUE, data = joint)

# Knee flexion percentage of gait cycle
t.test(KneeFlex. ~ Group, var.equal = TRUE, data = joint)

# ------------------------------------------------------------
# Boxplot example (Hip ROM)
# ------------------------------------------------------------
bxp <- ggboxplot(joint,
                 y = "HipROM",
                 color = "Group",
                 palette = "jco")
bxp
