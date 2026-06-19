# ============================================================
# Hierarchical Linear Modeling: Step Time and Step Length
# (Combined Final Models)
#
# Condensed final-model comparison for both step time and step
# length. The full bottom-up model selection process is in the
# separate files hlm_step_time.R and hlm_step_length.R.
#
# This script focuses on the comparison between intercept-only
# random effects and richer random-effects structures for the
# selected fixed-effects models.
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

# ------------------------------------------------------------
# Step time
# ------------------------------------------------------------

# Load step time dataset
ST <- read.csv("Data_HLM_ST.csv", header = TRUE)

# Set categorical predictors as factors with explicit levels
ST$loading <- factor(ST$loading,
                     levels = c("0", "1", "2", "3"))
ST$time <- factor(ST$time,
                  levels = c("1", "2", "3"))

# Random-intercept-only model
model1a <- lmer(steptime ~ 1 + time + side + loading + (1 | subjectid),
                data = ST)

# Random intercept and random slope for side, plus height as covariate
model1b <- lmer(steptime ~ 1 + time + side + loading + height
                + (1 + side | subjectid),
                data = ST)

# Display the more complex model
tab_model(model1b,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Compare via likelihood ratio test (REML)
anova(model1a, model1b, refit = FALSE)

# ------------------------------------------------------------
# Step length
# ------------------------------------------------------------

# Load step length dataset
SL <- read.csv("Data_HLM_SL.csv", header = TRUE)

# Set categorical predictors as factors with explicit levels
SL$loading <- factor(SL$loading,
                     levels = c("0", "1", "2", "3"))
SL$time <- factor(SL$time,
                  levels = c("1", "2", "3"))

# Random intercept and random slope for time
model2a <- lmer(steplength ~ 1 + time + loading
                + (1 + time | subjectid),
                data = SL)

# Random intercept and random slopes for both time and loading
model2b <- lmer(steplength ~ 1 + time + loading
                + (1 + time + loading | subjectid),
                data = SL)

# Display the more complex model
tab_model(model2b,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Compare via likelihood ratio test (REML)
anova(model2a, model2b, refit = FALSE)
