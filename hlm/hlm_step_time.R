# ============================================================
# Hierarchical Linear Modeling: Step Time
#
# Analysis of step time during treadmill walking with unilateral
# ankle loading. Step time data come from my Python spatiotemporal
# analysis pipeline (treadmill_spatiotemporal_analysis.ipynb).
#
# Bottom-up model selection workflow:
#   1. Null model (random intercept only)
#   2. Single fixed effects (interval, loading, side)
#   3. Combined fixed effects
#   4. Interaction terms (loading * side gives the best fit)
#   5. Random slopes (side and loading both improve fit)
#   6. Final model: model8d
#        steptime ~ 1 + loading * side + interval
#                   + (1 + side + loading | subjectid)
#
# Model comparison is done via likelihood ratio tests on REML
# fits (anova(..., refit = FALSE)), which is the appropriate
# choice when comparing models that differ only in random
# effects structure.
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

# Load step time dataset
# Columns: subjectid, steptime, interval (time), loading, side
ST <- read.csv("Data_HLM_ST.csv", header = TRUE)

# ------------------------------------------------------------
# Step 1: Null model (random intercept only)
#
# This baseline lets me see how much step time variance is
# between subjects vs within subjects across repeated measures.
# ------------------------------------------------------------
model0 <- lmer(steptime ~ 1 + (1 | subjectid), data = ST)
tab_model(model0,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# ------------------------------------------------------------
# Step 2: Single fixed effects
#
# Test each predictor on its own before combining them.
# ------------------------------------------------------------

# Time effect alone
model1 <- lmer(steptime ~ 1 + interval + (1 | subjectid), data = ST)
tab_model(model1,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Loading effect alone
model2 <- lmer(steptime ~ 1 + loading + (1 | subjectid), data = ST)
tab_model(model2,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Side (left vs right) effect alone
model3 <- lmer(steptime ~ 1 + side + (1 | subjectid), data = ST)
tab_model(model3,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# ------------------------------------------------------------
# Step 3: Combined fixed effects
# ------------------------------------------------------------

# Interval + loading
model4 <- lmer(steptime ~ 1 + interval + loading + (1 | subjectid), data = ST)
tab_model(model4,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# All three main effects
model5 <- lmer(steptime ~ 1 + interval + loading + side + (1 | subjectid), data = ST)
tab_model(model5,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# ------------------------------------------------------------
# Step 4: Interaction terms
#
# Test which two-way interaction (if any) improves the fit.
# ------------------------------------------------------------

# interval * loading
model6 <- lmer(steptime ~ 1 + interval * loading + side + (1 | subjectid), data = ST)
# not significant
tab_model(model6,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# interval * side
model7 <- lmer(steptime ~ 1 + interval * side + loading + (1 | subjectid), data = ST)
# not significant
tab_model(model7,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# loading * side  <-- this interaction is the meaningful one
model8 <- lmer(steptime ~ 1 + loading * side + interval + (1 | subjectid), data = ST)
tab_model(model8,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Compare model5 (no interaction) vs model8 (loading * side)
anova(model5, model8)
# Data: ST
# Models:
# model5: steptime ~ 1 + interval + loading + side + (1 | subjectid)
# model8: steptime ~ 1 + loading * side + interval + (1 | subjectid)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model5    9  808.75 840.08 -395.38   790.75
# model8   12  672.77 714.54 -324.39   648.77 141.98  3  < 2.2e-16 ***

# ------------------------------------------------------------
# Step 5: Random effects structure
#
# With the best fixed-effects structure identified (model8),
# test whether adding random slopes for each predictor improves
# model fit. Use REML (refit = FALSE) for these comparisons.
# ------------------------------------------------------------

# Random slope for interval (time)
model8a <- lmer(steptime ~ 1 + loading * side + interval
                + (1 + interval | subjectid), data = ST)
tab_model(model8a,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model8, model8a, refit = FALSE)
# not significant -> stay with model8
# Data: ST
# Models:
# model8:  steptime ~ 1 + loading * side + interval + (1 | subjectid)
# model8a: steptime ~ 1 + loading * side + interval + (1 + interval | subjectid)
#          npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# model8    12  687.06 728.82 -331.53   663.06
# model8a   17  696.76 755.93 -331.38   662.76 0.2926  5    0.9978

# Random slope for side
model8b <- lmer(steptime ~ 1 + loading * side + interval
                + (1 + side | subjectid), data = ST)
tab_model(model8b,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model8, model8b, refit = FALSE)
# significant -> random slope for side improves fit
# Data: ST
# Models:
# model8:  steptime ~ 1 + loading * side + interval + (1 | subjectid)
# model8b: steptime ~ 1 + loading * side + interval + (1 + side | subjectid)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model8    12  687.06 728.82 -331.53   663.06
# model8b   14  562.67 611.40 -267.34   534.67 128.38  2  < 2.2e-16 ***

# Random slope for loading
model8c <- lmer(steptime ~ 1 + loading * side + interval
                + (1 + loading | subjectid), data = ST)
tab_model(model8c,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model8, model8c, refit = FALSE)
# significant -> random slope for loading improves fit
# Data: ST
# Models:
# model8:  steptime ~ 1 + loading * side + interval + (1 | subjectid)
# model8c: steptime ~ 1 + loading * side + interval + (1 + loading | subjectid)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model8    12  687.06 728.82 -331.53   663.06
# model8c   21  680.38 753.47 -319.19   638.38 24.677  9   0.003351 **

# Random slopes for both side AND loading
# Uses bobyqa optimizer to handle the more complex covariance structure
model8d <- lmer(steptime ~ 1 + loading * side + interval
                + (1 + side + loading | subjectid), data = ST,
                control = lmerControl(optimizer = "bobyqa"))
tab_model(model8d,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model8c, model8d, refit = FALSE)
# significant -> stay with model8d (final model)
# Data: ST
# Models:
# model8c: steptime ~ 1 + loading * side + interval + (1 + loading | subjectid)
# model8d: steptime ~ 1 + loading * side + interval + (1 + side + loading | subjectid)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model8c   21  680.38 753.47 -319.19   638.38
# model8d   26  462.07 552.56 -205.03   410.07 228.31  5  < 2.2e-16 ***

# ------------------------------------------------------------
# Step 6: Diagnostics and prediction plots for the final model
# ------------------------------------------------------------

# Interaction plot
plot_model(model8d, type = "int")

# Multicollinearity check on the fixed-effects design matrix
round(cor(model.matrix(~ loading * side + interval, data = ST)[, -1]),
      digits = 2)
# No multicollinearity

# Diagnostic plots
plot_model(model8d, type = "diag")[3]  # residuals look somewhat normal
plot_model(model8d, type = "diag")[1]  # residuals look somewhat normal
plot_model(model8d, type = "diag")[4]  # residuals not constant across predictors
plot_model(model8d, type = "diag")[2]  # level-2 intercept residuals reasonable

# Predicted effects
plot_model(model8d, type = "pred", terms = c("loading [all]", "side")) +
  theme_sjplot2()

plot_model(model8d, type = "pred", terms = c("interval [all]", "side")) +
  theme_sjplot2()
