# ============================================================
# Hierarchical Linear Modeling: Step Length
#
# Analysis of step length during treadmill walking with
# unilateral ankle loading. Step length data come from my
# Python spatiotemporal analysis pipeline
# (treadmill_spatiotemporal_analysis.ipynb).
#
# Bottom-up model selection workflow:
#   1. Null model (random intercept only)
#   2. Single fixed effects (time, loading, side)
#   3. Combined fixed effects (time + loading is the best)
#   4. Interaction terms (none improve fit)
#   5. Random slopes (loading slope improves fit)
#   6. Final model: model4b
#        steplength ~ 1 + interval + loading
#                     + (1 + loading | subjectid)
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

# Load step length dataset
# Columns: subjectid, steplength, time, loading, side
SL <- read.csv("Data_HLM_SL.csv", header = TRUE)

# ------------------------------------------------------------
# Step 1: Null model (random intercept only)
# ------------------------------------------------------------
model0 <- lmer(steplength ~ 1 + (1 | subjectid), data = SL)
tab_model(model0,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# ------------------------------------------------------------
# Step 2: Single fixed effects
# ------------------------------------------------------------

# Time effect alone
model1 <- lmer(steplength ~ 1 + time + (1 | subjectid), data = SL)
tab_model(model1,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Loading effect alone
model2 <- lmer(steplength ~ 1 + loading + (1 | subjectid), data = SL)
tab_model(model2,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Side effect alone
model3 <- lmer(steplength ~ 1 + side + (1 | subjectid), data = SL)
# Not significant
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

# Time + loading
model4 <- lmer(steplength ~ 1 + time + loading + (1 | subjectid), data = SL)
tab_model(model4,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# ------------------------------------------------------------
# Step 4: Interaction terms
# ------------------------------------------------------------

# time * loading interaction
model5 <- lmer(steplength ~ 1 + time * loading + (1 | subjectid), data = SL)
# Not significant
tab_model(model5,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# side * time interaction
model6 <- lmer(steplength ~ 1 + side * time + loading + (1 | subjectid), data = SL)
# Not significant
tab_model(model6,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# side * loading interaction
model7 <- lmer(steplength ~ 1 + time + side * loading + (1 | subjectid), data = SL)
# Not significant
tab_model(model7,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

# Decision: no interaction terms improve fit. Keep model4 as the
# fixed-effects structure and explore random slopes from here.

# ------------------------------------------------------------
# Step 5: Random effects structure
# ------------------------------------------------------------

# Random slope for time
model4a <- lmer(steplength ~ 1 + time + loading
                + (1 + time | subjectid), data = SL)
tab_model(model4a,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model4, model4a, refit = FALSE)
# Random slope for time significantly improves fit
# Data: SL
# Models:
# model4:  steplength ~ 1 + time + loading + (1 | subjectid)
# model4a: steplength ~ 1 + time + loading + (1 + time | subjectid)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model4     8 387.94 415.78 -185.97   371.94
# model4a   13 374.00 419.25 -174.00   348.00 23.938  5  0.0002232 ***

# Random slope for loading
model4b <- lmer(steplength ~ 1 + interval + loading
                + (1 + loading | subjectid), data = SL)
tab_model(model4b,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model4, model4b, refit = FALSE)
# Random slope for loading significantly improves fit even more
# Data: SL
# Models:
# model4:  steplength ~ 1 + time + loading + (1 | subjectid)
# model4b: steplength ~ 1 + time + loading + (1 + loading | subjectid)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model4     8 387.94 415.78 -185.97   371.94
# model4b   17 372.21 431.38 -169.10   338.21 33.726  9  9.976e-05 ***

# Random slopes for both time AND loading
model4c <- lmer(steplength ~ 1 + time + loading
                + (1 + time + loading | subjectid), data = SL,
                control = lmerControl(optimizer = "bobyqa"))
# Note: boundary (singular) fit warning - see help('isSingular')
tab_model(model4c,
          dv.labels = "",
          show.se = TRUE,
          show.ci = FALSE,
          show.r2 = FALSE,
          string.pred = "Fixed Effects",
          string.se = "SE")

anova(model4b, model4c, refit = FALSE)
# model4c is significantly better, but the singular fit warning
# suggests the random-effects structure is overparameterized.
# Decision: choose model4b as the final model (parsimonious + clean fit).
# Data: SL
# Models:
# model4b: steplength ~ 1 + time + loading + (1 + loading | subjectid)
# model4c: steplength ~ 1 + time + loading + (1 + time + loading | subjectid)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# model4b   17 372.21 431.38 -169.10   338.21
# model4c   28 353.65 451.10 -148.82   297.65 40.564 11  2.864e-05 ***

# ------------------------------------------------------------
# Step 6: Diagnostics and prediction plots for the final model
# ------------------------------------------------------------

plot_model(model4b, type = "diag")[3]  # residuals look somewhat normal
plot_model(model4b, type = "diag")[1]  # residuals look somewhat normal
plot_model(model4b, type = "diag")[4]  # residuals not constant across predictors
plot_model(model4b, type = "diag")[2]  # level-2 intercept residuals reasonable

# Predicted step length across interval, faceted by loading
plot_model(model4b, type = "pred", terms = c("interval [all]", "loading")) +
  theme_sjplot2()

# Predicted step length across loading, faceted by interval
plot_model(model4b, type = "pred", terms = c("loading [all]", "interval")) +
  theme_sjplot2()
