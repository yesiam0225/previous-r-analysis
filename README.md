# previous-r-analysis

This repository documents the **statistical analysis workflow** I wrote in R during my PhD research in Kinesiology at Georgia State University (2021–2024). It is the R-language counterpart to my [`previous-analysis`](https://github.com/yesiam0225/previous-analysis) Python notebook archive — together the two repositories cover the full pipeline from raw motion-capture data through inferential statistics.

These scripts are the analyses I wrote by hand during dissertation work. They reflect my actual model selection and statistical reasoning at the time, rather than production-ready code. A separate repository will host AI-assisted, more automated versions of these analyses.

## Background

- **Author**: Yeon-Joo Kang
- **Institution**: Georgia State University, Department of Kinesiology and Health
- **Period**: 2021–2024
- **Companion repository**: [`previous-analysis`](https://github.com/yesiam0225/previous-analysis) (Python notebooks producing the data analyzed here)

## Repository structure

```
previous-r-analysis/
├── README.md
│
├── hlm/                                    # Hierarchical Linear Modeling
│   ├── hlm_step_time.R                     # full model selection workflow
│   ├── hlm_step_length.R                   # full model selection workflow
│   └── hlm_overview.R                      # condensed final models
│
├── anova/                                  # Mixed and repeated-measures ANOVA
│   └── anova_all.R                         # COP + joint kinematics + spatiotemporal
│
├── ancova/                                 # Analysis of Covariance
│   └── ancova_age_adjusted.R               # age-adjusted group comparisons
│
└── t-test/                                 # Independent-samples t-tests
    └── group_comparison_t_test.R           # joint kinematic group comparison
```

## Statistical methods used

### Hierarchical Linear Modeling (`hlm/`)

Mixed-effects models for nested, repeated-measures gait data. Fit with `lme4::lmer`, displayed with `sjPlot::tab_model`, compared via likelihood ratio tests on REML fits (`anova(..., refit = FALSE)`, which is the correct comparison when models differ only in random-effects structure).

The two main HLM files (`hlm_step_time.R`, `hlm_step_length.R`) walk through the bottom-up model selection process — null model, single fixed effects, combined fixed effects, interaction terms, random slope structures, diagnostics — ending with the selected final model.

### Mixed ANOVA (`anova/`)

Three-way mixed ANOVA with `afex::aov_car` using `Error(ID/(intervention*time))` notation to specify the within-subjects structure. Two-way mixed ANOVAs with `rstatix::anova_test`. Post-hoc contrasts via `emmeans` with `pairs()`.

This file covers three analyses bundled in one script (as in my original working file):
1. COP balance metrics from a Rocker Board / Wobble Board intervention study (adults and children)
2. Joint kinematic peaks and percentages for ankle, hip, and knee
3. Spatiotemporal parameter comparisons

### ANCOVA (`ancova/`)

Age-adjusted group comparisons using `rstatix::anova_test` with the covariate listed first in the formula (`Age + Group * Time`). Post-hoc comparisons via `emmeans_test` with Bonferroni correction. Linearity assumption checked with `ggpubr::ggscatter` and a loess smoother across covariate values within each Group × Time cell.

Age control is essential in my work because I compare children and adults whose absolute gait parameters differ substantially due to body size and developmental stage.

### Independent t-tests (`t-test/`)

Simple two-group comparisons on joint kinematic peaks and ROM with `t.test(..., var.equal = TRUE)`. Used as a quick first pass before more elaborate ANOVAs.

## Relationship to the Python pipeline

Each R script consumes data produced by a notebook in the companion `previous-analysis` repository:

| R script | Python notebook that produces the input |
|---|---|
| `hlm_step_time.R`, `hlm_step_length.R` | `treadmill_spatiotemporal_analysis.ipynb` |
| `anova_all.R` (COP portion) | `cop_analysis.ipynb` |
| `anova_all.R` (joint portion) | `vicon_joint_kinematics_peak_rom.ipynb` + `joint_kinematics_ensemble_avg.ipynb` |
| `anova_all.R` (spatiotemporal portion) | `treadmill_spatiotemporal_analysis.ipynb` and `overground_spatiotemporal_from_events.ipynb` |
| `ancova_age_adjusted.R` | joint kinematics pipeline + spatiotemporal pipeline |
| `group_comparison_t_test.R` | joint kinematics pipeline |

The split between Python (custom data extraction, ensemble averaging, descriptive analysis, time-series statistical analysis via SPM) and R (mixed-effects modeling, mixed ANOVA, ANCOVA, post-hoc contrasts) reflects how each language's ecosystem best supports those steps.

## Notes on the data

The CSV inputs referenced in these scripts are not included in the repository to protect participant privacy. The scripts document the analytical approach; the data themselves remain with the IRB-approved data custodianship at GSU.

## Tooling

- **Language**: R
- **Mixed-effects**: `lme4`, `Matrix`
- **ANOVA / mixed ANOVA**: `afex`, `rstatix`
- **Post-hoc**: `emmeans`
- **Display and diagnostics**: `sjPlot`, `ggpubr`, `ggplot2`, `broom`
- **Data wrangling**: `tidyverse`, `tidyr`

## License

This repository is published as a record of past analytical work. The scripts are shared as-is for transparency and are not intended as a redistributable tool.
