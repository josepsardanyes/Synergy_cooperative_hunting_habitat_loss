README
======

Overview
--------
This folder contains an R script and two data files used to analyse how habitat loss (D) and the fraction of cooperative hunters (sigma) affect the probability of alternative ecological outcomes.

The analysis considers two scenarios:

- RCP: coexistence of resource (R), consumer (C), and predator (P).
- RC: coexistence of resource (R) and consumer (C), with predator extinction.

The script transforms the estimated probabilities to logits and fits third-order bivariate polynomial models as functions of habitat loss (D) and cooperative hunting (sigma). Model simplification is performed using backward elimination.

Files
-----
D_sigma_fit.txt
    R script containing the statistical analysis for the RCP and RC scenarios.

D_sigma_RCP_ar.txt
    Data used for the RCP (resource-consumer-predator coexistence) analysis.

D_sigma_RC_ar.txt
    Data used for the RC (resource-consumer coexistence with predator extinction) analysis.

Data format
-----------
Both data files are plain-text tables without column headers. Each row contains three columns:

    Column 1: D      - fraction of habitat loss
    Column 2: sigma  - fraction of cooperative hunters
    Column 3: prob   - estimated probability (fraction of sampled parameter and initial-condition combinations) of the corresponding ecological scenario

The R script reads these columns as V1, V2, and V3, respectively.

Analysis
--------
For each ecological scenario, the script:

1. Reads the corresponding data file using read.table().
2. Assigns the three columns to D, sigma, and prob.
3. Replaces probabilities equal to zero with 0.0001 to allow the logit transformation.
4. Computes the logit-transformed response:

       y = log(prob / (1 - prob))

5. Constructs polynomial and interaction terms involving D and sigma.
6. Fits a third-order bivariate polynomial using lm().
7. Applies backward model selection using step(..., direction = "backward").
8. Reports the final model using summary().

For the RC scenario, a remaining non-significant sigma^3 term is removed manually after backward selection, and the reduced model is fitted again.

Requirements
------------
The analysis requires R. No additional R packages are needed because all functions used in the script are part of base R and the standard stats package.

Running the analysis
--------------------
Place all three files in the same working directory. In R, set the working directory to this folder if necessary and run:

    source("D_sigma_fit.txt")

The script will print the model-selection procedure and summaries of the fitted models to the R console.

Notes
-----
The data files contain probabilities obtained from the numerical analyses associated with the study. Values of D and sigma are sampled on a discrete grid.

Zero probabilities are replaced by 0.0001 only for the purpose of computing finite logits. The original data files are not modified.

The script also defines some fourth-order polynomial terms that are not included in the fitted third-order models.

Anonymity
---------
This README and the accompanying analysis files contain no author-identifying information and are intended for anonymous peer review.
