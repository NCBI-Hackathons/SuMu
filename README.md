# SuMu: Estimating Association of Genetic Features with Survival

SuMu is an R package that aims to make it easier to estimate the association of genetic features (such as mutations, gene expression, copy number, or methylation) with survival (time from diagnosis/treatment to clinical event) or benefit from therapy, in the context of clinical data.

## Installation

```
library(devtools)
install_github("NCBI-Hackathons/SuMu")
```

## Usage

## Description

This package aims to make it easier to fit [Stan](https://mc-stan.org) models in R to estimate the association of genetic variants with survival (time from diagnosis/treatment to clinical event) or benefit from therapy, in the context of clinical data.

Under the hood, it uses [rstan](https://cran.r-project.org/package=rstan) and [rstanarm](https://cran.r-project.org/package=rstanarm) to fit a variety of models. Given a dataframe of clinical data & a matrix of genetic features, this package helps prepare the data for input to the fit functions, providing an interface for these models which is tailored for the use case. It also provides several utilities appropriate for summarizing the output in a manner bioinformatics professionals may be accustomed to.

The analysis workflow would look something like:

1. Prepare your clinical data.frame & genetic features for analysis (we won't help with this step)
2. Input the data to the `fit_???` function (e.g. fit_glm or fit_surv)
3. Check convergence
4. Summarize the posterior estimates of parameters
5. Generate posterior-predicted distributions for observed data

## Background

[stan](https://mc-stan.org) has powerful interfaces to R through [rstan](https://cran.r-project.org/package=rstan) and its companion package [rstanarm](https://cran.r-project.org/package=rstanarm). Together these enable fitting of a wide range of models using the most advanced algorithms for Bayesian inference (NUTS and HMC), including multilevel or "mixed-effects" models, GAMMs, and others. 

These models have potential to be widely applicable to the genetic and bioinformatics research community, enabling improved methods for discovery of prognostic biomarkers, predictive biomarkers, surrogate biomarkers, and/or clustering among other use cases. 

Since genetic data are often expensive to obtain, sample sizes are often limited and so these analyses typically fall under an **exploratory data analysis** use case. This is a perfect use case for the application of rich, informative models to moderately-sized data.

However, it is not clear whether evaluation of genetic features using Stan is feasible or practical. The analysis of genetic features often presents a "large-p" problem, in which there are many more features (p) than observations (n). This is often summarized as (p >> n). 

This poses two challenges for analysis: 

1. The computational burden of working with large matrices of features
2. The need for regularizing priors to limit the number of features associated with outcome

Given Stan's support for sparse matrix multiplication, we expect that this style of analysis will be feasible with Stan. However, this remains to be demonstrated. This package is in part a proof-of-concept to assess feasibility.


