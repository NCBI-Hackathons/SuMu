# SuMu: Estimating Association of Genetic Features with Survival

A pipeline for association of genetic features and phenotypes

Current methods for estimating the association of genetic variants with survival have limitations.

The goal is to estimate the association of genetic variants with survival (time from diagnosis/treatment to clinical event) or benefit from therapy.

In this project, we will implement a possibly improved method for analysis using Bayesian survival models estimated using [Stan](https://mc-stan.org). The resulting package will make it easier for researchers to fit Bayesian survival models to genetic data.

The workflow would look something like:

1. Prepare your data for analysis (we won't help with this step)
2. Input the data to the `fit-model` function
3. Check convergence
4. Summarize the posterior estimates of parameters
5. Generate posterior-predicted distributions for observed data
