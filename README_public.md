# Cooperative Breeding and Female Song Evolution in Passerines

This repository contains phylogenetic comparative analyses investigating the correlated evolution of cooperative breeding and female song in passerine birds (Order Passeriformes).

## Overview

This project tests whether cooperative breeding (where individuals other than parents help raise offspring) and female song have evolved together across evolutionary time in songbirds. We examine how different social systems influence the evolution of vocal communication using multiple phylogenetic comparative methods.

## Key Components

- **Phylogenetic comparative analyses** using multiple statistical frameworks
- **Bias correction methods** addressing geographic and taxonomic sampling biases
- **Robustness testing** through jackknifing and multi-tree analyses
- **Multi-trait evolution models** examining interactions between territoriality, cooperative breeding, and female song

## Methods

The analyses employ:
- Stochastic character mapping (simmap)
- Brownie algorithm for differential evolutionary rates
- Phylogenetic path analysis with bias correction
- Bayesian phylogenetic mixed models (MCMCglmm)
- Phylogenetic logistic regression (phyloglm)

## Data Sources

Analyses integrate data from multiple published sources on:
- Cooperative breeding classifications (Cockburn 2006, Jetz et al. 2011, Riehl 2013, Griesser et al. 2017/2023)
- Female song presence and characteristics (Odom et al. 2014, Webb et al. 2016, others)
- Territoriality classifications (Tobias et al. 2016, Birds of the World)
- Phylogenetic trees based on Hackett et al. and Ericson et al. backbones

## Repository Status

This repository is actively being updated as analyses are refined for manuscript resubmission. The codebase includes both established methods and novel approaches for addressing reviewer concerns about bias and robustness.

## Requirements

Key R packages include:
- `phytools`, `ape`, `geiger` (phylogenetic analyses)
- `phylopath`, `phylolm` (path analysis and regression)
- `MCMCglmm` (Bayesian models)
- `tidyverse`, `ggplot2` (data manipulation and visualization)

## Citation

Manuscript in revision. Please contact for citation information.

## Contact

Kate T. Snyder

---
*Note: This repository represents ongoing research. Scripts and results may be updated as analyses are refined.*