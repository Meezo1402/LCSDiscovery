# LPSDiscovery
Longitudinal Psychological Subtype Discovery in An Active-Duty Military Cohort with Combat-Related PTSD
## Overview

This repository contains the analysis code accompanying Istanbouli et al. 2026. We apply a custom negative 
binomial hidden Markov model (HMM) to longitudinal psychological 
symptom data to identify latent psychiatric subtypes, and integrate 
multi-omics data (metabolomics, DNA methylation, clinical labs) via 
the DIABLO framework to discover molecular signatures associated 
with HMM-derived subtypes.

## Repository Structure

- `code/01_HMM_fitting.R` — HMM model fitting and state assignment
- `code/02_DIABLO_classification.R` — Multi-omics integration and classification
- `data/README.md` — Data availability statement

## Dependencies

- R 4.1.3 (arm64)
- -R 4.4.1 (arm64)
- depmixS4 1.5-0
- gamlss 5.4-12
- gamlss.dist 6.0-5
- mixOmics 6.28
- mice 3.18

## License

This project is licensed under the MIT License.
