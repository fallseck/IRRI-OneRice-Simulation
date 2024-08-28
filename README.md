# IRRI-OneRice-Simulation
The IRRI breeding program has initiated an optimization process in its strategy to develop modern rice varieties by incorporating faster breeding techniques along with genomic selection to increase the rate of genetic gain while shortening the time it takes.
The main objective of this simulation was to evaluate the performance of the OneRice breeding schemes according to the parents' recycling length and the impact of genomic selection.

# Evaluated Breeding schemes

## 5Y: The 5-Year parent recycling scheme is based on rapid generation advance and genomic selection. It involves two stages of multi-location yield trial to select elite lines for population improvement and product development based on the estimated breeding values. The scheme is designed as a closed system, and the recycling of the elite lines as parents for the subsequent cycle occurs in year 5.

## 3Y: The 3-Year recycling scheme is derived from the 5-Years framework, with an early exit from the line fixation stage at F4 generation. 
The recycling of elite lines takes place in year 3 based on the GEBVs.

## 2YBP: The 2-Year between-cohort prediction is actually similar to the 3-Year scheme with a shift to between-cohort prediction of untested lines based on previous data. The inter-cohort prediction enables estimating the GEBVs in year 2 rather than waiting for data from the Stage 1 in year 3. 

## 2YWP: The 2-Year within-cohort is an upgrade of the previous 3-Year scheme by reducing the RGA exit at the F3. 

## Baseline: the baseline scheme is included as a reference, corresponding to the 5-Year scheme without GS. In this approach, parent selection and advancement are based solely on phenotypic data. 

# Software used:
R Core Team (2021). Copyright (C) 2021 The R Foundation for Statistical Computing Platform: R version 4.1.0 (2021-05-18) \
AlphaSimR program (Gaynor et al., 2023). Breeding Program Simulations (CRAN - Package AlphaSimR version 1.5.3). \
HPC MESO@LR-Platform at the University of Montpellier.
