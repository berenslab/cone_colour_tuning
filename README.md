# Ancestral circuits for vertebrate colour vision emerge at the first retinal synapse

This repository includes the code of the model and the statistical analysis for the experimental recordings of the paper: "Ancestral circuits for vertebrate colour vision emerge at the first retinal synapse" (biorxiv: https://www.biorxiv.org/content/10.1101/2020.10.26.356089v2.full)

A list of the folders:
- data: some data used for the analysis in this repo. The full dataset of this paper can be found in the dryad archive. 
- EM_Confocal_Fig3_S3: contains the analysis of the numeric tables for the EM and the confocal microscopy data. 
- log_opsins_Fig1f_f_S4: analysis of the opsin curves: first fitted a log transform to the HC block condition (Fig.1f and Fig.2b), then the other way round for calculating the data distribution (Fig. S4)
- differences_Fig2d-h: Analysis of significant differences between the cone recordings under different conditions. Also includes calculation of the zero crossings for green and blue cone recordings and comparisons between retinal regions (cf. Fig.1g).
- HC_voltage_analysis_Fig4h-j: Analysis and clustering of the HC voltage recordings.
- modelling: includes the cone-HC interaction model and the tools for the SBI method and all analysis. See the readme file in this folder for more information.
