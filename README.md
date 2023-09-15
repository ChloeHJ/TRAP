# TRAP: A robust deep learning workflow to predict CD8+ T-cell epitopes

TRAP allows for the robust prediction of CD8+ T-cell epitopes from MHC-I ligands. It employs a deep learning-based platform for predicting immunogenicity and a DecisionTree classifier to estimate the degree of correctness.

The dash web application can be found here: http://35.176.114.129:8050/ 

All required data, models and supplementary notebooks are stored in: https://drive.google.com/drive/folders/15A2P5xP2c-q48vVGPRB7h7uHEMycPYoX?usp=drive_link 

We provide gbm_example_test_data.csv as an example test data. The input data should have the same structure as the example dataset, and contain peptides having 9-10 amino acids in length. Minimal requirments for the input data are Peptide sequence and -log2(NetMHCpan rank score). Please select the pathogenic or self-antigen (cancer, autoantigens, etc.) model and insert your list of peptides on the application. The model may take 5-10 minutes to generate the output. We highly recommend the users to have <100 test peptides at a time due to limited server space. The "Tutorial for TRAP web application" document contains a detailed tutorial on how to use the web application.

Interpretation of the result: 
* TRAP: Immunogenicity (TRAP > 0.5: Positive, TRAP < 0.5: Negative)
* MCDropout: MaxProb value from averaged 100 Monte Carlo dropouts 
* Confidence: Degree of correctness predicted using MCDropout (TRAP > ~0.75: High confidence) 

This GitHub repository includes documentation for:
* TRAP-CNN model structure (model.py) 
* RSAT (RSAT.R) 

The data/ includes training data for these modules:
* pathogenic_db.csv and selfantigen_data_forMODEL.csv are the training datasets for TRAP 
* ood_dropout_pathogenic.csv and ood_dropout_selfantigen.csv are training data for OOD MCDropout-based OOD classifier
* Autoantigen_5data.csv are autoantigen and cancer-associated antigens gathered from 5 databases. hsProtemeDF.RData is reference human proteome sequence for RSAT calculation. Autoantigen_peptide_vector.RData and iedb_dissimilar_proteome_blosum62.RData are processed files from RSAT calculations. 

Please cite TRAP as following: Lee, C. H., Huh, J., Buckley, P. R., Jang, M., Pinho, M. P., Fernandes, R. A., Antanaviciute, A., Simmons, A., & Koohy, H. (2023). A robust deep learning workflow to predict CD8 + T-cell epitopes. Genome Medicine, 15(1), 70. https://doi.org/10.1186/s13073-023-01225-z

Available from: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01225-z



Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
