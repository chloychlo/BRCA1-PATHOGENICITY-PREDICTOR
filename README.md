# BRCA1 Variant Pathogenicity Prediction

## Project Overview
This repository contains a machine learning pipeline designed to predict the pathogenicity of variants in the **BRCA1** gene in the context of the MeetEU project carried at Sorbonne Université. Using a frozen snapshot of ClinVar data, the model classifies variants to help distinguish between pathogenic mutations and benign variants, addressing the challenge of Variants of Uncertain Significance (VUS).

Note: Project in progress - beyond predicting the pathogenicity of variants, the goal of the project is to also ultimately identify what underlying mechanism is into play. 
---

## Dataset Description
The project is built upon a frozen ClinVar dataset specifically filtered for single amino-acid substitutions in BRCA1 (reference transcript: `NM_007294.4`, protein: `NP_009225.1`, GRCh38 assembly).

* **`brca1_clinvar_2026-10-05_reference_labels.tsv`**: High-confidence ClinVar variants (Pathogenic/Likely Pathogenic vs. Benign/Likely Benign) used as the training and cross-validation set.
* **`brca1_clinvar_2026-10-05_predict.tsv`**: The target prediction set consisting of VUS and conflicting classifications.
* **`brca1_clinvar_2026-10-05_all.tsv`**: Complete unfiltered dataset.

___

## Model
The classification pipeline is built using scikit-learn:

* **Model**: Random Forest classifier (`RandomForestClassifier`)
* **Library**: `scikit-learn` v1.2.0
* **Hyperparameters**:
  * `n_estimators=100`: number of trees
  * `random_state=42`
  * Stratified train-test split to maintain class balance between pathogenic and benign variants during training.
