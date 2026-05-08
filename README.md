# Acetylcholinesterase (AChE) Inhibitor Discovery Pipeline
## Computational Drug Discovery using ChEMBL, QSAR Modelling, and Molecular Docking

**Course:** Pharmacoinformatics — Final Project
**Institution:** National University of Sciences and Technology (NUST), School of Interdisciplinary Engineering and Sciences (SINES)
**Programme:** BS Bioinformatics, UG-1
**Date:** 8th May 2026

---

## Authors

| Name | Roll No |
|------|---------|
| Namra Basharat | 476203 |
| Ghania Munir | 460673 |
| Hania Fahad | 455132 |
| Alishba Saleem | 479300 |

---

## Project Overview

This repository contains the complete computational drug discovery pipeline developed for the Pharmacoinformatics final project. The pipeline targets **Acetylcholinesterase (AChE; EC 3.1.1.7; ChEMBL Target: CHEMBL220)**, a clinically validated therapeutic target for Alzheimer's disease (AD). The work spans nine interconnected tasks: bioactivity data acquisition and curation, molecular descriptor and fingerprint generation, QSAR modelling, structure-activity relationship (SAR) analysis, protein structure prediction and preparation, multi-parameter virtual screening, molecular docking, ADMET/toxicity filtering, and multi-parameter lead prioritisation.

Alzheimer's disease is the most prevalent neurodegenerative disorder worldwide, affecting approximately 55 million people globally. It is characterised by progressive cognitive decline and memory loss, driven in part by the degeneration of cholinergic neurons and the resulting deficit in acetylcholine (ACh) at synapses. AChE is the enzyme responsible for hydrolysing ACh, thereby terminating cholinergic neurotransmission. Inhibiting AChE prolongs synaptic ACh availability and temporarily compensates for the cholinergic deficit, representing the most clinically validated pharmacological strategy for symptomatic AD treatment. Approved AChE inhibitors include Donepezil (1996), Rivastigmine (1997), and Galantamine (2001).

The pipeline integrates cheminformatics, machine learning, structural bioinformatics, molecular docking, and ADMET profiling to identify and characterise potent, drug-like AChE inhibitors from the ChEMBL bioactivity database. All computational steps are documented in this repository with associated data files, scripts, figures, and reports.

---

## Key Statistics

| Metric | Value |
|--------|-------|
| Raw bioactivity records retrieved | 9,731 |
| Final curated dataset | 5,604 unique compounds |
| Data reduction | 42.5% (4,127 records removed) |
| Potent inhibitors (pIC50 >= 7; IC50 <= 100 nM) | 1,481 |
| Very potent inhibitors (pIC50 >= 8; IC50 <= 10 nM) | 584 |
| Sub-nanomolar inhibitors (pIC50 >= 9; IC50 < 1 nM) | 170 |
| Most potent compound | CHEMBL4205425 (pIC50 = 10.96; IC50 = 0.011 nM) |
| AlphaFold2 model confidence (mean pLDDT) | 91.2 (Very High) |
| Virtual screening library (after Ro5 filtering) | 660 drug-like compounds |
| Best docking binding energy | -13.67 kcal/mol (CHEMBL4209803) |
| Primary lead compound | CHEMBL4209803 (Ligand 07) |
| Secondary lead compound | CHEMBL4205954 (Ligand 06) |

---

## Repository Structure

```
AChE-Drug-Discovery-Pipeline/
│
├── task-1-data-curation/
│   ├── raw_ache_data.csv
│   ├── cleaned_ache_data.csv
│   ├── ache_with_descriptors.csv
│   ├── fingerprints.npy
│   ├── feature_matrix.csv
│   ├── ACHE_QSAR_dataset.csv
│   └── potent_ache_inhibitors.csv
│
├── task-3-qsar/
│   ├── QSAR_Modeling.Rmd
│   ├── QSAR_Modeling.html
│   ├── Figure1_Predicted_vs_Observed.png
│   ├── Figure2_Residual_Error_Plot.png
│   ├── qsar_train_set.csv
│   ├── qsar_validation_set.csv
│   ├── qsar_test_set.csv
│   └── random_forest_model.rds
│
├── task-4-sar/
│   ├── Figure3_SAR_Boxplots.png
│   ├── Figure4_Chemical_Space_PCA.png
│   ├── Figure5_RF_Variable_Importance.png
│   ├── Figure6_pIC50_Distribution.png
│   ├── sar_analysis_results.csv
│   ├── pca_loadings.csv
│   └── wilcoxon_test_results.csv
│
├── task-5-protein-structure-prep/
│   ├── alphafold2-model/
│   │   ├── ACHE_rank001.pdb
│   │   ├── ACHE_rank001_prepared.pdbqt
│   │   ├── predicted_aligned_error.json
│   │   ├── plddt_scores.json
│   │   └── msa_coverage_plot.png
│   │
│   ├── pymol-preparation/
│   │   ├── 4EY7_apo_chainA.pdb
│   │   ├── donepezil_reference.pdb
│   │   └── protein_preparation_commands.pml
│   │
│   ├── docking-configuration/
│   │   ├── vina_config.txt
│   │   ├── grid_box_parameters.txt
│   │   └── active_site_residues.txt
│   │
│   ├── Figure5-1_MSA_Coverage.png
│   ├── Figure5-2_PAE_Matrix.png
│   ├── Figure5-3_All_Model_pLDDT.png
│   ├── Figure5-4_pLDDT_Profile.png
│   ├── Figure5-5_pLDDT_Histogram.png
│   ├── Figure5-6_Model_Comparison.png
│   ├── Figure5-7_Active_Site_pLDDT.png
│   ├── Figure5-8_PyMOL_Commands.png
│   ├── Figure5-9_Protein_Preparation.png
│   └── Figure5-10_AutoDock_Preparation.png
│
├── task-6-virtual-screening/
│   ├── screening_library.csv
│   ├── composite_scoring_results.csv
│   ├── top10_selected_compounds.csv
│   ├── docking_candidates/
│   │   ├── ACHE_Cand_01.pdbqt
│   │   ├── ACHE_Cand_02.pdbqt
│   │   ├── ACHE_Cand_03.pdbqt
│   │   ├── ACHE_Cand_04.pdbqt
│   │   ├── ACHE_Cand_05.pdbqt
│   │   ├── ACHE_Cand_06.pdbqt
│   │   ├── ACHE_Cand_07.pdbqt
│   │   ├── ACHE_Cand_08.pdbqt
│   │   ├── ACHE_Cand_09.pdbqt
│   │   └── ACHE_Cand_10.pdbqt
│   │
│   ├── Figure6-1_VS_Funnel.png
│   ├── Figure6-2_Composite_Score_Distribution.png
│   ├── Figure6-3_Chemical_Space_Map.png
│   ├── Figure6-4_Ligand_Efficiency.png
│   ├── Figure6-5_Potency_Distribution.png
│   ├── Figure6-6_pIC50_vs_CompositeScore.png
│   └── Figure6-7_Multi_Property_Profile.png
│
├── task-7-molecular-docking/
│   ├── docking_results/
│   │   ├── vina_scores.csv
│   │   ├── best_binding_poses/
│   │   └── interaction_analysis/
│   │
│   ├── docking_logs/
│   └── docking_visualizations/
│
├── docs/
│   ├── Pharmacoinformatics_Final_Report.docx
│   ├── Pharmacoinformatics_Final_Report.pdf
│   └── presentation_slides.pptx
│
└── README.md
```

---

## Pipeline Description

### Task 1 — Target and Ligand Set Definition

The therapeutic target selected for this study is Human Acetylcholinesterase (AChE), encoded by the ACHE gene (UniProt: P22303), with ChEMBL target identifier CHEMBL220. AChE was selected on the basis of its clinical and pharmacological validation as the primary target for Alzheimer's disease therapy. The enzyme presents a deep catalytic gorge of approximately 20 Angstroms in depth, containing two key pharmacophoric subsites: the Catalytic Anionic Site (CAS), which houses the catalytic triad (Ser203–His447–Glu334), and the Peripheral Anionic Site (PAS) at the gorge entrance, lined by residues Trp286 and Tyr72. Both sites are exploited by dual-site inhibitors such as Donepezil.

Bioactivity data were retrieved from the ChEMBL database (version 33) using the IC50 activity type, organism filter Homo sapiens, and target type Single Protein. The initial retrieval returned 9,731 bioactivity records.

---

### Task 2 — Bioactivity Data Acquisition and Curation

The raw dataset of 9,731 records was subjected to a four-step curation pipeline to ensure data quality and consistency prior to modelling. The curation steps are described below.

**Step 1 — Deduplication and Record Filtering:** Duplicate records arising from the same compound measured in different assays were removed. Records with ambiguous or missing activity values and those with non-standard relations (e.g., >, <) were also excluded. Approximately 2,213 records were removed at this stage.

**Step 2 — Activity Type Standardisation:** Only IC50 measurements were retained to ensure a consistent activity scale across all compounds. Records with non-standard activity types (Ki, Kd, percentage inhibition) and those with pChEMBL values below 4.0 (corresponding to IC50 > 100 micromolar, biologically insignificant) were removed. Approximately 1,914 records were excluded.

**Step 3 — Structure Standardisation:** Valid SMILES strings were verified for all retained compounds. Salt stripping was performed using RDKit's SaltRemover module, and canonical SMILES were generated using RDKit's MolStandardize function.

**Step 4 — Activity Scale Conversion:** IC50 values in nanomolar units were converted to the negative logarithmic pIC50 scale using the formula pIC50 = -log10(IC50 x 10^-9). This transformation produces a linear scale suitable for regression-based QSAR modelling.

The final curated dataset comprises 5,604 unique compounds with IC50 measurements against human AChE.

| Curation Stage | Records | Notes |
|----------------|---------|-------|
| Raw ChEMBL download | 9,731 | All IC50 records for CHEMBL220 |
| After deduplication | ~7,518 | Duplicates removed |
| After activity filtering | ~5,604 | IC50 only, valid pChEMBL >= 4.0 |
| After structure standardisation | 5,604 | Final curated dataset |

**Potency Distribution of the Curated Dataset:**

| Category | pIC50 Range | IC50 Equivalent | Count | Percentage |
|----------|-------------|-----------------|-------|------------|
| Weak | < 6 | > 1 uM | 2,028 | 36.2% |
| Moderate | 6 - 7 | 100 nM - 1 uM | 2,095 | 37.4% |
| Potent | 7 - 8 | 10 - 100 nM | 897 | 16.0% |
| Very Potent | 8 - 9 | 1 - 10 nM | 414 | 7.4% |
| Sub-nanomolar | >= 9 | < 1 nM | 170 | 3.0% |
| Total | | | 5,604 | 100% |

---

### Task 3 — Molecular Descriptor Calculation and QSAR Modelling

**Descriptor Computation:** Six physicochemical descriptors were computed using RDKit for all 5,604 compounds: Molecular Weight (MW), Lipophilicity (AlogP), Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA), Rotatable Bonds, and Topological Polar Surface Area (TPSA). Lipinski's Rule of Five compliance was evaluated to assess drug-likeness; 59.5% of compounds were found to be fully RO5-compliant, indicating good oral bioavailability potential.

**Morgan Fingerprint Generation:** 2048-bit Morgan circular fingerprints (radius 2, equivalent to ECFP4) were generated for all compounds using RDKit's GetMorganFingerprintAsBitVect function. These fingerprints encode circular atomic neighbourhoods and are widely used in virtual screening and QSAR modelling. Fingerprints are stored in binary NumPy format (fingerprints.npy).

**QSAR Feature Matrix:** A QSAR-ready feature matrix of dimensions 5,604 x 7 (six physicochemical descriptors plus pIC50) was constructed. All descriptors were scaled to zero mean and unit variance using StandardScaler prior to model training. This matrix is provided as feature_matrix.csv.

**Dataset Splitting:** A stratified scaffold split (70/15/15) was employed: training set (3,897 compounds), validation set (835 compounds), and test set (834 compounds). Stratification was performed by binning pIC50 values into five equal-frequency groups prior to splitting, ensuring proportional representation of the full biological activity range (pIC50: 4.00–10.96) in all subsets.

**Models Trained:**

*Linear Regression (Baseline):* A multiple linear regression model was trained using R's lm() function on the six-descriptor feature matrix without regularisation. This model serves as an interpretable baseline.

*Random Forest (Non-linear):* A Random Forest regression model was configured with 300 trees (ntree = 300) and default mtry parameter. Random Forest models complex non-linear interactions between descriptors without requiring explicit feature engineering.

| Model | Split | R2 | RMSE | MAE |
|-------|-------|----|------|-----|
| Linear Regression | Train | 0.073 | 1.266 | 1.020 |
| Linear Regression | Test | 0.074 | 1.247 | 1.002 |
| Random Forest | Train | 0.789 | 0.689 | 0.541 |
| Random Forest | Test | 0.267 | 1.110 | 0.864 |

The near-identical training and test metrics for Linear Regression confirm underfitting — the model explains only 7% of variance in pIC50 and lacks the capacity to capture non-linear structure-activity relationships. The Random Forest substantially outperforms the baseline (test R2: 0.267 vs 0.074), confirming non-linearity in the data, but shows a large gap between training and test R2 consistent with overfitting. Both models would benefit from richer fingerprint representations and hyperparameter optimisation in a production pipeline.

**Output Files:**

| File | Description |
|------|-------------|
| cleaned_ache_data.csv | Curated ChEMBL dataset (5,604 compounds) |
| ache_with_descriptors.csv | Dataset with computed RDKit descriptors |
| feature_matrix.csv | Standardised QSAR feature matrix (5,604 x 7) |
| ACHE_QSAR_dataset.csv | Full QSAR-ready dataset |
| Person2_QSAR_Modeling.Rmd | R Markdown script for model training and evaluation |
| Person2_QSAR_Modeling.html | Rendered HTML report |
| Figure1_Predicted_vs_Observed.png | Predicted vs Observed pIC50 diagnostic plot |
| Figure2_Residual_Error_Plot.png | Residual error diagnostic plot |

---

### Task 4 — Structure-Activity Relationship (SAR) Analysis

SAR analysis was performed to identify which physicochemical properties are most strongly associated with AChE inhibitory potency. The curated dataset was divided at the median pIC50 (6.00) into two groups: High Activity (pIC50 >= 6.00; n = 2,807) and Low Activity (pIC50 < 6.00; n = 2,759). Four descriptors were compared between groups using the non-parametric Wilcoxon rank-sum test, as descriptor distributions are non-normal.

| Descriptor | High Activity Mean | Low Activity Mean | p-value | Significant |
|------------|-------------------|------------------|---------|-------------|
| Molecular Weight (Da) | 219.70 | 198.74 | < 0.0001 | Yes |
| logP | 4.929 | 4.247 | < 0.0001 | Yes |
| H-Bond Donors | 4.228 | 4.145 | 0.0048 | Yes |
| H-Bond Acceptors | 4.304 | 4.263 | 0.0426 | Yes |

All four descriptors show statistically significant differences between activity groups. High-activity compounds are notably larger (MW +21 Da) and more lipophilic (logP +0.68), consistent with the hydrophobic, gorge-shaped AChE binding site lined with aromatic residues (Trp84, Trp279, Phe330). MW and logP are the most strongly discriminating descriptors (p < 0.0001), while HBD and HBA show smaller but statistically significant differences, suggesting that H-bonding capacity to catalytic triad residues (Ser200, His440, Glu327) contributes modestly to potency.

**Chemical Space Visualisation (PCA):** PCA was performed on a five-descriptor matrix (MW, logP, HBD, HBA, RO5 violations), standardised prior to analysis. The first two principal components together capture 88.4% of total descriptor variance (PC1: 55.9%, PC2: 32.5%), providing a faithful 2D representation of the chemical space. A weak but visible gradient from low to high pIC50 is observed along the PC1 axis, which is dominated by MW and HBA loadings. The three activity ellipses (low, medium, high) show substantial overlap, explaining the difficulty of both QSAR models in achieving high R2.

**Random Forest Variable Importance:** The variable importance analysis confirms that MW-related descriptors and lipophilicity-associated features are the most influential predictors of pIC50, consistent with the known pharmacology of AChE inhibitors where molecular size and hydrophobicity are key determinants of binding at the dual CAS/PAS sites.

**Output Files:**

| File | Description |
|------|-------------|
| Figure3_SAR_Boxplots.png | Descriptor comparison by activity group |
| Figure4_Chemical_Space_PCA.png | PCA chemical space visualisation |
| Figure5_RF_Variable_Importance.png | Random Forest feature importance plot |

---

### Task 5 — Protein Structure Prediction and Preparation

**AlphaFold2 Structure Prediction:** The three-dimensional structure of Human AChE (UniProt P22303, 614 residues) was predicted de novo using ColabFold v1.5 (AlphaFold2-ptm implementation) with MMseqs2-based multiple sequence alignment searching UniRef and Environmental databases. A total of 19,679 homologous sequences were retrieved, providing a deep co-evolutionary signal for the catalytic domain. Five independent models were generated; the rank 1 model (model_3) was selected based on the highest mean pLDDT of 91.2 and pTM of 0.881.

| Rank | Model | Mean pLDDT | pTM | Status |
|------|-------|-----------|-----|--------|
| 1 | model_3 | 91.2 | 0.881 | Selected |
| 2 | model_4 | 90.9 | 0.875 | Not selected |
| 3 | model_2 | 90.4 | 0.875 | Not selected |
| 4 | model_5 | 90.2 | 0.876 | Not selected |
| 5 | model_1 | 90.1 | 0.868 | Not selected |

The catalytic core (residues 151–450) achieves a mean pLDDT of 95.9, with all three catalytic triad residues showing exceptional confidence: Ser203 (pLDDT 98.75), His447 (96.06), and Glu334 (97.62). The primary pharmacophoric anchors Trp86 (98.00) and Trp286 (97.69) are also predicted with very high confidence. The mean Predicted Aligned Error (PAE) for the catalytic core is 3.28 Angstroms, confirming that inter-residue spatial relationships are reliably predicted and the structure is suitable for molecular docking.

**Protein Preparation Protocol:** The rank 1 AlphaFold2 model was prepared for AutoDock Vina using a standardised 10-step protocol implemented in PyMOL and AutoDockTools (MGLTools). The experimental crystal structure PDB 4EY7 (human AChE co-crystallised with Donepezil, 2.35 Angstrom resolution) was used as a reference for binding site coordinate verification. Preparation steps included removal of non-standard residues, water molecules, and ligands; retention of chain A only; addition of hydrogen atoms at physiological pH 7.4; assignment of Gasteiger partial charges; and saving in PDBQT format (5,056 ATOM records).

**Docking Box Configuration:** The AutoDock Vina search box (40 x 40 x 40 Angstroms, exhaustiveness = 32) was centred on the CAS region at coordinates (x = -2.895, y = -40.109, z = 30.761), positioned 3.0 Angstroms from the Ser203 Calpha. The box was verified to encompass all six primary pharmacophoric residues: Ser203, His447, Glu334, Trp86 (CAS) and Trp286, Tyr72 (PAS).

---

### Task 6 — Virtual Screening and Lead Identification

**Screening Strategy:** A QSAR-based multi-parameter composite scoring approach was employed to select 10 compounds from a pool of 1,475 potent ChEMBL AChE inhibitors (pIC50 >= 7.0) for molecular docking. This approach integrates multiple orthogonal drug discovery criteria simultaneously and mirrors multi-parameter optimisation (MPO) scoring used in pharmaceutical industry lead identification.

**Drug-Likeness Filtering:** The 1,475-compound potent inhibitor library was filtered using Lipinski's Rule of Five criteria (MW <= 500 Da, logP <= 5.0, RO5 violations <= 1) and a fragment removal filter (MW >= 150 Da). This yielded a final screening library of 660 drug-like compounds (44.7% retention).

**Virtual Screening Funnel:**

| Step | Filter Applied | Compounds Retained | Removed | Retention % |
|------|---------------|-------------------|---------|-------------|
| 0 | Starting pool: pIC50 >= 7.0 | 1,475 | 0 | 100% |
| 1 | Lipinski Ro5: MW <= 500 Da, logP <= 5.0, violations <= 1 | 672 | 803 | 45.6% |
| 2 | Fragment removal: MW >= 150 Da | 660 | 12 | 44.7% |
| 3 | Hydrophilicity filter: logP >= -2.0 | 660 | 0 | 44.7% |
| 4 | Deduplication by ChEMBL ID | 660 | 0 | 44.7% |
| **Final** | **Final screening library** | **660** | — | **44.7%** |

**Composite Scoring Function:**

| Component | Weight | Rationale |
|-----------|--------|-----------|
| Normalised pIC50 | 50% | Primary driver of potency |
| Ligand Efficiency (LE = pIC50 / heavy atom count) | 25% | Rewards potency per atom; prevents selection of oversized molecules |
| Drug-likeness score (Ro5 violations penalty) | 15% | Promotes oral bioavailability |
| logP optimality score (Gaussian, optimal range 1.0-3.5) | 10% | Favours CNS-appropriate lipophilicity for BBB penetration |

**Top 10 Selected Docking Candidates:**

| Rank | ChEMBL ID | Docking ID | pIC50 | IC50 (nM) | MW (Da) | logP | LE | Comp. Score | Scaffold |
|------|-----------|------------|-------|-----------|---------|------|----|-------------|----------|
| 1 | CHEMBL3585780 | ACHE_Cand_01 | 10.57 | 0.027 | 230.3 | 0.55 | 0.622 | 0.9059 | o-MeO-phenyl sulfamide |
| 2 | CHEMBL3585781 | ACHE_Cand_02 | 10.57 | 0.027 | 230.3 | 0.55 | 0.622 | 0.9059 | p-MeO-phenyl sulfamide |
| 3 | CHEMBL3585784 | ACHE_Cand_03 | 10.57 | 0.027 | 260.3 | 0.56 | 0.556 | 0.8719 | 2,6-diMeO-phenyl sulfamide |
| 4 | CHEMBL3585783 | ACHE_Cand_04 | 10.57 | 0.027 | 260.3 | 0.56 | 0.556 | 0.8719 | 3,4-diMeO-phenyl sulfamide |
| 5 | CHEMBL3585782 | ACHE_Cand_05 | 10.57 | 0.027 | 260.3 | 0.56 | 0.556 | 0.8719 | 3,5-diMeO-phenyl sulfamide |
| 6 | CHEMBL4205954 | ACHE_Cand_06 | 10.88 | 0.013 | 404.9 | 4.32 | 0.375 | 0.8169 | THIQone-chlorobenzamide |
| 7 | CHEMBL4209803 | ACHE_Cand_07 | 10.87 | 0.014 | 398.5 | 4.28 | 0.375 | 0.8155 | THIQone-methylbenzamide |
| 8 | CHEMBL4205425 | ACHE_Cand_08 | 10.96 | 0.011 | 434.9 | 4.33 | 0.342 | 0.8100 | THIQone-MeO-chlorobenzamide |
| 9 | CHEMBL3585776 | ACHE_Cand_09 | 10.57 | 0.027 | 364.4 | 2.52 | 0.406 | 0.8091 | p-MeO-sulfonyl-carbamate |
| 10 | CHEMBL3585775 | ACHE_Cand_10 | 10.57 | 0.027 | 364.4 | 2.52 | 0.406 | 0.8091 | o-MeO-sulfonyl-carbamate |

The top 10 compounds represent two distinct chemical scaffolds: the **sulfamide series** (Ranks 1–5 and 9–10; MW 230–365 Da, LE 0.406–0.622, low logP 0.55–2.52) and the **THIQone-benzamide series** (Ranks 6–8; MW 399–435 Da, LE 0.342–0.375, logP 4.28–4.33). Scaffold redundancy among the sulfamide series is acknowledged as a limitation; a Tanimoto similarity filter (Tc < 0.7) was not applied within the scope of this project.

**Output Files:**

| File | Description |
|------|-------------|
| potent_ache_inhibitors.csv | 1,475 potent inhibitors (pIC50 >= 7.0) used as starting pool |
| drug_like_library_660.csv | Final 660-compound drug-like screening library |
| composite_scores.csv | Composite scores for all 660 screened compounds |
| top10_docking_candidates.csv | Final 10 candidates selected for docking with all properties |
| Figure6_Screening_Funnel.png | Virtual screening funnel visualisation |
| Figure7_Composite_Score_Distribution.png | Distribution of composite scores across 660 compounds |
| Figure8_Chemical_Space_Top10.png | MW vs logP chemical space map of screening library with top 10 highlighted |
| Figure9_Ligand_Efficiency_Comparison.png | LE comparison between full library and top 10 candidates |

---

### Task 7 — Molecular Docking

**Ligand Preparation:** SMILES strings for all 10 candidate compounds were sourced from ChEMBL and converted to 3D PDB format using Open Babel / RDKit with hydrogen addition and 3D coordinate generation. PDB structures were subsequently converted to PDBQT format using AutoDockTools (MGLTools), assigning Gasteiger partial charges and defining rotatable bonds for flexible ligand docking.

**Docking Execution:** Molecular docking was performed using AutoDock Vina against the prepared AlphaFold2 AChE receptor (ACHE_rank001_prepared.pdbqt) with the 40 x 40 x 40 Angstrom search box centred on the CAS gorge (exhaustiveness = 32). Binding energies (kcal/mol) for all 10 compounds are reported below.

**Docking Results:**

| Ligand ID | ChEMBL ID | Binding Energy (kcal/mol) |
|-----------|-----------|--------------------------|
| Ligand 1 | CHEMBL3585780 | -7.50 |
| Ligand 2 | CHEMBL3585781 | -7.45 |
| Ligand 3 | CHEMBL3585784 | -7.47 |
| Ligand 4 | CHEMBL3585783 | -7.55 |
| Ligand 5 | CHEMBL3585782 | -8.06 |
| Ligand 6 | CHEMBL4205954 | -13.07 |
| Ligand 7 | CHEMBL4209803 | **-13.67** |
| Ligand 8 | CHEMBL4205425 | -11.43 |
| Ligand 9 | CHEMBL3585776 | -10.01 |
| Ligand 10 | CHEMBL3585775 | -9.87 |

A clear bimodal distribution is observed in the docking scores. Ligands 1–5 (sulfamide series) cluster between -7.45 and -8.06 kcal/mol, representing moderate binding affinity. Ligands 6–8 (THIQone-benzamide series) exhibit markedly superior affinities (-11.43 to -13.67 kcal/mol), indicating a structurally distinct and more complementary interaction with the AChE gorge. Ligands 9 and 10 form an intermediate group (-9.87 to -10.01 kcal/mol).

**Lead Compounds Selected for Interaction Analysis:**

| Rank | ChEMBL ID | Binding Energy (kcal/mol) | Reason for Selection |
|------|-----------|--------------------------|----------------------|
| 1 | CHEMBL4209803 | -13.67 | Best binding affinity overall |
| 2 | CHEMBL4205954 | -13.07 | Strong binding; broad surface complementarity |
| 3 | CHEMBL4205425 | -11.43 | Good binding energy with balanced drug-like profile |

**Protein–Ligand Interaction Analysis (Discovery Studio Visualizer):**

*CHEMBL4209803 (Ligand 7; -13.67 kcal/mol):* The richest directional interaction profile of all candidates. Dual Pi-Pi stacking with TRP A:286 (3.88 Angstroms) and TYR A:341 engages both PAS and gorge mid-section simultaneously. Additional Pi-Sigma and Pi-Alkyl interactions with TRP A:286 and TYR A:72, respectively, together with van der Waals contacts with VAL A:294 and PHE A:295, establish a multi-mode aromatic network that mirrors the pharmacophore of clinically approved donepezil.

*CHEMBL4205954 (Ligand 6; -13.07 kcal/mol):* The interaction profile is composed entirely of van der Waals contacts spanning seven residues (PHE A:297, PHE A:295, TRP A:286, TYR A:341, PHE A:338, ASP A:74, TYR A:124), reflecting broad shape complementarity with the gorge. The absence of directional non-covalent bonds (no hydrogen bonds, no pi-stacking) raises concerns about binding specificity, though the large contact surface accounts for the strong binding energy.

*CHEMBL4205425 (Ligand 8; -11.43 kcal/mol):* A compact but directionally specific profile. Pi-Sigma interaction with TRP A:286 confirms PAS engagement; a carbon hydrogen bond with TYR A:72 provides an additional directional contact at the gorge entrance. Fewer total contacts than Ligands 6 and 7 account for the lower binding energy, but the interactions are specific and well-defined.

**Output Files:**

| File | Description |
|------|-------------|
| docking_results_all10.csv | Binding energies and top pose data for all 10 ligands |
| Figure10_Docking_Score_Comparison.png | Bar chart of binding energies for all 10 compounds |
| Figure11_Ligand07_3D_Pose.png | 3D docking pose of CHEMBL4209803 in AChE gorge |
| Figure12_Ligand07_2D_Interactions.png | 2D interaction diagram for CHEMBL4209803 |
| Figure13_Ligand06_3D_Pose.png | 3D docking pose of CHEMBL4205954 in AChE gorge |
| Figure14_Ligand06_2D_Interactions.png | 2D interaction diagram for CHEMBL4205954 |
| Figure15_Ligand08_3D_Pose.png | 3D docking pose of CHEMBL4205425 in AChE gorge |
| Figure16_Ligand08_2D_Interactions.png | 2D interaction diagram for CHEMBL4205425 |

---

### Task 8 — ADMET and Toxicity Filtering

ADMET profiling was performed on the three lead compounds using SwissADME (pharmacokinetics and drug-likeness), pkCSM (ADMET predictions), and ProTox-3.0 (toxicity prediction).

**SwissADME / Pharmacokinetic Profile:**

| Parameter | Threshold | CHEMBL4205954 (L06) | CHEMBL4209803 (L07) | CHEMBL4205425 (L08) |
|-----------|-----------|---------------------|---------------------|---------------------|
| MW (Da) | < 500 Da | 404.85 | 398.45 | 434.87 |
| Consensus LogP | < 5 | 3.56 | 3.68 | 3.60 |
| H-Bond Donors | < 5 | 1 | 1 | 1 |
| H-Bond Acceptors | < 10 | 3 | 3 | 4 |
| TPSA (A2) | < 90 A2 | 66.48 | 66.48 | 75.71 |
| BBB permeant | Yes | Yes | Yes | Yes |
| Pgp substrate | No | No | No | No |
| CYP2D6 inhibitor | No | No | No | No |
| CYP3A4 inhibitor | No | No | Yes | Yes |
| Lipinski compliance | Pass | Full pass | Full pass | Full pass |

All three compounds pass the Lipinski Rule of Five without a single violation, and critically all three are BBB permeant with TPSA values of 66.48–75.71 A2, well below the 90 A2 CNS penetration threshold. None are Pgp substrates, meaning P-glycoprotein efflux will not reduce CNS drug exposure. The shared CYP3A4 inhibition across Ligands 7 and 8 is a scaffold-level liability that will require structural optimisation in the hit-to-lead phase.

**Toxicity Profile (pkCSM and ProTox-3.0):**

| Parameter | Tool | Criterion | L06 (CHEMBL4205954) | L07 (CHEMBL4209803) | L08 (CHEMBL4205425) |
|-----------|------|-----------|---------------------|---------------------|---------------------|
| AMES mutagenicity | pkCSM | Non-mutagenic | Negative | Negative | Negative |
| hERG I inhibitor | pkCSM | Non-inhibitor | No | No | No |
| hERG II inhibitor | pkCSM | Non-inhibitor | Yes | Yes | Yes |
| Hepatotoxicity | ProTox | Non-hepatotoxic | Yes | Yes | Yes |
| Skin sensitization | pkCSM | Non-sensitizer | No | No | No |
| Max. tolerated dose (log mg/kg/day) | pkCSM | Higher preferred | 0.457 | 0.069 | 0.036 |
| Predicted LD50 rat oral (mg/kg) | ProTox-3.0 | > 1000 mg/kg | 1600 | 1600 | 2000 |
| WHO toxicity class | ProTox-3.0 | Class IV–V | Class IV | Class IV | Class IV |

The shared hERG II inhibition and hepatotoxicity flags across all three compounds reflect a class effect of the THIQone-benzamide scaffold core rather than compound-specific liabilities. The key differentiator is the maximum tolerated dose (MTD): Ligand 06 has the highest MTD (0.457 log mg/kg/day), Ligand 07 is intermediate (0.069), and Ligand 08 is the most restricted (0.036) — a critical liability for chronic daily dosing in elderly Alzheimer's patients.

**Output Files:**

| File | Description |
|------|-------------|
| swissadme_results.csv | Full SwissADME pharmacokinetic output for three leads |
| pkcsm_results.csv | pkCSM ADMET prediction results |
| protox3_results.csv | ProTox-3.0 toxicity prediction results |
| Figure17_ADMET_Summary.png | ADMET property comparison across three lead compounds |

---

### Task 9 — Multi-Parameter Lead Prioritisation

A weighted multi-parameter scoring system was applied to rank the three lead compounds across seven independent criteria, integrating docking affinity, pharmacokinetics, drug-likeness, and safety.

**Multi-Parameter Scoring Matrix:**

| Parameter | Weight | CHEMBL4205954 (L06) | CHEMBL4209803 (L07) | CHEMBL4205425 (L08) |
|-----------|--------|---------------------|---------------------|---------------------|
| Binding energy | 25% | 2 → 0.50 | 3 → 0.75 | 2 → 0.50 |
| BBB penetration | 20% | 3 → 0.60 | 3 → 0.60 | 3 → 0.60 |
| Oral bioavailability | 15% | 3 → 0.45 | 3 → 0.45 | 3 → 0.45 |
| Lipinski compliance | 10% | 3 → 0.30 | 3 → 0.30 | 3 → 0.30 |
| hERG safety | 15% | 1 → 0.15 | 1 → 0.15 | 1 → 0.15 |
| AMES mutagenicity | 10% | 3 → 0.30 | 3 → 0.30 | 3 → 0.30 |
| Hepatotoxicity (MTD-differentiated) | 5% | 1 → 0.05 | 1 → 0.05 | 0 → 0.00 |
| **Weighted Total Score** | **100%** | **2.35 — Rank 2** | **2.60 — Rank 1** | **2.30 — Rank 3** |

**CHEMBL4209803 — Primary Lead (Rank 1, Score: 2.60/3.00)**

CHEMBL4209803 is selected as the primary lead compound on the basis of its superior binding affinity (-13.67 kcal/mol), the lowest molecular weight of the series (398.45 Da, providing maximum lead-optimisation headroom), confirmed BBB permeability (TPSA 66.48 A2), full Lipinski compliance, clean skin sensitisation profile, and AMES-negative status. Its dual-site binding mode — engaging TRP A:286 at the PAS and aromatic gorge residues through a precise Pi-Pi stacking network — mirrors the pharmacophore of clinically approved donepezil, conferring strong mechanistic confidence. The shared hERG II inhibition and hepatotoxicity flags are scaffold-level class effects and do not differentiate this compound from the others; its MTD of 0.069 log mg/kg/day is intermediate. CHEMBL4209803 is the recommended candidate for next-stage hit-to-lead optimisation, including in vitro AChE IC50 determination and hepatocyte toxicity assays.

**CHEMBL4205954 — Secondary Lead (Rank 2, Score: 2.35/3.00)**

CHEMBL4205954 is retained as the secondary lead. It ranks second due to a lower binding energy score relative to the primary lead, but compensates with the highest maximum tolerated dose of the series (MTD: 0.457 log mg/kg/day) and an identical TPSA/BBB profile (TPSA 66.48 A2). Its PAS-biased binding mode — broad van der Waals contacts across seven gorge residues — represents a structurally distinct backup to the primary lead, supporting portfolio diversification. If CHEMBL4209803 encounters attrition in subsequent in vitro hepatotoxicity or cardiac safety assays, CHEMBL4205954 is the immediate next candidate.

**CHEMBL4205425 — Deprioritised (Rank 3, Score: 2.30/3.00)**

CHEMBL4205425 was ranked third and deprioritised based on the convergence of multiple liabilities: the highest molecular weight of the series (434.87 Da), the highest TPSA (75.71 A2), the lowest maximum tolerated dose (MTD: 0.036 log mg/kg/day) — less than one-tenth that of the secondary lead — and the lowest weighted total score (2.30/3.00). The exceptionally narrow human tolerability window is a critical concern for a drug intended for chronic daily administration in elderly patients. CHEMBL4205425 is retained for scaffold analysis to guide future bioisosteric replacement strategies targeting its chloro-substituted aromatic core.

**Output Files:**

| File | Description |
|------|-------------|
| multiparameter_scoring.csv | Full multi-parameter score matrix for three lead compounds |
| Figure18_Lead_Prioritisation_Radar.png | Radar chart comparing three leads across all seven parameters |

---

## Final Lead Summary

| Rank | ChEMBL ID | Binding Energy (kcal/mol) | MTD (log mg/kg/day) | BBB | Lipinski | Recommendation |
|------|-----------|--------------------------|---------------------|-----|----------|----------------|
| 1 | CHEMBL4209803 | -13.67 | 0.069 | Yes | Pass | **Primary lead — advance to hit-to-lead** |
| 2 | CHEMBL4205954 | -13.07 | 0.457 | Yes | Pass | **Secondary lead — backup candidate** |
| 3 | CHEMBL4205425 | -11.43 | 0.036 | Yes | Pass | Deprioritised — scaffold analysis only |

---

## Tools and Technologies

| Tool / Library | Version | Purpose |
|----------------|---------|---------|
| Python | 3.x | Data processing and cheminformatics |
| RDKit | Latest | Descriptor calculation, fingerprint generation, structure standardisation |
| R | 4.x | QSAR model training and statistical analysis |
| randomForest (R) | Latest | Random Forest regression model |
| ChEMBL Database | v33 | Bioactivity data source |
| ColabFold | v1.5 | AlphaFold2 protein structure prediction |
| PyMOL | 3.x | Protein structure preparation and visualisation |
| AutoDock Vina | Latest | Molecular docking |
| MGLTools / AutoDockTools | Latest | Receptor preparation (PDBQT format, charge assignment) |
| Open Babel | Latest | Ligand 3D coordinate generation and format conversion |
| Discovery Studio Visualizer | Latest | Docking pose visualisation and interaction analysis |
| SwissADME | Web | Pharmacokinetic and drug-likeness profiling |
| pkCSM | Web | ADMET property prediction |
| ProTox-3.0 | Web | Toxicity and LD50 prediction |

---

## Limitations

Several methodological limitations are acknowledged across the pipeline:

- **QSAR models:** Both Linear Regression and Random Forest models fall below the accepted R2 threshold of 0.60 on the test set. The descriptor set (six physicochemical descriptors) is sparse relative to the chemical diversity of 5,604 compounds. Extended fingerprint representations (e.g., full 2048-bit Morgan fingerprints as model features) and cross-validated hyperparameter optimisation would likely improve performance substantially.

- **Static receptor docking:** AlphaFold2 predicts a static ground-state structure and does not capture receptor flexibility or induced-fit effects. For the conformationally flexible AChE gorge, ensemble docking or flexible-residue docking (particularly for Phe295) would be preferable in a production pipeline.

- **Scaffold redundancy:** Seven of the ten selected docking candidates belong to the same sulfamide scaffold series. A Tanimoto similarity filter (Tc < 0.7) was not applied due to project scope constraints. In a real pipeline, scaffold clustering would be performed before final candidate selection.

- **Assay heterogeneity:** IC50 values originate from multiple laboratories and assay formats in ChEMBL, which may introduce systematic variability not fully addressed by the curation pipeline. The sub-picomolar pIC50 values (10.57) reported for the sulfamide series likely originate from a single assay publication and may reflect assay-specific artefacts.

- **ADMET scope:** Toxicity, metabolic stability, BBB penetration, and hERG liability profiling were performed in silico only, using predictive web tools. All ADMET predictions require experimental validation before advancing candidates to in vitro studies.

- **Docking box placement:** The search box centre was defined visually using AutoDockTools rather than by formal structural superimposition of the AlphaFold model onto PDB 4EY7. Formal alignment would yield more precise binding site transfer coordinates.

---

## Citation

If referencing this work, please cite as:

> Basharat N., Munir G., Fahad H., Saleem A. (2026). *AChE Inhibitor Discovery Pipeline: A Computational Drug Discovery Study using ChEMBL, QSAR Modelling, and Molecular Docking.* Pharmacoinformatics Final Project, BS Bioinformatics UG-1, NUST SEECS.

---

## References

- AutoDock Vina. https://vina.scripps.edu/
- MGLTools / AutoDockTools. https://ccsb.scripps.edu/mgltools/
- Open Babel. https://openbabel.org/
- Discovery Studio — Dassault Systèmes. https://www.3ds.com/products/biovia/discovery-studio
- RCSB Protein Data Bank. https://www.rcsb.org/
- ChEMBL Database. https://www.ebi.ac.uk/chembl/
- SwissADME. https://www.swissadme.ch/
- pkCSM. https://biosig.lab.uq.edu.au/pkcsm/
- ProTox-3.0. https://tox.charite.de/protox3/
- ColabFold / AlphaFold2. https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
- NAMRA-5/AChE-Drug-Discovery-Pipeline. GitHub. https://github.com/Namra-5/AChE-Drug-Discovery-Pipeline