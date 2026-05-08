# Acetylcholinesterase (AChE) Inhibitor Discovery Pipeline
## Computational Drug Discovery using ChEMBL, QSAR Modelling, and Molecular Docking

**Course:** Pharmacoinformatics — Final Project
**Institution:** National University of Sciences and Technology (NUST)
**School:** School of Interdisciplinary Engineering and Sciences (SINES)
**Programme:** BS Bioinformatics, UG-1
**Date:** 8th May 2026

---

## Authors

| Name | CMS ID |
|------|--------|
| Namra Basharat | 476203 |
| Ghania | 460673 |
| Hania Fahad | 455132 |
| Alishba Saleem | 479300 |

---

## Project Overview

This repository contains the complete computational drug discovery pipeline developed for the Pharmacoinformatics final project. The pipeline targets **Acetylcholinesterase (AChE; EC 3.1.1.7; ChEMBL Target: CHEMBL220)**, a clinically validated therapeutic target for Alzheimer's disease (AD).

The work spans six interconnected tasks: bioactivity data acquisition and curation, molecular descriptor and fingerprint generation, QSAR modelling, structure-activity relationship (SAR) analysis, protein structure prediction and preparation, and multi-parameter virtual screening with molecular docking.

Alzheimer's disease is the most prevalent neurodegenerative disorder worldwide, affecting approximately 55 million people globally. It is characterised by progressive cognitive decline and memory loss, driven in part by the degeneration of cholinergic neurons and the resulting deficit in acetylcholine (ACh) at synapses. AChE is the enzyme responsible for hydrolysing ACh, thereby terminating cholinergic neurotransmission. Inhibiting AChE prolongs synaptic ACh availability and temporarily compensates for the cholinergic deficit — representing the most clinically validated pharmacological strategy for symptomatic AD treatment. Approved AChE inhibitors include Donepezil (1996), Rivastigmine (1997), and Galantamine (2001).

The pipeline integrates cheminformatics, machine learning, structural bioinformatics, and molecular docking to identify and characterise potent, drug-like AChE inhibitors from the ChEMBL bioactivity database. All computational steps are documented in this repository with associated data files, scripts, figures, and reports.

---

## Key Statistics

| Metric | Value |
|--------|-------|
| Raw bioactivity records retrieved | 9,731 |
| Final curated dataset | 5,604 unique compounds |
| Data reduction | 42.5% (4,127 records removed) |
| Potent inhibitors (pIC50 ≥ 7; IC50 ≤ 100 nM) | 1,481 |
| Very potent inhibitors (pIC50 ≥ 8; IC50 ≤ 10 nM) | 584 |
| Sub-nanomolar inhibitors (pIC50 ≥ 9; IC50 < 1 nM) | 170 |
| Most potent compound | CHEMBL4205425 (pIC50 = 10.96; IC50 = 0.011 nM) |
| AlphaFold2 model confidence (mean pLDDT) | 91.2 (Very High) |
| Virtual screening library (after Ro5 filtering) | 660 drug-like compounds |
| Best docking binding energy | -13.67 kcal/mol (CHEMBL4209803) |

---

## Repository Structure

```
AChE-Drug-Discovery-Pipeline/
|
|-- task-1-target-selection/
|   |-- target_selection.R
|   |-- target_summary.csv
|   |-- LITERATURE_REVIEW.txt
|   `-- ache_structure.png
|
|-- task-2-bioactivity-curation/
|   |-- bioactivity_curation.R
|   |-- raw_ache_data.csv
|   `-- ACHE_QSAR_dataset.csv
|
|-- task-3-qsar/
|   |-- cleaned_ache_data.csv
|   |-- ache_with_descriptors.csv
|   |-- feature_matrix.csv
|   |-- ACHE_QSAR_dataset.csv
|   |-- QSAR_Modeling.Rmd
|   |-- QSAR_Modeling.html
|   |-- Figure1_Predicted_vs_Observed.png
|   `-- Figure2_Residual_Error_Plot.png
|
|-- task-4-sar/
|   |-- Figure3_SAR_Boxplots.png
|   |-- Figure4_Chemical_Space_PCA.png
|   `-- Figure5_RF_Variable_Importance.png
|
|-- task-5-protein-structure-prep/
|   `-- [AlphaFold2 model, preparation files, docking configuration]
|
|-- task-6-virtual-screening/
|   `-- [Composite scoring results, docking outputs, lead compound files]
|
`-- README.md
```

---

## Pipeline Description

### Task 1 — Target and Ligand Set Definition

The therapeutic target selected for this study is Human Acetylcholinesterase (AChE), encoded by the ACHE gene (UniProt: P22303), with ChEMBL target identifier CHEMBL220. AChE was selected on the basis of its clinical and pharmacological validation as the primary target for Alzheimer's disease therapy. The enzyme presents a deep catalytic gorge of approximately 20 Ångstroms in depth, containing two key pharmacophoric subsites: the Catalytic Anionic Site (CAS), which houses the catalytic triad (Ser203–His447–Glu334), and the Peripheral Anionic Site (PAS) at the gorge entrance, lined by residues Trp286 and Tyr72. Both sites are exploited by dual-site inhibitors such as Donepezil.

Bioactivity data were retrieved from the ChEMBL database (version 33) using the IC50 activity type, organism filter Homo sapiens, and target type Single Protein. The initial retrieval returned 9,731 bioactivity records.

---

### Task 2 — Bioactivity Data Acquisition and Curation

The raw dataset of 9,731 records was subjected to a four-step curation pipeline to ensure data quality and consistency prior to modelling.

**Step 1 — Deduplication and Record Filtering:** Duplicate records arising from the same compound measured in different assays were removed. Records with ambiguous or missing activity values and those with non-standard relations (e.g., >, <) were also excluded. Approximately 2,213 records were removed at this stage.

**Step 2 — Activity Type Standardisation:** Only IC50 measurements were retained to ensure a consistent activity scale across all compounds. Records with non-standard activity types (Ki, Kd, percentage inhibition) and those with pChEMBL values below 4.0 (corresponding to IC50 > 100 micromolar, biologically insignificant) were removed. Approximately 1,914 records were excluded.

**Step 3 — Structure Standardisation:** Valid SMILES strings were verified for all retained compounds. Salt stripping was performed using RDKit's SaltRemover module, and canonical SMILES were generated using RDKit's MolStandardize function.

**Step 4 — Activity Scale Conversion:** IC50 values in nanomolar units were converted to the negative logarithmic pIC50 scale using the formula pIC50 = -log10(IC50 × 10⁻⁹). This transformation produces a linear scale suitable for regression-based QSAR modelling.

The final curated dataset comprises 5,604 unique compounds with IC50 measurements against human AChE.

| Curation Stage | Records | Notes |
|----------------|---------|-------|
| Raw ChEMBL download | 9,731 | All IC50 records for CHEMBL220 |
| After deduplication | ~7,518 | Duplicates removed |
| After activity filtering | ~5,604 | IC50 only, valid pChEMBL ≥ 4.0 |
| After structure standardisation | 5,604 | Final curated dataset |

**Potency Distribution of the Curated Dataset:**

| Category | pIC50 Range | IC50 Equivalent | Count | Percentage |
|----------|-------------|-----------------|-------|------------|
| Weak | < 6 | > 1 µM | 2,028 | 36.2% |
| Moderate | 6 – 7 | 100 nM – 1 µM | 2,095 | 37.4% |
| Potent | 7 – 8 | 10 – 100 nM | 897 | 16.0% |
| Very Potent | 8 – 9 | 1 – 10 nM | 414 | 7.4% |
| Sub-nanomolar | ≥ 9 | < 1 nM | 170 | 3.0% |
| **Total** | | | **5,604** | **100%** |

---

### Task 3 — Molecular Descriptor Calculation and QSAR Modelling

**Descriptor Computation:** Six physicochemical descriptors were computed using RDKit for all 5,604 compounds: Molecular Weight (MW), Lipophilicity (AlogP), Hydrogen Bond Donors (HBD), Hydrogen Bond Acceptors (HBA), Rotatable Bonds, and Topological Polar Surface Area (TPSA). Lipinski's Rule of Five compliance was evaluated to assess drug-likeness; 59.5% of compounds were found to be fully RO5-compliant, indicating good oral bioavailability potential.

**Morgan Fingerprint Generation:** 2048-bit Morgan circular fingerprints (radius 2, equivalent to ECFP4) were generated for all compounds using RDKit's GetMorganFingerprintAsBitVect function. These fingerprints encode circular atomic neighbourhoods and are widely used in virtual screening and QSAR modelling. Fingerprints are stored in binary NumPy format (fingerprints.npy).

**QSAR Feature Matrix:** A QSAR-ready feature matrix of dimensions 5,604 × 7 (six physicochemical descriptors plus pIC50) was constructed. All descriptors were scaled to zero mean and unit variance using StandardScaler prior to model training. This matrix is provided as feature_matrix.csv.

**Dataset Splitting:** A stratified scaffold split (70/15/15) was employed: training set (3,897 compounds), validation set (835 compounds), and test set (834 compounds). Stratification was performed by binning pIC50 values into five equal-frequency groups prior to splitting, ensuring proportional representation of the full biological activity range (pIC50: 4.00–10.96) in all subsets.

**Models Trained:**

*Linear Regression (Baseline):* A multiple linear regression model was trained using R's lm() function on the six-descriptor feature matrix without regularisation. This model serves as an interpretable baseline.

*Random Forest (Non-linear):* A Random Forest regression model was configured with 300 trees (ntree = 300) and default mtry parameter. Random Forest models complex non-linear interactions between descriptors without requiring explicit feature engineering.

| Model | Split | R² | RMSE | MAE |
|-------|-------|----|------|-----|
| Linear Regression | Train | 0.073 | 1.266 | 1.020 |
| Linear Regression | Test | 0.074 | 1.247 | 1.002 |
| Random Forest | Train | 0.789 | 0.689 | 0.541 |
| Random Forest | Test | 0.267 | 1.110 | 0.864 |

The near-identical training and test metrics for Linear Regression confirm underfitting — the model explains only 7% of variance in pIC50 and lacks the capacity to capture non-linear structure-activity relationships. The Random Forest substantially outperforms the baseline (test R²: 0.267 vs 0.074), confirming non-linearity in the data, but shows a large gap between training and test R² consistent with overfitting. Both models would benefit from richer fingerprint representations and hyperparameter optimisation in a production pipeline.

**Output Files:**

| File | Description |
|------|-------------|
| cleaned_ache_data.csv | Curated ChEMBL dataset (5,604 compounds) |
| ache_with_descriptors.csv | Dataset with computed RDKit descriptors |
| feature_matrix.csv | Standardised QSAR feature matrix (5,604 × 7) |
| ACHE_QSAR_dataset.csv | Full QSAR-ready dataset |
| QSAR_Modeling.Rmd | R Markdown script for model training and evaluation |
| QSAR_Modeling.html | Rendered HTML report |
| Figure1_Predicted_vs_Observed.png | Predicted vs Observed pIC50 diagnostic plot |
| Figure2_Residual_Error_Plot.png | Residual error diagnostic plot |

---

### Task 4 — Structure-Activity Relationship (SAR) Analysis

SAR analysis was performed to identify which physicochemical properties are most strongly associated with AChE inhibitory potency. The curated dataset was divided at the median pIC50 (6.00) into two groups: High Activity (pIC50 ≥ 6.00; n = 2,807) and Low Activity (pIC50 < 6.00; n = 2,759). Four descriptors were compared between groups using the non-parametric Wilcoxon rank-sum test, as descriptor distributions are non-normal.

| Descriptor | High Activity Mean | Low Activity Mean | p-value | Significant |
|------------|:-----------------:|:----------------:|:-------:|:-----------:|
| Molecular Weight (Da) | 219.70 | 198.74 | < 0.0001 | Yes |
| logP | 4.929 | 4.247 | < 0.0001 | Yes |
| H-Bond Donors | 4.228 | 4.145 | 0.0048 | Yes |
| H-Bond Acceptors | 4.304 | 4.263 | 0.0426 | Yes |

All four descriptors show statistically significant differences between activity groups. High-activity compounds are notably larger (MW +21 Da) and more lipophilic (logP +0.68), consistent with the hydrophobic, gorge-shaped AChE binding site lined with aromatic residues (Trp84, Trp279, Phe330). MW and logP are the most strongly discriminating descriptors (p < 0.0001), while HBD and HBA show smaller but statistically significant differences, suggesting that H-bonding capacity to catalytic triad residues (Ser200, His440, Glu327) contributes modestly to potency.

**Chemical Space Visualisation (PCA):** PCA was performed on a five-descriptor matrix (MW, logP, HBD, HBA, RO5 violations), standardised prior to analysis. The first two principal components together capture 88.4% of total descriptor variance (PC1: 55.9%, PC2: 32.5%), providing a faithful 2D representation of the chemical space. A weak but visible gradient from low to high pIC50 is observed along the PC1 axis, which is dominated by MW and HBA loadings. The three activity ellipses (low, medium, high) show substantial overlap, explaining the difficulty of both QSAR models in achieving high R².

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
|:----:|-------|:----------:|:---:|--------|
| 1 | model_3 | 91.2 | 0.881 | ✅ Selected |
| 2 | model_4 | 90.9 | 0.875 | Not selected |
| 3 | model_2 | 90.4 | 0.875 | Not selected |
| 4 | model_5 | 90.2 | 0.876 | Not selected |
| 5 | model_1 | 90.1 | 0.868 | Not selected |

The catalytic core (residues 151–450) achieves a mean pLDDT of 95.9, with all three catalytic triad residues showing exceptional confidence: Ser203 (pLDDT 98.75), His447 (96.06), and Glu334 (97.62). The primary pharmacophoric anchors Trp86 (98.00) and Trp286 (97.69) are also predicted with very high confidence. The mean Predicted Aligned Error (PAE) for the catalytic core is 3.28 Ångstroms, confirming that inter-residue spatial relationships are reliably predicted and the structure is suitable for molecular docking.

**Protein Preparation Protocol:** The rank 1 AlphaFold2 model was prepared for AutoDock Vina using a standardised 10-step protocol implemented in PyMOL and AutoDockTools (MGLTools). The experimental crystal structure PDB 4EY7 (human AChE co-crystallised with Donepezil, 2.35 Ångstrom resolution) was used as a reference for binding site coordinate verification. Preparation steps included removal of non-standard residues, water molecules, and ligands; retention of chain A only; addition of hydrogen atoms at physiological pH 7.4; assignment of Gasteiger partial charges; and saving in PDBQT format (5,056 ATOM records).

**Docking Box Configuration:** The AutoDock Vina search box (40 × 40 × 40 Ångstroms, exhaustiveness = 32) was centred on the CAS region at coordinates (x = -2.895, y = -40.109, z = 30.761), positioned 3.0 Ångstroms from the Ser203 Cα. The box was verified to encompass all six primary pharmacophoric residues: Ser203, His447, Glu334, Trp86 (CAS) and Trp286, Tyr72 (PAS).

---

### Task 6 — Virtual Screening and Lead Identification

**Screening Strategy:** A QSAR-based multi-parameter composite scoring approach was employed to select 10 compounds from a pool of 1,475 potent ChEMBL AChE inhibitors (pIC50 ≥ 7.0) for molecular docking. This approach integrates multiple orthogonal drug discovery criteria simultaneously and mirrors multi-parameter optimisation (MPO) scoring used in pharmaceutical industry lead identification.

**Drug-Likeness Filtering:** The 1,475-compound potent inhibitor library was filtered using Lipinski's Rule of Five criteria (MW ≤ 500 Da, logP ≤ 5.0, RO5 violations ≤ 1) and a fragment removal filter (MW ≥ 150 Da). This yielded a final screening library of 660 drug-like compounds (44.7% retention).

**Composite Scoring Function:**

| Component | Weight | Rationale |
|-----------|:------:|-----------|
| Normalised pIC50 | 50% | Primary driver of potency |
| Ligand Efficiency (LE = pIC50 / heavy atom count) | 25% | Rewards potency per atom; prevents selection of oversized molecules |
| Drug-likeness score (Ro5 violations penalty) | 15% | Promotes oral bioavailability |
| logP optimality score (Gaussian, optimal range 1.0–3.5) | 10% | Favours CNS-appropriate lipophilicity for BBB penetration |

**Top 10 Selected Docking Candidates:**

| Rank | ChEMBL ID | Docking ID | pIC50 | IC50 (nM) | MW (Da) | logP | LE | Comp. Score | Scaffold |
|:----:|-----------|------------|:-----:|:---------:|:-------:|:----:|:--:|:-----------:|----------|
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

**Molecular Docking Results:**

| Ligand ID | ChEMBL ID | Binding Energy (kcal/mol) |
|-----------|-----------|:------------------------:|
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

**Lead Compounds Selected for Further Analysis:**

| Rank | ChEMBL ID | Binding Energy (kcal/mol) | Reason for Selection |
|:----:|-----------|:------------------------:|----------------------|
| 1 | CHEMBL4209803 | -13.67 | Best binding affinity overall |
| 2 | CHEMBL4205954 | -13.07 | Strong binding with consistent interactions |
| 3 | CHEMBL4205425 | -11.43 | Good binding energy with balanced drug-like profile |

Compounds CHEMBL4209803, CHEMBL4205954, and CHEMBL4205425, all belonging to the THIQone-benzamide scaffold series, were selected as final lead candidates on the basis of their substantially lower binding energies compared to all other screened compounds. These three compounds were further subjected to binding interaction analysis and ADMET profiling.

---

## Tools and Technologies

| Tool / Library | Version | Purpose |
|----------------|:-------:|---------|
| Python | 3.x | Data processing and cheminformatics |
| RDKit | Latest | Descriptor calculation, fingerprint generation, structure standardisation |
| R | 4.x | QSAR model training and statistical analysis |
| randomForest (R) | Latest | Random Forest regression model |
| ChEMBL Database | v33 | Bioactivity data source |
| ColabFold | v1.5 | AlphaFold2 protein structure prediction |
| PyMOL | 3.x | Protein structure preparation and visualisation |
| AutoDock Vina | Latest | Molecular docking |
| MGLTools / AutoDockTools | Latest | Receptor preparation (PDBQT format, charge assignment) |

---

## Limitations

Several methodological limitations are acknowledged across the pipeline:

- **QSAR models:** Both Linear Regression and Random Forest models fall below the accepted R² threshold of 0.60 on the test set. The descriptor set (six physicochemical descriptors) is sparse relative to the chemical diversity of 5,604 compounds. Extended fingerprint representations (e.g., full 2048-bit Morgan fingerprints as model features) and cross-validated hyperparameter optimisation would likely improve performance substantially.

- **Static receptor docking:** AlphaFold2 predicts a static ground-state structure and does not capture receptor flexibility or induced-fit effects. For the conformationally flexible AChE gorge, ensemble docking or flexible-residue docking would be preferable in a production pipeline.

- **Scaffold redundancy:** Seven of the ten selected docking candidates belong to the same sulfamide scaffold series. A Tanimoto similarity filter (Tc < 0.7) was not applied due to project scope constraints. In a real pipeline, scaffold clustering would be performed before final candidate selection.

- **Assay heterogeneity:** IC50 values originate from multiple laboratories and assay formats in ChEMBL, which may introduce systematic variability not fully addressed by the curation pipeline.

- **No ADMET filtering at screening stage:** Toxicity, metabolic stability, BBB penetration, and hERG liability profiling were not applied at the virtual screening stage and are deferred to downstream analysis.

---

## Citation

If referencing this work, please cite as:

> Basharat N., Ghania, Fahad H., Saleem A. (2026). *AChE Inhibitor Discovery Pipeline: A Computational Drug Discovery Study using ChEMBL, QSAR Modelling, and Molecular Docking.* Pharmacoinformatics Final Project, BS Bioinformatics UG-1, NUST SINES.