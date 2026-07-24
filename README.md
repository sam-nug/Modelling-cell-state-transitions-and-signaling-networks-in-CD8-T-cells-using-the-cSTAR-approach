# Modelling-cell-state-transitions-and-signaling-networks-in-CD8-T-cells-using-the-cSTAR-approach
Modelling cell state transitions and singaling networks in CD8+ T cells using the cSTAR approach
# Network Inference of CD8+ T-Cell Drug Perturbations

Reconstruction of a signalling network from CD8+ T-cell drug-perturbation
single-cell RNA-seq data using a three-stage computational pipeline:
**cSTAR → BMRA → MRA**.

## Overview

The pipeline infers direct connections between seven signalling modules that are perturbed by 12 targeted drugs. Each of the signallings are as a response of CD8+ T cells to the drugs. Modules are selected based on their dynamic phenotype descriptor (DPD), per-drug activity is estimated from IC50 values. Each of the pathway-activity coefficients are learned from cSTAR, and the module 
network is reconstructed using Modular Response Analysis (MRA) across three estimators (Bayesian regression, OLS, ODR).
Two phenotype descriptors are used: the DPD score dot product of stim vs unstim and a BTLA score. Both are regressed onto the module network in the final MRA stage via FPTU

The seven modules selected:
- **Modules:** ARFGAP, BCL, IGF1R, JAK, MTOR, NAE, SERCA

The 12 drugs taken:
- **Drugs (12):** QS-11, navitoclax, BMS-536924, AT9283, TG-101348, CYT-387,
  ruxolitinib, Deforolimus, Temsirolimus, Sapanisertib, Pevonedistat,
  Thapsigargin

## Data

Source: GEO accession **GSE306429** (CD8+ T cells only), DeMeo et al. 2025.
Download the series and place the derived per-compound expression / IC50 files
in the working directory.
Two other file are included BTLAvsTCR_24h.csv and cd8_limma_merged_filtered_targets_ic50_dpd.csv

## Pipeline stages

1. **Preparation** — Compute the DPD (dot product of the state-transition
   vector between stimulated and unstimulated cells); take the top-4 and
   bottom-4 modules (plus co-targeted drugs) in this case 7 modules, 12 drugs. Build the
   activity (g) matrix from dose and IC50:
   - Inhibitor: `g = 1 / (1 + I/K)`
   - Activator: `g = (1 + γ·I/K) / (1 + I/K)`

2. **Phenotype Descriptors** - DPD and BTLA DPD are calculated by the logfc against the Szabo reference vector using log2foldchange for DPD and stat for BTLA. BTLA value is only used during MRA. The top 4 and bottom 4 compounds based on DPD were selected.

3. **cSTAR pathway activity** — `predict_coeffs(...)` with sub-sampling set to
   `200,000` learns the α coefficients; the global response matrix is then
   computed by cSTAR integration.

4. **BMRA** — `run_bmra` (Bayesian regression) estimates the local response

5. **MRA** — `run_mra`
   matrix `r` with the diagonal normalised to `−1`; `run_mra_rp` (OLS/ODR) is
   used for cross-estimator comparison. The connection matrix and its
   minus-inverse (`−r⁻¹`) give network topology and the predicted
   systems-level drug response. The connection and minus-inverse both included modules + DPD + BTLA fitted onto the network. 
   Edges thresholds were exported to Cytoscape edge lists.

## Notebook

### Manual assignments to fill in before running

Module/drug selection, DPD, and BTLA rows are built automatically. Two
hardcoded lookups still need to be kept up to date by hand when compounds
change:

- **`compound_to_module`** (dict, prep section) — maps each compound to the
  module it perturbs. Add an entry here for any new drug.
- **`IC50_OVERRIDES_NM`** (dict, activity-matrix cell) — manual IC50 (nM)
  overrides for compounds whose `IC50_nM` field can't be parsed automatically
  (currently `Dorsomorphin`, `Apicidin`, `Serdemetan`).


## Outputs (`02_outputs/`)

- `dpd_sum_per_compound_raw_btla.csv`,`btla_deg_per_compound_btla.csv`—phenotype scores
- `top4_bottom4_selected_modules_drugs_ic50_btla.csv` — selected module/drug/IC50 table
- `a_coeffs_df_top4_bottom4_btla.csv`, `pathway_activity_top4_bottom4_btla.csv`
- `R_global_core_top4_bottom4_btla.csv`, `y_true_top4_bottom4_btla.csv`, `Data_norm.csv`
- `pert_matrix_top4_bottom4_btla.csv`
- `r_mean_ols_rp_fptu_prefinal_btla.csv`, `r_minv_ols_rp_fptu_prefinal_btla.csv` — combined modules+DPD+BTLA connection matrix and minus-inverse
- `CD8_r_ols_btla_net.txt` — Cytoscape-style edge list (From/To/Strength)
- `module_network_ols.png` (and`_br`/`_odr`variants) — network diagrams


## Author

Samuel Nugent
