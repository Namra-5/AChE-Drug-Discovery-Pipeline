#!/bin/bash
# ============================================================
# AChE Batch Docking Script -- Top 10 Candidates
# Requires: AutoDock Vina, OpenBabel (obabel)
# Usage:  bash Task6_run_docking.sh
# Output: docking_results/ folder
# NOTE: Top 5 compounds share a methoxyphenyl-sulfamide scaffold.
# In a production pipeline, a Tanimoto diversity filter (cutoff < 0.7)
# would be applied before this step. All 10 are retained per project scope.
# ============================================================

RECEPTOR="ACHE_rank001_prepared.pdbqt"
CONFIG="vina_config.txt"        # from Task 5 outputs
SMILES_FILE="Task6_top10_smiles.smi"
RESULTS_DIR="docking_results"

mkdir -p $RESULTS_DIR

echo "======================================================"
echo " AChE Docking Pipeline -- Top 10 Candidates"
echo " Receptor: $RECEPTOR"
echo " Config:   $CONFIG"
echo "======================================================"

while IFS=$\'\\t\' read -r smiles name; do
    echo ""
    echo "Processing: $name"

    # Step 1: Generate 3D conformer from SMILES (OpenBabel)
    obabel -:"$smiles" --gen3d -O "${RESULTS_DIR}/${name}.sdf" \
           -h --minimize --ff MMFF94 2>/dev/null

    # Step 2: Convert SDF -> PDBQT (adds rotatable bonds)
    obabel "${RESULTS_DIR}/${name}.sdf" \
           -O "${RESULTS_DIR}/${name}.pdbqt" -xh 2>/dev/null

    # Step 3: Dock with AutoDock Vina
    vina --config  $CONFIG \
         --ligand   "${RESULTS_DIR}/${name}.pdbqt" \
         --out      "${RESULTS_DIR}/${name}_poses.pdbqt" \
         --log      "${RESULTS_DIR}/${name}_log.txt"

    echo "  v Done: $name"
done < "$SMILES_FILE"

echo ""
echo "======================================================"
echo " All docking complete. Extracting best scores..."
echo "======================================================"
echo ""
grep -h "^   1 " ${RESULTS_DIR}/*_log.txt 2>/dev/null | \
  awk -v OFS="\\t" "{print \$2, \$3, \$4}" | \
  sort -k2 -n

