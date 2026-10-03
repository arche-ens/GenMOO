#!/bin/bash
set -euo pipefail

export CUDA_VISIBLE_DEVICES=1,2
echo CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}

echo;echo generating molecules... 
python step3_GenerateMolecules.py --model model/network_final.pth --output results/new_molecules_final.txt --num 20000

echo;echo calculate properties...
python step4_CalculateProperties.py results/new_molecules_final_active.csv

echo;echo non-dominated sorting...
python step5_NondominatedSorting.py results/new_molecules_final_active_properties.csv -p avg SAscore QED_w -o min min max --top 10000

cp results/new_molecules_final_active_properties_pareto_results/pareto_front_001.csv results/final_candidates.csv

echo done!!
echo find candidates at: results/final_candidates.csv
