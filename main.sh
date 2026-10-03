#!/bin/bash
set -euo pipefail

export CUDA_VISIBLE_DEVICES=1,3
echo CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}

for path in "results/" "train/"; do {
    if [[ ! -d $path ]]; then {
        mkdir $path
    }; fi
}; done

CODENAME=1

read -p "This is the first turn, do you want to continue? (yes/[no]): " TAG
TAG=${TAG:-no}   # default to "no"

while [[ "$TAG" != "yes" && "$TAG" != "no" ]]; do
    echo "[WARNING] unrecognized input: ${TAG}"
    read -p "This is the first turn, do you want to continue? (yes/[no]): " TAG
    TAG=${TAG:-no}   # default to "no"
done

if [[ "$TAG" == "yes" ]]; then {
    echo "You are running the first turn! Prepared to run:"
    echo "python step1_PrepareData.py"
    read -p "Press Enter to continue..."

    echo;echo preparing data...
    python step1_PrepareData.py
    echo;echo training...
    python step2_Train.py

    echo;echo generating molecules... 
    python step3_GenerateMolecules.py --num 20000

    echo;echo calculate properties...
    python step4_CalculateProperties.py results/new_molecules_active.csv

    echo;echo non-dominated sorting...
    python step5_NondominatedSorting.py results/new_molecules_active_properties.csv -p avg SAscore QED_w -o min min max --top 10000

    echo;echo refining the model...
    python step6_Refine.py --model train/network.pth --output train/network_refined${CODENAME}.pth --mols results/new_molecules_active_properties_pareto_results/top_10000_smiles.txt    

    echo first turn is over
    read -p "Press Enter to continue..."
}
elif [[ "$TAG" == "no" ]]; then
    echo "Skipping the first turn..."
fi


while [[ $CODENAME <= 5 ]]; do {
    NEXT=${NEXT:-$((CODENAME + 1))}
    
    echo;echo generating molecules... 
    python step3_GenerateMolecules.py --model train/network_refined${CODENAME}.pth --output results/new_molecules${CODENAME}.txt --num 20000

    echo;echo calculate properties...
    python step4_CalculateProperties.py results/new_molecules${CODENAME}_active.csv

    echo;echo non-dominated sorting...
    python step5_NondominatedSorting.py results/new_molecules${CODENAME}_active_properties.csv -p avg SAscore QED_w -o min min max --top 10000

    echo;echo refining the model...
    python step6_Refine.py --model train/network_refined${CODENAME}.pth --output train/network_refined${NEXT}.pth --mols results/new_molecules${CODENAME}_active_properties_pareto_results/top_10000_smiles.txt    

    CODENAME=${NEXT}

} done

echo done!!
