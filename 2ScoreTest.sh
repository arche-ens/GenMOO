#!/bin/bash
set -euo pipefail

export CUDA_VISIBLE_DEVICES=1,3
echo CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES}

for path in "2ScoreExample/" "2ScoreExample/results/" "2ScoreExample/train/"; do {
    if [[ ! -d $path ]]; then {
        mkdir $path
    }; fi
}; done

CODENAME=1
echo CODENAME: ${CODENAME}

read -p "This is the first turn, do you want to continue? (yes/[no]): " TAG
TAG=${TAG:-no}   # default to "no"

while [[ "$TAG" != "yes" && "$TAG" != "no" ]]; do
    echo "[WARNING] unrecognized input: ${TAG}"
    read -p "This is the first turn, do you want to continue? (yes/[no]): " TAG
    TAG=${TAG:-no}   # default to "no"
done

if [[ "$TAG" == "yes" ]]; then {
    echo "You are running the first turn! Prepared to run:"
    echo "python step3_GenerateMolecules.py --model train/network.pth --output 2ScoreExample/results/2SE_molecules${CODENAME}.txt --num 20000 --temp 0.7"
    read -p "Press Enter to continue..."

    echo;echo generating molecules... 
    python step3_GenerateMolecules.py --model train/network.pth --output 2ScoreExample/results/2SE_molecules${CODENAME}.txt --num 20000 --temp 0.7

    echo;echo calculate properties...
    python step4_CalculateProperties.py 2ScoreExample/results/2SE_molecules${CODENAME}.csv

    echo;echo non-dominated sorting...
    python step5_NondominatedSorting.py 2ScoreExample/results/2SE_molecules${CODENAME}_properties.csv -p SAscore QED_w -o min max --top 10000

    echo;echo refining the model...
    python step6_Refine.py --model train/network.pth --output 2ScoreExample/train/2SE_network_refined${CODENAME}.pth --mols 2ScoreExample/results/2SE_molecules${CODENAME}_properties_pareto_results/top_10000_smiles.txt    

    echo first turn is over
    read -p "Press Enter to continue..."
}
elif [[ "$TAG" == "no" ]]; then
    echo "Skipping the first turn..."
fi

TAG="yes"
while [[ "$TAG" == "yes" ]]; do {
    NEXT=$((CODENAME + 1))
    read -p "Codename for the next turn [default: ${NEXT}]: " CODENAME
    CODENAME=${CODENAME:-${NEXT}}
    
    read -p "Using codename: $CODENAME, do you want to continue? (yes/[no]): " TAG
    TAG=${TAG:-no}   # default to "no"
    while [[ "$TAG" != "yes" && "$TAG" != "no" ]]; do {
        echo "[WARNING] unrecognized input :( Please enter 'yes' or 'no'"
        read -p "Using codename: ${CODENAME}, do you want to continue? (yes/[no]): " TAG
        TAG=${TAG:-no}   # default to "no"
    } done

    if [[ "$TAG" == "yes" ]]; then {
        echo "Prepared to run:"
        echo "python step3_GenerateMolecules.py --model 2ScoreExample/train/2SE_network_refined$((CODENAME - 1)).pth --output 2ScoreExample/results/2SE_molecules${CODENAME}.txt --num 20000 --temp 0.7"
        read -p "Press Enter to continue..."

        echo;echo generating molecules... 
        python step3_GenerateMolecules.py --model 2ScoreExample/train/2SE_network_refined$((CODENAME - 1)).pth --output 2ScoreExample/results/2SE_molecules${CODENAME}.txt --num 20000 --temp 0.7

        echo;echo calculate properties...
        python step4_CalculateProperties.py 2ScoreExample/results/2SE_molecules${CODENAME}.csv

        echo;echo non-dominated sorting...
        python step5_NondominatedSorting.py 2ScoreExample/results/2SE_molecules${CODENAME}_properties.csv -p SAscore QED_w -o min max --top 10000

        echo;echo refining the model...
        python step6_Refine.py --model 2ScoreExample/train/2SE_network_refined$((CODENAME - 1)).pth --output 2ScoreExample/train/2SE_network_refined${CODENAME}.pth --mols 2ScoreExample/results/2SE_molecules${CODENAME}_properties_pareto_results/top_10000_smiles.txt    
    }
    elif [[ "$TAG" == "no" ]]; then
        echo "no codename used"
    fi
} done




