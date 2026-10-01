## Breif manual for several tool scripts
### `clean.py`

> If you are not sure whether your ligand library is clean (i.e. no metal ions, no salts, no duplicates, etc.), use this script to clean up the data.

Each row contains *one* SMILES string in `test_data.txt`. Run:

```bash
python clean.py test_data.txt
```

### `add_smiles_to_docking_data.py`

> The docking result file from Schrödinger may lack the SMILES column. This script helps to add SMILES metadata to the docking result file based on ID matching.

Suppose `SP_HTVS_OUT.csv` is the docking result file, `molecules.csv` is the original file before docking (contains the ID-Smiles information). 

Run as:
```bash
python add_smiles_to_docking_data.py -i SP_HTVS_OUT.csv -v molecules.csv
```
`SP_HTVS_OUT.csv` is replaced by the integrated file.

### `extract_docking_data.py`

> This script is used to merge the docking results with multiple (three) conformations into one CSV file. Active state and inactive state is dealed separately.

Prepare the input folder:
```
docking_result
|-active_state_1
| |-SP_HTVS_OUT.csv
|-active_state_2  
|-active_state_3  
...
|-inactive_state_1  
| |-SP_HTVS_OUT.csv
|-inactive_state_2  
|-inactive_state_3
...
```
Modify script arguments:
```python
DOCKING_DIR = os.path.join(SCRIPT_DIR, "docking_result")
...
CLASSES = {
    "active": "Generate0_active.csv",
    "inactive": "Generate0_inactive.csv",
}
```
Run the script:
```bash
python extract_docking_data.py
```
Generate:
```
Generate0_active.csv
Generate0_inactive.csv
```
```
ID,Smiles,score_1,score_2,score_3,avg,best,n_valid
...
```

<details open><summary> only one state? </summary>

```
docking_result
|-active_state_1
| |-
...
|-inactive_state_1  
| |-SP_HTVS_OUT.csv 
...
```

If you only have docking result for one state (active/inactive), an easy way is to move your target state forward and neglect the `FileNotFoundError` error.

```diff
CLASSES = {
+   "inactive": "Generate0_inactive.csv",
    "active": "Generate0_active.csv",
-   "inactive": "Generate0_inactive.csv",
}
```

</details>