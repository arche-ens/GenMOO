import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting


def non_dominated_sort(df: pd.DataFrame, objectives: dict):
    """
    Perform non-dominated sorting.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing objective columns.
    objectives : dict
        Mapping from column name to "min" or "max".
        Example:
            {
                "c1": "min",
                "c2": "max",
                "c3": "min",
                "c4": "max",
                "c5": "min"
            }

    Returns
    -------
    ranks : numpy.ndarray
        Pareto rank for every row. First front = rank 1.
    fronts : list of numpy.ndarray
        Each array contains the row indices belonging to that front.
    """
    columns = list(objectives.keys())

    # Extract objective values
    values = df[columns].to_numpy(dtype=float)

    # pymoo's non-dominated sorting assumes MINIMIZATION.
    # Therefore, multiply maximization objectives by -1.
    for i, column in enumerate(columns):
        if objectives[column].lower() == "max":
            values[:, i] *= -1
        elif objectives[column].lower() != "min":
            raise ValueError(
                f"Objective for {column} must be 'min' or 'max'"
            )

    # Non-dominated sorting
    fronts = NonDominatedSorting().do(values)
    # NonDominatedSorting(method='fast_non_dominated_sort')
            # fast_non_dominated_sort
            # efficient_non_dominated_sort
            # tree_based_non_dominated_sort
            # dominance_degree_non_dominated_sort

    # Assign Pareto rank
    ranks = np.empty(len(df), dtype=int)
    for rank, front in enumerate(fronts, start=1):
        ranks[front] = rank

    return ranks, fronts


def main():
    parser = argparse.ArgumentParser(description="Pareto non-dominated sorting")
    parser.add_argument("input", type=Path, help="Input CSV file")
    parser.add_argument("-p", "--properties", nargs="+", required=True, help="Column name of properties that participate in NDS.")
    parser.add_argument("-o", "--objectives", nargs="+", choices=["min", "max"], required=True, help="Optimization direction for each objective column")
    parser.add_argument("--top", nargs="+", type=int, default=None, help="Select the top-N molecules by Pareto rank (add whole fronts until exceeding N). Multiple values allowed, e.g. --top 50 5000. If omitted, keep all.")
    parser.add_argument("--top_fronts", type=int, default=None, help="Merge the top-N Pareto fronts into a single CSV file.")
    args = parser.parse_args()

    input_csv = args.input
    df = pd.read_csv(input_csv)
    objective_columns = args.properties # ["avg", "SAscore", "My_QED"]

    if "ID" not in df.columns:
        print("[WARNING] 'ID' column is not found in the file. please check the file.\n")
    if len(args.objectives) != len(objective_columns):
        raise ValueError(
            f"Found {len(objective_columns)} objective columns, but received {len(args.objectives)} objective directions."
        )

    objectives = dict(zip(objective_columns, args.objectives))

    print("Objectives:")
    for column, direction in objectives.items():
        print(f"  {column}: {direction}")

    print(f"\nNumber of solutions: {len(df):,}")

    ########## Non-dominated sorting ##########
    ranks, fronts = non_dominated_sort(df, objectives)

    df["pareto_rank"] = ranks

    output_dir = Path(input_csv).parent / (input_csv.stem + "_pareto_results")
    output_dir.mkdir(parents=True, exist_ok=True)
    combined_file = output_dir / "all_solutions.csv"
    df.to_csv(combined_file, index=False)

    for rank, front in enumerate(fronts, start=1):
        front_df = df.iloc[front].copy()
        front_file = (output_dir / f"pareto_front_{rank:03d}.csv")
        front_df.to_csv(front_file, index=False)
        print(
            f"Front {rank:3d}: "
            f"{len(front):6d} solutions -> "
            f"{front_file}"
        )

    if args.top is not None:
        for top in args.top:
            selected = []
            n_fronts = 0
            for front in fronts:
                if len(selected) + len(front) > top:
                    break
                selected.extend(front)
                n_fronts += 1

            if "Smiles" not in df.columns:
                print("\n[WARNING] 'Smiles' column not found. Skipping --top output.")
                break
            elif len(selected) == 0:
                print(
                    f"\n[WARNING] First front already exceeds --top {top}; "
                    "no molecules selected."
                )
            else:
                smiles_file = output_dir / f"top_{top}_smiles.txt"
                df.iloc[selected]["Smiles"].to_csv(
                    smiles_file, index=False, header=False
                )
                print(
                    f"\nSelected {len(selected):,} molecules "
                    f"(<= {top:,}) across {n_fronts} fronts -> "
                    f"{smiles_file}"
                )

    if args.top_fronts is not None:
        n = min(args.top_fronts, len(fronts))
        merged = df[df["pareto_rank"] <= n]
        merged_file = output_dir / f"top_fronts_{n:03d}.csv"
        merged.to_csv(merged_file, index=False)
        print(
            f"\nMerged top {n} fronts "
            f"({len(merged):,} solutions) -> {merged_file}"
        )

    print("\nDone.")
    print(f"Total Pareto fronts: {len(fronts)}")
    print(f"Results saved to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()