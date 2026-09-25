def calc_novelty(mols_file, library):
    with open(mols_file, "r") as fnew, open(library, "r") as flib:
        new_mols = [s.strip() for s in fnew.readlines()]
        lib_mols = [s.strip() for s in flib.readlines()]

    novelty = len(list(set(new_mols + lib_mols))) / (len(new_mols + lib_mols))

    return novelty
