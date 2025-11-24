import h5py
import numpy as np
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

# fonctions de comparaison générées avec chatGPT, fonctionne bien


def compare_datasets(dset1, dset2, path=""):
    """Compare deux datasets HDF5."""
    if dset1.shape != dset2.shape:
        print(f"Shape différent pour {path}: {dset1.shape} vs {dset2.shape}")
        return False
    if dset1.dtype != dset2.dtype:
        print(f"Type différent pour {path}: {dset1.dtype} vs {dset2.dtype}")
        return False
    if not np.array_equal(dset1[()], dset2[()]):
        print(f"Valeurs différentes dans {path}")
        return False
    return True


def compare_attrs(attrs1, attrs2, path=""):
    """Compare les attributs d'un groupe ou dataset."""
    keys1 = set(attrs1.keys())
    keys2 = set(attrs2.keys())
    if keys1 != keys2:
        print(f"Attributs différents dans {path}: {keys1 ^ keys2}")
        return False
    for key in keys1:
        if not np.array_equal(attrs1[key], attrs2[key]):
            print(f"Attribut '{key}' différent dans {path}")
            return False
    return True


def compare_groups(group1, group2, path=""):
    """Compare récursivement deux groupes HDF5."""
    # Vérifie les datasets et groupes
    names1 = set(group1.keys())
    names2 = set(group2.keys())
    if names1 != names2:
        print(f"Noms différents dans {path}: {names1 ^ names2}")

    all_ok = True
    for name in names1 & names2:
        item1 = group1[name]
        item2 = group2[name]
        current_path = f"{path}/{name}" if path else name

        if isinstance(item1, h5py.Dataset) and isinstance(item2, h5py.Dataset):
            if not compare_datasets(item1, item2, current_path):
                all_ok = False
            if not compare_attrs(item1.attrs, item2.attrs, current_path):
                all_ok = False
        elif isinstance(item1, h5py.Group) and isinstance(item2, h5py.Group):
            if not compare_attrs(item1.attrs, item2.attrs, current_path):
                all_ok = False
            if not compare_groups(item1, item2, current_path):
                all_ok = False
        else:
            print(f"Type différent pour {current_path}")
            all_ok = False
    return all_ok


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Compare two hdf5 files", formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--files",
        nargs=2,
        required=True,
        help="List of files to compare : reference, prediction",
    )
    args = parser.parse_args()

    with h5py.File(args.files[0], "r") as f1, h5py.File(args.files[1], "r") as f2:
        all_ok = compare_groups(f1, f2)
        if not all_ok:
            raise ValueError("Differences detected")
