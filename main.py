from typing import *
from atom import Atom
from cyclohexane import Cyclohexane
from molecule import Molecule, MoleculeType


def load_file(file_name: str) -> list[str]:
    with open(file_name) as file:
        return file.readlines()


def load_names(file_name: str) -> dict[str, list[list[str]]]:
    lines = load_file(file_name)
    out: dict[str, list[list[str]]] = dict()
    for line in lines:
        split_line = line.split()
        ligand = split_line[0]
        if ligand in out:
            for i, name in enumerate(split_line[1:]):
                out[ligand][i].append(name)
        else:
            out[ligand] = []
            for name in split_line[1:]:
                out[ligand].append([name])
    return out


def load_file_list(file_name: str) -> list[str]:
    return load_file(file_name)


def print_dict(dct):
    for name in dct:
        print(f"Ligand: {name}")
        for lst in dct[name]:
            out = ", ".join(x for x in lst)
            print(f"   -> {out}")


def run():
    molecule_type = MoleculeType.Cyclohexane
    files = load_file("./filtered_ligands/paths_to_pdbs.txt")
    # files = load_file("./dataset/sing_test_file_path")
    names = load_names("./filtered_ligands/atom_names.txt")
    # print(f"Loaded:")
    # print_dict(names)

    Molecule.initialize(names)

    for file in files:
        filename = file.replace("\n", "").replace("\r", "")
        data = load_file(filename)
        molecule: Optional[Molecule, Cyclohexane] = None  # TODO: Remove later when we got certainty matching will work

        match molecule_type:
            case MoleculeType.Cyclohexane:
                molecule = Cyclohexane(data)
            case MoleculeType.Cyclopentane:
                ...
            case MoleculeType.Benzene:
                ...
            case MoleculeType.Oxane:
                ...

        if molecule.is_valid:
            Molecule.molecules.append(molecule)

        else:
            ligand = filename.split('/')[-1].split('_')[0]
            print(f"{filename}: ommited")
            # print(f"URL: https://www.rcsb.org/ligand/{ligand}")
            # print(f"URL: https://cdn.rcsb.org/images/ccd/labeled/{ligand[0]}/{ligand}.svg\n")
        # print(f"{filename.split('/')[-1]}: {molecule.conformation.name.upper()}")

    # print_dict(Molecule.names)
    Molecule.molecules[0].print_statistics()


if __name__ == "__main__":
    run()
