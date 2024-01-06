from typing import *
from atom import Atom
from benzene import Benzene
from cyclohexane import Cyclohexane
from cyclopentane import Cyclopentane
from oxane import Oxane
from molecule import Molecule, MoleculeType, Conformation
import time
from config import Config

from multiprocessing import Pool



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


parallel = True


def work_file(resources) -> Optional[Molecule]:
    file = resources[0]
    molecule_type = resources[3]
    if parallel:
        Molecule.names = resources[1]
        Molecule.config = resources[2]
    filename = file.replace("\n", "").replace("\r", "")
    data = load_file(filename)
    molecule: Optional[Molecule, Cyclohexane] = None  # TODO: Remove later when we got certainty matching will work

    match molecule_type:
        case MoleculeType.Oxane:
            molecule = Oxane(data, file)
        case MoleculeType.Cyclohexane:
            molecule = Cyclohexane(data)
        case MoleculeType.Cyclopentane:
            molecule = Cyclopentane(data)
        case MoleculeType.Benzene:
            molecule = Benzene(data)

    if molecule.is_valid:
        return molecule

    else:
        # ligand = filename.split('/')[-1].split('_')[0]
        print(f"{filename}: ommited")
        return None
        # print(f"URL: https://www.rcsb.org/ligand/{ligand}")
        # print(f"URL: https://cdn.rcsb.org/images/ccd/labeled/{ligand[0]}/{ligand}.svg\n")


def run():
    molecule_type = MoleculeType.Oxane
    mol_type = ""
    match molecule_type:
        case MoleculeType.Cyclohexane:
            mol_type = "cyclohexanes"
        case MoleculeType.Cyclopentane:
            mol_type = "cyclopentanes"
        case MoleculeType.Oxane:
            mol_type = "oxanes"
        case MoleculeType.Benzene:
            mol_type = "benzenes"

    files = load_file(f"./{mol_type}/paths_to_pdbs.txt")
    # files = load_file(f"./{mol_type}/ommited.txt")
    # files = ["C:/oxanes/6UZ/patterns/6UZ_5klu_0.pdb"]
    # files = load_file("./dataset/sing_test_file_path")
    names = load_names(f"./{mol_type}/atom_names.txt")
    cfg = Config()
    # print(f"Loaded:")
    # print_dict(names)

    Molecule.initialize(names)

    data = [(file, names, cfg, molecule_type) for file in files]

    if parallel:
        with Pool() as p:
            Molecule.molecules = list(p.map(work_file, data))
    else:
        for file in data:
            Molecule.molecules.append(work_file(file))
            # print(f"{filename.split('/')[-1]}: {molecule.conformation.name.upper()}")
    # print_dict(Molecule.names)
    # for lst in names["18Y"]:
    #     out = ", ".join(x for x in lst)
    #     print(f"   -> {out}")
    Molecule.molecules = [x for x in Molecule.molecules if x is not None]
    if len(Molecule.molecules) == 0:
        print("No molecules detected!")
    else:
        Molecule.molecules[0].print_statistics()


if __name__ == "__main__":
    start_time = time.perf_counter()
    run()
    print(f"Program finished after {time.perf_counter() - start_time} seconds")
