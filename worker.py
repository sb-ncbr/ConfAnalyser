from typing import Optional
from benzene import Benzene
from cyclohexane import Cyclohexane
from cyclopentane import Cyclopentane
from molecule import Molecule, MoleculeType
from oxane import Oxane
from oxane_v2 import Oxane_v2

parallel = False


def load_file(file_name: str) -> list[str]:
    try:
        with open(file_name) as file:
            return file.readlines()
    except OSError:
        print(f"[ERROR] File `{file_name}` could not be opened!\nCheck the path and try again.")
        exit()


def load_names(file_name: str) -> dict[str, list[set[str]]]:
    lines = load_file(file_name)
    out: dict[str, list[set[str]]] = dict()
    for line in lines:
        split_line = line.split()
        ligand = split_line[0]
        if ligand[2] == "_":
            ligand = ligand[0:2]
        if ligand in out:
            for i, name in enumerate(split_line[1:]):
                out[ligand][i].add(name)
        else:
            out[ligand] = []
            for name in split_line[1:]:
                out[ligand].append({name})
    return out


def load_file_list(file_name: str) -> list[str]:
    return load_file(file_name)


def print_dict(dct):
    for name in dct:
        print(f"Ligand: {name}")
        for lst in dct[name]:
            out = ", ".join(x for x in lst)
            print(f"   -> {out}")

def work_file(resources) -> Optional[Molecule]:
    """
    Process a single file and analyze it.
    """
    file = resources[0]
    molecule_type = resources[3]
    if parallel:
        Molecule.names = resources[1]
        Molecule.config = resources[2]
    filename = file.replace("\n", "").replace("\r", "")
    data = load_file(filename)

    match molecule_type:
        case MoleculeType.Oxane:
            molecule = Oxane_v2(data)
            #molecule.updated_version = Oxane_v2(data)
        case MoleculeType.Cyclohexane:
            molecule = Cyclohexane(data)
        case MoleculeType.Cyclopentane:
            molecule = Cyclopentane(data)
        case MoleculeType.Benzene:
            molecule = Benzene(data)
        case _:
            molecule = None
    molecule.set_file_name(file)

    if molecule.is_valid:
        return molecule

    else:
        # ligand = filename.split('/')[-1].split('_')[0]
        print(f"{filename}: ommited")
        return None
        # print(f"URL: https://www.rcsb.org/ligand/{ligand}")
        # print(f"URL: https://cdn.rcsb.org/images/ccd/labeled/{ligand[0]}/{ligand}.svg\n")