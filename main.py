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
from argparse import ArgumentParser


# default option for parallel processing
parallel = True


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


def work_file(resources) -> Optional[Molecule]:
    file = resources[0]
    molecule_type = resources[3]
    if parallel:
        Molecule.names = resources[1]
        Molecule.config = resources[2]
    filename = file.replace("\n", "").replace("\r", "")
    data = load_file(filename)
    # TODO: Meh
    molecule: Optional[Molecule, Cyclohexane, Cyclopentane, Oxane, Benzene] = None

    match molecule_type:
        case MoleculeType.Oxane:
            molecule = Oxane(data)
        case MoleculeType.Cyclohexane:
            molecule = Cyclohexane(data)
        case MoleculeType.Cyclopentane:
            molecule = Cyclopentane(data)
        case MoleculeType.Benzene:
            molecule = Benzene(data)
    molecule.set_file_name(file)

    if molecule.is_valid:
        return molecule

    else:
        # ligand = filename.split('/')[-1].split('_')[0]
        print(f"{filename}: ommited")
        return None
        # print(f"URL: https://www.rcsb.org/ligand/{ligand}")
        # print(f"URL: https://cdn.rcsb.org/images/ccd/labeled/{ligand[0]}/{ligand}.svg\n")


def run():
    global parallel
    parallel = True
    molecule_type = MoleculeType.Oxane
    mol_type = ""
    print_summary = True
    print_list = False
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
        with Pool(1) as p:
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
    if print_list:
        for m in Molecule.molecules:
            print(m)
    if print_summary:
        if len(Molecule.molecules) == 0:
            print("No molecules detected!")
        else:
            Molecule.molecules[0].print_statistics()


def main():
    parser = ArgumentParser(prog="ConfAnalyser",
                            description="TODO")
    required = parser.add_argument_group('Required')
    required.add_argument('-i', '--input_list', required=True, type=str,
                          help=f'Read list of molecules to process from FILE - '
                               f'each line is treated as path to single PDB file')
    required.add_argument("-n", "--name_list", required=True, type=str, action="store",
                          help="Read list of names of an atom ring FILE. Each line represents one ligand where"
                               "first word on each line is the ligand's name and all the following words on the line"
                               "are trated as atom names. If ligand is not known or if name of the atom is not found"
                               "in this list, atom will not be processed and will be omitted. In case of multiple"
                               "name variatons, more lines with the same ligand name need to be present. Order of"
                               "atoms decides the order of atoms within the ring.")

    optional = parser.add_argument_group('Optional')
    optional.add_argument("-l", "--list", required=False, action='store_true',
                          help="Display results only as a list of molecules and their conformations.")
    optional.add_argument("-s", "--summary", required=False, action="store_true",
                          help="Display results only as a short summary of relative occurances "
                               "of conformations among tested molecules.")
    optional.add_argument("-a", "--all", required=False, action="store_true",
                          help="Display both list and summary. Default option when neither -s or -l is used.")
    args = parser.parse_args()


if __name__ == "__main__":
    start_time = time.perf_counter()
    # run()
    main()
    print(f"Program finished after {time.perf_counter() - start_time} seconds")
