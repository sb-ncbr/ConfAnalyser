from typing import *
from atom import Atom
from cyclohexane import Cyclohexane
from molecule import Molecule


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
            out[ligand].append(split_line[1:])
        else:
            out[ligand] = [split_line[1:]]
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
    files = load_file("./dataset/file_list.txt")
    names = load_names("./dataset/atom_names.txt")
    print(f"Loaded:")
    print_dict(names)
    Molecule.names = names

    for file in files:
        data = load_file(file.replace("\n", "").replace("\r", ""))
        molecule = Cyclohexane(data)
        print(molecule)


if __name__ == "__main__":
    run()
