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


def print_dict(dct):
    for name in dct:
        print(f"Ligand: {name}")
        for lst in dct[name]:
            out = ", ".join(x for x in lst)
            print(f"   -> {out}")


if __name__ == "__main__":
    data = load_file("./dataset/00K_1a46_1.pdb")
    Molecule.names = load_names("./dataset/atom_names.txt")
    hexane = Cyclohexane(data)
    print_dict(hexane.names)
    print(f"Got: {hexane}, expected: CHAIR\n")

    # data = load_file("./dataset/00P_1d9i_0.pdb")
    # molecule = Cyclohexane(data)
    # print(f"Got: {molecule}, expected: UNKNOWN\n")
