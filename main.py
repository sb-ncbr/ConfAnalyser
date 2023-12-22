from typing import *
from atom import Atom


def load_file(file_name: str) -> list[str]:
    with open(file_name) as file:
        return file.readlines()


if __name__ == "__main__":
    data = load_file("./dataset/sample_1.pdb")
    atom1 = Atom(data[1])
    print(atom1)
