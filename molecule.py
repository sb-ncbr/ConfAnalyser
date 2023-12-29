from atom import Atom
from typing import *


class Molecule:
    # A static variable containing data loaded from atom_names file on startup
    names: dict[str, list[list[str]]] = None

    def __init__(self):
        self.atoms: list[Atom] = []

