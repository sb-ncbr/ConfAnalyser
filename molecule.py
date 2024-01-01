from __future__ import annotations
from atom import Atom
from enum import Enum
from config import Config


class MoleculeType(Enum):
    Undefined = 0
    Cyclopentane = 1
    Cyclohexane = 2
    Benzene = 3
    Oxane = 4


class Conformation(Enum):
    Undefined = -1  # Fallback option when no conformation was found
    Unanalysed = 0  # Default state of new molecule, no conformation found yet
    Flat = 1  # Cyclopentane, Cyclohexane, Benzene, Oxane
    Half_Chair = 2  # Cyclohexane, Oxane
    Chair = 3  # Cyclohexane, Oxane
    Boat = 4  # Cyclohexane, Oxane
    Twisted_Boat = 5  # Cyclohexane, Benzene
    Envelope = 6  # Cyclopentane, Oxane
    Twist = 7  # Cyclopentane
    Skew = 8  # Oxane


class Molecule:
    # A static variable containing data loaded from atom_names file on startup
    names: dict[str, list[list[str]]] = None
    # Config from which to load tolerances later
    config: Config = None
    # A list of all saved molecules
    molecules: list[Molecule] = []

    @staticmethod
    def initialize(names: dict[str, list[list[str]]]):
        if Molecule.names is None:
            Molecule.names = names

        if Molecule.config is None:
            Molecule.config = Config()

    def __init__(self):
        # A list of all atoms within a given molecule
        self.atoms: list[Atom] = []
        # A list of possible conformation states the molecule can be in
        self.conformations: list[Conformation] = [Conformation.Unanalysed, Conformation.Undefined]
        self.conformation = Conformation.Unanalysed
        self.molecule_type = MoleculeType.Undefined

    def print_statistics(self) -> None:
        print("SUMMARY\n-------")
        total = len(Molecule.molecules)
        for conf in self.conformations:
            count = sum([1 if x.conformation == conf else 0 for x in Molecule.molecules])
            percentage = (count / total) * 100
            percentage = int(percentage) if percentage % 1 == 0 else percentage
            print(f"{(conf.name.upper() + ':'):14}{count} ({percentage}%)")
        print(f"{'TOTAL:':14}{total}")

    def set_conformations(self) -> None:
        match self.molecule_type:
            case MoleculeType.Cyclohexane:
                self.set_conf_cyclohexane()
            case MoleculeType.Cyclopentane:
                self.set_conf_cyclopentane()
            case MoleculeType.Benzene:
                self.set_conf_benzene()
            case MoleculeType.Oxane:
                self.set_conf_oxane()

    def set_conf_cyclohexane(self) -> None:
        self.conformations = [Conformation.Boat, Conformation.Chair,
                              Conformation.Flat, Conformation.Half_Chair,
                              Conformation.Twisted_Boat] + self.conformations

    def set_conf_cyclopentane(self) -> None:
        self.conformations = [Conformation.Twist, Conformation.Envelope] + self.conformations

    def set_conf_benzene(self) -> None:
        self.conformations = [Conformation.Twisted_Boat] + self.conformations

    def set_conf_oxane(self) -> None:
        self.conformations = [Conformation.Half_Chair, Conformation.Chair, Conformation.Boat,
                              Conformation.Envelope, Conformation.Skew] + self.conformations
