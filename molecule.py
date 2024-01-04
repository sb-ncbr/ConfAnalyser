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
        """
        Initializes the static parameters of the entire Molecule class
        on the first call.
        """
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

        # TODO: add validity check for when creating the molecule
        self.is_valid = True

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
        """
        Sets the correct order of possible conformations for each molecule
        which are later used for printing statistics.
        """
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

    def get_atom_count(self):
        match self.molecule_type:
            case MoleculeType.Cyclohexane:
                return 6
            case MoleculeType.Benzene:
                return 6
            case MoleculeType.Oxane:
                return 6
            case MoleculeType.Cyclopentane:
                return 5

    def create_from_source(self, source: list[str]) -> None:
        """
        Creates atoms from the source file
        """
        for i in range(self.get_atom_count()):
            self.atoms.append(Atom(source[i + 1]))

    def validate_atoms(self) -> None:
        """
        Validates whether all the atoms of a molecule have their name
        present in the names file and then creates the ordering
        of atoms based on the order from the names file. In case of
        duplicate names within the same index of atom position, error
        out and cancel processing of this molecule as a result.
        """
        atom_count = self.get_atom_count()
        new_lst = [None for _ in range(atom_count)]
        for atom in self.atoms:
            for i in range(atom_count):
                if atom.name in self.names[self.ligand][i]:
                    if new_lst[i] is not None:
                        print(f"{atom.name} atom found twice!")
                        self.is_valid = False
                        return
                    new_lst[i] = atom
                    break
        if None in new_lst:
            print("Not all atoms were found")
            self.is_valid = False
            return
        self.atoms = new_lst
