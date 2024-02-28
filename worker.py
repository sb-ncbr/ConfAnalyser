from typing import Optional
from benzene import Benzene
from cyclohexane import Cyclohexane
from cyclopentane import Cyclopentane
from molecule import Molecule, MoleculeType
from oxane_v2 import Oxane_v2
from utils import load_file

def work_file(resources) -> Optional[Molecule]:
    """
    Process a single file and analyze it.

    Note:
    This function has been extracted into own file to allow performance
    testing of the entire program.
    """
    file = resources[0]
    molecule_type = resources[3]
    if resources[4]:  # Is true if we're computing data multi-threaded
        Molecule.names = resources[1]
        Molecule.config = resources[2]
    filename = file.replace("\n", "").replace("\r", "")
    data = load_file(filename)

    match molecule_type:
        case MoleculeType.Oxane:
            molecule = Oxane_v2(data)
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
        print(f"{filename}: ommited")
        return None
