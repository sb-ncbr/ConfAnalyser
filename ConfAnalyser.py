from argparse import ArgumentParser
from os import name as os_name

from Molecules.Components.molecule import Molecule, MoleculeType
from Utils.config import Config
from Utils.worker import work_file
from Utils.utils import load_file, load_names

PERF_TEST = False
PARALLEL = True

def run(paths_file: str, names_file: str, molecule_type: MoleculeType,
        print_list: bool, print_summary: bool, print_all: bool) -> None:
    """
    Main runner function, loads all the data, creates dataset to run the
    program on, initializes config file and then runs the program
    either in single-thread mode or in multi-thread mode.
    """

    # Load data from the drive
    files = load_file(paths_file)
    names = load_names(names_file)
    cfg = Config()

    # Initialize the molecule structure with dictionary of atom names
    Molecule.initialize(names)

    # Prepare the dataset to work on
    data = [(file, names, cfg, molecule_type, PARALLEL) for file in files]

    # Run the analysis
    if PARALLEL:
        with Pool() as p:
            Molecule.molecules = list(p.map(work_file, data))
    else:
        for file in data:
            Molecule.molecules.append(work_file(file))

    # Remove all the invalid entries due to errors
    Molecule.molecules = [x for x in Molecule.molecules if x is not None]

    # Print out results based on input parameters
    if print_list or print_all:
        for m in Molecule.molecules:
            print(m)
    if print_summary or print_all:
        if len(Molecule.molecules) == 0:
            print("No molecules detected!")
        else:
            if print_list or print_all:
                print()

            # Call the print of statistics, the call is callable from any
            # molecule, using 1st one is safe option.
            Molecule.molecules[0].print_statistics()


def argument_parser():
    """
    Parse arguments from the command line.
    """
    parser = ArgumentParser(prog="python ConfAnalyser.py")
    if os_name == "posix":
        parser = ArgumentParser(prog="python3 ConfAnalyser.py")

    required = parser.add_argument_group('Required')
    required.add_argument('-i', '--input_list', required=True, type=str,
                          help=f'Read list of molecules to process from FILE - '
                               f'each line is treated as path to a single PDB file')
    required.add_argument("-n", "--name_list", required=True, type=str, action="store",
                          help="Read list of names of an atom ring FILE. Each line represents one ligand where"
                               "first word on each line is the ligand's name and all the following words on the line"
                               "are treated as atom names. If ligand is not known or if name of the atom is not found"
                               "in this list, atom will not be processed and will be omitted. In case of multiple"
                               "name variations, more lines with the same ligand name need to be present. Order of"
                               "atoms decides the order of atoms within the ring.")

    group = required.add_mutually_exclusive_group(required=True)
    group.add_argument("--cyclohexane", action="store_true", help="Set molecule type to cyclohexane.")
    group.add_argument("--cyclopentane", action="store_true", help="Set molecule type to cyclopentane.")
    group.add_argument("--benzene", action="store_true", help="Set molecule type to benzene.")

    optional = parser.add_argument_group('Optional').add_mutually_exclusive_group()
    optional.add_argument("-l", "--list", required=False, action='store_true',
                          help="Display results only as a list of molecules and their conformations.")
    optional.add_argument("-s", "--summary", required=False, action="store_true",
                          help="Display results only as a short summary of relative occurrences "
                               "of conformations among tested molecules.")
    optional.add_argument("-a", "--all", required=False, action="store_true",
                          help="Display both list and summary. Default option when neither -s or -l is used.")
    return parser.parse_args()

def main():
    args = argument_parser()

    paths_file = args.input_list
    names_file = args.name_list

    if args.cyclohexane:
        molecule_type = MoleculeType.Cyclohexane
    elif args.cyclopentane:
        molecule_type = MoleculeType.Cyclopentane
    elif args.benzene:
        molecule_type = MoleculeType.Benzene
    else:
        molecule_type = MoleculeType.Undefined

    print_list = args.list
    print_summary = args.summary
    print_all = args.all or (not print_list and not print_summary)

    run(paths_file, names_file, molecule_type, print_list, print_summary, print_all)


if __name__ == "__main__":
    # Additional imports based on what is needed.
    if PERF_TEST:
        from time import perf_counter
    if PARALLEL:
        from multiprocessing import Pool

    if PERF_TEST:
        start_time = perf_counter()
        main()
        print(f"Program finished after {perf_counter() - start_time} seconds")
    else:
        main()
