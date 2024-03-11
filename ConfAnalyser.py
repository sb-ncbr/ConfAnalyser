from argparse import ArgumentParser
from os import name as os_name

from Molecules.Components.molecule import Molecule, MoleculeType
from Utils.config import Config
from Utils.worker import work_file
from Utils.utils import load_file, load_names

PERF_TEST = True
PARALLEL = True

class ConfAnalyser:
    # Exists so users can use this property when calling CA without further imports
    MoleculeType = MoleculeType

    def __init__(self, paths_file: str, names_file: str,
                 molecule_type: MoleculeType, print_list: bool = False,
                 print_summary: bool = False, print_all: bool = True,
                 parallel: bool = False, perf_test: bool = False) -> None:
        """
        Driver class for ConfAnalyser.

        This class exists to be imported by outside programs.

        Example:

        from ConfAnalyser import ConfAnalyser
        ConfAnalyser(paths_file=paths_to_pdbs_path, names_file=atom_names_path,
                     molecule_type=ConfAnalyser.MoleculeType.Cyclohexane)


        :param paths_file: a path to the file with all pdb file paths
        :param names_file:  a path to the file with atom names
        :param molecule_type: a molecule type of processed files
        :param print_list: after analysis is done, print list of files and
                           their determined conformations
        :param print_summary: after analysis is done, print summary result
        :param print_all:  after analysis is done, print both list and summary
        :param parallel:   enable parallel processing of data
        :param perf_test:  enable performance test, displaying runtime in seconds
        """

        self.paths_file = paths_file
        self.names_file = names_file
        self.molecule_type = molecule_type
        self.print_list = print_list
        self.print_summary = print_summary
        self.print_all = print_all
        self.parallel = parallel or PARALLEL and False  # temporarily disabled
        self.perf_test = perf_test or PERF_TEST


        if self.parallel:
            from multiprocessing import Pool


        if self.perf_test:
            from time import perf_counter

            start_time = perf_counter()
            self.run()
            print(f"Program finished after {perf_counter() - start_time} seconds")

        else:
            self.run()

    def run(self) -> None:
        """
        Main runner function, loads all the data, creates dataset to run the
        program on, initializes config file and then runs the program
        either in single-thread mode or in multi-thread mode.
        """

        # Load data from the drive
        files = load_file(self.paths_file)
        names = load_names(self.names_file)
        cfg = Config()

        # Initialize the molecule structure with dictionary of atom names
        Molecule.initialize(names)

        # Prepare the dataset to work on
        data = [(file, names, cfg, self.molecule_type, PARALLEL) for file in files]

        # Run the analysis
        if self.parallel:
            with Pool() as p:
                Molecule.molecules = list(p.map(work_file, data))
        else:
            for file in data:
                Molecule.molecules.append(work_file(file))

        # Remove all the invalid entries due to errors
        Molecule.molecules = [x for x in Molecule.molecules if x is not None]

        # Print out results based on input parameters
        if self.print_list or self.print_all:
            for m in Molecule.molecules:
                print(m)
        if self.print_summary or self.print_all:
            if len(Molecule.molecules) == 0:
                print("No molecules detected!")
            else:
                if self.print_list or self.print_all:
                    print()

                # Call the print of statistics, the call is callable from any
                # molecule, using 1st one is safe option.
                Molecule.molecules[0].print_statistics()


def argument_parser():
    """
    Parse arguments from the command line.

    This function only gets called if ConfAnalyser is called as a program
    directly from the command line.
    """
    parser = ArgumentParser(prog="python ConfAnalyser.py")
    if os_name == "posix":
        parser = ArgumentParser(prog="python3 ConfAnalyser.py")

    required = parser.add_argument_group('Required')
    required.add_argument('-i', '--input_list', required=True, type=str,
                          help=f'Read list of molecules to process from FILE - '
                               f'each line is treated as path to a single PDB file')
    required.add_argument("-n", "--name_list", required=True, type=str, action="store",
                          help="Read list of names of an atom ring FILE. Each line represents one ligand where "
                               "first word on each line is the ligand's name and all the following words on the line "
                               "are treated as atom names. If ligand is not known or if name of the atom is not found "
                               "in this list, atom will not be processed and will be omitted. In case of multiple "
                               "name variations, more lines with the same ligand name need to be present. Order of "
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

    ConfAnalyser(paths_file, names_file, molecule_type, print_list, print_summary, print_all)


if __name__ == "__main__":
    # Additional imports based on what is needed.
    main()
