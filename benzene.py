from geometries import Plane
from rings import SixAtomRing
from molecule import MoleculeType, Conformation


class Benzene(SixAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Benzene)

        # Set the needed parameters
        self.source_file: list[str] = source_line

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(e)

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """

        self.find_plane(self.config.b.t_flat_in)

        if self.is_flat():
            self.conformation = Conformation.Flat
        else:
            self.conformation = Conformation.Undefined

    def is_flat(self):
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.

        Is shared by all atoms
        """
        if not self.has_plane:
            return False

        right_plane = Plane(self.index(0), self.index(1), self.index(3))
        left_plane = Plane(self.index(0), self.index(1), self.index(4))
        return (right_plane.is_on_plane(self.index(2), self.config.b.t_flat_in) and
                right_plane.is_on_plane(self.index(5), self.config.b.t_flat_in) and
                left_plane.is_on_plane(self.index(2), self.config.b.t_flat_in) and
                left_plane.is_on_plane(self.index(5), self.config.b.t_flat_in))
