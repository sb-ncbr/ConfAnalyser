from geometries import Plane
from rings import SixAtomRing
from molecule import MoleculeType, Conformation


class Benzene(SixAtomRing):
    def __init__(self, source_file: list[str]):
        # Initialize the parent structure
        super().__init__()
        self.molecule_type = MoleculeType.Benzene
        self.set_conformations()
        self.conformation = Conformation.Undefined

        # Set the needed parameters
        self.source_file: list[str] = source_file

        try:
            self.create_from_source(source_file)
            self.ligand = self.atoms[0].residue_name if self.atoms else "Ligand not recognized!"  # TODO: raise error?

            self.validate_atoms()
            self.analyze()
        except Exception as e:
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
