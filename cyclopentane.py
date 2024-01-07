from geometries import Plane
from rings import FiveAtomRing
from molecule import MoleculeType, Conformation
from angle import dihedral_angle


class Cyclopentane(FiveAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Cyclopentane)

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

        self.find_plane(self.config.cp.t_in)

        if self.is_flat():
            self.conformation = Conformation.Flat
        elif self.is_envelope():
            self.conformation = Conformation.Envelope
        elif self.is_twist():
            self.conformation = Conformation.Twist
        else:
            self.conformation = Conformation.Undefined

    def is_flat(self) -> bool:
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.
        """
        if not self.has_plane:
            return False

        plane = Plane(self.index(0), self.index(1), self.index(2))
        return (plane.is_on_plane(self.index(3), self.config.cp.t_in)
                and plane.is_on_plane(self.index(4), self.config.cp.t_in))

    def is_envelope(self) -> bool:
        """
        Decides whether this molecule's conformation is envelope.
        Envelope conformation of molecule is determined by having all but one atom on a same plane.
        """
        if not self.has_plane:
            return False

        plane = Plane(self.index(0), self.index(1), self.index(2))
        return abs(plane.distance_from(self.index(4))) > self.config.cp.t_out

    def is_twist(self):
        """
        Decides whether this molecule's conformation is envelope.
        Twist conformation of molecule is determined by having no plane within the ring.
        """

        left_plane = Plane(self.index(0), self.index(1), self.index(3))
        right_plane = Plane(self.index(0), self.index(2), self.index(3))
        left_distance = left_plane.distance_from(self.index(4))
        right_distance = right_plane.distance_from(self.index(4))
        twist_angle = dihedral_angle(self.index(0), self.index(1), self.index(2), self.index(3))

        return ((abs(abs(twist_angle) - self.config.cp.t_tw_boat_angle) < self.config.cp.t_angle)
                and (abs(right_distance) > self.config.cp.t_tw_out)
                and (abs(left_distance) > self.config.cp.t_tw_out)
                and (right_distance * left_distance > 0))
