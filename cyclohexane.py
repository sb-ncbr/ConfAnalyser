from rings import SixAtomRing
from atom import Atom
from geometries import Plane
from angle import dihedral_angle
from molecule import MoleculeType, Conformation


class Cyclohexane(SixAtomRing):
    def __init__(self, source_file: list[str]):
        # Initialize the parent structure
        super().__init__()
        self.molecule_type = MoleculeType.Cyclohexane
        self.set_conformations()
        self.conformation = Conformation.Undefined

        # Set the needed parameters
        self.source_file: list[str] = source_file
        # TODO: add validity check for when creating the molecule
        self.is_valid = True

        # TODO: Move tolerances into separate class to be loadable from file
        # To be replaced with Molecule's `config` parameter
        self.tolerance_in = 0.1
        self.tolerance_flat_in = 0.1
        self.tolerance_out = 0.6
        self.tolerance_tw_out = 0.4
        self.angle_tw_boat = 17.1
        self.angle_tolerance = 1.0

        try:
            self.create_from_source(source_file)
            self.ligand = self.atoms[0].residue_name if self.atoms else "Ligand not recognized!"  # TODO: raise error?

            self.validate_atoms()
            self.analyze()
        except Exception as e:
            print(e)

    def validate_atoms(self) -> None:
        """
        Does magic and shoots out sparks
        """
        new_lst = [None for _ in range(6)]
        for atom in self.atoms:
            for i in range(6):
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

    def create_from_source(self, source: list[str]) -> None:
        """
        Creates atoms from the source file
        """
        for i in range(6):
            self.atoms.append(Atom(source[i + 1]))

    def __str__(self) -> str:
        out = "Molecule of Cyclohexane:\n"
        out += self.conformation.name
        return out

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """

        self.find_plane(self.tolerance_in)

        if self.is_flat():
            self.conformation = Conformation.Flat
        elif self.is_half_chair():
            self.conformation = Conformation.Half_Chair
        elif self.is_chair():
            self.conformation = Conformation.Chair
        elif self.is_boat():
            self.conformation = Conformation.Boat
        elif self.is_twisted_boat():
            self.conformation = Conformation.Twisted_Boat
        else:
            self.conformation = Conformation.Undefined

    def is_flat(self) -> bool:
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.
        """
        if not self.has_plane:
            return False
        
        right_plane = Plane(self.index(0), self.index(1), self.index(3))
        left_plane = Plane(self.index(0), self.index(1), self.index(4))
        return (right_plane.is_on_plane(self.index(2), self.tolerance_flat_in) and
                right_plane.is_on_plane(self.index(5), self.tolerance_flat_in) and
                left_plane.is_on_plane(self.index(2), self.tolerance_flat_in) and
                left_plane.is_on_plane(self.index(5), self.tolerance_flat_in))

    def is_half_chair(self) -> bool:
        """
        Decides whether this molecule's conformation is half chair.
        Half chair conformation of molecule is determined by having all but one atom in one plane.
        """
        if not self.has_plane:
            return False

        plane = Plane(self.index(0), self.index(1), self.index(3))
        right_dist = abs(plane.distance_from(self.index(2)))
        left_dist = abs(plane.distance_from(self.index(5)))
        return (plane.is_on_plane(self.index(2), self.tolerance_flat_in) !=
                plane.is_on_plane(self.index(5), self.tolerance_flat_in)) and \
            plane.is_on_plane(self.index(4), self.tolerance_flat_in) and \
            ((right_dist > self.tolerance_out) != (left_dist > self.tolerance_out))

    def is_chair(self) -> bool:
        """
        Decides whether the molecule's conformation is a chair.
        Chair conformation is determined by having two opposite atoms on opposite sides of
        the main molecule plane.
        """
        if not self.has_plane:
            return False

        plane = Plane(self.index(0), self.index(1), self.index(3))
        right_dist = plane.distance_from(self.index(2))
        left_dist = plane.distance_from(self.index(5))
        return abs(right_dist) > self.tolerance_out and \
            abs(left_dist) > self.tolerance_out and \
            (right_dist * left_dist < 0)

    def is_boat(self) -> bool:
        """
        Decides whether the molecule's conformation is a boat.
        Boat conformation is determined by having two opposite atoms on the same side of a plane.
        """
        if not self.has_plane:
            return False

        plane = Plane(self.index(0), self.index(1), self.index(3))
        right_dist = plane.distance_from(self.index(2))
        left_dist = plane.distance_from(self.index(5))
        return abs(right_dist) > self.tolerance_out and \
            abs(left_dist) > self.tolerance_out and \
            (right_dist * left_dist > 0)

    def is_twisted_boat(self) -> bool:
        """
        Determine whether the molecule's conformation is a twisted boat.
        Twisted boat conformation is determined by having no plane within the ring.
        """
        right_plane = Plane(self.index(0), self.index(1), self.index(3))
        left_plane = Plane(self.index(0), self.index(1), self.index(4))
        right_dist = right_plane.distance_from(self.index(2))
        left_dist = left_plane.distance_from(self.index(5))
        tw_angle = abs(dihedral_angle(self.index(1), self.index(3), self.index(4), self.index(0)))
        return ((self.angle_tw_boat - self.angle_tolerance) < tw_angle < (self.angle_tw_boat + self.angle_tolerance)
                and abs(right_dist) > self.tolerance_tw_out
                and abs(left_dist) > self.tolerance_tw_out
                and (right_dist * left_dist) > 0)

