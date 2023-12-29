from rings import SixAtomRing
from atom import Atom
from geometries import Plane
from angle import dihedral_angle


class Cyclohexane(SixAtomRing):
    def __init__(self, source_file: list[str]):
        super().__init__()
        self.source_file: list[str] = source_file
        # TODO: add validity check for when creating the molecule
        self.is_valid = True

        # TODO: Move tolerances into separate class to be loadable from file
        self.tolerance_in = 0.1
        self.tolerance_flat_in = 0.1
        self.tolerance_out = 0.6
        self.tolerance_tw_out = 0.4
        self.angle_tw_boat = 17.1
        self.angle_tolerance = 1.0

        self.create_from_source(source_file)
        self.ligand = self.atoms[0].residue_name if self.atoms else "Ligand not recognized!"  # TODO: raise error?

        self.sort_atoms()
        self.analyze()

    def sort_atoms(self) -> None:
        new_lst = []
        for atom_name in self.names[self.ligand][0]:
            for atom in self.atoms:
                if atom.name == atom_name:
                    print(f"To atom name {atom_name} assigning atom {atom}")
                    new_lst.append(atom)
        if len(new_lst) != 6:
            # TODO: throw some exception, handle, idk
            print("Not all atoms were found!")
        self.atoms = new_lst

    def create_from_source(self, source: list[str]) -> None:
        for i in range(6):
            self.atoms.append(Atom(source[i + 1]))

    def __str__(self) -> str:
        out = "Molecule of Cyclohexane:\n"
        out += self.conformation
        return out

    def analyze(self):
        print("Analyzing..")
        for atom in self.atoms:
            print(atom)
        self.find_plane(self.tolerance_in)
        print(f"Begin at: {self.begin}")
        if self.is_flat():
            self.conformation = "FLAT"
        elif self.is_half_chair():
            self.conformation = "HALF CHAIR"
        elif self.is_chair():
            self.conformation = "CHAIR"
        elif self.is_boat():
            self.conformation = "BOAT"
        elif self.is_twisted_boat():
            self.conformation = "TWISTED BOAT"
        else:
            self.conformation = "UNKNOWN"

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
            (right_dist > self.tolerance_out != left_dist > self.tolerance_out)

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
        print(right_dist)
        print(left_dist)
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

