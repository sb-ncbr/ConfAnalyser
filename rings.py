from molecule import Molecule, MoleculeType
from geometries import Plane
from atom import Atom


class Ring(Molecule):
    def __init__(self, molecule_type: MoleculeType):
        super().__init__(molecule_type)
        # TODO: Rework conformation into own enum? A class perhaps?
        self.conformation = "UNDEFINED"
        self.has_plane = False
        self.begin = 0

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.atoms[(self.begin + item) % self.get_atom_count()]
        return None


class SixAtomRing(Ring):
    def __init__(self, molecule_type: MoleculeType):
        super().__init__(molecule_type)

    def find_plane(self, tolerance: float, dist1: int = 1, dist2: int = 3, dist3: int = 4) -> bool:
        """
        Find the best suitable plane to start working with and set the `begin` parameter
        of the molecule to be used in other conformation validators.
        Also decide whether the molecule even has any valid plane at all.
        FindIdealPlane(ToleranceIn, Distance1 = 1, Distance 2 = 3, Distance 3 = 4):
            distance = infinity
            Molecule has plane = FALSE
            FOR i = 1 to i == 6
                Atom1 = Molecule's Atom i modulo 6
                Atom2 = Molecule's Atom (i + Distance1) modulo 6
                Atom3 = Molecule's Atom (i + Distance2) modulo 6
                Atom4 = Molecule's Atom (i + Distance3) modulo 6
                plane = Plane between Atom1, Atom2 and Atom3
                distanceFromPlane = distance between plane and Atom3
                IF distanceFromPlane < distance
                THEN
                    Molecule's begin = 1
                    distance = distanceFromPlane
                    IF distance between plane and Atom4 < ToleranceIn THEN
                        Molecule has plane = TRUE
                    FI
                FI
            END FOR
            RETURN Molecule has plane
        """
        # print(f"[find_plane@SixAtomRing] CALLED")
        distance = float("inf")
        self.has_plane = False

        for i in range(6):
            # print(f"\n[find_plane@SixAtomRing] ITERATION {i}")
            # print(f"[find_plane@SixAtomRing] Checking atoms: {self.atoms[i % 6], self.atoms[(i + dist1) % 6], self.atoms[(i + dist2) % 6]}")
            plane = Plane(self.atoms[i % 6], self.atoms[(i + dist1) % 6], self.atoms[(i + dist2) % 6])
            dist_from_plane = abs(plane.distance_from(self.atoms[(i + dist3) % 6]))
            # print(f"[find_plane@SixAtomRing] Atom [{i}]: {self.atoms[i]}")
            # print(f"[find_plane@SixAtomRing] Distance {i}: {dist_from_plane}")
            # print(f"[find_plane@SixAtomRing] Plane: < {plane.a}, {plane.b}, {plane.c}, {plane.d} >")
            if dist_from_plane < distance:
                # print(f"[find_plane@SixAtomRing] Smaller distance found, setting min dist to {dist_from_plane}")
                self.begin = i
                distance = dist_from_plane

                if plane.is_on_plane(self.atoms[(i + dist3) % 6], tolerance):
                    self.has_plane = True
        # print(f"[find_plane@SixAtomRing] - returning {self.has_plane}")
        return self.has_plane


class FiveAtomRing(Ring):
    def __init__(self, molecule_type: MoleculeType):
        super().__init__(molecule_type)

    def find_plane(self, tolerance: float, dist1: int = 1, dist2: int = 2, dist3: int = 3) -> None:
        """
        Find the best suitable plane to start working with and set the `begin` parameter
        of the molecule to be used in other conformation validators.
        Also decide whether the molecule even has any valid plane at all.
        """
        # print(f"[find_plane@FiveAtomRing] CALLED")
        distance = float("inf")

        for i in range(5):
            # print(f"\n[find_plane@FiveAtomRing] ITERATION {i}")
            # print(f"[find_plane@FiveAtomRing] Checking atoms: {self.atoms[i % 5], self.atoms[(i + dist1) % 5], self.atoms[(i + dist2) % 5]}")
            plane = Plane(self.atoms[i % 5], self.atoms[(i + dist1) % 5], self.atoms[(i + dist2) % 5])
            dist_from_plane = abs(plane.distance_from(self.atoms[(i + dist3) % 5]))
            # print(f"[find_plane@FiveAtomRing] Atom [{i}]: {self.atoms[i]}")
            # print(f"[find_plane@FiveAtomRing] Distance {i}: {dist_from_plane}")
            # print(f"[find_plane@FiveAtomRing] Plane: < {plane.a}, {plane.b}, {plane.c}, {plane.d} >")
            if dist_from_plane < distance:
                self.begin = i
                distance = dist_from_plane

                if plane.is_on_plane(self.atoms[(i + dist3) % 5], tolerance):
                    self.has_plane = True
