from molecule import Molecule
from geometries import Plane
from atom import Atom


class Ring(Molecule):
    def __init__(self):
        super().__init__()
        # TODO: Rework conformation into own enum? A class perhaps?
        self.conformation = "UNDEFINED"
        self.has_plane = False
        self.begin = 0

    def index(self, i: int) -> Atom:
        # TODO: Think of a better name? index kinda sus
        """
        Returns ring's atom at a given index while including the `begin` index determined
        by the find_plane function.
        """
        return self.atoms[(self.begin + i) % self.get_atom_count()]


class SixAtomRing(Ring):
    def __init__(self):
        super().__init__()

    def find_plane(self, tolerance: float, dist1: int = 1, dist2: int = 3, dist3: int = 4) -> None:
        """
        Find the best suitable plane to start working with and set the `begin` parameter
        of the molecule to be used in other conformation validators.
        Also decide whether the molecule even has any valid plane at all.
        """
        # print(f"[find_plane@SixAtomRing] CALLED")
        distance = float("inf")

        for i in range(6):
            # print(f"\n[find_plane@SixAtomRing] ITERATION {i}")
            # print(f"[find_plane@SixAtomRing] Checking atoms: {self.atoms[i % 6], self.atoms[(i + dist1) % 6], self.atoms[(i + dist2) % 6]}")
            plane = Plane(self.atoms[i % 6], self.atoms[(i + dist1) % 6], self.atoms[(i + dist2) % 6])
            dist_from_plane = abs(plane.distance_from(self.atoms[(i + dist3) % 6]))
            # print(f"[find_plane@SixAtomRing] Atom [{i}]: {self.atoms[i]}")
            # print(f"[find_plane@SixAtomRing] Distance {i}: {dist_from_plane}")
            # print(f"[find_plane@SixAtomRing] Plane: < {plane.a}, {plane.b}, {plane.c}, {plane.d} >")
            if dist_from_plane < distance:
                self.begin = i
                distance = dist_from_plane

                if plane.is_on_plane(self.atoms[(i + dist3) % 6], tolerance):
                    self.has_plane = True


class FiveAtomRing(Ring):
    def __init__(self):
        super().__init__()

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
