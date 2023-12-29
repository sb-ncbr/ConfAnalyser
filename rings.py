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


class SixAtomRing(Ring):
    def __init__(self):
        super().__init__()

    def find_plane(self, tolerance: float, dist1: int = 1, dist2: int = 3, dist3: int = 4) -> None:
        distance = float("inf")

        for i in range(6):
            plane = Plane(self.atoms[i % 6], self.atoms[(i + dist1) % 6], self.atoms[(i + dist2) % 6])
            dist_from_plane = abs(plane.distance_from(self.atoms[(i + dist3) % 6]))
            print(f"Atom [{i}]: {self.atoms[i]}")
            print(f"Distance {i}: {dist_from_plane}")
            print(f"Plane: < {plane.a}, {plane.b}, {plane.c}, {plane.d} >")
            if dist_from_plane < distance:
                self.begin = i
                distance = dist_from_plane

                if plane.is_on_plane(self.atoms[(i + dist3) % 6], tolerance):
                    self.has_plane = True

    def index(self, i: int) -> Atom:
        # TODO: Think of a better name? index kinda sus
        """
        Returns ring's atom at a given index while including the `begin` index determined
        by the find_plane function.
        """
        return self.atoms[(self.begin + i) % 6]
