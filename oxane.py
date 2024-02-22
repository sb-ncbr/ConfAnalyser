from typing import Optional
from rings import SixAtomRing
from molecule import MoleculeType, Conformation
from geometries import Plane
from atom import Atom
import traceback

# constants for denoting position of out-of-plane atoms
ABOVE = 0
UNDER = 1


class OutOfPlaneAtom:
    def __init__(self):
        self.atom: Optional[Atom] = None
        self.position: Optional[int] = None
        self.presence: Optional[bool] = None


class Oxane(SixAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Oxane)

        # Set the needed parameters
        self.source_line: list[str] = source_line
        self.oxygen_position: int = -1
        # Currently useless??
        self.out_of_plane_atoms: list[OutOfPlaneAtom] = [OutOfPlaneAtom(), OutOfPlaneAtom()]

        self.updated_version = None

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(f"Error 3: {e}")

    def __str__(self):
        #angles = ", ".join(str(x) for x in self.updated_version.angles)
        return (f"{self.file_name}: {self.conformation.name.upper()}")

    def validate_atoms(self) -> None:
        """
        Oxane-specific validator for atoms
        """
        atom_count = self.get_atom_count()
        new_lst = [None for _ in range(atom_count)]
        found_oxygen = False
        for atom in self.atoms:
            for i in range(atom_count):
                if atom.name in self.names[self.ligand][i]:
                    if new_lst[i] is not None:
                        print(f"{atom.name} atom found twice!")
                        self.is_valid = False
                        return
                    new_lst[i] = atom
                    if atom.element_symbol == "O":
                        if found_oxygen:
                            self.is_valid = False
                            print("Oxygen found twice!")
                            return
                        found_oxygen = True
                        self.oxygen_position = i
                    # break  # TODO: this break is here in other 3 molecules but according to C++ code, here it's not
        if None in new_lst:
            print("Not all atoms were found")
            self.is_valid = False
            return
        if not found_oxygen:
            print("Unable to find oxygen atom position using element_name PDB field.")
        self.atoms = new_lst

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """

        if self.find_plane(self.config.o.t_in):
            ... # flat, chair, boat, envelope
            if self.is_flat():
                self.conformation = Conformation.Flat
            elif self.is_chair():
                self.conformation = Conformation.Chair
            elif self.is_boat():
                self.conformation = Conformation.Boat
            elif self.is_envelope():
                self.conformation = Conformation.Envelope
            else:
                self.conformation = Conformation.Undefined

        elif self.find_plane(self.config.o.t_in, 1, 2, 3):
            ... # half chair
            if self.is_half_chair():  # 1, 2, 3
                self.conformation = Conformation.Half_Chair
            else:
                self.conformation = Conformation.Undefined
        elif self.find_plane(self.config.o.t_in, 1, 2, 4):
            ... # skew
            if self.is_skew():  # 1, 2, 4
                self.conformation = Conformation.Skew
            else:
                self.conformation = Conformation.Undefined
        else:
            self.conformation = Conformation.Undefined


    def is_flat(self) -> bool:
        left_plane = Plane(self[0], self[1], self[4])
        right_plane = Plane(self[0], self[1], self[3])

        return (left_plane.is_on_plane(self[2], self.config.o.t_in) and
                left_plane.is_on_plane(self[5], self.config.o.t_in) and
                right_plane.is_on_plane(self[2], self.config.o.t_in) and
                right_plane.is_on_plane(self[5], self.config.o.t_in))

    def is_chair(self) -> bool:
        left_plane = Plane(self[0], self[1], self[4])
        right_plane = Plane(self[0], self[1], self[3])

        distance_1 = right_plane.distance_from(self[2])
        distance_2 = left_plane.distance_from(self[2])
        right_distance = distance_1 if abs(distance_1) < abs(distance_2) else distance_2

        distance_3 = right_plane.distance_from(self[5])
        distance_4 = left_plane.distance_from(self[5])
        left_distance = distance_3 if abs(distance_3) < abs(distance_4) else distance_4

        is_chair = (abs(right_distance) > self.config.o.t_out and
                    abs(left_distance) > self.config.o.t_out and
                    (right_distance * left_distance < 0))

        if is_chair:
            self.set_oopa(right_distance, left_distance, 2, 5)

        return is_chair

    def is_half_chair(self) -> bool:
        left_plane = Plane(self[0], self[1], self[3])
        right_plane = Plane(self[0], self[1], self[2])

        right_dist = min(right_plane.signed_distance_from(self[4]),
                         left_plane.signed_distance_from(self[4]),
                         key=abs)
        left_dist = min(right_plane.signed_distance_from(self[5]),
                        left_plane.signed_distance_from(self[5]),
                        key=abs)

        is_half_chair = (abs(right_dist) > self.config.o.t_out and
                         abs(left_dist) > self.config.o.t_out and
                         (right_dist * left_dist < 0))

        if is_half_chair:
            self.set_oopa(right_dist, left_dist, 4, 5)

        return is_half_chair

    def is_boat(self) -> bool:
        left_plane = Plane(self[0], self[1], self[4])
        right_plane = Plane(self[0], self[1], self[3])

        dist1 = right_plane.distance_from(self[2])
        dist2 = left_plane.distance_from(self[2])
        dist3 = right_plane.distance_from(self[5])
        dist4 = left_plane.distance_from(self[5])
        right_distance = dist1 if abs(dist1) < abs(dist2) else dist2
        left_distance = dist3 if abs(dist3) < abs(dist4) else dist4

        is_boat = (abs(right_distance) > self.config.o.t_out and
                   abs(left_distance) > self.config.o.t_out and
                   (right_distance * left_distance) > 0)

        if is_boat:
            self.set_oopa(right_distance, left_distance, 2, 5)

            # TODO: Sort of atoms? Why? What for?
            if self.out_of_plane_atoms[0].atom.name > self.out_of_plane_atoms[1].atom.name:
                self.out_of_plane_atoms[0], self.out_of_plane_atoms[1] = (
                    self.out_of_plane_atoms[1], self.out_of_plane_atoms[0])
        return is_boat

    def set_oopa(self, rd, lf, i1, i2):
        self.out_of_plane_atoms[0].presence = True
        self.out_of_plane_atoms[1].presence = True

        self.out_of_plane_atoms[0].position = ABOVE if rd > 0 else UNDER
        self.out_of_plane_atoms[1].position = ABOVE if lf > 0 else UNDER

        self.out_of_plane_atoms[0].atom = self[self.get_index_by_oxygen(i1)]
        self.out_of_plane_atoms[1].atom = self[self.get_index_by_oxygen(i2)]

    def is_envelope(self) -> bool:
        """
        IsEnvelope(Molecule, ToleranceIn):
            Atom 1, ..., Atom 6 = Molecule's atoms

            leftPlane = Plane of Atom 1, Atom 2 and Atom 5
            rightPlane = Plane of Atom 1, Atom 2 and Atom 4

            rightDistance = minimum of distances of Atom 3 from rightPlane and leftPlane
            leftDistance = minimum of distances of Atom 6 from rightPlane and leftPlane

            atom3IsOnMainPlane = distance of Atom 3 from leftPlane < ToleranceIn AND
                                 distance of Atom 3 from rightPlane < ToleranceIn

            atom6IsOnMainPlane = distance of Atom 6 from leftPlane < ToleranceIn AND
                                 distance of Atom 6 from rightPlane < ToleranceIn

            atom3IsOnBothPlanes = distance of Atom 3 from leftPlane < ToleranceIn ==
                                  distance of Atom 3 from rightPlane < ToleranceIn

            atom6IsOnBothPlanes = distance of Atom 6 from leftPlane < ToleranceIn ==
                                  distance of Atom 6 from rightPlane < ToleranceIn

            onlyOneAtomIsOnMainPlane = atom3IsOnMainPlane XOR atom6IsOnMainPlane
            bothAtomsAreEitherOnOrOffMainPlane = atom3IsOnBothPlanes AND atom6IsOnBothPlanes

            IF onlyOneAtomIsOnMainPlane AND bothAtomsAreEitherOnOrOffMainPlane
            THEN
                RETURN TRUE
            ELSE
                RETURN FALSE
            FI

        """
        left_plane = Plane(self[0], self[1], self[4])
        right_plane = Plane(self[0], self[1], self[3])
        right_distance = min(right_plane.true_distance_from(self[2]),
                             left_plane.true_distance_from(self[2]))
        left_distance = min(right_plane.true_distance_from(self[5]),
                            left_plane.true_distance_from(self[5]))

        left_on_plane_2 = left_plane.is_on_plane(self[2], self.config.o.t_in)
        right_on_plane_2 = right_plane.is_on_plane(self[2], self.config.o.t_in)
        left_on_plane_5 = left_plane.is_on_plane(self[5], self.config.o.t_in)
        right_on_plane_5 = right_plane.is_on_plane(self[5], self.config.o.t_in)

        is_envelope = (((left_on_plane_2 and right_on_plane_2) !=
                       (left_on_plane_5 and right_on_plane_5)) and
                       ((left_on_plane_2 ==right_on_plane_2) and
                        (left_on_plane_5 ==right_on_plane_5)))

        if is_envelope:
            self.set_oopa(right_distance, left_distance, 2, 5)
            self.out_of_plane_atoms[0].presence = not left_plane.is_on_plane(self[2], self.config.o.t_in)
            self.out_of_plane_atoms[1].presence = not left_plane.is_on_plane(self[5], self.config.o.t_in)

        return is_envelope

    def is_skew(self) -> bool:
        left_plane = Plane(self[0], self[1], self[4])
        right_plane = Plane(self[0], self[1], self[2])
        right_distance = min(right_plane.distance_from(self[3]),
                             left_plane.distance_from(self[3]),
                             key=abs)
        left_distance = min(right_plane.distance_from(self[5]),
                            left_plane.distance_from(self[5]),
                            key=abs)

        is_skew = ((abs(right_distance) > self.config.o.t_out and
                   abs(left_distance) > self.config.o.t_out) and
                   (right_distance * left_distance < 0))

        if is_skew:
            self.set_oopa(right_distance, left_distance, 3, 5)

        return is_skew

    def get_index_by_oxygen(self, delta_begin):
        """
        /* Convert to numbering relative to oxygen atom.
         * Oxygen atom is marked as 6, carbons will be marked 1 to 5 */
        short unsigned Oxane::get_index_by_oxygen(short int delta_begin)  {
            int delta_oxygen = 6 - (oxygen_position);  // real O index vs. desired mark 6
            return ((begin + delta_begin   // atom number relative to `begin`
                    + delta_oxygen         // shifted relative to oxygen position
                    - 1)                   // -1 so that 6 stays 6 after modulo 6
                    % 6 + 1);              // keep it in range 1 to 6
        }

        """
        delta_oxygen = 6 - self.oxygen_position
        return (self.begin + delta_begin + delta_oxygen - 1) % 6 + 1
