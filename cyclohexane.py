from rings import SixAtomRing
from atom import Atom
from geometries import Plane
from angle import dihedral_angle
from molecule import MoleculeType, Conformation


class Cyclohexane(SixAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Cyclohexane)

        # Set the needed parameters
        self.source_file: list[str] = source_line

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(f"Error1: {e}")

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.

        Analyze(Molecule, Tolerances...):
            FindIdealPlane(ToleranceIn)

            rightPlane = Plane between Atom 1, Atom 2 and Atom 4
            leftPlane = Plane between Atom 1, Atom 2 and Atom 5

            IF Molecule has ideal plane THEN
                IF IsFlat(Molecule, Tolerances) THEN
                    conformation = Flat
                ELSE IF IsHalfChair(Molecule, Tolerances) THEN
                    conformation = Half Chair
                ELSE IF IsChair(Molecule, Tolerances) THEN
                    conformation = Chair
                ELSE IF IsBoat(Molecule, Tolerances) THEN
                    conformation = Boat
                ELSE
                    conformation = Undefined
                FI
            ELSE
                IF IsTwistedBoat(Molecule, Tolerances) THEN
                    conformation = Twisted Boat
                ELSE
                    conformation = Undefined
                FI
            FI
        """

        self.find_plane(self.config.ch.t_in)

        self.right_plane = Plane(self[0], self[1], self[3])
        self.left_plane = Plane(self[0], self[1], self[4])

        if self.has_plane:
            if self.is_flat():
                self.conformation = Conformation.Flat
            elif self.is_half_chair():
                self.conformation = Conformation.Half_Chair
            elif self.is_chair():
                self.conformation = Conformation.Chair
            elif self.is_boat():
                self.conformation = Conformation.Boat
            else:
                self.conformation = Conformation.Undefined
        else:
            if self.is_twisted_boat():
                self.conformation = Conformation.Twisted_Boat
            else:
                self.conformation = Conformation.Undefined

    def is_flat(self):
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.

        Is shared by all atoms

        IsFlat(Molecule, FlatTolerance):
            Atom 1, ..., Atom 6 = Molecule Atoms

            IF distance between rightPlane and Atom 3 < FlatTolerance AND
               distance between rightPlane and Atom 6 < FlatTolerance AND
               distance between leftPlane and Atom 3 < FlatTolerance AND
               distance between leftPlane and Atom 6 < FlatTolerance
            THEN
                RETURN TRUE
            ELSE
                RETURN FALSE
            FI
        """
        return (self.right_plane.is_on_plane(self[2], self.config.ch.t_flat_in) and
                self.right_plane.is_on_plane(self[5], self.config.ch.t_flat_in) and
                self.left_plane.is_on_plane(self[2], self.config.ch.t_flat_in) and
                self.left_plane.is_on_plane(self[5], self.config.ch.t_flat_in))

    def is_half_chair(self) -> bool:
        """
        Decides whether this molecule's conformation is half chair.
        Half chair conformation of molecule is determined by having all but one atom in one plane.

        IsHalfChair(Molecule, FlatTolerance, XTolerance):
            Atom 1, ..., Atom 6 = Molecule Atoms

            mainPlane = Molecule right plane
            atom3Distance = distance between mainPlane and Atom 3
            atom6Distance = distance between mainPlane and Atom 6

            atom5IsOnTheMainPlane = distance between mainPlane and Atom 5 < FlatTolerance

            onlyOneAtomIsOutOfMainPlane = atom3Distance < FlatTolerance XOR
                                          atom6Distance < FlatTolerance

            onlyOneAtomIsFarFromPlane = atom3Distance > XTolerance XOR
                                        atom5Distance > XTolerance

            IF onlyOneAtomIsOutOfMainPlane AND
               onlyOneAtomIsFarFromPlane AND
               atom5IsOnTheMainPlane
            THEN
                RETURN TRUE
            ELSE
                RETURN FALSE
            FI
        """
        right_dist = self.right_plane.true_distance_from(self[2])
        left_dist = self.right_plane.true_distance_from(self[5])
        return (self.right_plane.is_on_plane(self[2], self.config.ch.t_flat_in) !=
                self.right_plane.is_on_plane(self[5], self.config.ch.t_flat_in)) and \
                self.right_plane.is_on_plane(self[4], self.config.ch.t_flat_in) and \
                ((right_dist > self.config.ch.t_out) != (left_dist > self.config.ch.t_out))

    def is_chair(self) -> bool:
        """
        Decides whether the molecule's conformation is a chair.
        Chair conformation is determined by having two opposite atoms on
        opposite sides of the main molecule plane.

        IsChair(Molecule, ToleranceOut):
            Atom 1, ..., Atom 6 = Molecule Atoms

            mainPlane = Molecule right plane

            IF distance between mainPlane and Atom 3 > ToleranceOut AND
               distance between mainPlane and Atom 6 > ToleranceOut AND
               Atom 3 and Atom 6 are on opposite sides of mainPlane
            THEN
                RETURN TRUE
            ELSE
                RETURN FALSE
            FI
        """
        right_dist = self.right_plane.true_distance_from(self[2])
        left_dist = self.right_plane.true_distance_from(self[5])
        return right_dist > self.config.ch.t_out and \
               left_dist > self.config.ch.t_out and  \
               self.right_plane.are_opposite_side(self[2], self[5])

    def is_boat(self) -> bool:
        """
        Decides whether the molecule's conformation is a boat.
        Boat conformation is determined by having two opposite atoms on the
        same side of a plane.

        IsBoat(Molecule, ToleranceOut):
            Atom 1, ..., Atom 6 = Molecule Atoms

            mainPlane = Molecule right plane

            IF distance between mainPlane and Atom 3 > ToleranceOut AND
               distance between mainPlane and Atom 6 > ToleranceOut AND
               Atom 3 and Atom 6 are the same side of mainPlane
            THEN
                RETURN TRUE
            ELSE
                RETURN FALSE
            FI
        """
        right_dist = self.right_plane.true_distance_from(self[2])
        left_dist = self.right_plane.true_distance_from(self[5])
        return right_dist > self.config.ch.t_out and \
               left_dist > self.config.ch.t_out and \
               self.right_plane.are_same_side(self[2], self[5])

    def is_twisted_boat(self) -> bool:
        """
        Determine whether the molecule's conformation is a twisted boat.
        Twisted boat conformation is determined by having no plane within the ring.

        IsTwistedBoat(Molecule, ToleranceAngle, ToleranceTwBoatAngle, ToleranceTwOut):
            Atom 1, ..., Atom 6 = Molecule Atoms

            atom3Distance = distance of Atom 3 from Molecule's right plane
            atom6Distance = distance of Atom 6 from Molecule's left plane
            twistAngle = absolute value of dihedral angle of Atom 2, Atom 4, Atom 5, Atom 1
            IF ToleranceTwBoatAngle - ToleranceAngle < twistAngle  AND
               twistAngle < ToleranceTwBoatAngle + ToleranceAngle AND
               absolute value of atom3Distance > ToleranceTwOut AND
               absolute value of atom6Distance > ToleranceTwOut AND
               atom3Distance * atom5Distance > 0
            THEN
                RETURN TRUE
            ELSE
                RETURN FALSE
            FI
        """
        right_dist = self.right_plane.signed_distance_from(self[2])
        left_dist = self.left_plane.signed_distance_from(self[5])
        tw_angle = abs(dihedral_angle(self[1], self[3], self[4], self[0]))
        return ((self.config.ch.t_tw_boat_angle - self.config.ch.t_angle) < tw_angle < (self.config.ch.t_tw_boat_angle + self.config.ch.t_angle)
                and abs(right_dist) > self.config.ch.t_tw_out
                and abs(left_dist) > self.config.ch.t_tw_out
                and (right_dist * left_dist) > 0)

