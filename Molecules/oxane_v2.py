from Molecules.Components.rings import SixAtomRing
from Molecules.Components.molecule import MoleculeType, Conformation
from Molecules.Components.geometries import Plane, Vector
from math import pi, acos, degrees

class Oxane_v2(SixAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Oxane)

        # Set the needed parameters
        self.source_line: list[str] = source_line
        self.oxygen_position: int = -1
        self.angles = None

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(e)

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

    def calculate_properties(self):
        """
        This function calculates main properties of the ring based on which
        the conformations are later decided.

        CalculateProperties(Molecule):
            Atom 1, ..., Atom 6 = Molecule's atoms

            axes = empty list
            FOR i = 1 TO i == 3 DO
                axe = Vector from Atom 2 * (i + 1) to Atom 2 * i
                add axe vector to axes
            END FOR

            mainPlaneNormal = Normal vector to Axe 1 and Axe 2 made by crossing them

            bondVectors = empty list
            FOR i = 1 TO i == 6 DO
                bondVector = Vector from Atom (i+1) modulo 6 to Atom i
                add `bondVector` into bondVectors
            END FOR

            normals = empty list
            FOR i = 1 TO i == 3 DO
                normalVector = Vector made by crossing bond vector (i * 2 + 1) mod 6
                               with bond vector i * 2
                add `normalVector` into normals
            END FOR

            qVectors = empty list
            FOR i = 1 TO i == 3 DO
                qVector = Vector made by crossing axe i with normal vector i
                add `qVector` into qVectors
            END FOR

            angles = empty list
            FOR i = 1 TO i == 3 DO
                numerator = dot product of `mainPlaneNormal` and q vector i
                denominator = size of q vector i * size of `mainPlaneNormal`
                acos = arccosine of `numerator` / `denominator`
                theta = pi / 2 - acos
                degrees = theta converted into degrees value
                add `degrees` into angles
            END FOR

            set Molecule's angles variable to `angles`
        """
        # 1.) create a main plane from atoms indexed 0, 2, 4 along with axes
        axes = [] # deonted as a_i
        for i in range(3):
            axes.append(Vector(self[2 * (i+1)], self[2*i]))
        main_plane = Plane(self[0], self[2], self[4])
        # 2.) in this plane we then create a normal vector - n
        n = main_plane.normal  # denoted as n
        # 3.) we create (bond) vectors for all the bonds of atoms, 6 in total
        bond_vectors = []  # denoted as r_i
        for i in range(6):
            v = Vector(self[(i+1)%6], self[i])
            bond_vectors.append(v)
        # 4.) using bond vectors on either side we create normal (orientation)
        # vectors by crossing them
        normals = [] # denoted as p_i
        for i in range(3):
            nv = bond_vectors[(i*2+1)%6].cross(bond_vectors[i*2])
            normals.append(nv)
        # 5.) we create angle of puckering by crossing the normal vector of the
        # far atom of a flap with the vector going between the atoms forming the
        # hinge of a flap to create q_i vector
        q_vectors = []
        for i in range(3):
            q = axes[i].cross(normals[i])
            q_vectors.append(q)
        # 6.) now, the angle of intersection between vector n and q_i is equal
        # to pi/2 - theta_i which gives positive theta_i when the flap is above
        # the plane and pi/2 + theta_i giving a negative theta_i when flap is
        # under the plane.
        # theta = pi/2 - cos^(-1)[(q_i * n) * (||q_i|| * ||n||)^(-1)]
        angles = []
        for i in range(3):
            inside_acos = n.dot(q_vectors[i]) / (q_vectors[i].length() * n.length())
            theta = pi/2 - acos(inside_acos)
            angles.append(degrees(theta))

        self.angles = angles
        # 7.) At this point we've produced the angles theta_i of the three flaps
        # that define shape of the molecule. Based on those angles, whether or
        # not they fall within some limits and (I assume) what is the
        # relationship between their signs, we can determine shape of the molecule.

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """
        self.calculate_properties()
        if self.is_chair():
            self.conformation = Conformation.Chair
        elif self.is_boat():
            self.conformation = Conformation.Boat
        elif self.is_half_chair():
            self.conformation = Conformation.Half_Chair
        elif self.is_skew():
            self.conformation = Conformation.Skew
        elif self.is_envelope():
            self.conformation = Conformation.Envelope
        elif self.is_flat():
            self.conformation = Conformation.Flat
        else:
            self.conformation = Conformation.Undefined

    def is_flat(self) -> bool:
        """
        Validate whether the angles correspond to the values expected for
        the flat conformation.

        IsFlat(Molecule, Tolerances):
            FOR EACH angle IN Molecule's Angles DO
                bottomLimit = Tolerance for flat angle * -1
                upperLimit = Tolerance for flat angle
                IF angle is smaller than bottomLimit or angle is bigger than upperLimit
                THEN
                    RETURN FALSE
                FI
            END FOR
            RETURN TRUE
        :return:
        """
        for angle in self.angles:
            if not -self.config.o2.t_flat_angle < angle < self.config.o2.t_flat_angle:
                return False
        return True

    def is_chair(self) -> bool:
        """
        Validate whether the angles and their properties correspond to the
        expected values for chair conformation.

        IsChair(Molecule, Tolerances):
            IF all three angles don't have same sign
            THEN
                RETURN FALSE
            FI

            lowerLimit = expected value of chair angle - chair angle tolerance
            upperLimit = expected value of chair angle + chair angle tolerance

            FOR EACH angle IN Molecule's angles DO
                IF angle is lower than lowerLimit OR angle is bigger than upperLimit
                THEN
                    RETURN FALSE
                FI
            END FOR
            RETURN TRUE
        """
        all_same_side = (self.angles[0] < 0 and self.angles[1] < 0 and self.angles[2] < 0) or \
                        (self.angles[0] > 0 and self.angles[1] > 0 and self.angles[2] > 0)
        if not all_same_side:
            return False

        upper_limit = self.config.o2.t_chair_degree + self.config.o2.t_chair_angle
        lower_limit = self.config.o2.t_chair_degree - self.config.o2.t_chair_angle
        for angle in self.angles:
            if not lower_limit < abs(angle) < upper_limit:
                return False
        return True


    def is_half_chair(self) -> bool:
        """
        Validate whether the angles and their properties correspond to the
        expected values for half chair conformation.

        IsHalfChair(Molecule, Tolerances):
            IF all three angles have same sign OR any of angles is 0
            THEN
                RETURN FALSE
            FI

            sortedAngles = sorted absolute values of Molecule's angles from
                           smalles to biggest
            small = angle from sortedAngles at index 1
            medium = angle from sortedAngles at index 2
            big = angle from sortedAngles at index 3

            IF small is smaller than lower bound for small angle OR
               small is bigger than upper bound for small angle  OR
               medium is smaller than lower bound for medium angle OR
               medium is bigger than upper bound for medium angle OR
               big is smaller than lower bound for big angle OR
               big is bigger than upper bound for big angle
            THEN
                RETURN FALSE
            FI
            RETURN TRUE
        """
        # we check if it's true that result have two negatives and one positive
        # or two positives and one negative result. Otherwise, not half chair
        signs = sum([-1 if x < 0 else 1 for x in self.angles if x != 0])
        if signs != -1 and signs != 1:
            return False

        sorted_atoms = sorted([abs(x) for x in self.angles])
        small = sorted_atoms[0]
        medium = sorted_atoms[1]
        big = sorted_atoms[2]
        if not self.config.o2.t_half_s_degree - self.config.o2.t_half_angle < small <  \
            self.config.o2.t_half_s_degree + self.config.o2.t_half_angle:
            return False
        if not self.config.o2.t_half_m_degree - self.config.o2.t_half_angle < medium < \
            self.config.o2.t_half_m_degree + self.config.o2.t_half_angle:
            return False
        if not self.config.o2.t_half_b_degree - self.config.o2.t_half_angle < big < \
            self.config.o2.t_half_b_degree + self.config.o2.t_half_angle:
            return False
        return True

    def is_boat(self) -> bool:
        """
        Validate whether the angles and their properties correspond to the
        expected values for the boat conformation.

        IsBoat(Molecule, Tolerances):
            IF all three angles have same sign OR any of angles is 0
            THEN
                RETURN FALSE
            FI

            sortedAngles = sorted absolute values of Molecule's angles from
                           smallest to biggest
            small = angle from sortedAngles at index 1
            medium = angle from sortedAngles at index 2
            big = angle from sortedAngles at index 3

            IF small is smaller than lower bound for smaller angle OR
               small is bigger than upper bound for smaller angle  OR
               medium is smaller than lower bound for smaller angle OR
               medium is bigger than upper bound for smaller angle OR
               big is smaller than lower bound for big angle OR
               big is bigger than upper bound for big angle
            THEN
                RETURN FALSE
            FI
            RETURN TRUE
        """
        # we check if it's true that result have two negatives and one positive
        # or two positives and one negative result. Otherwise, not half chair
        signs = sum([-1 if x < 0 else 1 for x in self.angles if x != 0])
        if signs != -1 and signs != 1:
            return False

        sorted_atoms = sorted([abs(x) for x in self.angles])
        small = sorted_atoms[0]
        medium = sorted_atoms[1]
        big = sorted_atoms[2]
        if not self.config.o2.t_boat_s_degree - self.config.o2.t_boat_angle < small < \
               self.config.o2.t_boat_s_degree + self.config.o2.t_boat_angle:
            return False
        if not self.config.o2.t_boat_s_degree - self.config.o2.t_boat_angle < medium < \
               self.config.o2.t_boat_s_degree + self.config.o2.t_boat_angle:
            return False
        if not self.config.o2.t_boat_b_degree - self.config.o2.t_boat_angle < big < \
               self.config.o2.t_boat_b_degree + self.config.o2.t_boat_angle:
            return False
        return True

    def is_envelope(self) -> bool:
        """
        Validate whether the angles and their properties correspond to the
        expected values for the envelope conformation.

        IsEnvelope(Molecule, Tolerances):
            sortedAngles = sorted values of Molecule's angles from
                           smallest to biggest
            middle = angle from sortedAngles at index 2
            IF middle is bigger than angle tolerance for envelope * -1 AND
               middle is smaller than angle tolerance for envelope
            THEN
                small = smaller value of absolute values of angles sortedAngles 1 and sortedAngles 3
                big = bigger value of absolute values of angles sortedAngles 1 and sortedAngles 3

                IF small is smaller than angle tolerance for envelope * -1 OR
                   small is bigger than angle tolerance for envelope
                THEN
                    RETURN FALSE
                FI

                IF big is smaller than lower limit of envelope tolerance for big degree OR
                   big is bigger than upper limit of envelope's tolerance for big degree
                THEN
                    RETURN FALSE
                FI
            ELSE
                mid = absolute value of angle from sortedAngles at index 2

                IF first two angles in sortedAngles are negative and third is positive
                THEN
                    small = absolute value of angle from sortedAngles at index 3
                    big = absolute value of angle from sortedAngles at index 1
                ELSE IF first angle in sortedAngles is negative and last two are positive
                THEN
                    small = absolute value of angle from sortedAngles at index 1
                    big = absolute value of angle from sortedAngles at index 3
                ELSE
                    RETURN FALSE
                FI

                IF small is smaller than lower limit of envelope's tolerance for small angles OR
                   small is bigger than upper limit of envelope's tolerance for small angles OR
                   mid or big are smaller than lower limit for envelope's tolerance for middle angles OR
                   mid or big are bigger than upper limit for envelope's tolerance for middle angles
                THEN
                    RETURN FALSE
                FI
            FI
            RETURN TRUE
        """
        sorted_atoms = sorted(self.angles)

        # double zero case
        if -self.config.o2.t_envelope_angle < sorted_atoms[1] < self.config.o2.t_envelope_angle:
            bigger = max(abs(sorted_atoms[0]), abs(sorted_atoms[2]))
            smaller = min(abs(sorted_atoms[0]), abs(sorted_atoms[2]))
            if not -self.config.o2.t_envelope_angle < smaller < self.config.o2.t_envelope_angle:
                return False
            if not self.config.o2.t_envelope_b_degree - self.config.o2.t_envelope_angle < bigger < \
                self.config.o2.t_envelope_b_degree + self.config.o2.t_envelope_angle:
                return False
        else: # non-zero case
            medium = abs(sorted_atoms[1])
            if sorted_atoms[0] < 0 and sorted_atoms[1] < 0 and sorted_atoms[2] > 0:
                # --+
                small = abs(sorted_atoms[2])
                big = abs(sorted_atoms[0])
            elif sorted_atoms[0] < 0 and sorted_atoms[1] > 0 and sorted_atoms[2] > 0:
                # -++
                small = abs(sorted_atoms[0])
                big = abs(sorted_atoms[2])
            else:
                return False

            if not self.config.o2.t_envelope_s_degree - self.config.o2.t_envelope_angle < small < \
                self.config.o2.t_envelope_s_degree + self.config.o2.t_envelope_angle:
                return False
            if not self.config.o2.t_envelope_m_degree - self.config.o2.t_envelope_angle < medium < \
                self.config.o2.t_envelope_m_degree + self.config.o2.t_envelope_angle:
                return False
            if not self.config.o2.t_envelope_m_degree - self.config.o2.t_envelope_angle < big < \
                self.config.o2.t_envelope_m_degree + self.config.o2.t_envelope_angle:
                return False
        return True


    def is_skew(self) -> bool:
        """
        Validate whether the angles and their properties correspond to the
        expected values for the skew conformation.

        IsSkew(Molecule, Tolerances):
            sortedAngles = sorted values of Molecule's angles from
                           smallest to biggest
            small = angle from sortedAngles at index 1
            mid = angle from sortedAngles at index 2
            big = angle from sortedAngles at index 3

            IF small is smaller than lower bound for skew's angle tolerance OR
               small is bigger than upper bound for skew's angle tolerance OR
               mid is smaller than skew's angle tolerance * -1 OR
               mid is bigger than skew's angle tolerance OR
               big is smaller than lower bound for skew's angle tolerance OR
               big is bigger than upper bound for skew's angle tolerance
            THEN
                RETURN FALSE
            FI
            RETURN TRUE
        """
        sorted_atoms = sorted(self.angles)
        small = sorted_atoms[0]
        medium = sorted_atoms[1]
        big = sorted_atoms[2]
        if small >= 0 or \
            self.config.o2.t_skew_degree - self.config.o2.t_skew_angle < abs(small) < \
            self.config.o2.t_skew_degree + self.config.o2.t_skew_angle:
            return False
        if -self.config.o2.t_skew_angle < medium < self.config.o2.t_skew_angle:
            return False
        if big <= 0 or \
            self.config.o2.t_skew_degree - self.config.o2.t_skew_angle < big < \
            self.config.o2.t_skew_degree + self.config.o2.t_skew_angle:
            return False
        return True