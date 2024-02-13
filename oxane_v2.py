from typing import Optional
from rings import SixAtomRing
from molecule import MoleculeType, Conformation
from geometries import Plane, Vector
from atom import Atom
import traceback
from math import pi, acos, degrees

# TODO: Temporary constant
INVALID = 10000000.12345

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
        # 1.) create a main plane from atoms indexed 0, 2, 4 along with axes
        # print(f"Stage 1: main plane")
        axes = [] # deonted as a_i
        for i in range(3):
            axes.append(Vector(self[2 * (i+1)], self[2*i]))
        main_plane = Plane(self[0], self[2], self[4])
        # print(main_plane)
        # 2.) in this plane we then create a normal vector - n
        # print("\nStage 2: normal vector:")
        n = main_plane.normal  # denoted as n
        # print(n)
        # 3.) we create (bond) vectors for all the bonds of atoms, 6 in total
        # print("\nStage 3: bond vectors")
        bond_vectors = []  # denoted as r_i
        for i in range(6):
            v = Vector(self[(i+1)%6], self[i])
            bond_vectors.append(v)
            # print(v)
        # 4.) using bond vectors on either side we create normal (orientation)
        # vectors by crossing them
        # print("\nStage 4: normal vectors")
        normals = [] # denoted as p_i
        for i in range(3):
            nv = bond_vectors[(i*2+1)%6].cross(bond_vectors[i*2])
            normals.append(nv)
            # print(nv)
        # 5.) we create angle of puckering by crossing the normal vector of the
        # far atom of a flap with the vector going between the atoms forming the
        # hinge of a flap to create q_i vector
        # print("\nStage 5: q-vectors")
        q_vectors = []
        for i in range(3):
            q = axes[i].cross(normals[i])
            q_vectors.append(q)
            # print(q)
        # 6.) now, the angle of intersection between vector n and q_i is equal
        # to pi/2 - theta_i which gives positive theta_i when the flap is above
        # the plane and pi/2 + theta_i giving a negative theta_i when flap is
        # under the plane.
        # theta = pi/2 - cos^(-1)[(q_i * n) * (||q_i|| * ||n||)^(-1)]
        # print("\nStage 6: thetas")
        angles = []
        for i in range(3):
            numerator = n.dot(q_vectors[i])
            denominator = (q_vectors[i].length() * n.length())
            inside = numerator / denominator
            ac = acos(inside)
            theta = pi/2 - ac
            angles.append(degrees(theta))
            # print(f"Inside: {numerator}/{denominator} = {inside}, acos: {ac}, theta: {theta}, result: {degrees(theta)}")

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
        if self.is_flat():
            self.conformation = Conformation.Flat
        elif self.is_chair():
            self.conformation = Conformation.Chair
        elif self.is_boat():
            self.conformation = Conformation.Boat
        elif self.is_half_chair():
            self.conformation = Conformation.Half_Chair
        elif self.is_skew():
            self.conformation = Conformation.Skew
        elif self.is_envelope():
            self.conformation = Conformation.Envelope
        else:
            self.conformation = Conformation.Undefined

    def is_flat(self) -> bool:
        distance = 0
        for angle in self.angles:
            if -self.config.o2.t_flat_angle < angle < self.config.o2.t_flat_angle:
                distance += abs(angle)
            else:
                return False
        return True

    def is_chair(self) -> bool:
        distance = 0
        all_same_side = (self.angles[0] < 0 and self.angles[1] < 0 and self.angles[2] < 0) or \
                        (self.angles[0] > 0 and self.angles[1] > 0 and self.angles[2] > 0)
        if not all_same_side:
            return False

        upper_limit = self.config.o2.t_chair_degree + self.config.o2.t_chair_angle
        lower_limit = self.config.o2.t_chair_degree - self.config.o2.t_chair_angle
        for angle in self.angles:
            if lower_limit < abs(angle) < upper_limit:
                distance += abs(angle - self.config.o2.t_chair_degree)
            else:
                return False
        return True


    def is_half_chair(self) -> bool:
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
        sorted_atoms = sorted(self.angles)

        if -self.config.o2.t_envelope_angle < sorted_atoms[1] < self.config.o2.t_envelope_angle:
            ... #deal with double zero case
            bigger = max(abs(sorted_atoms[0]), abs(sorted_atoms[2]))
            smaller = min(abs(sorted_atoms[0]), abs(sorted_atoms[2]))
            if not -self.config.o2.t_envelope_angle < smaller < self.config.o2.t_envelope_angle:
                return False
            if not self.config.o2.t_envelope_b_degree - self.config.o2.t_envelope_angle < bigger < \
                self.config.o2.t_envelope_b_degree + self.config.o2.t_envelope_angle:
                return False
        else:
            signs = sum([-1 if x < 0 else 1 for x in self.angles if x != 0])
            if signs != -1 and signs != 1:
                return False
            sorted_atoms = sorted([abs(x) for x in self.angles])
            small = sorted_atoms[0]
            medium = sorted_atoms[1]
            big = sorted_atoms[2]
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