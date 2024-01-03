from geometries import Vector, Point
from math import degrees, acos, atan2


def radians_to_degrees(value: float) -> float:
    return degrees(value)


def vector_angle(vector1: Vector, vector2: Vector) -> float:
    return degrees(acos((vector1.normalize().dot(vector2.normalize()))))


def dihedral_angle(p1: Point, p2: Point, p3: Point, p4: Point) -> float:
    u = Vector(p2, p1).cross(Vector(p2, p3))
    v = Vector(p3, p2).cross(Vector(p3, p4))
    unit_vector = Vector(p2, p3).normalize()

    return degrees(atan2(u.cross(v).dot(unit_vector), u.dot(v)))

