from __future__ import annotations
from math import sqrt
from typing import Union


class Point:
    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
        """
        A class representing a single point in a 3-dimensional space, consisting of
        three float numbers.
        """
        self.x = x
        self.y = y
        self.z = z

    def distance_from(self, other: Point) -> float:
        """
        Calculates the distance of this point from the other point
        Formula: sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)

        :param other: a point to which we calculate distance
        :return: absolute value of a distance between points represented as a floating point number
        """
        return abs(sqrt((other.x - self.x) ** 2 + (other.y - self.y) ** 2 + (other.z - self.z) ** 2))


class Vector(Point):
    def __init__(self, x: Union[float, Point] = 0.0, y: Union[float, Point] = 0.0, z: float = 0.0):
        # Creating vector from a single point
        if isinstance(x, Point) and isinstance(y, float) and isinstance(z, float):
            p1 = x
            super().__init__(p1.x, p1.y, p1.z)

        # Creating vector from two points
        elif isinstance(x, Point) and isinstance(y, Point) and isinstance(z, float):
            p1, p2 = x, y
            super().__init__(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)

        # Creatubg vectir from three coordinates
        elif isinstance(x, float) and isinstance(y, float) and isinstance(z, float):
            super().__init__(x, y, z)

    def from_point(self, point: Point) -> Vector:
        self.x = point.x
        self.y = point.y
        self.z = point.z

        return self

    def from_points(self, point1: Point, point2: Point) -> Vector:
        self.x = point1.x - point2.x
        self.y = point1.y - point2.y
        self.z = point1.z - point2.z

        return self

    def __add__(self, other: Vector) -> Vector:
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vector) -> Vector:
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other: float) -> Vector:
        return Vector(self.x * other, self.y * other, self.z * other)

    # TODO: If required, implement __divmod__ too

    def length(self) -> float:
        """
        Calculates the length/size of this vector
        """
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self) -> Vector:
        """
        Calculates a normalized, unit vector of this vector and returns it
        """
        length = self.length()
        return Vector(self.x / length, self.y / length, self.z / length)

    def dot(self, other) -> float:
        """
        Calculates a dot product between this and the other vector
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other) -> Vector:
        """
        Calculates a cross product between this and the other vector
        """
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)


class Plane:
    def __init__(self, point1: Point, point2: Point, point3: Point):
        self.u = Vector().from_points(point2, point1)
        self.v = Vector().from_points(point2, point3)
        self.normal = self.u.cross(self.v)
        self.d = -(self.normal.x * point1.x + self.normal.y * point1.y + self.normal.z * point1.z)
        self.a = self.normal.x
        self.b = self.normal.y
        self.c = self.normal.z

    def distance_from(self, point: Point) -> float:
        """
        Calculates a distance of a given point from this plane
        """
        return ((self.a * point.x + self.b * point.y + self.c * point.z + self.d) /
                sqrt(self.a ** 2 + self.b ** 2 + self.c ** 2))

    # TODO: tolerance may6 be globalizes as an input parameter from file
    def is_on_plane(self, point: Point, tolerance: float) -> bool:
        """
        Decides whether a given point is located on this plane.
        Threshold determined by passed `tolerance` parameter
        """
        return abs(self.distance_from(point)) <= tolerance
