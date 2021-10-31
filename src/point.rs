use crate::Vector;
use std::{
    fmt,
    ops::{Add, Div, Mul, Neg, Sub},
};
#[derive(Debug, Copy, Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}
impl Add for Point {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Add<Vector> for Point {
    type Output = Self;

    fn add(self, other: Vector) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl Sub for Point {
    type Output = Vector;

    fn sub(self, rhs: Self) -> Vector {
        Vector {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Sub<Vector> for Point {
    type Output = Self;

    fn sub(self, rhs: Vector) -> Self::Output {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Neg for Point {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Mul<Point> for Point {
    type Output = Self;

    fn mul(self, rhs: Point) -> Self {
        Self {
            x: self.x * rhs.x,
            y: self.y * rhs.y,
            z: self.z * rhs.z,
        }
    }
}

impl Mul<f64> for Point {
    // Multiply by scalar
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f64> for Point {
    type Output = Self;

    fn div(self, other: f64) -> Self::Output {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "x: {0:>10} y: {1:>10} z: {2:>10}",
            format!("{0:.5}", self.x),
            format!("{0:.5}", self.y),
            format!("{0:.5}", self.z),
        )
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        ((self.x - other.x).abs() < f64::EPSILON)
            & ((self.y - other.y).abs() < f64::EPSILON)
            & ((self.z - other.z).abs() < f64::EPSILON)
    }
}

#[cfg(test)]
mod point_tests {
    use super::*;
    use crate::eq_f64;

    #[test]
    fn create_point() {
        let point = Point::new(3.0, 4.0, 5.0);
        assert!(eq_f64(point.x, 3.0));
        assert!(eq_f64(point.y, 4.0));
        assert!(eq_f64(point.z, 5.0));
    }
    #[test]
    fn add_two_points() {
        let tuple1 = Point::new(3.0, -2.0, 5.0);
        let tuple2 = Point::new(-2.0, 3.0, 1.0);
        assert_eq!(tuple1 + tuple2, Point::new(1.0, 1.0, 6.0));
    }
    #[test]
    fn add_vector_to_point_return_vector() {
        let p = Point::new(3.0, -2.0, 5.0);
        let v = Vector::new(-2.0, 3.0, 1.0);
        assert_eq!(p + v, Point::new(1.0, 1.0, 6.0));
    }
    #[test]
    fn subtracting_two_points() {
        let p1 = Point::new(3.0, 2.0, 1.0);
        let p2 = Point::new(5.0, 6.0, 7.0);
        assert_eq!(p1 - p2, Vector::new(-2.0, -4.0, -6.0));
    }
    #[test]
    fn negating_point() {
        let p = Point::new(1.0, -2.0, 3.0);
        assert_eq!(-p, Point::new(-1.0, 2.0, -3.0));
    }
    #[test]
    fn multiply_point_by_scalar() {
        let p = Point::new(1.0, -2.0, 3.0);
        assert_eq!(p * 3.5, Point::new(3.5, -7.0, 10.5));
    }
    #[test]
    fn multiply_point_by_fraction() {
        let p = Point::new(1.0, -2.0, 3.0);
        assert_eq!(p * 0.5, Point::new(0.5, -1.0, 1.5));
    }
    #[test]
    fn dividing_a_point_by_scalar() {
        let p = Point::new(1.0, -2.0, 3.0);
        assert_eq!(p / 2.0, Point::new(0.5, -1.0, 1.5));
    }
}
