use std::{
    fmt,
    ops::{Add, Div, Mul, Neg, Sub},
};

#[derive(Debug, Copy, Clone)]
pub struct Vector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn magnitude(&self) -> f64 {
        (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt()
    }
    pub fn normalize(&self) -> Self {
        let magnitude = self.magnitude();
        Vector {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
        }
    }
    pub fn dot(&self, other: &Self) -> f64 {
        // Only applicable to Vectors
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross(&self, other: &Self) -> Self {
        // Only applicable to vectors
        Vector::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}
impl Add for Vector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Vector {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Neg for Vector {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Mul<f64> for Vector {
    // Multiply by scalar
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl Div<f64> for Vector {
    type Output = Self;

    fn div(self, other: f64) -> Self::Output {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl fmt::Display for Vector {
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

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        ((self.x - other.x).abs() < f64::EPSILON)
            & ((self.y - other.y).abs() < f64::EPSILON)
            & ((self.z - other.z).abs() < f64::EPSILON)
    }
}

#[cfg(test)]
mod vector_tests {
    use crate::eq_f64;

    use super::*;

    #[test]
    fn create_vector() {
        let v = Vector::new(3.0, 4.0, 5.0);
        assert!(eq_f64(v.x, 3.0));
        assert!(eq_f64(v.y, 4.0));
        assert!(eq_f64(v.z, 5.0));
    }
    #[test]
    fn add_two_vectors_returns_a_vector() {
        let v1 = Vector::new(3.0, -2.0, 5.0);
        let v2 = Vector::new(-2.0, 3.0, 1.0);
        assert_eq!(v1 + v2, Vector::new(1.0, 1.0, 6.0));
    }
    #[test]
    fn subtract_two_vectors() {
        let v1 = Vector::new(3.0, 2.0, 1.0);
        let v2 = Vector::new(5.0, 6.0, 7.0);
        assert_eq!(v1 - v2, Vector::new(-2.0, -4.0, -6.0));
    }
    #[test]
    fn negate_vector() {
        let v = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(-v, Vector::new(-1.0, 2.0, -3.0));
    }
    #[test]
    fn multiply_vector_by_scalar() {
        let p = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(p * 3.5, Vector::new(3.5, -7.0, 10.5));
    }
    #[test]
    fn multiply_vector_by_fraction() {
        let p = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(p * 0.5, Vector::new(0.5, -1.0, 1.5));
    }
    #[test]
    fn dividing_a_vector_by_scalar() {
        let p = Vector::new(1.0, -2.0, 3.0);
        assert_eq!(p / 2.0, Vector::new(0.5, -1.0, 1.5));
    }
    #[test]
    fn compute_magnitude() {
        let v = Vector::new(1.0, 0.0, 0.0);
        assert!(eq_f64(v.magnitude(), 1.0));
    }
    #[test]
    fn compute_magnitude_123() {
        let v = Vector::new(-1.0, -2.0, -3.0);
        assert!(eq_f64(v.magnitude(), 14_f64.sqrt()));
    }
    #[test]
    fn normalize_123() {
        let v = Vector::new(1.0, 2.0, 3.0);
        assert_eq!(
            v.normalize(),
            Vector::new(
                1.0 / 14.0_f64.sqrt(),
                2.0 / 14.0_f64.sqrt(),
                3.0 / 14.0_f64.sqrt()
            )
        );
    }
    #[test]
    fn magn_norm() {
        let v = Vector::new(1.0, 2.0, 3.0);
        let norm = v.normalize();
        assert!(eq_f64(norm.magnitude(), 1.0));
    }
    #[test]
    fn dot_product() {
        let a = Vector::new(1.0, 2.0, 3.0);
        let b = Vector::new(2.0, 3.0, 4.0);
        assert!(eq_f64(a.dot(&b), 20.0));
    }
    #[test]
    fn cross_product() {
        let a = Vector::new(1.0, 2.0, 3.0);
        let b = Vector::new(2.0, 3.0, 4.0);
        assert_eq!(a.cross(&b), Vector::new(-1.0, 2.0, -1.0));
        assert_eq!(b.cross(&a), Vector::new(1.0, -2.0, 1.0));
    }
}
