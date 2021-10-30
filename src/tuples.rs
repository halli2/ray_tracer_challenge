use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

#[allow(unused_macros)]
#[macro_export]
macro_rules! point {
    {$x:expr, $y:expr, $z:expr} => (Tuple::new_point($x, $y, $z))
}
#[allow(unused_macros)]
#[macro_export]
macro_rules! vector {
    {$x:expr, $y:expr, $z:expr} => (Tuple::new_vector($x, $y, $z))
}

#[allow(dead_code)]
impl Tuple {
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Tuple { x, y, z, w }
    }
    pub fn new_point(x: f64, y: f64, z: f64) -> Self {
        Tuple { x, y, z, w: 1.0 }
    }
    pub fn new_vector(x: f64, y: f64, z: f64) -> Self {
        Tuple { x, y, z, w: 0.0 }
    }
    pub fn is_point(&self) -> bool {
        (self.w - 1.0).abs() < f64::EPSILON
    }
    pub fn is_vector(&self) -> bool {
        self.w == 0.0
    }
    pub fn magnitude(&self) -> f64 {
        (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0) + self.w.powf(2.0)).sqrt()
    }
    pub fn normalize(&self) -> Self {
        let magnitude = self.magnitude();
        Tuple {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
            w: self.w / magnitude,
        }
    }
    pub fn dot(&self, other: &Self) -> f64 {
        // Only applicable to Vectors
        self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w
    }
    pub fn cross(&self, other: &Self) -> Self {
        // Only applicable to vectors
        Tuple::new_vector(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

impl Add for Tuple {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            w: self.w + other.w,
        }
    }
}

impl Sub for Tuple {
    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            w: self.w - other.w,
        }
    }
}

impl Neg for Tuple {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: -self.w,
        }
    }
}

impl Mul<f64> for Tuple {
    // Multiply by scalar
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Self {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
            w: self.w * other,
        }
    }
}

impl Div<f64> for Tuple {
    type Output = Self;

    fn div(self, other: f64) -> Self::Output {
        Self {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
            w: self.w / other,
        }
    }
}

trait ApproxEq {
    fn approx_eq(self, other: Self) -> bool;
}

impl ApproxEq for Tuple {
    fn approx_eq(self, other: Self) -> bool {
        ((self.x - other.x).abs() < f64::EPSILON)
            & ((self.y - other.y).abs() < f64::EPSILON)
            & ((self.z - other.z).abs() < f64::EPSILON)
            & ((self.w - other.w).abs() < f64::EPSILON)
    }
}

#[allow(dead_code)]
pub fn compare_f64(x: f64, y: f64) -> bool {
    (x - y).abs() < f64::EPSILON
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuple_as_point() {
        let point = Tuple::new(4.3, -4.2, 3.1, 1.0);
        assert!(point.approx_eq(Tuple::new_point(4.3, -4.2, 3.1)));
        assert!(point.is_point());
        assert!(!point.is_vector());
    }
    #[test]
    fn tuple_as_vector() {
        let vector = Tuple::new(4.3, -4.2, 3.1, 0.0);
        assert!(!vector.is_point());
        assert!(vector.is_vector());
    }
    #[test]
    fn create_point() {
        let point = Tuple::new_point(3.0, 4.0, 5.0);
        assert!(point.is_point());
    }
    #[test]
    fn create_vector() {
        let vector = Tuple::new_vector(3.0, 4.0, 5.0);
        assert!(vector.is_vector());
    }
    #[test]
    fn add_two_tuples() {
        let tuple1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let tuple2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);
        assert_eq!(tuple1 + tuple2, Tuple::new(1.0, 1.0, 6.0, 1.0));
    }
    #[test]
    fn subtract_two_tuples() {
        let p1 = Tuple::new_point(3.0, 2.0, 1.0);
        let p2 = Tuple::new_point(5.0, 6.0, 7.0);
        assert_eq!(p1 - p2, Tuple::new_vector(-2.0, -4.0, -6.0));
    }
    #[test]
    fn subtract_vector_from_point() {
        let p = Tuple::new_point(3.0, 2.0, 1.0);
        let v = Tuple::new_vector(5.0, 6.0, 7.0);
        assert_eq!(p - v, Tuple::new_point(-2.0, -4.0, -6.0));
    }
    #[test]
    fn subtract_two_vectors() {
        let v1 = Tuple::new_vector(3.0, 2.0, 1.0);
        let v2 = Tuple::new_vector(5.0, 6.0, 7.0);
        assert_eq!(v1 - v2, Tuple::new_vector(-2.0, -4.0, -6.0));
    }
    #[test]
    fn negate_tuple() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(-a, Tuple::new(-1.0, 2.0, -3.0, 4.0));
    }
    #[test]
    fn multiply_by_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a * 3.5, Tuple::new(3.5, -7.0, 10.5, -14.0));
    }
    #[test]
    fn multiply_by_fraction() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a * 0.5, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }
    #[test]
    fn divide_by_scalar() {
        let a = Tuple::new(1.0, -2.0, 3.0, -4.0);
        assert_eq!(a / 2.0, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }
    #[test]
    fn compute_magnitude_vector() {
        let v = Tuple::new_vector(1.0, 0.0, 0.0);
        assert!(compare_f64(v.magnitude(), 1.0));
        let v = Tuple::new_vector(0.0, 1.0, 0.0);
        assert!(compare_f64(v.magnitude(), 1.0));
        let v = Tuple::new_vector(0.0, 0.0, 1.0);
        assert!(compare_f64(v.magnitude(), 1.0));
        let v = Tuple::new_vector(-1.0, -2.0, -3.0);
        assert!(compare_f64(v.magnitude(), 14.0_f64.sqrt()));
    }
    #[test]
    fn normalize_vector() {
        let v = Tuple::new_vector(4.0, 0.0, 0.0);
        assert_eq!(v.normalize(), Tuple::new_vector(1.0, 0.0, 0.0));
        let v = Tuple::new_vector(1.0, 2.0, 3.0);
        assert_eq!(
            v.normalize(),
            Tuple::new_vector(
                1.0 / 14.0_f64.sqrt(),
                2.0 / 14.0_f64.sqrt(),
                3.0 / 14.0_f64.sqrt()
            )
        );
        assert!(compare_f64(v.normalize().magnitude(), 1.0));
    }
    #[test]
    fn macro_point() {
        let v = point!(1.0, 2.0, 3.0);
        assert_eq!(v, Tuple::new_point(1.0, 2.0, 3.0));
    }
    #[test]
    fn dot() {
        let a = vector!(1.0, 2.0, 3.0);
        let b = vector!(2.0, 3.0, 4.0);
        assert!(compare_f64(a.dot(&b), 20.0));
    }
    #[test]
    fn cross() {
        let a = vector!(1.0, 2.0, 3.0);
        let b = vector!(2.0, 3.0, 4.0);
        assert_eq!(a.cross(&b), vector!(-1.0, 2.0, -1.0));
        assert_eq!(b.cross(&a), vector!(1.0, -2.0, 1.0));
    }
}
