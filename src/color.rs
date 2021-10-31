use std::{
    fmt,
    ops::{Add, Div, Mul, Sub},
};

#[derive(Debug, Copy, Clone)]
pub struct Color {
    pub red: f64,
    pub green: f64,
    pub blue: f64,
}

impl Color {
    pub fn new(red: f64, green: f64, blue: f64) -> Self {
        Color { red, green, blue }
    }

    pub fn new_from_u8(red: u8, green: u8, blue: u8) -> Self {
        Color {
            red: red as f64 / 255.0,
            green: green as f64 / 255.0,
            blue: blue as f64 / 255.0,
        }
    }

    pub fn to_rbg_string(&self) -> [String; 3] {
        [
            format!("{} ", (self.red.clamp(0.0, 1.0) * 256.0) as u8),
            format!("{} ", (self.green.clamp(0.0, 1.0) * 256.0) as u8),
            format!("{} ", (self.blue.clamp(0.0, 1.0) * 256.0) as u8),
        ]
    }
}

impl Add for Color {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            red: self.red + rhs.red,
            green: self.green + rhs.green,
            blue: self.blue + rhs.blue,
        }
    }
}

impl Mul for Color {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Self {
            red: self.red * other.red,
            green: self.green * other.green,
            blue: self.blue * other.blue,
        }
    }
}

impl Mul<f64> for Color {
    // Multiply by scalar
    type Output = Self;

    fn mul(self, other: f64) -> Self::Output {
        Self {
            red: self.red * other,
            green: self.green * other,
            blue: self.blue * other,
        }
    }
}

impl Div<f64> for Color {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            red: self.red / rhs,
            green: self.green / rhs,
            blue: self.blue / rhs,
        }
    }
}

impl Sub for Color {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            red: self.red - rhs.red,
            green: self.green - rhs.green,
            blue: self.blue - rhs.blue,
        }
    }
}

impl fmt::Display for Color {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Red: {0:>10} Green: {1:>10} Blue: {2:>10}",
            format!("{0:.5}", self.red),
            format!("{0:.5}", self.green),
            format!("{0:.5}", self.blue),
        )
    }
}

impl PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        ((self.red - other.red).abs() < f64::EPSILON)
            & ((self.green - other.green).abs() < f64::EPSILON)
            & ((self.blue - other.blue).abs() < f64::EPSILON)
    }
}

#[cfg(test)]
mod color_tests {
    use super::*;
    use crate::eq_f64;

    #[test]
    fn create_color() {
        let c = Color::new(0.5, 0.4, 1.7);
        assert!(eq_f64(c.red, 0.5));
        assert!(eq_f64(c.green, 0.4));
        assert!(eq_f64(c.blue, 1.7));
    }
    #[test]
    fn color_math() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 + c2, Color::new(1.6, 0.7, 1.0));
        let c = Color::new(0.2, 0.3, 0.4);
        assert_eq!(c1 - c2, Color::new(0.2, 0.5, 0.5));
        assert_eq!(c * 2.0, Color::new(0.4, 0.6, 0.8));
        let c1 = Color::new(1.0, 0.2, 0.4);
        let c2 = Color::new(0.9, 1.0, 0.1);
        assert_eq!(c1 * c2, Color::new(0.9, 0.2, 0.04));
    }
}
