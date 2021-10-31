pub mod canvas;
pub mod color;
pub mod point;
pub mod vector;

pub use crate::point::Point;
pub use crate::vector::Vector;

#[allow(dead_code)]
pub fn eq_f64(x: f64, y: f64) -> bool {
    (x - y).abs() < f64::EPSILON
}
