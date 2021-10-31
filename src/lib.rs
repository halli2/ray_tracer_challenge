pub mod canvas;
pub mod color;
pub mod matrix;
pub mod point;
pub mod vector;

pub use crate::point::Point;
pub use crate::vector::Vector;

#[allow(dead_code)]
pub fn eq_f64(x: f64, y: f64) -> bool {
    // Matrix multiplication needed a bigger epsilon..
    if (x - y).abs() < f64::EPSILON * 100.0 {
        return true;
    }
    eprintln!(
        "
eq_f64: `(left == right)`
left: `{}`,
right: `{}`
",
        x, y
    );
    false
}
