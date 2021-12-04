pub mod canvas;
pub mod color;
pub mod matrix;
pub mod point;
pub mod rays;
pub mod vector;

pub use crate::canvas::Canvas;
pub use crate::matrix::Matrix;
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
pub fn float_eq(left: &f64, right: &f64) -> bool {
    (left - right).abs() < f64::EPSILON * 100.0
}

#[macro_export]
macro_rules! assert_eq_f64 {
    ($left:expr, $right:expr $(,)?) => {{
        match (&$left, &$right) {
            (left_val, right_val) => {
                if !$crate::float_eq(&$left, &$right) {
                    panic!(
                        "assertion failed: `(left == right)`
left: `{}`,
rigth: `{}`",
                        &*left_val, &*right_val
                    );
                }
            }
        }
    }};
}
