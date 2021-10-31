use std::{
    fmt,
    ops::{Index, IndexMut, Mul},
};

use crate::{eq_f64, Point, Vector};

#[derive(Debug, Copy, Clone)]
pub struct Matrix {
    data: [[f64; 4]; 4],
}

pub const IDENTITY: Matrix = Matrix {
    data: [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
};

#[allow(dead_code)]
impl Matrix {
    pub fn new(data: [[f64; 4]; 4]) -> Self {
        Self { data }
    }

    fn transpose(&self) -> Self {
        let mut result = [[0.0; 4]; 4];

        for i in 0..4 {
            (0..4).for_each(|j| {
                result[j][i] = self[i][j];
            });
        }
        Self { data: result }
    }

    fn determinant(self, s: usize) -> f64 {
        // Where s is size of matrix
        if s == 2 {
            self[0][0] * self[1][1] - self[0][1] * self[1][0]
        } else {
            let mut det = 0.;
            for i in 0..s {
                det += self[0][i] * self.cofactor(0, i, s - 1);
            }
            det
        }
    }
    fn submatrix(self, row: usize, col: usize) -> Self {
        let mut m = [[0.; 4]; 4];

        for (r, ri) in (0..4).into_iter().filter(|&x| x != row).enumerate() {
            for (c, ci) in (0..4).into_iter().filter(|&x| x != col).enumerate() {
                m[r][c] = self[ri][ci]
            }
        }
        Self::new(m)
    }
    fn minor(&self, row: usize, col: usize, s: usize) -> f64 {
        // Takes current size of matrix, removes row and col, then gives determinant
        self.submatrix(row, col).determinant(s)
    }
    fn cofactor(self, row: usize, col: usize, s: usize) -> f64 {
        if (row + col) & 1 != 0 {
            // Finds odd numbers on integers
            return -self.minor(row, col, s);
        }
        self.minor(row, col, s)
    }
    fn invertible(self, s: usize) -> bool {
        if self.determinant(s) == 0.0 {
            return false;
        }
        true
    }
    fn inverse(self, s: usize) -> Option<Self> {
        if !self.invertible(s) {
            return None;
        }
        let mut inverse_matrix = Matrix::new([[0.; 4]; 4]);
        let determinant = self.determinant(s);
        for row in 0..s {
            for col in 0..s {
                inverse_matrix[col][row] = self.cofactor(row, col, s - 1) / determinant;
            }
        }
        Some(inverse_matrix)
    }

    // TRANSFORMATION SPECIFIC METHODS
    pub fn translation(x: f64, y: f64, z: f64) -> Self {
        Matrix::new([
            [1.0, 0.0, 0.0, x],
            [0.0, 1.0, 0.0, y],
            [0.0, 0.0, 1.0, z],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
    pub fn scaling(x: f64, y: f64, z: f64) -> Self {
        Matrix::new([
            [x, 0.0, 0.0, 0.0],
            [0.0, y, 0.0, 0.0],
            [0.0, 0.0, z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
    pub fn rotate_x(r: f64) -> Self {
        Self::new([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, r.cos(), -r.sin(), 0.0],
            [0.0, r.sin(), r.cos(), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
    pub fn rotate_y(r: f64) -> Self {
        Self::new([
            [r.cos(), 0.0, r.sin(), 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [-r.sin(), 0.0, r.cos(), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
    pub fn rotate_z(r: f64) -> Self {
        Self::new([
            [r.cos(), -r.sin(), 0.0, 0.0],
            [r.sin(), r.cos(), 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
    fn shearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) -> Self {
        Self::new([
            [1.0, xy, xz, 0.0],
            [yx, 1.0, yz, 0.0],
            [zx, zy, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
}

impl Index<usize> for Matrix {
    type Output = [f64; 4];

    fn index(&self, i: usize) -> &Self::Output {
        &self.data[i]
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, i: usize) -> &mut Self::Output {
        &mut self.data[i]
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = [[0.0; 4]; 4];
        for row in 0..4 {
            for col in 0..4 {
                result[row][col] = self[row][0] * rhs[0][col]
                    + self[row][1] * rhs[1][col]
                    + self[row][2] * rhs[2][col]
                    + self[row][3] * rhs[3][col]
            }
        }
        Matrix::new(result)
    }
}

impl Mul<Point> for Matrix {
    type Output = Point;

    fn mul(self, rhs: Point) -> Point {
        let x =
            (self[0][0] * rhs.x) + (self[0][1] * rhs.y) + (self[0][2] * rhs.z) + (self[0][3] * 1.0);
        let y =
            (self[1][0] * rhs.x) + (self[1][1] * rhs.y) + (self[1][2] * rhs.z) + (self[1][3] * 1.0);
        let z =
            (self[2][0] * rhs.x) + (self[2][1] * rhs.y) + (self[2][2] * rhs.z) + (self[2][3] * 1.0);
        Point::new(x, y, z)
    }
}
impl Mul<Vector> for Matrix {
    type Output = Vector;

    fn mul(self, rhs: Vector) -> Vector {
        let x =
            (self[0][0] * rhs.x) + (self[0][1] * rhs.y) + (self[0][2] * rhs.z) + (self[0][3] * 0.0);
        let y =
            (self[1][0] * rhs.x) + (self[1][1] * rhs.y) + (self[1][2] * rhs.z) + (self[1][3] * 0.0);
        let z =
            (self[2][0] * rhs.x) + (self[2][1] * rhs.y) + (self[2][2] * rhs.z) + (self[2][3] * 0.0);
        Vector::new(x, y, z)
    }
}

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{0:>10} {1:>10} {2:>10} {3:>10}
            {4:>10} {5:>10} {6:>10} {7:>10}
            {8:>10} {9:>10} {10:>10} {11:>10}
            {12:>10} {13:>10} {14:>10} {15:>10}",
            format!("{0:.5}", self[0][0]),
            format!("{0:.5}", self[0][1]),
            format!("{0:.5}", self[0][2]),
            format!("{0:.5}", self[0][3]),
            format!("{0:.5}", self[1][0]),
            format!("{0:.5}", self[1][1]),
            format!("{0:.5}", self[1][2]),
            format!("{0:.5}", self[1][3]),
            format!("{0:.5}", self[2][0]),
            format!("{0:.5}", self[2][1]),
            format!("{0:.5}", self[2][2]),
            format!("{0:.5}", self[2][3]),
            format!("{0:.5}", self[3][0]),
            format!("{0:.5}", self[3][1]),
            format!("{0:.5}", self[3][2]),
            format!("{0:.5}", self[3][3]),
        )
    }
}

impl PartialEq for Matrix {
    fn eq(&self, other: &Self) -> bool {
        for r in 0..4 {
            for c in 0..4 {
                if !eq_f64(self[r][c], other[r][c]) {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod matrix_tests {
    use super::*;
    use crate::eq_f64;

    #[test]
    fn construct() {
        let matrix = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);
        assert!(eq_f64(matrix[0][0], 1.0));
        assert!(eq_f64(matrix[0][3], 4.0));
        assert!(eq_f64(matrix[1][0], 5.5));
        assert!(eq_f64(matrix[1][2], 7.5));
        assert!(eq_f64(matrix[2][2], 11.0));
        assert!(eq_f64(matrix[3][0], 13.5));
        assert!(eq_f64(matrix[3][2], 15.5));
    }

    #[test]
    fn matrix_2x2() {
        let m = Matrix::new([
            [-3.0, 5.0, 0.0, 0.0],
            [1.0, -2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);

        assert!(eq_f64(m[0][0], -3.0));
        assert!(eq_f64(m[1][0], 1.0));
        assert!(eq_f64(m[1][1], -2.0));
    }

    #[test]
    fn compare_matrix() {
        let m1 = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let m2 = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let m3 = Matrix::new([
            [2.0, 2.0, 3.0, 4.0],
            [4.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        println!("{}", m1);
        assert_eq!(m1, m2);
        assert_ne!(m1, m3);
    }

    #[test]
    fn matr_mul() {
        let m1 = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let m2 = Matrix::new([
            [-2.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, -1.0],
            [4.0, 3.0, 6.0, 5.0],
            [1.0, 2.0, 7.0, 8.0],
        ]);
        assert_eq!(
            m1 * m2,
            Matrix::new([
                [20.0, 22.0, 50.0, 48.0],
                [44.0, 54.0, 114.0, 108.0],
                [40.0, 58.0, 110.0, 102.0],
                [16.0, 26.0, 46.0, 42.0],
            ])
        );
    }
    #[test]
    fn tuple_mul() {
        let m = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let v = Vector::new(1.0, 2.0, 3.0);
        let p = Point::new(1.0, 2.0, 3.0);

        assert_eq!(m * v, Vector::new(18.0, 24.0, 33.0));
        assert_eq!(m * p, Point::new(18.0, 24.0, 33.0));
    }

    #[test]
    fn m_times_identity() {
        let m = Matrix::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        assert_eq!(m * IDENTITY, m);
    }

    #[test]
    fn transpose_matrix() {
        let a = Matrix::new([
            [0.0, 9.0, 3.0, 0.0],
            [9.0, 8.0, 0.0, 8.0],
            [1.0, 8.0, 5.0, 3.0],
            [0.0, 0.0, 5.0, 8.0],
        ]);
        let a_trans = Matrix::new([
            [0.0, 9.0, 1.0, 0.0],
            [9.0, 8.0, 8.0, 0.0],
            [3.0, 0.0, 5.0, 5.0],
            [0.0, 8.0, 3.0, 8.0],
        ]);
        assert_eq!(a.transpose(), a_trans);
    }
    #[test]
    fn transpose_of_identity() {
        assert_eq!(IDENTITY.transpose(), IDENTITY);
    }
    #[test]
    fn determinant() {
        let a = Matrix::new([
            [1.0, 5.0, 0.0, 0.0],
            [-3.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        assert!(eq_f64(a.determinant(2), 17.0));
    }
    #[test]
    fn submatrix() {
        let a = Matrix::new([
            [1.0, 5.0, 0.0, 0.0],
            [-3.0, 2.0, 7.0, 0.0],
            [0.0, 6.0, -3.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        assert_eq!(
            a.submatrix(0, 2),
            Matrix::new([
                [-3.0, 2.0, 0.0, 0.0],
                [0.0, 6.0, -0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ])
        );
        let a = Matrix::new([
            [-6.0, 1.0, 1.0, 6.0],
            [-8.0, 5.0, 8.0, 6.0],
            [-1.0, 0.0, 8.0, 2.0],
            [-7.0, 1.0, -1.0, 1.0],
        ]);
        assert_eq!(
            a.submatrix(2, 1),
            Matrix::new([
                [-6.0, 1.0, 6.0, 0.0],
                [-8.0, 8.0, 6.0, 0.0],
                [-7.0, -1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ])
        );
    }
    #[test]
    fn calculate_minor() {
        let a = Matrix::new([
            [3.0, 5.0, 0.0, 0.0],
            [2.0, -1.0, -7.0, 0.0],
            [6.0, -1.0, 5.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        let b = a.submatrix(1, 0);
        assert!(eq_f64(b.determinant(2), 25.0));
        assert!(eq_f64(a.minor(1, 0, 2), 25.0));
    }
    #[test]
    fn calculate_cofactor() {
        let a = Matrix::new([
            [3.0, 5.0, 0.0, 0.0],
            [2.0, -1.0, -7.0, 0.0],
            [6.0, -1.0, 5.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        assert!(eq_f64(a.minor(0, 0, 2), -12.0));
        assert!(eq_f64(a.cofactor(0, 0, 2), -12.0));
        assert!(eq_f64(a.minor(1, 0, 2), 25.0));
        assert!(eq_f64(a.cofactor(1, 0, 2), -25.0));
    }
    #[test]
    fn calculate_determinant_3_3_4_4() {
        let a = Matrix::new([
            [1.0, 2.0, 6.0, 0.0],
            [-5.0, 8.0, -4.0, 0.0],
            [2.0, 6.0, 4.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        assert!(eq_f64(a.cofactor(0, 0, 2), 56.0));
        assert!(eq_f64(a.cofactor(0, 1, 2), 12.0));
        assert!(eq_f64(a.cofactor(0, 2, 2), -46.0));
        assert!(eq_f64(a.determinant(3), -196.0));
        let b = Matrix::new([
            [-2.0, -8.0, 3.0, 5.0],
            [-3.0, 1.0, 7.0, 3.0],
            [1.0, 2.0, -9.0, 6.0],
            [-6.0, 7.0, 7.0, -9.0],
        ]);
        assert!(eq_f64(b.cofactor(0, 0, 3), 690.0));
        assert!(eq_f64(b.cofactor(0, 1, 3), 447.0));
        assert!(eq_f64(b.cofactor(0, 2, 3), 210.0));
        assert!(eq_f64(b.cofactor(0, 3, 3), 51.0));
        assert!(eq_f64(b.determinant(4), -4071.0));
    }
    #[test]
    fn is_invertible() {
        let a = Matrix::new([
            [6.0, 4.0, 4.0, 4.0],
            [5.0, 5.0, 7.0, 6.0],
            [4.0, -9.0, 3.0, -7.0],
            [9.0, 1.0, 7.0, -6.0],
        ]);
        let b = Matrix::new([
            [-4.0, 2.0, -2.0, -3.0],
            [9.0, 6.0, 2.0, 6.0],
            [0.0, -5.0, 1.0, -5.0],
            [0.0, 0.0, 0.0, 0.0],
        ]);
        assert!(eq_f64(a.determinant(4), -2120.0));
        assert!(a.invertible(4));
        assert!(!b.invertible(4));
    }
    #[test]
    fn inverse_matrix() {
        let a = Matrix::new([
            [-5.0, 2.0, 6.0, -8.0],
            [1.0, -5.0, 1.0, 8.0],
            [7.0, 7.0, -6.0, -7.0],
            [1.0, -3.0, 7.0, 4.0],
        ]);
        let b = a.inverse(4).unwrap();
        println!("{}", b);
        assert!(eq_f64(a.determinant(4), 532.0));
        assert!(eq_f64(a.cofactor(2, 3, 3), -160.0));
        assert!(eq_f64(b[3][2], -160.0 / 532.0));
        assert!(eq_f64(a.cofactor(3, 2, 3), 105.0));
        assert!(eq_f64(b[2][3], 105.0 / 532.0));
        assert_eq!(
            b,
            Matrix::new([
                [
                    0.21804511278195488,
                    0.45112781954887216,
                    0.24060150375939848,
                    -0.045112781954887216
                ],
                [
                    -0.8082706766917294,
                    -1.4567669172932332,
                    -0.44360902255639095,
                    0.5206766917293233
                ],
                [
                    -0.07894736842105263,
                    -0.2236842105263158,
                    -0.05263157894736842,
                    0.19736842105263158
                ],
                [
                    -0.5225563909774437,
                    -0.8139097744360902,
                    -0.3007518796992481,
                    0.30639097744360905
                ]
            ])
        );
    }
    #[test]
    fn multiply_by_inverse() {
        let a = Matrix::new([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);
        let b = Matrix::new([
            [8.0, 2.0, 2.0, 2.0],
            [3.0, -1.0, 7.0, 0.0],
            [7.0, 0.0, 5.0, 4.0],
            [6.0, -2.0, 0.0, 5.0],
        ]);
        let c = a * b;
        assert_eq!(c * b.inverse(4).unwrap(), a);
    }
}
#[cfg(test)]
mod transformation_matrix_tests {
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    use super::*;
    #[test]
    fn multiple_by_translation_matrix() {
        let transform = Matrix::translation(5.0, -3.0, 2.0);
        let p = Point::new(-3.0, 4.0, 5.0);
        assert_eq!(transform * p, Point::new(2.0, 1.0, 7.0));
    }
    #[test]
    fn multiple_by_inverse_translation_matrix() {
        let transform = Matrix::translation(5.0, -3.0, 2.0);
        let p = Point::new(-3.0, 4.0, 5.0);
        assert_eq!(
            transform.inverse(4).unwrap() * p,
            Point::new(-8.0, 7.0, 3.0)
        );
    }
    #[test]
    fn translation_does_not_affect_vectors() {
        let transform = Matrix::translation(5.0, -3.0, 2.0);
        let v = Vector::new(-3.0, 4.0, 5.0);
        assert_eq!(transform * v, v);
    }
    #[test]
    fn scaling_matrix_applied_point() {
        let transform = Matrix::scaling(2.0, 3.0, 4.0);
        let p = Point::new(-4.0, 6.0, 8.0);
        assert_eq!(transform * p, Point::new(-8.0, 18.0, 32.0));
    }
    #[test]
    fn scaling_matrix_applied_vector() {
        let transform = Matrix::scaling(2.0, 3.0, 4.0);
        let v = Vector::new(-4.0, 6.0, 8.0);
        assert_eq!(transform * v, Vector::new(-8.0, 18.0, 32.0));
    }
    #[test]
    fn multiplying_inverse_scaling() {
        let transform = Matrix::scaling(2.0, 3.0, 4.0);
        let v = Vector::new(-4.0, 6.0, 8.0);
        assert_eq!(
            transform.inverse(4).unwrap() * v,
            Vector::new(-2.0, 2.0, 2.0)
        );
    }
    #[test]
    fn reflection_is_scaling_negative() {
        let transform = Matrix::scaling(-1.0, 1.0, 1.0);
        let p = Point::new(2.0, 3.0, 4.0);
        assert_eq!(transform * p, Point::new(-2.0, 3.0, 4.0));
    }
    #[test]
    fn rotate_around_x() {
        let p = Point::new(0.0, 1.0, 0.0);
        let half_quarter = Matrix::rotate_x(FRAC_PI_4);
        let full_quarter = Matrix::rotate_x(FRAC_PI_2);
        assert_eq!(
            half_quarter * p,
            Point::new(0.0, 2_f64.sqrt() / 2.0, 2_f64.sqrt() / 2.0)
        );
        assert_eq!(full_quarter * p, Point::new(0.0, 0.0, 1.0));
    }
    #[test]
    fn rotate_around_x_inverse() {
        let p = Point::new(0.0, 1.0, 0.0);
        let half_quarter = Matrix::rotate_x(FRAC_PI_4);
        assert_eq!(
            half_quarter.inverse(4).unwrap() * p,
            Point::new(0.0, 2_f64.sqrt() / 2.0, -(2_f64.sqrt() / 2.0))
        );
    }
    #[test]
    fn rotate_around_y() {
        let p = Point::new(0.0, 0.0, 1.0);
        let half_quarter = Matrix::rotate_y(FRAC_PI_4);
        let full_quarter = Matrix::rotate_y(FRAC_PI_2);
        assert_eq!(
            half_quarter * p,
            Point::new(2_f64.sqrt() / 2.0, 0.0, 2_f64.sqrt() / 2.0)
        );
        assert_eq!(full_quarter * p, Point::new(1.0, 0.0, 0.0));
    }
    #[test]
    fn rotate_around_z() {
        let p = Point::new(0.0, 1.0, 0.0);
        let half_quarter = Matrix::rotate_z(FRAC_PI_4);
        let full_quarter = Matrix::rotate_z(FRAC_PI_2);
        assert_eq!(
            half_quarter * p,
            Point::new(-(2_f64.sqrt() / 2.0), 2_f64.sqrt() / 2.0, 0.0)
        );
        assert_eq!(full_quarter * p, Point::new(-1.0, 0.0, 0.0));
    }
    #[test]
    fn shearing() {
        let p = Point::new(2.0, 3.0, 4.0);
        let t1 = Matrix::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let t2 = Matrix::shearing(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        let t3 = Matrix::shearing(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        let t4 = Matrix::shearing(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        let t5 = Matrix::shearing(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        let t6 = Matrix::shearing(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        assert_eq!(t1 * p, Point::new(5.0, 3.0, 4.0));
        assert_eq!(t2 * p, Point::new(6.0, 3.0, 4.0));
        assert_eq!(t3 * p, Point::new(2.0, 5.0, 4.0));
        assert_eq!(t4 * p, Point::new(2.0, 7.0, 4.0));
        assert_eq!(t5 * p, Point::new(2.0, 3.0, 6.0));
        assert_eq!(t6 * p, Point::new(2.0, 3.0, 7.0));
    }
    #[test]
    fn transformations_in_sequence() {
        let p = Point::new(1.0, 0.0, 1.0);
        let a = Matrix::rotate_x(FRAC_PI_2);
        let b = Matrix::scaling(5.0, 5.0, 5.0);
        let c = Matrix::translation(10.0, 5.0, 7.0);
        let p2 = a * p;
        assert_eq!(p2, Point::new(1.0, -1.0, 0.0));
        let p3 = b * p2;
        assert_eq!(p3, Point::new(5.0, -5.0, 0.0));
        let p4 = c * p3;
        assert_eq!(p4, Point::new(15.0, 0.0, 7.0));

        assert_eq!(c * b * a * p, Point::new(15.0, 0.0, 7.0));
    }
}
