#[cfg(test)]
mod tests {
    use crate::{Point, Vector};

    use super::*;

    #[test]
    fn create_and_query_ray() {
        let origin = Point::new(1.0, 2.0, 3.0);
        let direction = Vector::new(4.0, 5.0, 6.0);
    }
}
