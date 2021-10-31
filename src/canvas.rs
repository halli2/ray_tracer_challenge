use crate::{color::Color, point::Point, vector::Vector};

pub struct Canvas {
    width: usize,
    height: usize,
    pixels: Vec<Color>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Self {
        Canvas {
            width,
            height,
            pixels: vec![Color::new(0.0, 0.0, 0.0); height * width],
        }
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> Color {
        self.pixels[x + y * self.width]
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, color: Color) {
        self.pixels[x + y * self.width] = color;
    }
}

impl IntoIterator for Canvas {
    // Iterates over all pixels
    type Item = Color;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.pixels.into_iter()
    }
}

#[cfg(test)]
mod canvas_tests {
    use super::*;

    #[test]
    fn create_canvas() {
        let c = Canvas::new(10, 20);
        assert_eq!(c.width, 10);
        assert_eq!(c.height, 20);
        let black = Color::new(0.0, 0.0, 0.0);
        for pixel in c {
            assert_eq!(pixel, black);
        }
    }
    #[test]
    fn write_pixels_to_canvas() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1.0, 0.0, 0.0);
        c.write_pixel(2, 3, red);
        assert_eq!(c.pixel_at(2, 3), red);
    }
}
