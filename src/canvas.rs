use crate::{color::Color, point::Point, vector::Vector};

pub struct Canvas {
    pub width: usize,
    pub height: usize,
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
        if (x + y * self.width) >= self.pixels.len() {
            return;
        }
        self.pixels[x + y * self.width] = color;
    }
    pub fn to_ppm(self) -> String {
        let mut buffer = ["P3", &format!("{} {}", self.width, self.height), "255"].join("\n");
        buffer.push('\n');

        let mut line_len = 0;
        for (index, color) in self.pixels.iter().enumerate() {
            if index % self.width == 0 {
                buffer.pop();
                buffer.push('\n');
                line_len = 0;
            }
            for color in &color.to_rbg_string() {
                line_len += color.len();
                if line_len > 70 {
                    line_len = color.len();
                    buffer.pop();
                    buffer.push('\n');
                }
                buffer.push_str(color);
            }
        }
        buffer.pop(); // Remove last space
        buffer.push('\n');
        buffer
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
    #[test]
    fn construct_ppm_header() {
        let c = Canvas::new(5, 3);
        let ppm = c.to_ppm();
        let split_at_header = ppm.split_at(10); // Split at index 10
        assert_eq!(split_at_header.0, "P3\n5 3\n255");
    }
    #[test]
    fn construct_ppm_pix_data() {
        let mut c = Canvas::new(5, 3);
        let c1 = Color::new(1.5, 0.0, 0.0);
        let c2 = Color::new(0.0, 0.5, 0.0);
        let c3 = Color::new(0.0, 0.0, 1.0);
        c.write_pixel(0, 0, c1);
        c.write_pixel(2, 1, c2);
        c.write_pixel(4, 2, c3);
        let ppm = c.to_ppm();
        let split_at_header = ppm.split_at(10); // Split at index 10
        assert_eq!(split_at_header.0, "P3\n5 3\n255");
        assert_eq!(
            split_at_header.1,
            "
255 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 128 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 255
" // SHOULD END WITH NEWLINE!
        );
    }
    #[test]
    fn split_long_lines() {
        let mut c = Canvas::new(10, 2);
        let c1 = Color::new(1.0, 0.8, 0.6);
        for x in 0..10 {
            for y in 0..2 {
                c.write_pixel(x, y, c1);
            }
        }
        let ppm = c.to_ppm();
        let split = ppm.split('\n').collect::<Vec<_>>();

        assert_eq!(
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204",
            split[3]
        );
        assert_eq!(
            "153 255 204 153 255 204 153 255 204 153 255 204 153",
            split[4]
        );
        assert_eq!(
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204",
            split[5]
        );
        assert_eq!(
            "153 255 204 153 255 204 153 255 204 153 255 204 153",
            split[6]
        );
    }
}
