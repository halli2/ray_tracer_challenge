use ray_tracer_challenge::{color::Color, Canvas, Matrix, Point};
use std::{f64::consts::PI, fs::File, io::Write, path::Path};

fn main() {
    let mut canvas = Canvas::new(200, 200);
    let twelve_o_clock = Point::new(0.0, 0.0, 1.0);
    let red = Color::new(1.0, 0.0, 0.0);
    for i in 0..1000 {
        let r = Matrix::rotate_y(i as f64 * PI / 500.0);
        let next_hour = r * twelve_o_clock;
        let radius = (3.0 / 8.0) * canvas.height as f64;
        let x = (canvas.width as f64 / 2.0 + (next_hour.x * radius)) as usize;
        let y = (canvas.height as f64 / 2.0 + (next_hour.z * radius)) as usize;
        canvas.write_pixel(x, y, red);
    }

    let ppm = canvas.to_ppm();
    let path = Path::new("clock.ppm");
    let mut file = match File::create(path) {
        Err(err) => panic!("couldn't create {}: {}", path.display(), err),
        Ok(file) => file,
    };
    match file.write_all(ppm.as_bytes()) {
        Err(err) => panic!("couldn't write to file {}: {}", path.display(), err),
        Ok(_) => println!("Successfully wrote to {}", path.display()),
    };
}
