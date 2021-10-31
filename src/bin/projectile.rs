use ray_tracer_challenge::{canvas::Canvas, color::Color, point::Point, vector::Vector};
use std::{fs::File, io::Write, thread::sleep, time::Duration};

#[derive(Debug)]
struct Projectile {
    position: Point,
    velocity: Vector,
}

struct Environment {
    gravity: Vector,
    wind: Vector,
}

fn tick(env: &Environment, proj: &mut Projectile, can: &mut Canvas, color: Color) {
    let position = proj.position + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    *proj = Projectile { position, velocity };
    let x = position.x as usize;
    let y = can.height - position.y as usize;
    // let y = (can.height - position.y as usize).clamp(0, can.height - 1);
    can.write_pixel(x, y, color)
}

fn main() -> std::io::Result<()> {
    let start = Point::new(0.0, 1.0, 0.0);
    let vel = Vector::new(1.0, 1.8, 0.0).normalize() * 11.25;
    let grav = Vector::new(0.0, -0.1, 0.0);
    let wind_vec = Vector::new(-0.01, 0.0, 0.0);
    let mut proj = Projectile {
        position: start,
        velocity: vel,
    };
    let env = Environment {
        gravity: grav,
        wind: wind_vec,
    };

    let mut canvas = Canvas::new(900, 550);
    let red = Color::new(1.0, 0.0, 0.0);

    while proj.position.y >= 0.0 {
        tick(&env, &mut proj, &mut canvas, red);
    }
    let ppm = canvas.to_ppm();
    let mut file = File::create("projectile.ppm")?;
    file.write_all(ppm.as_bytes())?;
    Ok(())
}
