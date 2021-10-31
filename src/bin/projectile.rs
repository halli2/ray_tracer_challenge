use ray_tracer_challenge::{point::Point, vector::Vector};
use std::{thread::sleep, time::Duration};

#[derive(Debug)]
struct Projectile {
    position: Point,
    velocity: Vector,
}

struct Environment {
    gravity: Vector,
    wind: Vector,
}

fn tick(env: &Environment, proj: &mut Projectile) {
    let position = proj.position + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    *proj = Projectile { position, velocity };
}

fn main() {
    let mut proj = Projectile {
        position: Point::new(0.0, 1.0, 0.0),
        velocity: Vector::new(1.0, 1.0, 0.0).normalize(),
    };
    let env = Environment {
        gravity: Vector::new(0.0, -0.01, 0.0),
        wind: Vector::new(-0.01, 0.0, 0.0),
    };

    while proj.position.y >= 0.0 {
        tick(&env, &mut proj);
        println!("{:?}", proj);
        sleep(Duration::from_millis(200));
    }
}
