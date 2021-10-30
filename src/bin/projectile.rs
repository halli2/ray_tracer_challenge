#[macro_use(vector, point)]
extern crate ray_tracer_challenge;

use ray_tracer_challenge::tuples::*;
use std::{thread::sleep, time::Duration};

#[derive(Debug)]
struct Projectile {
    position: Tuple,
    velocity: Tuple,
}

struct Environment {
    gravity: Tuple,
    wind: Tuple,
}

fn tick(env: &Environment, proj: &mut Projectile) {
    let position = proj.position + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    *proj = Projectile { position, velocity };
}

fn main() {
    let mut proj = Projectile {
        position: point!(0.0, 1.0, 0.0),
        velocity: vector!(1.0, 1.0, 0.0).normalize(),
    };
    let env = Environment {
        gravity: vector!(0.0, -0.01, 0.0),
        wind: vector!(-0.01, 0.0, 0.0),
    };

    while proj.position.y >= 0.0 {
        tick(&env, &mut proj);
        println!("{:?}", proj);
        sleep(Duration::from_millis(200));
    }
}
