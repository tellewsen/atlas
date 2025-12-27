// Gravitational constant (m^3 kg^-1 s^-2)
pub const G: f64 = 6.67430e-11;
// Astronomical Unit in meters
pub const AU: f64 = 1.496e11;

#[derive(Debug, Clone)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&self) -> Vec3 {
        let mag = self.magnitude();
        Vec3 {
            x: self.x / mag,
            y: self.y / mag,
            z: self.z / mag,
        }
    }

    pub fn add(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    pub fn sub(&self, other: &Vec3) -> Vec3 {
        Vec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    pub fn scale(&self, scalar: f64) -> Vec3 {
        Vec3 {
            x: self.x * scalar,
            y: self.y * scalar,
            z: self.z * scalar,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Body {
    pub name: String,
    pub mass: f64,      // kg
    pub position: Vec3, // m
    pub velocity: Vec3, // m/s
}

impl Body {
    pub fn new(name: &str, mass: f64, position: Vec3, velocity: Vec3) -> Self {
        Body {
            name: name.to_string(),
            mass,
            position,
            velocity,
        }
    }
}

pub fn print_body_state(body: &Body) {
    let r = body.position.magnitude();
    let v = body.velocity.magnitude();
    println!(
        "  {}: distance = {:.3e} m ({:.4} AU), velocity = {:.3e} m/s",
        body.name,
        r,
        r / AU,
        v
    );
}

/// Calculate gravitational acceleration on a body due to all other bodies
pub fn calculate_acceleration(body: &Body, other_bodies: &[Body]) -> Vec3 {
    let mut acceleration = Vec3::new(0.0, 0.0, 0.0);

    for other in other_bodies {
        if body.name == other.name {
            continue;
        }

        let r_vec = other.position.sub(&body.position);
        let r = r_vec.magnitude();

        // F = G * m1 * m2 / r^2
        // a = G * m2 / r^2
        let a_magnitude = G * other.mass / (r * r);
        let a_vec = r_vec.normalize().scale(a_magnitude);

        acceleration = acceleration.add(&a_vec);
    }

    acceleration
}

/// Calculate the total kinetic energy of the system
pub fn total_kinetic_energy(bodies: &[Body]) -> f64 {
    bodies
        .iter()
        .map(|body| 0.5 * body.mass * body.velocity.magnitude().powi(2))
        .sum()
}

/// Calculate the total potential energy of the system
pub fn total_potential_energy(bodies: &[Body]) -> f64 {
    let mut potential_energy = 0.0;

    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let r_vec = bodies[j].position.sub(&bodies[i].position);
            let r = r_vec.magnitude();
            // U = -G * m1 * m2 / r
            potential_energy -= G * bodies[i].mass * bodies[j].mass / r;
        }
    }

    potential_energy
}

/// Calculate the total mechanical energy of the system
pub fn total_energy(bodies: &[Body]) -> f64 {
    total_kinetic_energy(bodies) + total_potential_energy(bodies)
}

/// Calculate the total angular momentum of the system (returns magnitude)
pub fn total_angular_momentum(bodies: &[Body]) -> f64 {
    let mut angular_momentum = Vec3::new(0.0, 0.0, 0.0);

    for body in bodies {
        // L = r × p = r × (m * v)
        let momentum = body.velocity.scale(body.mass);
        let l_x = body.position.y * momentum.z - body.position.z * momentum.y;
        let l_y = body.position.z * momentum.x - body.position.x * momentum.z;
        let l_z = body.position.x * momentum.y - body.position.y * momentum.x;

        angular_momentum = angular_momentum.add(&Vec3::new(l_x, l_y, l_z));
    }

    angular_momentum.magnitude()
}

/// Calculate the theoretical orbital period for a two-body circular orbit
/// T = 2π * sqrt(a³ / (G * (m1 + m2)))
/// where a is the semi-major axis (total separation / 2 for equal masses)
pub fn orbital_period(bodies: &[Body]) -> f64 {
    if bodies.len() != 2 {
        panic!("Orbital period calculation only valid for two-body systems");
    }

    let total_mass = bodies[0].mass + bodies[1].mass;
    let separation = bodies[0].position.sub(&bodies[1].position).magnitude();

    // For the two-body problem, the orbital period is:
    // T = 2π * sqrt(a³ / (G * M_total))
    // where a is the semi-major axis (separation for circular orbits)
    2.0 * std::f64::consts::PI * (separation.powi(3) / (G * total_mass)).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrators::euler_step;
    use crate::scenarios::two_body_circular;

    #[test]
    fn test_two_body_energy_conservation() {
        let mut bodies = two_body_circular().bodies;
        let initial_energy = total_energy(&bodies);

        let dt = 100.0; // 100 seconds
        let steps = 10000; // Simulate for a while

        for _ in 0..steps {
            euler_step(&mut bodies, dt);
        }

        let final_energy = total_energy(&bodies);
        let energy_drift = (final_energy - initial_energy).abs() / initial_energy.abs();

        // Euler method will have some energy drift, but it should be bounded
        // For a perfect integrator, this would be < 1e-10
        // For Euler, we expect some drift but not catastrophic
        assert!(
            energy_drift < 0.5,
            "Energy drift too large: {:.2}% (initial: {:.3e}, final: {:.3e})",
            energy_drift * 100.0,
            initial_energy,
            final_energy
        );
    }

    #[test]
    fn test_two_body_angular_momentum_conservation() {
        let mut bodies = two_body_circular().bodies;
        let initial_angular_momentum = total_angular_momentum(&bodies);

        let dt = 100.0;
        let steps = 10000;

        for _ in 0..steps {
            euler_step(&mut bodies, dt);
        }

        let final_angular_momentum = total_angular_momentum(&bodies);
        let angular_momentum_drift = (final_angular_momentum - initial_angular_momentum).abs()
            / initial_angular_momentum.abs();

        // Angular momentum should be better conserved than energy with Euler
        assert!(
            angular_momentum_drift < 0.1,
            "Angular momentum drift too large: {:.2}% (initial: {:.3e}, final: {:.3e})",
            angular_momentum_drift * 100.0,
            initial_angular_momentum,
            final_angular_momentum
        );
    }

    #[test]
    fn test_two_body_center_of_mass() {
        let mut bodies = two_body_circular().bodies;

        // For equal masses starting at equal distances from origin with opposite velocities,
        // the center of mass should remain at the origin
        let initial_com = Vec3::new(
            (bodies[0].position.x * bodies[0].mass + bodies[1].position.x * bodies[1].mass)
                / (bodies[0].mass + bodies[1].mass),
            (bodies[0].position.y * bodies[0].mass + bodies[1].position.y * bodies[1].mass)
                / (bodies[0].mass + bodies[1].mass),
            (bodies[0].position.z * bodies[0].mass + bodies[1].position.z * bodies[1].mass)
                / (bodies[0].mass + bodies[1].mass),
        );

        assert!(
            initial_com.magnitude() < 1e-6,
            "Initial center of mass should be at origin"
        );

        let dt = 100.0;
        let steps = 5000;

        for _ in 0..steps {
            euler_step(&mut bodies, dt);
        }

        let final_com = Vec3::new(
            (bodies[0].position.x * bodies[0].mass + bodies[1].position.x * bodies[1].mass)
                / (bodies[0].mass + bodies[1].mass),
            (bodies[0].position.y * bodies[0].mass + bodies[1].position.y * bodies[1].mass)
                / (bodies[0].mass + bodies[1].mass),
            (bodies[0].position.z * bodies[0].mass + bodies[1].position.z * bodies[1].mass)
                / (bodies[0].mass + bodies[1].mass),
        );

        // Center of mass should remain near origin (some drift expected with Euler)
        assert!(
            final_com.magnitude() < 1e8,
            "Center of mass drifted too far: {:.3e} m",
            final_com.magnitude()
        );
    }

    #[test]
    fn test_two_body_orbital_stability() {
        let scenario = two_body_circular();
        let mut bodies = scenario.bodies;
        let initial_separation = bodies[0].position.sub(&bodies[1].position).magnitude();
        let expected_period = orbital_period(&bodies);

        println!("Initial separation: {:.3e} m", initial_separation);
        println!(
            "Expected orbital period: {:.3e} seconds ({:.2} days)",
            expected_period,
            expected_period / 86400.0
        );

        // Use small time step for stability
        let dt = 10.0; // 10 seconds
        let mut time = 0.0;

        // Simulate for just 1% of one orbital period to avoid excessive drift with Euler
        let simulation_time = expected_period * 0.01;

        while time < simulation_time {
            euler_step(&mut bodies, dt);
            time += dt;
        }

        let final_separation = bodies[0].position.sub(&bodies[1].position).magnitude();
        let separation_change = (final_separation - initial_separation).abs() / initial_separation;

        println!("Final separation: {:.3e} m", final_separation);
        println!("Separation change: {:.2}%", separation_change * 100.0);
        println!(
            "Simulation time: {:.2e} s ({:.2} days, {:.2}% of orbit)",
            simulation_time,
            simulation_time / 86400.0,
            100.0 * simulation_time / expected_period
        );

        // For a circular orbit, separation should remain roughly constant over short periods
        // Allow 5% variation even with Euler over 1% of an orbit
        assert!(
            separation_change < 0.05,
            "Orbital separation changed too much: {:.2}% (initial: {:.3e} m, final: {:.3e} m)",
            separation_change * 100.0,
            initial_separation,
            final_separation
        );
    }
}
