use crate::physics::{Body, Vec3};

pub struct Scenario {
    pub name: String,
    pub bodies: Vec<Body>,
    pub time_step: f64,        // seconds
    pub duration_days: f64,
    pub print_interval_days: f64,
}

impl Scenario {
    pub fn total_steps(&self) -> usize {
        ((self.duration_days * 24.0 * 3600.0) / self.time_step) as usize
    }

    pub fn print_interval_steps(&self) -> usize {
        ((self.print_interval_days * 24.0 * 3600.0) / self.time_step) as usize
    }
}

/// Create a solar system scenario with Sun, Earth, and Mars
pub fn solar_system() -> Scenario {
    let bodies =
    vec![
        Body::new(
            "Sun",
            1.989e30,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
        ),
        Body::new(
            "Earth",
            5.972e24,
            Vec3::new(1.496e11, 0.0, 0.0),
            Vec3::new(0.0, 29780.0, 0.0),
        ),
        Body::new(
            "Mars",
            6.4171e23,
            Vec3::new(2.279e11, 0.0, 0.0),
            Vec3::new(0.0, 24077.0, 0.0),
        ),
    ];

    Scenario {
        name: "Solar System (Earth and Mars)".to_string(),
        bodies,
        time_step: 3600.0,           // 1 hour
        duration_days: 365.25,       // 1 year
        print_interval_days: 30.0,   // Every 30 days
    }
}

/// Create a perfect two-body circular orbit scenario
///
/// This creates two equal-mass bodies orbiting their common center of mass (barycenter).
/// The system is set up in the center-of-mass frame where:
/// - Total momentum is zero
/// - Bodies orbit in a circular pattern
/// - Angular momentum and total energy should be conserved
///
/// For circular orbits: v = sqrt(G * M_total / (2 * r))
/// where r is the distance from each body to the center of mass
/// and M_total is the sum of both masses.
///
/// Orbital period: T = 2π * sqrt(r³ / (G * M_total / 2))
pub fn two_body_circular() -> Scenario {
    const G: f64 = 6.67430e-11;

    // Two equal masses (e.g., binary star system)
    let mass = 1.0e30; // kg (approximately half a solar mass each)

    // Separation distance: distance from each mass to the center of mass
    let separation = 1.0e11; // m (100 million km from barycenter each)

    // For circular orbit with equal masses:
    // v = sqrt(G * M_total / (2 * r)) = sqrt(G * 2m / (2 * r)) = sqrt(G * m / r)
    let orbital_velocity = (G * mass / separation).sqrt();

    let bodies = vec![
        Body::new(
            "Body1",
            mass,
            Vec3::new(separation, 0.0, 0.0), // On the positive x-axis
            Vec3::new(0.0, orbital_velocity, 0.0), // Moving in +y direction
        ),
        Body::new(
            "Body2",
            mass,
            Vec3::new(-separation, 0.0, 0.0), // On the negative x-axis
            Vec3::new(0.0, -orbital_velocity, 0.0), // Moving in -y direction
        ),
    ];

    Scenario {
        name: "Two-Body Circular Orbit".to_string(),
        bodies,
        time_step: 10.0,           // 10 seconds (small for accuracy)
        duration_days: 50.0,       // Simulate ~50 days
        print_interval_days: 5.0,  // Print every 5 days
    }
}
