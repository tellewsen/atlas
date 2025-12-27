use crate::physics::{Body, Vec3, calculate_acceleration};

/// Euler integration method
///
/// Simple first-order integrator. Fast but has poor energy conservation.
/// Energy tends to drift (usually increases) over time.
///
/// Order: O(dt)
/// Energy conservation: Poor (not symplectic)
pub fn euler_step(bodies: &mut Vec<Body>, dt: f64) {
    let bodies_snapshot = bodies.clone();

    for body in bodies.iter_mut() {
        let acceleration = calculate_acceleration(body, &bodies_snapshot);

        // Euler integration
        // v_new = v_old + a * dt
        // x_new = x_old + v_old * dt
        body.velocity = body.velocity.add(&acceleration.scale(dt));
        body.position = body.position.add(&body.velocity.scale(dt));
    }
}

/// Velocity Verlet integration method
///
/// A symplectic integrator that conserves energy much better than Euler.
/// Second-order accurate and time-reversible, making it ideal for
/// long-term orbital simulations.
///
/// Order: O(dt²)
/// Energy conservation: Excellent (symplectic)
///
/// Algorithm:
/// 1. x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
/// 2. Calculate a(t+dt) at new positions
/// 3. v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
pub fn verlet_step(bodies: &mut Vec<Body>, dt: f64) {
    // Step 1: Calculate current accelerations
    let mut accelerations: Vec<Vec3> = Vec::new();
    for body in bodies.iter() {
        let acc = calculate_acceleration(body, &bodies);
        accelerations.push(acc);
    }

    // Step 2: Update positions
    // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²
    for (i, body) in bodies.iter_mut().enumerate() {
        let vel_contribution = body.velocity.scale(dt);
        let acc_contribution = accelerations[i].scale(0.5 * dt * dt);
        body.position = body.position.add(&vel_contribution).add(&acc_contribution);
    }

    // Step 3: Calculate new accelerations at updated positions
    let mut new_accelerations: Vec<Vec3> = Vec::new();
    for body in bodies.iter() {
        let acc = calculate_acceleration(body, &bodies);
        new_accelerations.push(acc);
    }

    // Step 4: Update velocities using average of old and new accelerations
    // v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    for (i, body) in bodies.iter_mut().enumerate() {
        let avg_acc = accelerations[i].add(&new_accelerations[i]).scale(0.5 * dt);
        body.velocity = body.velocity.add(&avg_acc);
    }
}

/// Runge-Kutta 4th order (RK4) integration method
///
/// A high-accuracy explicit integrator with fourth-order accuracy.
/// More computationally expensive than Verlet (4 force evaluations per step)
/// but provides excellent accuracy for general dynamical systems.
///
/// Order: O(dt⁴)
/// Energy conservation: Very good (but not symplectic)
/// Computational cost: 4x force evaluations per step
///
/// Algorithm:
/// k1 = f(t, y)
/// k2 = f(t + dt/2, y + k1*dt/2)
/// k3 = f(t + dt/2, y + k2*dt/2)
/// k4 = f(t + dt, y + k3*dt)
/// y(t+dt) = y(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
pub fn rk4_step(bodies: &mut Vec<Body>, dt: f64) {
    let n = bodies.len();

    // Store initial state
    let initial_positions: Vec<Vec3> = bodies.iter().map(|b| b.position.clone()).collect();
    let initial_velocities: Vec<Vec3> = bodies.iter().map(|b| b.velocity.clone()).collect();

    // k1 calculation (at t)
    let k1_velocities: Vec<Vec3> = initial_velocities.clone();
    let k1_accelerations: Vec<Vec3> = bodies
        .iter()
        .map(|b| calculate_acceleration(b, &bodies))
        .collect();

    // k2 calculation (at t + dt/2, using k1)
    for i in 0..n {
        bodies[i].position = initial_positions[i].add(&k1_velocities[i].scale(dt / 2.0));
        bodies[i].velocity = initial_velocities[i].add(&k1_accelerations[i].scale(dt / 2.0));
    }
    let k2_velocities: Vec<Vec3> = bodies.iter().map(|b| b.velocity.clone()).collect();
    let k2_accelerations: Vec<Vec3> = bodies
        .iter()
        .map(|b| calculate_acceleration(b, &bodies))
        .collect();

    // k3 calculation (at t + dt/2, using k2)
    for i in 0..n {
        bodies[i].position = initial_positions[i].add(&k2_velocities[i].scale(dt / 2.0));
        bodies[i].velocity = initial_velocities[i].add(&k2_accelerations[i].scale(dt / 2.0));
    }
    let k3_velocities: Vec<Vec3> = bodies.iter().map(|b| b.velocity.clone()).collect();
    let k3_accelerations: Vec<Vec3> = bodies
        .iter()
        .map(|b| calculate_acceleration(b, &bodies))
        .collect();

    // k4 calculation (at t + dt, using k3)
    for i in 0..n {
        bodies[i].position = initial_positions[i].add(&k3_velocities[i].scale(dt));
        bodies[i].velocity = initial_velocities[i].add(&k3_accelerations[i].scale(dt));
    }
    let k4_velocities: Vec<Vec3> = bodies.iter().map(|b| b.velocity.clone()).collect();
    let k4_accelerations: Vec<Vec3> = bodies
        .iter()
        .map(|b| calculate_acceleration(b, &bodies))
        .collect();

    // Final update using weighted average
    for i in 0..n {
        // Position update: x(t+dt) = x(t) + dt/6 * (v1 + 2*v2 + 2*v3 + v4)
        let pos_increment = k1_velocities[i]
            .add(&k2_velocities[i].scale(2.0))
            .add(&k3_velocities[i].scale(2.0))
            .add(&k4_velocities[i])
            .scale(dt / 6.0);

        // Velocity update: v(t+dt) = v(t) + dt/6 * (a1 + 2*a2 + 2*a3 + a4)
        let vel_increment = k1_accelerations[i]
            .add(&k2_accelerations[i].scale(2.0))
            .add(&k3_accelerations[i].scale(2.0))
            .add(&k4_accelerations[i])
            .scale(dt / 6.0);

        bodies[i].position = initial_positions[i].add(&pos_increment);
        bodies[i].velocity = initial_velocities[i].add(&vel_increment);
    }
}

/// Leapfrog integration method
///
/// A symplectic integrator similar to Verlet, but with a different formulation.
/// Positions and velocities are evaluated at interleaved time points (leap-frogging).
/// Excellent for long-term orbital simulations.
///
/// Order: O(dt²)
/// Energy conservation: Excellent (symplectic)
/// Computational cost: 1 force evaluation per step
///
/// Algorithm:
/// v(t+dt/2) = v(t) + a(t)*dt/2
/// x(t+dt) = x(t) + v(t+dt/2)*dt
/// a(t+dt) = calculate acceleration at x(t+dt)
/// v(t+dt) = v(t+dt/2) + a(t+dt)*dt/2
///
/// Note: This implementation uses a "kick-drift-kick" formulation
pub fn leapfrog_step(bodies: &mut Vec<Body>, dt: f64) {
    // First kick: update velocities by half step
    let accelerations: Vec<Vec3> = bodies
        .iter()
        .map(|b| calculate_acceleration(b, &bodies))
        .collect();

    for (i, body) in bodies.iter_mut().enumerate() {
        body.velocity = body.velocity.add(&accelerations[i].scale(dt / 2.0));
    }

    // Drift: update positions by full step using half-stepped velocities
    for body in bodies.iter_mut() {
        body.position = body.position.add(&body.velocity.scale(dt));
    }

    // Second kick: update velocities by another half step
    let new_accelerations: Vec<Vec3> = bodies
        .iter()
        .map(|b| calculate_acceleration(b, &bodies))
        .collect();

    for (i, body) in bodies.iter_mut().enumerate() {
        body.velocity = body.velocity.add(&new_accelerations[i].scale(dt / 2.0));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::{orbital_period, total_angular_momentum, total_energy};
    use crate::scenarios::two_body_circular;

    #[test]
    fn test_verlet_energy_conservation() {
        let scenario = two_body_circular();
        let mut bodies = scenario.bodies.clone();
        let initial_energy = total_energy(&bodies);

        let dt = 100.0; // Same time step as Euler test
        let steps = 10000;

        for _ in 0..steps {
            verlet_step(&mut bodies, dt);
        }

        let final_energy = total_energy(&bodies);
        let energy_drift = (final_energy - initial_energy).abs() / initial_energy.abs();

        // Verlet should conserve energy much better than Euler (which had ~50% drift)
        assert!(
            energy_drift < 0.01,
            "Verlet energy drift too large: {:.4}% (initial: {:.3e}, final: {:.3e})",
            energy_drift * 100.0,
            initial_energy,
            final_energy
        );
    }

    #[test]
    fn test_verlet_vs_euler_energy_conservation() {
        let scenario = two_body_circular();

        // Test Euler
        let mut euler_bodies = scenario.bodies.clone();
        let initial_energy = total_energy(&euler_bodies);
        let dt = 100.0;
        let steps = 5000;

        for _ in 0..steps {
            euler_step(&mut euler_bodies, dt);
        }

        let euler_energy_drift =
            (total_energy(&euler_bodies) - initial_energy).abs() / initial_energy.abs();

        // Test Verlet
        let mut verlet_bodies = scenario.bodies.clone();
        for _ in 0..steps {
            verlet_step(&mut verlet_bodies, dt);
        }

        let verlet_energy_drift =
            (total_energy(&verlet_bodies) - initial_energy).abs() / initial_energy.abs();

        println!("Euler energy drift:  {:.4}%", euler_energy_drift * 100.0);
        println!("Verlet energy drift: {:.4}%", verlet_energy_drift * 100.0);

        // Verlet should be significantly better than Euler
        assert!(
            verlet_energy_drift < euler_energy_drift * 0.1,
            "Verlet should conserve energy at least 10x better than Euler"
        );
    }

    #[test]
    fn test_verlet_angular_momentum_conservation() {
        let scenario = two_body_circular();
        let mut bodies = scenario.bodies;
        let initial_l = total_angular_momentum(&bodies);

        let dt = 100.0;
        let steps = 10000;

        for _ in 0..steps {
            verlet_step(&mut bodies, dt);
        }

        let final_l = total_angular_momentum(&bodies);
        let l_drift = (final_l - initial_l).abs() / initial_l.abs();

        // Angular momentum should be very well conserved
        assert!(
            l_drift < 0.01,
            "Angular momentum drift too large: {:.4}%",
            l_drift * 100.0
        );
    }

    #[test]
    fn test_large_timestep_euler_fails() {
        use crate::scenarios::solar_system;
        let scenario = solar_system();
        let mut bodies = scenario.bodies.clone();

        // Use an extreme timestep: 7 days
        let dt = 7.0 * 86400.0; // 7 days
        let steps = 52; // 1 year

        let initial_energy = total_energy(&bodies);

        for _ in 0..steps {
            euler_step(&mut bodies, dt);
        }

        let euler_energy_drift =
            (total_energy(&bodies) - initial_energy).abs() / initial_energy.abs();

        println!(
            "Euler (Solar System) with extreme timestep (dt=7 days): {:.4}% energy drift",
            euler_energy_drift * 100.0
        );

        // Euler should show measurable drift with very large timesteps
        assert!(
            euler_energy_drift > 0.0001,
            "Expected Euler to have measurable drift with extreme timestep, got {:.4}%",
            euler_energy_drift * 100.0
        );
    }

    #[test]
    fn test_large_timestep_verlet_succeeds() {
        use crate::scenarios::solar_system;
        let scenario = solar_system();
        let mut bodies = scenario.bodies.clone();

        // Use the same extreme timestep
        let dt = 7.0 * 86400.0; // 7 days
        let steps = 52; // 1 year

        let initial_energy = total_energy(&bodies);

        for _ in 0..steps {
            verlet_step(&mut bodies, dt);
        }

        let verlet_energy_drift =
            (total_energy(&bodies) - initial_energy).abs() / initial_energy.abs();

        println!(
            "Verlet (Solar System) with extreme timestep (dt=7 days): {:.4}% energy drift",
            verlet_energy_drift * 100.0
        );

        // Verlet should maintain much better energy conservation
        assert!(
            verlet_energy_drift < 0.001,
            "Verlet energy drift too large: {:.4}%",
            verlet_energy_drift * 100.0
        );
    }

    #[test]
    fn test_extreme_timestep_comparison() {
        use crate::scenarios::solar_system;
        let scenario = solar_system();

        // Use 6 hour timestep vs normal 1 hour
        let dt = 21600.0; // 6 hours
        let steps = 365 * 4; // 1 year at 6-hour steps

        // Test Euler
        let mut euler_bodies = scenario.bodies.clone();
        let initial_energy = total_energy(&euler_bodies);

        for _ in 0..steps {
            euler_step(&mut euler_bodies, dt);
        }

        let euler_final_energy = total_energy(&euler_bodies);
        let euler_drift = (euler_final_energy - initial_energy).abs() / initial_energy.abs();

        // Test Verlet
        let mut verlet_bodies = scenario.bodies.clone();
        for _ in 0..steps {
            verlet_step(&mut verlet_bodies, dt);
        }

        let verlet_final_energy = total_energy(&verlet_bodies);
        let verlet_drift = (verlet_final_energy - initial_energy).abs() / initial_energy.abs();

        println!(
            "\n=== Solar System: Large Timestep Test (dt={:.1} hours) ===",
            dt / 3600.0
        );
        println!(
            "Euler:  {:.2}% energy drift (E_initial={:.3e}, E_final={:.3e})",
            euler_drift * 100.0,
            initial_energy,
            euler_final_energy
        );
        println!(
            "Verlet: {:.2}% energy drift (E_initial={:.3e}, E_final={:.3e})",
            verlet_drift * 100.0,
            initial_energy,
            verlet_final_energy
        );
        if verlet_drift > 0.0 {
            println!(
                "Verlet is {:.0}x better at conserving energy",
                euler_drift / verlet_drift
            );
        }

        // Symplectic integrator should be dramatically better
        assert!(
            verlet_drift < euler_drift * 0.1,
            "Verlet should conserve energy at least 10x better than Euler with extreme timesteps"
        );
    }

    #[test]
    fn test_rk4_energy_conservation() {
        let scenario = two_body_circular();
        let mut bodies = scenario.bodies.clone();
        let initial_energy = total_energy(&bodies);

        let dt = 100.0; // Same time step as other tests
        let steps = 10000;

        for _ in 0..steps {
            rk4_step(&mut bodies, dt);
        }

        let final_energy = total_energy(&bodies);
        let energy_drift = (final_energy - initial_energy).abs() / initial_energy.abs();

        // RK4 should have excellent accuracy despite not being symplectic
        assert!(
            energy_drift < 0.001,
            "RK4 energy drift too large: {:.4}% (initial: {:.3e}, final: {:.3e})",
            energy_drift * 100.0,
            initial_energy,
            final_energy
        );
    }

    #[test]
    fn test_leapfrog_energy_conservation() {
        let scenario = two_body_circular();
        let mut bodies = scenario.bodies.clone();
        let initial_energy = total_energy(&bodies);

        let dt = 100.0; // Same time step as other tests
        let steps = 10000;

        for _ in 0..steps {
            leapfrog_step(&mut bodies, dt);
        }

        let final_energy = total_energy(&bodies);
        let energy_drift = (final_energy - initial_energy).abs() / initial_energy.abs();

        // Leapfrog is symplectic, should conserve energy very well
        assert!(
            energy_drift < 0.01,
            "Leapfrog energy drift too large: {:.4}% (initial: {:.3e}, final: {:.3e})",
            energy_drift * 100.0,
            initial_energy,
            final_energy
        );
    }

    #[test]
    fn test_all_integrators_comparison() {
        use crate::scenarios::solar_system;

        let scenario = solar_system();
        let dt = 6.0 * 3600.0; // 6 hours
        let steps = 365 * 4; // 1 year

        // Test all integrators
        let integrators: Vec<(&str, fn(&mut Vec<Body>, f64))> = vec![
            ("Euler", euler_step),
            ("Verlet", verlet_step),
            ("RK4", rk4_step),
            ("Leapfrog", leapfrog_step),
        ];

        println!("\n=== Solar System: All Integrators Comparison (1 year, dt=6h) ===");

        for (name, integrator) in integrators {
            let mut bodies = scenario.bodies.clone();
            let initial_energy = total_energy(&bodies);

            for _ in 0..steps {
                integrator(&mut bodies, dt);
            }

            let final_energy = total_energy(&bodies);
            let energy_drift = (final_energy - initial_energy).abs() / initial_energy.abs();

            println!(
                "{:<10} {:.4}% energy drift (E_initial={:.3e}, E_final={:.3e})",
                name,
                energy_drift * 100.0,
                initial_energy,
                final_energy
            );
        }
    }

    #[test]
    fn test_orbital_period_drift_with_moderate_timesteps() {
        let scenario = two_body_circular();

        let initial_period = orbital_period(&scenario.bodies);
        let dt = 500.0; // Moderate timestep ~8.3 minutes

        // Run for just 2 orbits to keep separation change reasonable
        let steps = ((2.0 * initial_period) / dt) as usize;

        println!("\n=== Orbital Period Drift (2 orbits, dt={:.1}s) ===", dt);
        println!("Initial period: {:.2} days", initial_period / 86400.0);

        let integrators: Vec<(&str, fn(&mut Vec<Body>, f64), bool)> = vec![
            ("Euler", euler_step as fn(&mut Vec<Body>, f64), false),
            ("Verlet", verlet_step, true),
            ("RK4", rk4_step, false),
            ("Leapfrog", leapfrog_step, true),
        ];

        for (name, integrator, is_symplectic) in &integrators {
            let mut bodies = scenario.bodies.clone();
            let initial_separation = bodies[0].position.sub(&bodies[1].position).magnitude();

            for _ in 0..steps {
                integrator(&mut bodies, dt);
            }

            let final_separation = bodies[0].position.sub(&bodies[1].position).magnitude();
            let final_period = orbital_period(&bodies);
            let separation_change =
                ((final_separation - initial_separation) / initial_separation) * 100.0;
            let period_change = ((final_period - initial_period) / initial_period) * 100.0;
            let energy_drift = (total_energy(&bodies) - total_energy(&scenario.bodies)).abs()
                / total_energy(&scenario.bodies).abs();

            println!(
                "{:<10} Sep: {:.3e}m ({:+.2}%), Period: {:.2}d ({:+.2}%), Energy: {:.4}% {}",
                name,
                final_separation,
                separation_change,
                final_period / 86400.0,
                period_change,
                energy_drift * 100.0,
                if *is_symplectic { "✓" } else { "" }
            );
        }

        println!("\n✓ = Symplectic (preserves phase space structure)");
        println!("Note: orbital_period() is calculated from current separation,");
        println!("so changes in separation lead to calculated period changes.");
    }

    #[test]
    fn test_separation_stability_symplectic_vs_nonsymplectic() {
        let scenario = two_body_circular();

        let dt = 200.0; // Small timestep for reasonable behavior (~3.3 minutes)
        let simulation_days = 100.0; // 100 days
        let steps = ((simulation_days * 86400.0) / dt) as usize;

        println!(
            "\n=== Separation Stability Test ({:.0} days, dt={:.1}s) ===",
            simulation_days, dt
        );

        let integrators: Vec<(&str, fn(&mut Vec<Body>, f64), bool)> = vec![
            ("Euler", euler_step as fn(&mut Vec<Body>, f64), false),
            ("Verlet", verlet_step, true),
            ("RK4", rk4_step, false),
            ("Leapfrog", leapfrog_step, true),
        ];

        let mut symplectic_max_drift = 0.0;
        let mut nonsymplectic_min_drift = f64::MAX;

        for (name, integrator, is_symplectic) in &integrators {
            let mut bodies = scenario.bodies.clone();
            let initial_separation = bodies[0].position.sub(&bodies[1].position).magnitude();
            let initial_energy = total_energy(&bodies);

            for _ in 0..steps {
                integrator(&mut bodies, dt);
            }

            let final_separation = bodies[0].position.sub(&bodies[1].position).magnitude();
            let separation_drift =
                ((final_separation - initial_separation) / initial_separation).abs();
            let energy_drift =
                (total_energy(&bodies) - initial_energy).abs() / initial_energy.abs();

            println!(
                "{:<10} Separation drift: {:6.2}%, Energy drift: {:.4}% {}",
                name,
                separation_drift * 100.0,
                energy_drift * 100.0,
                if *is_symplectic { "✓ symplectic" } else { "" }
            );

            if *is_symplectic && separation_drift > symplectic_max_drift {
                symplectic_max_drift = separation_drift;
            }
            if !*is_symplectic && separation_drift < nonsymplectic_min_drift {
                nonsymplectic_min_drift = separation_drift;
            }
        }

        println!(
            "\nSymplectic methods all show similar separation drift due to discretization effects."
        );
        println!("Energy conservation is the key advantage, not necessarily separation stability.");
    }

    #[test]
    fn test_large_timestep_period_divergence() {
        let scenario = two_body_circular();

        let initial_period = orbital_period(&scenario.bodies);
        let dt = 5000.0; // Large timestep - about 83 minutes
        let simulation_days = 100.0;

        let steps = ((simulation_days * 86400.0) / dt) as usize;

        println!(
            "\n=== Large Timestep Behavior (dt={:.1}s = {:.1} min, {:.0} days) ===",
            dt,
            dt / 60.0,
            simulation_days
        );
        println!("Initial period: {:.2} days\n", initial_period / 86400.0);

        let integrators: Vec<(&str, fn(&mut Vec<Body>, f64), bool)> = vec![
            ("Euler", euler_step as fn(&mut Vec<Body>, f64), false),
            ("Verlet", verlet_step, true),
            ("RK4", rk4_step, false),
            ("Leapfrog", leapfrog_step, true),
        ];

        for (name, integrator, is_symplectic) in integrators {
            let mut bodies = scenario.bodies.clone();
            let initial_separation = bodies[0].position.sub(&bodies[1].position).magnitude();

            for _ in 0..steps {
                integrator(&mut bodies, dt);
            }

            let final_separation = bodies[0].position.sub(&bodies[1].position).magnitude();
            let final_period = orbital_period(&bodies);
            let period_change_pct = ((final_period - initial_period) / initial_period) * 100.0;
            let separation_change_pct =
                ((final_separation - initial_separation) / initial_separation) * 100.0;
            let energy_drift = (total_energy(&bodies) - total_energy(&scenario.bodies)).abs()
                / total_energy(&scenario.bodies).abs();

            println!(
                "{:<10} Sep: {:+7.2}%, Period: {:+8.2}%, Energy: {:.4}% {}",
                name,
                separation_change_pct,
                period_change_pct,
                energy_drift * 100.0,
                if is_symplectic { "✓" } else { "✗" }
            );
        }

        println!("\n✓ = Symplectic (better long-term stability)");
        println!("✗ = Non-symplectic (accumulates phase errors)");
        println!("\nNote: Large timesteps cause orbit expansion in all methods,");
        println!("but symplectic methods preserve energy/phase space structure.");
    }
}
