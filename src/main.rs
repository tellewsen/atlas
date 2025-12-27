mod integrators;
mod physics;
mod scenarios;

use integrators::{euler_step, leapfrog_step, rk4_step, verlet_step};
use physics::{print_body_state, total_angular_momentum, total_energy};
use scenarios::{solar_system, two_body_circular};
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    // Parse command line arguments
    let mut integrator_name = "euler";
    let mut scenario_name = "solar";

    for arg in args.iter().skip(1) {
        match arg.as_str() {
            "two-body" => scenario_name = "two-body",
            "solar" => scenario_name = "solar",
            "euler" => integrator_name = "euler",
            "verlet" => integrator_name = "verlet",
            "rk4" => integrator_name = "rk4",
            "leapfrog" => integrator_name = "leapfrog",
            _ => {}
        }
    }

    let scenario = if scenario_name == "two-body" {
        two_body_circular()
    } else {
        solar_system()
    };

    let total_steps = scenario.total_steps();
    let print_interval = scenario.print_interval_steps();
    let mut bodies = scenario.bodies;

    println!("Starting gravity simulation: {}", scenario.name);
    println!("Integration method: {}", integrator_name);
    println!(
        "Time step: {} seconds ({} hour)",
        scenario.time_step,
        scenario.time_step / 3600.0
    );
    println!("Total duration: {} days\n", scenario.duration_days);

    // Print initial state
    println!("Initial state:");
    for body in &bodies {
        print_body_state(body);
    }
    println!("  Total energy: {:.3e} J", total_energy(&bodies));
    println!(
        "  Total angular momentum: {:.3e} kg·m²/s",
        total_angular_momentum(&bodies)
    );

    let initial_energy = total_energy(&bodies);
    let initial_angular_momentum = total_angular_momentum(&bodies);

    // Run simulation
    for step in 0..total_steps {
        match integrator_name {
            "euler" => euler_step(&mut bodies, scenario.time_step),
            "verlet" => verlet_step(&mut bodies, scenario.time_step),
            "rk4" => rk4_step(&mut bodies, scenario.time_step),
            "leapfrog" => leapfrog_step(&mut bodies, scenario.time_step),
            _ => panic!("Unknown integrator: {}", integrator_name),
        }

        // Print state periodically
        if step % print_interval == 0 && step > 0 {
            let days = (step as f64 * scenario.time_step) / 86400.0;
            println!("\nDay {:.1}:", days);
            for body in &bodies {
                print_body_state(body);
            }

            let current_energy = total_energy(&bodies);
            let energy_drift = ((current_energy - initial_energy) / initial_energy).abs();
            println!("  Energy drift: {:.2}%", energy_drift * 100.0);
        }
    }

    // Print final state
    println!("\nFinal state (after {} days):", scenario.duration_days);
    for body in &bodies {
        print_body_state(body);
    }

    let final_energy = total_energy(&bodies);
    let final_angular_momentum = total_angular_momentum(&bodies);
    let energy_drift = ((final_energy - initial_energy) / initial_energy).abs();
    let angular_momentum_drift =
        ((final_angular_momentum - initial_angular_momentum) / initial_angular_momentum).abs();

    println!("\nConservation analysis:");
    println!("  Energy drift: {:.2}%", energy_drift * 100.0);
    println!(
        "  Angular momentum drift: {:.2}%",
        angular_momentum_drift * 100.0
    );
}
