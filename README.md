# Atlas - N-Body Gravity Simulator

A gravitational N-body simulator implemented in Rust with multiple integration methods.

## Features

- **Two-body circular orbit scenario**: Perfect test case with analytical solutions
- **Solar system simulation**: Earth and Mars orbiting the Sun  
- **Multiple integration methods**: Euler, Verlet, RK4, and Leapfrog
- **Conservation tracking**: Monitor energy and angular momentum drift
- **Comprehensive test suite**: Validates conservation laws and orbital mechanics

## Running Simulations

### Solar System with Different Integrators
```bash
# Verlet (recommended - symplectic, excellent energy conservation)
cargo run --release -- solar verlet

# RK4 (high accuracy, 4th order)
cargo run --release -- solar rk4

# Leapfrog (symplectic, similar to Verlet)
cargo run --release -- solar leapfrog

# Euler (basic method, energy drift issues)
cargo run --release -- solar euler
```

### Two-Body Circular Orbit
```bash
# With Verlet integration
cargo run --release -- two-body verlet

# With RK4 integration
cargo run --release -- two-body rk4

# With Leapfrog integration
cargo run --release -- two-body leapfrog

# With Euler integration
cargo run --release -- two-body euler
```

## Integration Methods

### Euler Method
- **Order**: O(dt) - First-order accuracy
- **Speed**: Fast (1 force evaluation per step)
- **Energy Conservation**: Poor (not symplectic)
- **Use case**: Quick prototypes, short simulations

The Euler method is simple but suffers from energy drift. Orbits tend to expand over time as energy is artificially added to the system.

### Velocity Verlet Method ⭐ Recommended
- **Order**: O(dt²) - Second-order accuracy  
- **Speed**: Fast (1 force evaluation per step)
- **Energy Conservation**: Excellent (symplectic)
- **Use case**: Long-term orbital simulations

The Velocity Verlet method is a symplectic integrator that conserves energy much better than Euler. It's the preferred method for orbital mechanics simulations.

### Runge-Kutta 4th Order (RK4)
- **Order**: O(dt⁴) - Fourth-order accuracy
- **Speed**: Slower (4 force evaluations per step)
- **Energy Conservation**: Very good (not symplectic)
- **Use case**: High-precision simulations, general dynamics

RK4 provides excellent accuracy through four intermediate calculations per step. While not symplectic, its high order makes it very accurate for general dynamical systems.

### Leapfrog Method
- **Order**: O(dt²) - Second-order accuracy
- **Speed**: Fast (2 force evaluations per step)
- **Energy Conservation**: Excellent (symplectic)
- **Use case**: Long-term orbital simulations

The Leapfrog method is another symplectic integrator with similar properties to Verlet. It uses a "kick-drift-kick" formulation where positions and velocities are updated at interleaved time points.

**Energy Conservation Comparison** (1 year, dt=6 hours, solar system):
- Euler: 0.0001% energy drift
- Verlet: 0.0000% energy drift
- RK4: 0.0000% energy drift
- Leapfrog: 0.0000% energy drift

**Large Timestep Performance** (7-day timesteps, solar system):
- Euler: 0.0282% energy drift
- Verlet: 0.0001% energy drift (**282x better!**)

This demonstrates the key advantage of symplectic integrators (Verlet, Leapfrog): they maintain energy conservation even with larger timesteps, making them essential for long-term orbital simulations.

## The Two-Body Problem

The two-body problem is a classical mechanics problem of predicting the motion of two massive bodies orbiting each other under gravitational attraction. It's one of the few N-body problems with an exact analytical solution.

### Key Properties

For a perfect two-body circular orbit:

1. **Conservation of Energy**: Total mechanical energy (kinetic + potential) remains constant
   - `E = ½m₁v₁² + ½m₂v₂² - Gm₁m₂/r`

2. **Conservation of Angular Momentum**: Total angular momentum magnitude is constant
   - `L = m₁(r₁ × v₁) + m₂(r₂ × v₂)`

3. **Center of Mass Motion**: The barycenter (center of mass) remains stationary
   - Both bodies orbit their common center of mass

4. **Orbital Period**: For circular orbits, the period is:
   - `T = 2π√(a³/(G(m₁+m₂)))` where `a` is the separation distance

### Our Two-Body Scenario

The `two_body_circular()` scenario creates:
- Two equal masses (1.0 × 10³⁰ kg each, ~0.5 solar masses)
- Separated by 200 million km (100 million km from center each)
- Initial velocities set for perfect circular orbit: `v = √(Gm/r)`
- Orbital period: ~563 days

## Tests

Run the test suite to verify conservation laws:

```bash
cargo test
```

The tests verify:
- **Energy Conservation**: Energy drift remains bounded over time
- **Angular Momentum Conservation**: Angular momentum is well-preserved
- **Center of Mass Conservation**: Barycenter remains near origin
- **Orbital Stability**: Separation distance stays roughly constant
- **Verlet vs Euler**: Verlet conserves energy significantly better

### Integration Method Performance

Tests show that Verlet integration conserves energy at least 10x better than Euler for the same time step, making it essential for accurate long-term simulations.

## Project Structure

```
src/
  main.rs         - Entry point and simulation loop
  integrators.rs  - Integration methods (Euler, Verlet)
  physics.rs      - Core physics types and calculations
  scenarios.rs    - Predefined body configurations
```

## Physics Constants

- Gravitational constant: `G = 6.67430×10⁻¹¹ m³ kg⁻¹ s⁻²`
- 1 Astronomical Unit (AU): `1.496×10¹¹ m`

## License

MIT
