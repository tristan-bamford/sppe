# SPPE
SPPE (Simple Particle Physics Engine) is a small C++ library designed to simulate various particle interactions. It provides a basic framework for modeling physical systems with particles that is easy setup and covers many simple use cases. The goal here is to provide a fast (>20k particles at 60 frames/s) and compact system that runs on a single thread without use of a GPU.

## Features
- Particle simulation with customizable forces.
- Modeling of particles with varying radius, mass, elasticity, etc.
- Spatial mapping for efficient particle interaction calculations.
- Boundary control and collision handling.

## Prerequisites
C++ compiler with C++20 support. SPPE uses [num_array](https://github.com/tristan-bamford/num_array) to handle vector calculations.

## Usage
### Constructing an `SPPE::Particle`
```cpp
// SPPE vector type
SPPE::Vector_type position = { 0.0, 0.0 };
SPPE::Vector_type velocity = { 0.0, 0.0 };

// SPPE float type
SPPE::Float_type radius = 1.0;
SPPE::Float_type density = 1.0;
SPPE::Float_type elasticity = 1.0;

// Construct a Particle with minimum arguments
SPPE::Particle particle1(position, radius);

// Construct a Particle with all arguments
SPPE::Particle particle2(position, radius, velocity, density, elasticity);
```
### Constructing an `SPPE::System` and adding `Particles`
```cpp
// Create an SPPE::System with 20,000 maximum particles
SPPE::System sppe(20'000);

// Add some previously constructed Particles
sppe.add_particle(particle1);
sppe.add_particle(particle2);
// ...or more directly
sppe.add_particle({ { 0.0, 0.0 }, 1.0, { 0.0, 0.0 } });
```
### Running the `SPPE::System`
```cpp
// Optionally set the system boundaries by passing a vector that points to a 
// corner of the bounding box
sppe.set_boundary({ 1000.0, 1000.0 });

// Boundary can be enabled/disabled, default = enabled
sppe.enable_boundary(true);

// Collisions can be enabled/disabled, default = enabled
sppe.enable_collisions(true);

while(true)
{
  // Step through the system in real time
  sppe.run();

  // Access an std::span of the Particles in the System through 
  // SPPE::System::particles()
  for (const SPPE::Particle& particle : sppe.particles())
  {
    // Do stuff with 'particle'
    // ...
  }
}
```
### Adding custom unary force functions
```cpp
// Unary_force_f = std::function<Vector_type(const Particle&)>;

// Define a gravity function
const auto gravity = [](const SPPE::Particle& particle) -> SPPE::Vector_type
{ return { 0.0, 9.80665 * particle.mass() }; };

// Define viscous drag
const auto viscosity = [](const SPPE::Particle& particle) -> SPPE::Vector_type
{ return particle.velocity() * -0.5; };

// Add the forces to the system
sppe.add_force_function(gravity);
sppe.add_force_function(viscosity);
```
## Demos
### Classic particles bouncing around
https://github.com/tristan-bamford/sppe/assets/120840025/e5a5d0b8-03ab-4949-a454-fde884f75c6e
### Add a gravity function
https://github.com/tristan-bamford/sppe/assets/120840025/378e6420-21bc-4fc3-8622-34d4db15e2b0
### 10k particles
https://github.com/tristan-bamford/sppe/assets/120840025/e1638630-f37b-4fb4-870b-7a87f9586242

## Contributing

Contributions to this library are welcome! If you find any issues or have ideas for improvements, please open an issue or create a pull request on the [GitHub repository](https://github.com/tristan-bamford/sppe).

## License

This project is licensed under the [MIT License](LICENSE).
