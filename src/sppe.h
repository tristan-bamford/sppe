#ifndef SPPE_H
#define SPPE_H

#include <functional>
#include <chrono>

#include "particle.h"
#include "spatial_map.h"

namespace SPPE {
  
  using Unary_force_f = std::function<Vector_type(const Particle&)>;
  
  class System {
  public:

    // Construct an SPPE particle system that can handle a max of n particles.
    System(std::size_t n = 20'000)
      : max_particles_(n), spatial_map_(1, n) { particles_.reserve(n); }
    
    // Advance the system state in real time.
    Float_type run();
    // Advance the system state by dt seconds.
    Float_type step(Float_type dt = 1/60.0);

    // Reset the system state.
    void reset();
    // Reset the simulation time.
    void reset_time()
    { simulation_time_ = std::chrono::high_resolution_clock::now(); }
    
    // Add a particle to the system. Returns its pointer, or nullptr on failure.
    Particle* add_particle(const Particle& particle);
    
    // Add a force function to the system.
    void add_force_function(const Unary_force_f& f) 
    { unary_forces_.push_back(f); }

    // Define a bounded region for particles to be contained.
    void set_boundary(Vector_type boundary) { boundary_ = boundary; }
    // Returns the boundary vector.
    const Vector_type& get_boundary() const { return boundary_; }

    // Enable/Disable the boundary.
    void enable_boundary(bool enable = true)
    { opt_boundary_enabled_ = enable; }
    // Enable/Disable collisions.
    void enable_collisions(bool enable = true) 
    { opt_collisions_enabled_ = enable; }

    // Returns a std::span<Particle> of the particles in the system.
    auto particles() const { return std::span{particles_}; }

    // Return the number of particles in the system.
    auto n_particles() const { return particles_.size(); }
    // Return the maximum number of particles that can be handled by the system.
    auto max_particles() const { return max_particles_; }

    bool is_boundary_enabled() const { return opt_boundary_enabled_; }
    bool is_collisions_enabled() const { return opt_collisions_enabled_; }    
  private:
    // Apply unary force functions to an individual particle.
    void apply_unary_forces(Particle& particle) 
    { for (const auto& f : unary_forces_) particle.apply_force(f(particle)); }
    // Update the system state.
    void update(Float_type dt);
    // Check and resolve boundary collisions.
    void check_boundary(Particle& particle);
    // Resolve particle-particle collisions.
    void resolve_collisions();
  private:
    const std::size_t max_particles_;
    Spatial_map spatial_map_;
    std::vector<Particle> particles_;           // particle data
    std::vector<Unary_force_f> unary_forces_;   // unary force functions
    
    // The simulation boundary can be defined as a rectangular region using a 
    // single vector that points to a corner of the region.
    Vector_type boundary_{1000.0};
    
    // The simulation time is used by run() to determine the dt step since last
    // call.
    std::chrono::time_point<std::chrono::high_resolution_clock> 
        simulation_time_{std::chrono::high_resolution_clock::now()};

    // The largest radius is tracked and used by spatial map queries.
    Float_type largest_radius_{0.0};

    bool opt_boundary_enabled_{true};
    bool opt_collisions_enabled_{true};
  };
} // end namespace SPPE
#endif//SPPE_H