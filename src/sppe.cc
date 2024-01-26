#include "sppe.h"

#include "spatial_map.h"
#include "vector.h"

namespace SPPE {

  // Step through the particle system, where dt = time, in seconds, 
  // elapsed since last run call or construction if called for the first time. 
  // reset_time() can be called to reset the time since last call.
  Float_type 
  System::run()
  {
    const auto t = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<Float_type> seconds_elapsed(t - simulation_time_);
    simulation_time_ = t;

    return step(seconds_elapsed.count());
  }

  // Step through the particle system, where dt = time step in seconds.
  Float_type 
  System::step(Float_type dt)
  {
    // dynamically create new spatial map each step; spacing is based on largest 
    // radius
    Spatial_map spatial_map(particles_.size(), largest_radius_ * 2);
    
    // First pass - Update particles, check against the boundary and initialize 
    // the index table of the spatial map.
    update(spatial_map, dt);

    // Second pass - finish initializing the index table and build the particle 
    // indices(pointers) in the spatial map.
    spatial_map.build();

    // Third pass - Now that the spatial map is set up, we can act on
    // relationships between particles, ie. spatial forces and collision
    // detection. Additionally, unary forces can be applied in this pass.
    resolve_collisions(spatial_map);
    // NOTE: if no spatial relationships exist apply_unary_forces() can by 
    // called in the first pass after(of before) the update step, and everything
    // can be done in one pass.
    return dt;
  }

  // Update the system state, and initialize the spatial map.
  void 
  System::update(Spatial_map& spatial_map, Float_type dt)
  {
    for (Particle& particle : particles_)
    {
      particle.update(dt); // update position, velocity, etc.
      if (is_boundary_enabled()) check_boundary(particle);
      spatial_map.index(particle);
      // if a particle radius is changed, largest_radius_ will need to be 
      // checked again.
      if (particle.radius() > largest_radius_) {
        largest_radius_ = particle.radius();
      }
    }
  }

  // Iterate through the particles and test for spatial relationships 
  // (collisions)
  void 
  System::resolve_collisions(const Spatial_map &spatial_map)
  {
    for (auto& particle : particles_)
    {
      // Not involved in resolving collisions, but unary forces can be applied 
      // in this pass.
      apply_unary_forces(particle);

      const int query_r = spatial_map.discretize(particle.radius() + largest_radius_);
      const int query_x = spatial_map.discretize(particle.position_[0]);
      const int query_y = spatial_map.discretize(particle.position_[1]);

      for (int x = query_x - query_r; x <= query_x + query_r; ++x)
      {
        for (int y = query_y - query_r; y <= query_y + query_r; ++y)
        {
          for (Particle* p2 : spatial_map.query(x, y))
          {
            // Since iteration occurs through a contiguous array with pointer, 
            // p1: ++p1 > p1. Duplicate collision tests can be avoided by 
            // testing p2 > p1. This will also skip same particles and nullptr.
            if (p2 > &particle)
            {
              particle.collide(*p2);
            }
          }
        }
      }
    }
  }

  // Reset and clear system state.
  void
  System::reset()
  {
    unary_forces_.clear();
    particles_.clear();
    reset_time();
  }

  // Add a particle to the system and return its pointer.
  Particle* 
  System::add_particle(const Particle& particle)
  {
    if (!(particles_.size() < particles_.capacity())) return nullptr;

    // Track the largest radius for use in spatial mapping
    if (particle.radius() > largest_radius_) largest_radius_ = particle.radius();

    particles_.push_back(particle);
    return &particles_.back();
  }

  // Test for and resolve particle collisions with the "world" boundary, works 
  // for any dimension.
  void 
  System::check_boundary(Particle &particle)
  {
    for (std::size_t i = 0; i < Vector_type::size(); ++i)
    {
      Float_type& v = particle.velocity_[i];
      Float_type& x = particle.position_[i];
      Float_type boundary = boundary_[i] - particle.radius();
      if (x < -boundary)
      {
        x = -boundary;
        v = -v * particle.elasticity();
      }
      if (x > boundary)
      {
        x = boundary;
        v = -v * particle.elasticity();
      }
    }
  }
} // namespace SPPE