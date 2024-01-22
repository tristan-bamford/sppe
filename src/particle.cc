#include "particle.h"

namespace SPPE {

  // Update the particle state and clear the force variable.
  void
  Particle::update(Float_type dt)
  {
    // Update using the mid-point method
    const Vector_type previous_velocity = velocity_;
    velocity_ += force_ * (dt / mass_);
    position_ += (velocity_ + previous_velocity) * (dt * 0.5);
    force_ = 0;
  }

  // Particle-Particle collision solver
  void
  Particle::collide(Particle& p2)
  {
    auto& p1 = *this; // convenient alias

    Vector_type normal = p1.position_ - p2.position_;
    const Float_type sum_of_radii = p1.radius_ + p2.radius_;

    //if (std::abs(normal[0]) > sum_of_radii || std::abs(normal[1] > sum_of_radii)) return;

    Float_type distance_sqrd = dot_product(normal, normal);
    if (distance_sqrd < (sum_of_radii * sum_of_radii))
    {
      if (distance_sqrd == 0)
      {
        normal = 1.0;
        distance_sqrd = 2;
      }
      const Float_type distance = std::sqrt(distance_sqrd);
      const Float_type penetration = sum_of_radii - distance;
      const Vector_type unit_normal = normal / distance;
      const Vector_type displacement = unit_normal * (penetration * 0.5);
      p2.position_ -= displacement;
      p1.position_ += displacement;

      // Apply force when particles overlap
      const auto f =  50.0 * std::pow(penetration+1,2);
      p2.apply_force(unit_normal * p2.mass_ * -f);
      p1.apply_force(unit_normal * p1.mass_ * f);

      // The following condition checks to see if the particles are
      // moving towards each other.
      const Vector_type delta_v = p1.velocity_ - p2.velocity_;
      const Float_type incidence = dot_product(delta_v, normal);
      if (incidence < 0)
      {
        const Float_type mf = 2.0 / (p1.mass_ + p2.mass_);
        const Vector_type proj = normal * (incidence / distance_sqrd);// projection(normal, delta_v);

        // FIXME: Need to derive a proper vector equation for inelastic
        // collisions. Using elasticity like this isn't exactly correct,
        // but it produces decent results under most circumstances.
        p1.velocity_ -= (p2.mass_ * mf * p2.elasticity_) * proj;
        p2.velocity_ += (p1.mass_ * mf * p1.elasticity_) * proj;
      }
    }
  }
} //end namespace SPPE