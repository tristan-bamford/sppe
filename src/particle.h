#ifndef SPPE_PARTICLE_H
#define SPPE_PARTICLE_H

#include "types.h"
#include "vector.h"

namespace SPPE {

  class Particle {
  public:
    //Particle() = delete;
    constexpr Particle(const Vector_type& position, 
                       const Float_type radius = 1.0,
                       const Vector_type& velocity = 0.0,
                       const Float_type density = 1.0,
                       const Float_type elasticity = 1.0);

    static constexpr Float_type compute_volume(Float_type r);

    // Update state variables
    void update(Float_type dt);

    // Physical properties
    Float_type radius()     const { return radius_; }
    Float_type mass()       const { return mass_; }
    Float_type density()    const { return mass_ / volume(); }
    Float_type elasticity() const { return elasticity_; }
    Float_type volume()     const { return compute_volume(radius_); }

    // State properties
    const Vector_type& position() const { return position_; }
    const Vector_type& velocity() const { return velocity_; }
    const Vector_type& force()    const { return force_; }

    Vector_type momentum() const { return velocity_ * mass_; }
    Float_type kinetic_energy() const 
    { return (mass_ * 0.5) * dot_product(velocity_, velocity_); }

    void apply_force(const Vector_type& f) { force_ += f; }
    
    // Resolve particle-particle collision
    void collide(Particle& other);
  public:
    // Public state variables
    Vector_type position_;
    Vector_type velocity_;
  private:
    Vector_type force_{0.0};
    Float_type radius_;
    Float_type mass_;
    Float_type elasticity_;
  };

  // Compute the Particle volume given the radius
  inline constexpr Float_type
  Particle::compute_volume(Float_type r)
  {
    if (Vector_type::size() == 3) return (4.0/3.0 * PI) * r * r * r;
    else return PI * r * r;
  }

  inline constexpr 
  Particle::Particle(const Vector_type& position, const Float_type radius,
                     const Vector_type& velocity, const Float_type density,
                     const Float_type elasticity)
    : position_(position), 
      velocity_(velocity),
      force_(0.0),
      radius_(radius), 
      mass_(compute_volume(radius) * density), 
      elasticity_(elasticity) {} 

} //end namespace SPPE

#endif//SPPE_PARTICLE_H