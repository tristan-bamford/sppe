#ifndef SPPE_SPATIAL_MAP_H
#define SPPE_SPATIAL_MAP_H

#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "particle.h"

namespace SPPE {
  
  // Spatial_map is a sparse spatial hash map.
  // On the surface, it is meant to mimic the usage of a typical map; insert(), 
  // query(), clear(), etc...Behind the scenes, the map needs to be rebuilt 
  // every time a particle is inserted. In practice, this is actually done when 
  // query() is first called. This means that inserting and querying is 
  // optimally efficient when done in two separate passes.
  
  class Spatial_map {
    static constexpr std::size_t n_cells = 200 * 200; // a good number for n_cells is n_particles
    static constexpr std::size_t count_array_sz = n_cells + 1;
  public:
    // Construct a spatial map add reserve space for n particles
    Spatial_map(Float_type spacing, std::size_t n = n_cells);

    // Reset the map.
    void clear();

    // Return a span of particle pointers that correspond to a grid position
    std::span<Particle* const> query(int x, int y);

    // "insert" a particle into the map
    void insert(Particle& particle);

    // Set the cell size/spacing. Note: this will delete the map.
    void set_spacing(Float_type spacing) { spacing_ = spacing; clear(); }

    int discretize(Float_type x) const { return std::ceil(x / spacing_); }
  private:
    static constexpr int primes[] = { 92837111, 689287499, 283923481 };
    std::size_t hash(int x, int y) const 
    { return (x * primes[0] ^ y * primes[1]) % n_cells; }
    // convenience overload for hash
    std::size_t hash(const Particle& p)
    { return hash(discretize(p.position()[0]), discretize(p.position()[1])); }
    // Build the map given the particles inserted so far
    void build();
  private:
    Float_type spacing_;                  // granularity of the spatial grid
    std::vector<Particle*> particles_;    // array of inserted particle pointers
    std::vector<Particle*> buckets_;      // above array sorted into "buckets"
    int count_array_[count_array_sz];     // count of particles in each bucket
    bool is_built_{false};
  };

  inline
  Spatial_map::Spatial_map(Float_type spacing, std::size_t n)
    : spacing_(spacing), buckets_(n, nullptr)
  {
    clear();
    particles_.reserve(n);
  }

  inline void
  Spatial_map::insert(Particle& particle) 
  { 
    ++count_array_[hash(particle)];
    particles_.push_back(&particle);
    is_built_ = false; // Adding or updating particles invalidates the map
  }

  inline std::span<Particle* const>
  Spatial_map::query(int x, int y)
  {
    if (!is_built_) build();

    const auto i = hash(x, y);
    const std::size_t n = count_array_[i + 1] - count_array_[i];
    // ASSERT: 0 < n <= size - offset
    return { buckets_.begin() + count_array_[i], n };
  }

  inline void 
  Spatial_map::clear() 
  {
    std::fill_n(count_array_, count_array_sz, 0); 
    particles_.clear();
    is_built_ = false;
  }

  inline void
  Spatial_map::build()
  {
    std::partial_sum(count_array_, count_array_ + count_array_sz, count_array_);
    
    buckets_.resize(particles_.size());
    for (auto particle_p : particles_)
    {
      const auto index = --count_array_[hash(*particle_p)];
      buckets_[index] = particle_p;
    }
    is_built_ = true;
  }
} // namespace SPPE

#endif//SPPE_SPATIAL_MAP_H