#ifndef SPPE_SPATIAL_MAP_H
#define SPPE_SPATIAL_MAP_H

#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "particle.h"

namespace SPPE {
  
  // Sparse spatial hash map.
  // This works in two passes. In the first, the number of particles that map to each cell are counted.index() is called on each particle
  // until the count_array_ is initialized to
  // the number of particles that map to each location in count_array_.
  // In the second pass, particle indices are arranged in a contiguous array 
  // such that "buckets" can then be queried as spans.
  // The map can be reset with reset();
  // Adding or updating particles invalidates the map
  class Spatial_map {
    static constexpr std::size_t n_cells = 200 * 200; // a good number for n_cells is n_particles
    static constexpr std::size_t table_size = n_cells + 1;
  public:
    // Construct a spatial map meant for n particles
    Spatial_map(std::size_t n, Float_type spacing)
      : spacing_(spacing) { reset(); particle_indices_.reserve(n); particles_.reserve(n); }

    // Reset the map.
    void reset() { std::fill_n(count_array_, table_size, 0); particles_.clear(); }
    
    // Return a span of particle pointers that correspond to the grid position {x,y}
    std::span<Particle* const> query(int x, int y) const;

    // "indexing" a particle consists of incrementing the count of the cell that 
    // it maps to.
    void index(Particle& p) 
    { 
      ++count_array_[hash(discretize(p.position()[0]), 
                          discretize(p.position()[1]))];
      particles_.push_back(&p);
      is_built_ = false;

    } // index a particle

    void build();
    int discretize(Float_type x) const { return std::ceil(x / spacing_); }
  
  private:
    static constexpr int primes[] = { 92837111, 689287499, 283923481 };
    std::size_t hash(int x, int y) const 
    { return (x * primes[0] ^ y * primes[1]) % n_cells; }
    std::size_t hash(const Particle& p)
    { return hash(discretize(p.position()[0]), discretize(p.position()[1])); }

  private:
    int count_array_[table_size];
    std::vector<Particle*> particle_indices_;
    std::vector<Particle*> particles_;
    const Float_type spacing_;
    bool is_built_{false};
  };
    
  inline std::span<Particle* const>
  Spatial_map::query(int x, int y) const
  {
    if (!is_built_) x = 0;

    const auto index = hash(x, y);
    const std::size_t n = count_array_[index + 1] - count_array_[index];
    // ASSERT: 0 < n <= size - offset, (size - offset) is positive when (offset < size)
    return { particle_indices_.begin() + count_array_[index], n };
  }

  inline void
  Spatial_map::build()

  {
    std::partial_sum(count_array_, count_array_ + table_size, count_array_);
    
    for (auto particle : particles_)
    {
      const auto index = hash(discretize(particle->position()[0]),
                              discretize(particle->position()[1]));

      particle_indices_[--count_array_[index]] = particle;
    }
    is_built_ = true;
  }
} // namespace SPPE

#endif//SPPE_SPATIAL_MAP_H