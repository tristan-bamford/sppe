#ifndef SPPE_SPATIAL_MAP_H
#define SPPE_SPATIAL_MAP_H

#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "sppe.h"

namespace SPPE {
  
  // Sparse spatial hash map.
  // This works in two passes. In the first, index() is called on each particle
  // until the count_array_ is initialized to
  // the number of particles that map to each location in count_array_.
  // In the second pass, particle indices are arranged in a contiguous array 
  // such that "buckets" can then be queried as spans.
  // The map can be reset with reset(); 
  class Spatial_map {
    static constexpr std::size_t n_cells = 200 * 200; // a good number for n_cells is n_particles
    static constexpr std::size_t table_size = n_cells + 1;
  public:
    // Construct a spatial map meant for n particles
    Spatial_map(std::size_t n, Float_type spacing)
      : particle_indices_(n, nullptr), spacing_(spacing) { reset(); }

    // Reset the map.
    void reset() { std::fill_n(count_array_, table_size, 0); }
    
    // Return a span of particle pointers that correspond to the grid position {x,y}
    std::span<Particle* const> query(int x, int y) const;

    // "indexing" a particle consists of incrementing the count of the cell that 
    // it maps to.
    void index(const Particle& p) 
    { ++count_array_[hash(discretize(p.position()[0]), 
                          discretize(p.position()[1]))]; } // index a particle
    
    void build(const std::span<Particle>&);
    int discretize(Float_type x) const { return std::ceil(x / spacing_); }
  
  private:
    static constexpr int primes[] = { 92837111, 689287499, 283923481 };
    std::size_t hash(int x, int y) const 
    { return (x * primes[0] ^ y * primes[1]) % n_cells; }
  private:
    int count_array_[table_size];
    std::vector<Particle*> particle_indices_;
    const Float_type spacing_;
  };
    
  inline std::span<Particle* const>
  Spatial_map::query(int x, int y) const
  {
    const auto index = hash(x, y);
    const std::size_t n = count_array_[index + 1] - count_array_[index];
    return { particle_indices_.begin() + count_array_[index], n };
  }

  inline void
  Spatial_map::build(const std::span<Particle>& particle_span)
  {
    std::partial_sum(count_array_, count_array_ + table_size, count_array_);
    
    for (auto& particle : particle_span)
    {
      const auto index = hash(discretize(particle.position()[0]),
                              discretize(particle.position()[1]));

      particle_indices_[--count_array_[index]] = &particle;
    }
  }
} // namespace SPPE

#endif//SPPE_SPATIAL_MAP_H