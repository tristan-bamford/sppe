#ifndef SPPE_SPATIAL_MAP_H
#define SPPE_SPATIAL_MAP_H

#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "particle.h"

namespace SPPE {
  
  // Sparse spatial hash map.
  // This works in two passes. In the first, index() is called on each particle
  // until the count_array_ is initialized to
  // the number of particles that map to each location in count_array_.
  // In the second pass, particle indices are arranged in a contiguous array 
  // such that "buckets" can then be queried as spans.
  // The map can be reset with reset(); 
  class Spatial_map {
    static constexpr std::size_t map_size = 200 * 200; // a good number for map_size is n particles
    static constexpr std::size_t count_array_sz = map_size + 1;
  public:
    // Construct a spatial map meant for n particles
    Spatial_map(Float_type spacing, std::size_t n = 20'000)
      : buckets_(n, nullptr), spacing_(spacing) { reset(); }

    // Reset the map.
    void reset() { std::fill_n(count_array_, count_array_sz, 0); }
    
    // Return a span of particle pointers that correspond to the grid position {x,y}
    std::span<Particle* const> query(int x, int y) const;

    // "indexing" a particle consists of incrementing the count of the cell that 
    // it maps to.
    void index(const Particle& particle) { ++count_array_[hash(particle)]; }
    
    void build(const std::span<Particle>&);

    int discretize(Float_type x) const { return std::ceil(x / spacing_); }
  
    void set_spacing(Float_type spacing) { spacing_ = spacing; reset(); }
  private:
    static constexpr int primes[] = { 92837111, 689287499, 283923481 };
    std::size_t hash(int x, int y) const 
    { return (x * primes[0] ^ y * primes[1]) % map_size; }
    std::size_t hash(const Particle& p) const 
    { return hash(discretize(p.position()[0]), discretize(p.position()[1])); }

  private:
    int count_array_[count_array_sz];
    std::vector<Particle*> buckets_;
    Float_type spacing_;
  };
    
  inline std::span<Particle* const>
  Spatial_map::query(int x, int y) const
  {
    const auto i = hash(x, y);
    const std::size_t n = count_array_[i + 1] - count_array_[i];
    return { buckets_.begin() + count_array_[i], n };
  }

  inline void
  Spatial_map::build(const std::span<Particle>& particle_span)
  {
    std::partial_sum(count_array_, count_array_ + count_array_sz, count_array_);
    // The number of particles might not change much, but if it does, buckets_
    // needs to be resized.
    const std::size_t n_particles = count_array_[map_size]; 
    buckets_.resize(n_particles);
    // ASSERT: particle_span.size() == n_particles
    for (auto& particle : particle_span)
    {
      const auto index = --count_array_[hash(particle)];
      buckets_[index] = &particle;
    }
  }
} // namespace SPPE

#endif//SPPE_SPATIAL_MAP_H