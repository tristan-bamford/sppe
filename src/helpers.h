#pragma once

#include "sppe.h"



auto random_number(SPPE::Float_type from = -1, SPPE::Float_type to = 1)
{
  return from + rand() * ((to - from) / RAND_MAX);
}

auto random_vector(SPPE::Vector_type from = -1, SPPE::Vector_type to = 1)
{
  SPPE::Vector_type result;
  for (std::size_t i = 0; i < result.size(); ++i) {
    result[i] = random_number(from[i], to[i]);
  }
  return result;
}

void random_particles(SPPE::System& sppe, int n)
{
  for (auto i = 0; i < n; ++i) {
    const auto position = random_vector(-sppe.get_boundary(), sppe.get_boundary());
    const auto velocity = random_vector(-100, 100);
    const auto radius = random_number(2,2);
    const auto density = 1;
    const auto elasticity = 1;//random_number(0.9,1.0);
    sppe.add_particle({ position, radius, velocity, density, elasticity }); 
  }
}

auto momentum(const SPPE::System& sppe)
{
  SPPE::Vector_type result = 0;
  for (const auto& p : sppe.particles()) result += p.momentum();
  return result;
}

auto kinetic_energy(const SPPE::System& sppe)
{
  SPPE::Float_type result = 0;
  for (const auto& p : sppe.particles()) result += p.kinetic_energy();
  return result;
}