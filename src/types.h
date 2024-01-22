#ifndef SPPE_TYPES_H
#define SPPE_TYPES_H

#include "num_array.h"

// Enable three-dimensional SPPE
//#define SPPE_3D;

namespace SPPE {

  using Float_type = double;

#ifdef SPPE_3D
  using Vector_type = tb::math::num_array<Float_type, 3>;
#else
  using Vector_type = tb::math::num_array<Float_type, 2>;
#endif

  constexpr Float_type PI = 3.14159265358979;
  constexpr Float_type SQRT_2 = 1.41421356237;

}
#endif//SPPE_TYPES_H
