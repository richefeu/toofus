#ifndef PERIODICTOOLS2D_HPP
#define PERIODICTOOLS2D_HPP

#include "mat4.hpp"

namespace periodictools2D {

/**
 * @brief Calculate the Green-Lagrange strain tensor.
 *
 * @param h0 Initial configuration (undeformed or reference state)
 * @param h Final configuration (deformed state)
 * @param E Strain tensor to be filled
 *
 * Formula: E = 0.5 * (F^T * F - I) = 0.5 * (h * h0^{-1})^T * (h * h0^{-1}) - I
 *
 * This function modifies @p E.
 */
void GreenLagrangeStrain(mat4r &h0, mat4r &h, mat4r &E) {
  mat4r F = h * h0.get_inverse();
  return 0.5 * (F.transposed() * F - mat4r::unit());
}

/**
 * @brief Adjusts the velocity vector based on periodic movement. ix and iy can be -1.0, 0.0, or 1.0
 *
 * @param vh Transformation matrix representing the velocity field.
 * @param ix Periodic displacement along the x-axis.
 * @param iy Periodic displacement along the y-axis.
 * @return Corrected velocity vector.
 */
vec2r correctVelocity(mat4r &vh, double ix, double iy) {
  return vec2r(vh.xx * ix + vh.xy * iy, vh.yx * ix + vh.yy * iy);
}

/**
 * @brief Adjusts the velocity vector based on periodic movement.
 *
 * @param vh Transformation matrix representing the velocity field.
 * @param numPerioMove Vector of periodic displacement along the x and y axes.
 * @return Corrected velocity vector.
 */
vec2r correctVelocity(mat4r &vh, vec2r &numPerioMove) {
  return vh * numPerioMove;
}

} // namespace periodictools2D

#endif /* end of include guard: PERIODICTOOLS_HPP */
