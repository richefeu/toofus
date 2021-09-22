#include "../Mth.hpp"
#include "../quat.hpp"

#include <iostream>

int main(int argc, char const *argv[]) {
  quat QA;
  quat QB;

  QA.set_axis_angle(vec3r::unit_z(), Mth::pi_2);
  QB.set_axis_angle(vec3r::unit_z(), 3.0 * Mth::pi_4);

  quat QAB = QA.get_conjugated() * QB;

  std::cout << "QA  = " << QA.get_angle() / Mth::pi << " x pi\n";
  std::cout << "QB  = " << QB.get_angle() / Mth::pi << " x pi\n";
  std::cout << "QAB = " << QAB.get_angle() / Mth::pi << " x pi\n\n";
  
  // ----
  
  quat Q;
  Q.set_axis_angle(vec3r::one(), Mth::pi_4);
  std::cout << "Q axis  = " << Q.get_axis() << "\n";
  std::cout << "Q angle = " << Q.get_angle() / Mth::pi << " x pi\n";
  std::cout << "Q Pitch = " << Q.get_Pitch() / Mth::pi << " x pi\n";
  std::cout << "Q Yaw   = " << Q.get_Yaw() / Mth::pi << " x pi\n";
  std::cout << "Q Roll  = " << Q.get_Roll() / Mth::pi << " x pi\n";
  
  return 0;
}