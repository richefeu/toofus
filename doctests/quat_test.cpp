#include "doctest.h"
#include "quat.hpp"

TEST_CASE("quat") {
  quat q;
  REQUIRE(q.v.x == 0.0);
  REQUIRE(q.s == 1.0);
}

#include <cmath>

TEST_CASE("quat default and identity") {
  quat q;
  CHECK(q.v == vec3r(0, 0, 0));
  CHECK(q.s == 1.0);
  CHECK(quat::identity() == q);

  // multiplying by the identity leaves a quaternion unchanged
  quat a(0.1, 0.2, 0.3, 0.4);
  CHECK(a * quat::identity() == a);
  CHECK(quat::identity() * a == a);
}

TEST_CASE("quat set_axis_angle rotates a vector") {
  quat q;
  q.set_axis_angle(vec3r(0, 0, 1), M_PI / 2.0); // 90 deg about z

  // both rotation paths agree and map x -> y
  vec3r byOperator = q * vec3r(1, 0, 0);
  CHECK(byOperator.x == doctest::Approx(0.0));
  CHECK(byOperator.y == doctest::Approx(1.0));
  CHECK(byOperator.z == doctest::Approx(0.0));

  vec3r byRotate = q.rotate(vec3r(1, 0, 0));
  CHECK(byRotate.x == doctest::Approx(0.0));
  CHECK(byRotate.y == doctest::Approx(1.0));
  CHECK(byRotate.z == doctest::Approx(0.0));

  CHECK(q.get_angle() == doctest::Approx(M_PI / 2.0));
  vec3r axis = q.get_axis();
  CHECK(axis.x == doctest::Approx(0.0));
  CHECK(axis.y == doctest::Approx(0.0));
  CHECK(axis.z == doctest::Approx(1.0));
}

TEST_CASE("quat rotate then unrotate is the identity") {
  quat q;
  q.set_axis_angle(vec3r(1, 1, 0), 0.7);
  q.normalize(); // rotate/unrotate round-trip requires a unit quaternion
  vec3r v(2, -1, 3);
  vec3r back = q.unrotate(q.rotate(v));
  CHECK(back.x == doctest::Approx(v.x));
  CHECK(back.y == doctest::Approx(v.y));
  CHECK(back.z == doctest::Approx(v.z));
}

TEST_CASE("quat conjugate and normalize") {
  quat q;
  q.set_axis_angle(vec3r(0, 0, 1), M_PI / 2.0); // already a unit quaternion

  // a unit quaternion times its conjugate gives the identity
  quat prod = q * q.get_conjugated();
  CHECK(prod.v.x == doctest::Approx(0.0));
  CHECK(prod.v.y == doctest::Approx(0.0));
  CHECK(prod.v.z == doctest::Approx(0.0));
  CHECK(prod.s == doctest::Approx(1.0));

  // normalize() scales to unit length and returns the previous norm
  quat r;
  r.set(0.0, 0.0, 0.0, 2.0);
  double n = r.normalize();
  CHECK(n == doctest::Approx(2.0));
  CHECK(r.s == doctest::Approx(1.0));
}