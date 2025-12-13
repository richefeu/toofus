#include "doctest.h"
#include "AABB.hpp"

TEST_CASE("Testing AABB class") {
    SUBCASE("Testing constructors and assignment operator") {
        AABB aabb1;
        AABB aabb2(aabb1);
        AABB aabb3;
        aabb3 = aabb2;

        CHECK(aabb1.min == aabb2.min);
        CHECK(aabb1.max == aabb2.max);

        CHECK(aabb2.min == aabb3.min);
        CHECK(aabb2.max == aabb3.max);
    }

    SUBCASE("Testing enlarge function") {
        AABB aabb;
        aabb.enlarge(1.0);

        CHECK(aabb.min.x == -1.0);
        CHECK(aabb.min.y == -1.0);
        CHECK(aabb.min.z == -1.0);
        CHECK(aabb.max.x == 1.0);
        CHECK(aabb.max.y == 1.0);
        CHECK(aabb.max.z == 1.0);
    }

    SUBCASE("Testing translate function") {
        AABB aabb;
        vec3r v(1.0, 2.0, 3.0);
        aabb.translate(v);

        CHECK(aabb.min.x == 1.0);
        CHECK(aabb.min.y == 2.0);
        CHECK(aabb.min.z == 3.0);
        CHECK(aabb.max.x == 1.0);
        CHECK(aabb.max.y == 2.0);
        CHECK(aabb.max.z == 3.0);
    }

    SUBCASE("Testing intersect function") {
        AABB aabb1, aabb2;
        aabb1.enlarge(1.0);
        aabb2.enlarge(1.0);

        CHECK(aabb1.intersect(aabb2));

        aabb2.translate(vec3r(2.2, 0.0, 0.0));
        CHECK(!aabb1.intersect(aabb2));
    }

    SUBCASE("Testing volume and surface area calculations") {
        AABB aabb(vec3r(0.0, 0.0, 0.0), vec3r(2.0, 2.0, 2.0));
        CHECK(aabb.volume() == 8.0);
    }

    SUBCASE("Testing point containment") {
        AABB aabb(vec3r(0.0, 0.0, 0.0), vec3r(2.0, 2.0, 2.0));
        vec3r pointInside(1.0, 1.0, 1.0);
        vec3r pointOutside(3.0, 3.0, 3.0);

        CHECK(aabb.contains(pointInside));
        CHECK(!aabb.contains(pointOutside));
    }

    SUBCASE("Testing merge function") {
        AABB aabb1(vec3r(0.0, 0.0, 0.0), vec3r(1.0, 1.0, 1.0));
        AABB aabb2(vec3r(1.0, 1.0, 1.0), vec3r(2.0, 2.0, 2.0));
        AABB merged = aabb1;
        merged.merge(aabb2);

        CHECK(merged.min == vec3r(0.0, 0.0, 0.0));
        CHECK(merged.max == vec3r(2.0, 2.0, 2.0));
    }
}