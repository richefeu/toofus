#include "doctest.h"
#include "OBB.hpp"

TEST_CASE("Testing OBB class") {
    SUBCASE("Testing constructors and assignment operator") {
        OBB obb1;
        OBB obb2(obb1);
        OBB obb3;
        obb3 = obb2;

        CHECK(obb1.center == obb2.center);
        CHECK(obb1.extent == obb2.extent);
        for(int i = 0; i < 3; ++i) {
            CHECK(obb1.e[i] == obb2.e[i]);
        }

        CHECK(obb2.center == obb3.center);
        CHECK(obb2.extent == obb3.extent);
        for(int i = 0; i < 3; ++i) {
            CHECK(obb2.e[i] == obb3.e[i]);
        }
    }

    SUBCASE("Testing enlarge function") {
        OBB obb;
        obb.enlarge(1.0);

        CHECK(obb.extent.x == 1.0);
        CHECK(obb.extent.y == 1.0);
        CHECK(obb.extent.z == 1.0);
    }

    SUBCASE("Testing translate function") {
        OBB obb;
        vec3r v(1.0, 2.0, 3.0);
        obb.translate(v);

        CHECK(obb.center.x == 1.0);
        CHECK(obb.center.y == 2.0);
        CHECK(obb.center.z == 3.0);
    }

    SUBCASE("Testing intersect function") {
        // This test assumes that you have a working intersect function
        OBB obb1, obb2;
        obb1.enlarge(1.0);
        obb2.enlarge(1.0);

        CHECK(obb1.intersect(obb2));

        obb2.translate(vec3r(2.2, 0.0, 0.0));
        CHECK(!obb1.intersect(obb2));
    }
}
