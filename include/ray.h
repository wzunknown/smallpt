#pragma once
#include "vector.h"


class Ray {
public:
    Vec origin;
    Vec direction;

    Ray() {}

    Ray(Vec o_, Vec d_) : origin(o_), direction(d_) {
        direction.normalize();
    }
};