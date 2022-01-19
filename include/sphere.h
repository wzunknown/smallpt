#pragma once
#include "vector.h"
#include "color.h"
#include "ray.h"
#include <string>


// enum ReflectionType {
//     NONE, // no reflection
//     DIFFUSE, // diffuse reflection
//     NORMAL, // use refractive index
// };


// enum TransmissionType {
//     NONE, // no transimission
//     DIFFUSE, // diffuse transimission
//     NORMAL, // use refractive index
// }

enum class SurfaceType {
    NONE, // no transmission or reflection
    DIFFUSE, // diffuse
    NORMAL, // use refractive index
};


class Sphere {
public:
    std::string name;
    double radius;
    Vec center; // center of the ball
    Color emission; // emission color
    Color color; // reflect color

    // double R_diffuse = 0.9;

    SurfaceType surf_refl;
    SurfaceType surf_tran;

    double n_refr; // refractive index, (<0 means specular reflection)

    Vec absorption;
    // Vec reflection;
    // Vec transmission;


    static constexpr double eps = 1e-7;

    Sphere() {
        name = "?????";
    }


    Sphere(std::string name_, double radius_, Vec center_, Color emission_, Color color_, SurfaceType refl_, SurfaceType tran_, double n_=0, Vec absorp_=Vec());
    double intersect(const Ray& ray) const;
};


