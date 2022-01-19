#pragma once
#include "vector.h"
#include "color.h"
#include "ray.h"


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


    Sphere(double radius_, Vec center_, Color emission_, Color color_, SurfaceType refl_, SurfaceType tran_, double n_=0, Vec absorp_=Vec());
    double intersect(const Ray& ray) const;
};


/**
 * @brief 
 * 
 * @param n relative refractive index
 * @param sin_i angle of incidence
 * @return double 
 */
double Fresnel_reflection(double n, double cos_i) {
    double sin_i = sqrt(1 - cos_i * cos_i);
    double sin_t = sin_i / n; // angle of refraction
    if (sin_t > 1 || n < 0) { // total reflection or specular reflection
        return 1; 
    }
    // double cos_i = sqrt(1 - sin_i * sin_i);
    double cos_t = sqrt(1 - sin_t * sin_t);
    double rs = (cos_i - n * cos_t) / (cos_i + n * cos_t);
    rs *= rs;
    double rp = (cos_t - n * cos_i) / (cos_t + n * cos_i);
    rp *= rp;
    return (rs + rp) / 2;
}


// double Fresnel_reflection(double n, double sin_i) {
//     double sin_t = sin_i / n; // angle of refraction
//     if (sin_t > 1 || n < 0) { // total reflection or specular reflection
//         return 1; 
//     }
//     double cos_i = sqrt(1 - sin_i * sin_i);
//     double cos_t = sqrt(1 - sin_t * sin_t);
//     double rs = (cos_i - n * cos_t) / (cos_i + n * cos_t);
//     rs *= rs;
//     double rp = (cos_t - n * cos_i) / (cos_t + n * cos_i);
//     rp *= rp;
//     return (rs + rp) / 2;
// }


double Fresnel_transmission(double n, double cos_i) {
    return 1 - Fresnel_reflection(n, cos_i);
}

// double Fresnel_transmission(double n, double sin_i) {
//     return 1 - Fresnel_reflection(n, sin_i);
// }


Sphere::Sphere(double radius_, Vec center_, Color emission_, Color color_, SurfaceType refl_, SurfaceType tran_, double n_, Vec absorp_) 
    : radius(radius_), center(center_), emission(emission_), color(color_), 
      surf_refl(refl_), surf_tran(tran_), n_refr(n_), absorption(absorp_)
{
    // if (surf_refl == SurfaceType::DIFFUSE && surf_tran == SurfaceType::DIFFUSE) {
    //     R_diffuse = 0.2;
    // }
}



double Sphere::intersect(const Ray& ray) const {
    // intersect point: ray.origin + t * ray.direction
    Vec oc = center - ray.origin;
    double p = ray.direction.dot(oc);
    double deter = p * p - oc.dot(oc) + radius * radius;
    if (deter < 0) {
        return std::nan("");
    }
    double delta = sqrt(deter);
    if (p + delta < eps) {
        return std::nan("");
    }
    return p - delta < eps ? p + delta : p - delta;
}