#include "sphere.h"
#include <string>



Sphere::Sphere(std::string name_, double radius_, Vec center_, Color emission_, Color color_, SurfaceType refl_, SurfaceType tran_, double n_, Vec absorp_) 
    : name(name_), radius(radius_), center(center_), emission(emission_), color(color_), 
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
    double delta = std::sqrt(deter);
    if (p + delta < eps) {
        return std::nan("");
    }
    return p - delta < eps ? p + delta : p - delta;
}