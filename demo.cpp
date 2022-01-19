#include "vector.h"
#include "color.h"
#include "sphere.h"
#include "ray.h"
#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include <fstream>
#include <omp.h>


const double eps = 1e-7;

// global configurations
// Vec air_absorp = Vec(0.0001, 0.0001, 0.0001);
Vec air_absorp = Vec();

// std::vector<Sphere> spheres;

double urand() {
    // static std::default_random_engine rng(time(0));
    // static std::uniform_real_distribution<double> udist(0.0, 1.0);
    static thread_local std::default_random_engine rng(time(0));
    static thread_local std::uniform_real_distribution<double> udist(0.0, 1.0);
    return udist(rng);
}




// sRGB
// std::vector<Sphere> spheres {
//     // Scene: radius, center, emission, color, refl type, trans type, n, absorp
//     Sphere("left", 1e5, Vec(1e5 + 1, 40.8, 81.6), Color(), Color(.75, .25, .25), SurfaceType::DIFFUSE, SurfaceType::NONE),   // Left
//     Sphere("right", 1e5, Vec(-1e5 + 99, 40.8, 81.6), Color(), Color(.25, .25, .75), SurfaceType::DIFFUSE, SurfaceType::NONE), // Rght
//     Sphere("back", 1e5, Vec(50, 40.8, 1e5), Color(), Color(.75, .75, .75), SurfaceType::DIFFUSE, SurfaceType::NONE),         // Back
//     Sphere("front", 1e5, Vec(50, 40.8, -1e5 + 170), Color(), Color(), SurfaceType::DIFFUSE, SurfaceType::NONE),               // Frnt
//     Sphere("bottom", 1e5, Vec(50, 1e5, 81.6), Color(), Color(.75, .75, .75), SurfaceType::DIFFUSE, SurfaceType::NONE),         // Botm
//     Sphere("top",  1e5, Vec(50, -1e5 + 81.6, 81.6), Color(), Color(.75, .75, .75), SurfaceType::DIFFUSE, SurfaceType::NONE), // Top
//     Sphere("mirror", 16.5, Vec(27, 16.5, 47), Color(), Color(1, 1, 1) * .999, SurfaceType::NORMAL, SurfaceType::NORMAL, -1.0),        // Mirr
//     Sphere("glass", 16.5, Vec(73, 16.5, 78), Color(), Color(1, 1, 1) * .999, SurfaceType::NORMAL, SurfaceType::NORMAL, 1.3, Vec(0.01, 0.09, 0.03)),        // Glas
//     Sphere("light", 600, Vec(50, 681.6 - .27, 81.6), Color(3,3,3), Color(), SurfaceType::DIFFUSE, SurfaceType::NONE)     // Lite
// };

// CIE XYZ
std::vector<Sphere> spheres {
    // Scene: name, radius, center, emission, color, refl type, trans type, n, absorp
    Sphere("left", 1e5, Vec(1e5 + 1, 40.8, 81.6), Color(), Color(0.623, 0.521, 0.438), SurfaceType::DIFFUSE, SurfaceType::NONE),   // Left
    Sphere("right", 1e5, Vec(-1e5 + 99, 40.8, 81.6), Color(), Color(0.572, 0.616, 1.035), SurfaceType::DIFFUSE, SurfaceType::NONE), // Rght
    Sphere("back", 1e5, Vec(50, 40.8, 1e5), Color(), Color(0.543, 0.571, 0.622), SurfaceType::DIFFUSE, SurfaceType::NONE),         // Back
    Sphere("front", 1e5, Vec(50, 40.8, -1e5 + 170), Color(), Color(0.543, 0.571, 0.622), SurfaceType::DIFFUSE, SurfaceType::NONE),               // Frnt
    Sphere("bottom", 1e5, Vec(50, 1e5, 81.6), Color(), Color(0.543, 0.571, 0.622), SurfaceType::DIFFUSE, SurfaceType::NONE),         // Botm
    Sphere("top", 1e5, Vec(50, -1e5 + 81.6, 81.6), Color(), Color(0.543, 0.571, 0.622), SurfaceType::DIFFUSE, SurfaceType::NONE), // Top
    Sphere("mirror", 16.5, Vec(27, 16.5, 47), Color(), Color(1, 1, 1) * .999, SurfaceType::NORMAL, SurfaceType::NORMAL, -0.7),        // Mirr
    Sphere("glass", 16.5, Vec(73, 16.5, 78), Color(), Color(1, 1, 1) * .999, SurfaceType::NORMAL, SurfaceType::NORMAL, 1.3, Vec(0.002, 0.01, 0.003)),        // Glas
    Sphere("light", 600, Vec(50, 681.6 - .27, 81.6), Color(4, 4, 4), Color(), SurfaceType::DIFFUSE, SurfaceType::NONE)     // Lite
};


/**
 * @brief 
 * 
 * @param n relative refractive index
 * @param sin_i angle of incidence
 * @return double 
 */
double Fresnel_reflection(double n, double cos_i) {
    double sin_i = std::sqrt(1 - cos_i * cos_i);
    double sin_t = sin_i / n; // angle of refraction
    if (sin_t > 1) { // total reflection or specular reflection
        return 1; 
    }
    if (n < 0) {
        return -n;
    }
    // double cos_i = std::sqrt(1 - sin_i * sin_i);
    double cos_t = std::sqrt(1 - sin_t * sin_t);
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
//     double cos_i = std::sqrt(1 - sin_i * sin_i);
//     double cos_t = std::sqrt(1 - sin_t * sin_t);
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



bool intersect_spheres(const std::vector<Sphere>& spheres, const Ray& ray, double& t, size_t& id) {
    t = std::numeric_limits<double>::infinity();
    id = -1;
    for (size_t i = 0; i < spheres.size(); ++i) {
        const auto& sph = spheres[i];
        double d = sph.intersect(ray);
        if (d < t) {
            id = i;
            t = d;
        }
    }
    return id != -1;
}

// TODO: ball in ball

Color radiance(const Ray& ray, int depth, Vec absorp) {
    double t; // distance to intersection
    size_t id;
    if (!intersect_spheres(spheres, ray, t, id)) { // no intersect
        return Color(); // return black
    }
    const Sphere& obj = spheres[id];

    // fprintf(stderr, "depth: %d, id: %zu, t: %f\n", depth, id, t);
    // ray.origin.show();
    // ray.direction.show();

    // if (obj.emission.norm() > eps) {
    //     fprintf(stderr, "\rfind");
    // }

    if (++depth > 5) {
        // if (urand() > 0.5)
        return obj.emission;
    }

    // normal vector
    Vec new_origin = ray.origin + t * ray.direction;
    Vec normal_vec = (new_origin - obj.center).normalize();
    double cos_i = normal_vec.dot(ray.direction);
    bool is_in = cos_i < 0;
    Vec eh = is_in ? normal_vec : -1 * normal_vec;
    cos_i = std::abs(cos_i);
    Color light = obj.emission;
    // cos_i = std::abs(cos_i);

    // reflection
    if (obj.surf_refl == SurfaceType::NORMAL) {
        double n_rel = is_in ? obj.n_refr : 1.0 / obj.n_refr;
        double R = Fresnel_reflection(n_rel, cos_i);
        Vec new_direction = ray.direction + 2 * cos_i * eh;
        light += R * obj.color.mult(radiance(Ray(new_origin, new_direction), depth, absorp));
    } else if (obj.surf_refl == SurfaceType::DIFFUSE) {
        // random output
        double phi = 2 * M_PI * urand();
        double sin_theta2 = urand();
        double sin_theta = std::sqrt(sin_theta2);
        Vec eu = ((std::abs(eh.x) > 0.1 ? Vec(0, 1) : Vec(1)) % eh).normalize();
        Vec ev = eh % eu;
        Vec new_direction = eu * std::cos(phi) * sin_theta + ev * std::sin(phi) * sin_theta + eh * std::sqrt(1 - sin_theta2);
        light += obj.color.mult(radiance(Ray(new_origin, new_direction), depth, absorp));
    }


    // transmission
    if (obj.surf_tran == SurfaceType::NORMAL && obj.n_refr > 0) {
        double n_rel = is_in ? obj.n_refr : 1.0 / obj.n_refr;
        double sin_t = std::sqrt(1 - cos_i * cos_i) / n_rel;
        if (sin_t < 1) {
            double T = Fresnel_transmission(n_rel, cos_i);
            Vec new_direction;
            if (sin_t < eps) {
                new_direction = -1 * eh;
            } else {
                Vec eu = (ray.direction + cos_i * eh).normalize();
                double cos_t = std::sqrt(1 - sin_t * sin_t);
                new_direction = sin_t * eu - cos_t * eh;
            }
            light += T * obj.color.mult(radiance(Ray(new_origin, new_direction), depth, is_in ? obj.absorption : air_absorp));
        } 
    }

    // absorption
    // for (int idx = 0; idx < 3; ++idx) {
    //     light[idx] *= exp(-t * absorp[idx]);
    // }
    light.x *= exp(-t * absorp.x);
    light.y *= exp(-t * absorp.y);
    light.z *= exp(-t * absorp.z);

    // light.show();
    return light;
}


int main(int argc, char* argv[]) {
    // canvas
    int width = 1024;
    int height = 768;
    std::vector<std::vector<Color>> canvas(width, std::vector<Color>(height, Color()));

    // samples
    std::vector<int> grid {2, 2};
    int total_grid = grid[0] * grid[1];
    int samples = argc == 2? atoi(argv[1]) / total_grid : 1;

    // default camera setup
    Ray camera(Vec(50, 15, 280.6), Vec(0, 0.15, -1));
    double camera_length = 140.;
    Vec cam_x = Vec(0.5135 / height);
    Vec cam_y = (cam_x % camera.direction).normalize() * 0.5135 / height;


    // TODO: load yaml config

    // Monte Carlo process
    // Color temp_color;
#pragma omp parallel for schedule(dynamic, 1)   
    for (int y = 0; y < height; ++y) {
        // if (y % (height / 50) == 0) {
        fprintf(stderr, "\rRendering (%d spp): %4.2f%%", samples * total_grid, 100.0 * y / height);
        // }
        for (int x = 0; x < width; ++x) {
            Color temp_color = Color();
            for (int sy = 0; sy < grid[1]; ++sy) {
                for (int sx = 0; sx < grid[0]; ++sx) {
                    for (int samp = 0; samp < samples; ++samp) {
                        // tent filter
                        double rand;
                        rand = urand() * 2;
                        double dx = rand < 1 ? std::sqrt(rand) - 1 : 1 - std::sqrt(2 - rand);
                        rand = urand() * 2;
                        double dy = rand < 1 ? std::sqrt(rand) - 1 : 1 - std::sqrt(2 - rand);

                        Vec dir = cam_x * (x + (sx + dx  - grid[0] / 2.) / 2. - width / 2.)
                                    + cam_y * (y + (sy + dy - grid[1] / 2.) / 2. - height / 2.)
                                    + camera.direction;
                        // dir.show();
                        
                        temp_color += radiance(Ray(camera.origin + 140 * dir, dir), 0, air_absorp) * (1.0 / samples);
                    }
                    // temp_color.show();
                    canvas[x][y] += temp_color / total_grid;
                }
            }
            // for (int idx = 0; idx < 3; ++idx) {
            //     canvas[x][y][idx] *= exp(-camera_length * air_absorp[idx]);
            // }
        }
    }

    // write to file
    FILE* fp = fopen("image.ppm", "w");
    fprintf(fp, "P3\n%d %d\n%d\n", width, height, 255);
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            Vec srgb = canvas[x][y].to_sRGB();
            fprintf(fp, "%d %d %d ", Color::to_8bit(srgb.x), Color::to_8bit(srgb.y), Color::to_8bit(srgb.z));
        }
    }
    fclose(fp);
    return 0;
}