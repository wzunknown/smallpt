#include "vector.h"
#include "color.h"
#include "sphere.h"
#include "ray.h"
#include <stdio.h>
#include <vector>
#include <random>
#include <math.h>
#include <fstream>



Vec air_absorp = Vec();
const double eps = 1e-7;

double urand() {
    static std::default_random_engine rng(time(0));
    static std::uniform_real_distribution<double> udist(0.0, 1.0);
    return udist(rng);
}



std::vector<Sphere> spheres {
    // Scene: radius, center, emission, color, refl type, trans type, n, absorp
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Color(), Color(.75, .25, .25), SurfaceType::DIFFUSE, SurfaceType::NONE),   // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Color(), Color(.25, .25, .75), SurfaceType::DIFFUSE, SurfaceType::NONE), // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Color(), Color(.75, .75, .75), SurfaceType::DIFFUSE, SurfaceType::NONE),         // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Color(), Color(), SurfaceType::DIFFUSE, SurfaceType::NONE),               // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Color(), Color(.75, .75, .75), SurfaceType::DIFFUSE, SurfaceType::NONE),         // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Color(), Color(.75, .75, .75), SurfaceType::DIFFUSE, SurfaceType::NONE), // Top
    Sphere(16.5, Vec(27, 16.5, 47), Color(), Color(1, 1, 1) * .999, SurfaceType::NORMAL, SurfaceType::NORMAL, -1.0),        // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Color(), Color(1, 1, 1) * .999, SurfaceType::NORMAL, SurfaceType::NORMAL, 1.3, Vec(0.01, 0.09, 0.03)),        // Glas
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Color(3,3,3), Color(), SurfaceType::DIFFUSE, SurfaceType::NONE)     // Lite
};



bool intersect_spheres(const std::vector<Sphere>& spheres, const Ray& ray, double& t, size_t& id) {
    t = std::numeric_limits<double>::infinity();
    id = -1;
    for (size_t i = 0; i < spheres.size(); ++i) {
        const auto& sph = spheres[i];
        double d = sph.intersect(ray);
        if (isnan(d)) {
            continue;
        } else if (d < t) {
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
    cos_i = abs(cos_i);
    Color light = obj.emission;
    // cos_i = abs(cos_i);

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
        double sin_theta = sqrt(sin_theta2);
        Vec eu = ((abs(eh.x) > 0.1 ? Vec(0, 1) : Vec(1)) % eh).normalize();
        Vec ev = eh % eu;
        Vec new_direction = eu * cos(phi) * sin_theta + ev * sin(phi) * sin_theta + eh * sqrt(1 - sin_theta2);
        light += obj.color.mult(radiance(Ray(new_origin, new_direction), depth, absorp));
    }


    // transmission
    if (obj.surf_tran == SurfaceType::NORMAL) {
        double n_rel = is_in ? obj.n_refr : 1.0 / obj.n_refr;
        double sin_t = sqrt(1 - cos_i * cos_i) / n_rel;
        if (sin_t < 1) {
            double T = Fresnel_transmission(n_rel, cos_i);
            Vec new_direction;
            if (sin_t < eps) {
                new_direction = -1 * eh;
            } else {
                Vec eu = (ray.direction + cos_i * eh).normalize();
                double cos_t = sqrt(1 - sin_t * sin_t);
                new_direction = sin_t * eu - cos_t * eh;
            }
            light += T * obj.color.mult(radiance(Ray(new_origin, new_direction), depth, is_in ? obj.absorption : air_absorp));
        } 
    }

    // absorption
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
    int samples = argc == 2? atoi(argv[1]) / (grid[0] * grid[1]) : 1;

    // camera setup
    Ray camera(Vec(50, 52, 295.6), Vec(0, -0.042612, -1));
    Vec cam_x = Vec(0.5135 / height);
    Vec cam_y = (cam_x % camera.direction).normalize() * 0.5135 / height;


    // Monte Carlo process
    Color temp_color;
#pragma omp parallel for schedule(dynamic, 1) private(temp_color)    
    for (int y = 0; y < height; ++y) {
        if (y % (height / 50) == 0) {
            fprintf(stderr, "\rRendering (%d spp): %4.2f%%", samples * grid[0] * grid[1], 100.0 * y / height);
        }
        for (int x = 0; x < width; ++x) {
            temp_color = Color();
            for (int sy = 0; sy < grid[1]; ++sy) {
                for (int sx = 0; sx < grid[0]; ++sx) {
                    for (int samp = 0; samp < samples; ++samp) {
                        // tent filter
                        double rand;
                        rand = urand() * 2;
                        double dx = rand < 1 ? sqrt(rand) - 1 : 1 - sqrt(2 - rand);
                        rand = urand() * 2;
                        double dy = rand < 1 ? sqrt(rand) - 1 : 1 - sqrt(2 - rand);

                        Vec dir = cam_x * (x + (sx + dx) / 2 - width / 2.0)
                                    + cam_y * (y + (sy + dy) / 2 - height / 2.0)
                                    + camera.direction;
                        // dir.show();
                        
                        temp_color += radiance(Ray(camera.origin + 140 * dir, dir), 0, air_absorp) * (1.0 / samples);
                    }
                    // temp_color.show();
                    canvas[x][y] += temp_color * 0.25;
                }
            }
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