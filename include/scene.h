#pragma once
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include "vector.h"
#include "sphere.h"
#include "ray.h"
#include "color.h"


class Scene {
private:
    
    double urand() {
        static thread_local std::default_random_engine rng(time(0));
        static thread_local std::uniform_real_distribution<double> udist(0.0, 1.0);
        return udist(rng);
    }


    Color radiance(const Ray& ray, int depth, Vec absorp);
    bool intersect_spheres(const Ray& ray, double& t, size_t& id);

public:
    static constexpr double eps = 1e-7;

    size_t width = 1024;
    size_t height = 768;


    std::vector<Sphere> spheres;
    std::vector<std::vector<Color>> canvas;

    std::string image_name = "image";

    Vec air_absorp = Vec();

    Ray camera = Ray(Vec(50, 15, 280.6), Vec(0, 0.15, -1));
    double camera_length = 140.;
    double fov = 0.5135;

    size_t samples_per_pixel=1;

    std::vector<size_t> grid {2, 2};

    // load default.yaml
    Scene(std::string yaml_file="default.yaml");

    // create new Scene
    Scene(const Ray& cam_, double cam_l_, size_t w_=1024, size_t h_=768, double fov_=0.5135);

    int render();
    int render_raw(int xmin, int xmax, int ymin, int ymax, std::vector<double>& data_raw);

    void load_yaml(std::string yaml_file="config.yaml");
    void save_ppm(std::string img_file="image.ppm");
    void save_png(std::string img_file="image.png");
    
    void add_spheres(const Sphere& sph) {
        spheres.push_back(sph);
    }

    void set_grid(int sx, int sy) {
        grid[0] = sx;
        grid[1] = sy;
    }
};