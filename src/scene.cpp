#include "scene.h"
#include <cmath>
#include <stdexcept>
#include "yaml-cpp/yaml.h"

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
    // double cos_i = std::sqrt(1 - sin_i * sin_i);
    double cos_t = std::sqrt(1 - sin_t * sin_t);
    double rs = (cos_i - n * cos_t) / (cos_i + n * cos_t);
    rs *= rs;
    double rp = (cos_t - n * cos_i) / (cos_t + n * cos_i);
    rp *= rp;
    return (rs + rp) / 2;
}


double Fresnel_transmission(double n, double cos_i) {
    return 1 - Fresnel_reflection(n, cos_i);
}



Scene::Scene(std::string yaml_file/* ="default.yaml" */) {
    load_yaml(yaml_file);
}


Scene::Scene(const Ray& cam_, double cam_l_, size_t w_/* =1024 */, size_t h_/* =768 */, double fov_/* =0.5135 */)
    : camera(cam_), camera_length(cam_l_), width(w_), height(h_), fov(fov_)
{

}


bool Scene::intersect_spheres(const Ray& ray, double& t, size_t& id) {
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

Color Scene::radiance(const Ray& ray, int depth, Vec absorp) {
    double t; // distance to intersection
    size_t id;
    if (!intersect_spheres(ray, t, id)) { // no intersect
        return Color(); // return black
    }
    const Sphere& obj = spheres[id];

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
        double R = obj.n_refr > 0 ? Fresnel_reflection(n_rel, cos_i) : -obj.n_refr;
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
    for (int idx = 0; idx < 3; ++idx) {
        light[idx] *= exp(-t * absorp[idx]);
    }
    // light.x *= exp(-t * absorp.x);
    // light.y *= exp(-t * absorp.y);
    // light.z *= exp(-t * absorp.z);

    // light.show();
    return light;
}


int Scene::render_raw(int xmin, int xmax, int ymin, int ymax, std::vector<double>& data_raw) {
    // canvas
    int X = xmax - xmin;
    int Y = ymax - ymin;
    int total_grid = grid[0] * grid[1];
    int samples = samples_per_pixel / total_grid;


    // default camera setup
    Vec cam_x = Vec(fov / height);
    Vec cam_y = (cam_x % camera.direction).normalize() * fov / height;

    // Monte Carlo process
#pragma omp parallel for ordered schedule(dynamic, 1)   
    for (int y = ymin; y < ymax; ++y) {
        for (int x = xmin; x < xmax; ++x) {
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
                        
                        temp_color += radiance(Ray(camera.origin + camera_length * dir, dir), 0, air_absorp) * (1.0 / samples_per_pixel);
                    }
                    // temp_color.show();
                }
            }
            int idx = (y - ymin) * X * 3 + (x - xmin) * 3;
            data_raw[idx] = temp_color.x;
            data_raw[idx + 1] = temp_color.y;
            data_raw[idx + 2] = temp_color.z;
        }
    }
    return 0;    
}


int Scene::render() {
    // canvas
    canvas = std::vector<std::vector<Color>>(width, std::vector<Color>(height, Color()));
    int total_grid = grid[0] * grid[1];
    int samples = samples_per_pixel / total_grid;

    // default camera setup
    Vec cam_x = Vec(fov / height);
    Vec cam_y = (cam_x % camera.direction).normalize() * fov / height;

    // Monte Carlo process
#pragma omp parallel for ordered schedule(dynamic, 1)   
    for (int y = 0; y < height; ++y) {
        // if (y % (height / 50) == 0) {
        fprintf(stderr, "\rRendering (%d spp): %4.2f%%", samples_per_pixel, 100.0 * y / height);
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
                        
                        temp_color += radiance(Ray(camera.origin + camera_length * dir, dir), 0, air_absorp) * (1.0 / samples_per_pixel);
                    }
                    // temp_color.show();
                }
            }
            canvas[x][y] = temp_color;
            // for (int idx = 0; idx < 3; ++idx) {
            //     canvas[x][y][idx] *= exp(-camera_length * air_absorp[idx]);
            // }
        }
    }
    return 0;    
}


// load Vec and Color
namespace YAML {
    template<>
    struct convert<Vec> {
        static Node encode(const Vec& rhs) {
            Node node;
            node.push_back(rhs.x);
            node.push_back(rhs.y);
            node.push_back(rhs.z);
            return node;
        }

        static bool decode(const Node& node, Vec& rhs) {
            if (!node.IsSequence() || node.size() != 3) {
                return false;
            }

            rhs.x = node[0].as<double>();
            rhs.y = node[1].as<double>();
            rhs.z = node[2].as<double>();
            return true;
        }
    };


    template<>
    struct convert<Color> {
        static Node encode(const Color& rhs) {
            Node node;
            node.push_back(rhs.x);
            node.push_back(rhs.y);
            node.push_back(rhs.z);
            return node;
        }

        static bool decode(const Node& node, Color& rhs) {
            if (!node.IsSequence() || node.size() != 3) {
                return false;
            }

            rhs.x = node[0].as<double>();
            rhs.y = node[1].as<double>();
            rhs.z = node[2].as<double>();
            return true;
        }
    };


    template<>
    struct convert<SurfaceType> {
        static Node encode(const SurfaceType& rhs) {
            Node node("unknown");
            if (rhs == SurfaceType::NONE) node = Node("none");
            else if (rhs == SurfaceType::DIFFUSE) node = Node("diffuse");
            else if (rhs == SurfaceType::NORMAL) node =Node("normal");                    
            
            return node;
        }

        static bool decode(const Node& node, SurfaceType& rhs) {
            if (!node.IsScalar()) {
                return false;
            }

            std::string surf_type = node.as<std::string>();
            
            if (surf_type == "none") rhs = SurfaceType::NONE;
            else if (surf_type == "diffuse") rhs = SurfaceType::DIFFUSE;
            else if (surf_type == "normal") rhs = SurfaceType::NORMAL;
            
            return true;
        }
    };

    template<>
    struct convert<Sphere> {
        static Node encode(const Sphere& rhs) {
            throw std::runtime_error("not implemented");
            return Node();
        }

        static bool decode(const Node& node, Sphere& rhs) {
            if (!node.IsMap()) {
                return false;
            }

            rhs.name = node["name"] ? node["name"].as<std::string>() : "";
            rhs.radius = node["radius"].as<double>();
            rhs.center = node["center"].as<Vec>();
            rhs.emission = node["emission"] ? node["emission"].as<Color>() : Color();
            rhs.color = node["color"] ? node["color"].as<Color>() : Color();
            rhs.surf_refl = node["reflection"] ? node["reflection"].as<SurfaceType>() : SurfaceType::NONE;
            rhs.surf_tran = node["transmission"] ? node["transmission"].as<SurfaceType>() : SurfaceType::NONE;
            rhs.n_refr = node["refractive_index"] ? node["refractive_index"].as<double>() : 1.5;
            rhs.absorption = node["absorption"] ? node["absorption"].as<Vec>() : Vec();

            return true;
        }
    };
}


void Scene::load_yaml(std::string yaml_file/* ="config.yaml" */) {
    YAML::Node config = YAML::LoadFile(yaml_file);

    // grid
    if (config["grid"]) {
        grid[0] = config["grid"][0].as<size_t>();
        grid[1] = config["grid"][1].as<size_t>();
    }
    
    // samples
    samples_per_pixel = config["samples_per_pixel"] ? config["samples_per_pixel"].as<size_t>() : 1;

    // canvas
    width = config["width"] ? config["width"].as<size_t>() : 1024;
    height = config["height"] ? config["height"].as<size_t>() : 768;

    // camera
    if (config["camera"]) {
        const YAML::Node camera_config = config["camera"];
        camera.origin = camera_config["position"].as<Vec>();
        camera.direction = camera_config["direction"].as<Vec>();
        camera_length = camera_config["length"] ? camera_config["length"].as<double>() : 140.;
        fov = camera_config["fov"] ? camera_config["fov"].as<double>() : 0.5135;
    }

    // spheres
    if (config["spheres"]) {
        for (const auto& sph_config : config["spheres"]) {
            add_spheres(sph_config.as<Sphere>());
        }
    }

    // air_absorp
    air_absorp = config["air_absorption"].as<Vec>();

}


void Scene::save_ppm(std::string img_file/* ="image.ppm" */) {
    FILE* fp = fopen(img_file.c_str(), "w");
    fprintf(fp, "P3\n%d %d\n%d\n", width, height, 255);
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            Vec srgb = canvas[x][y].to_sRGB();
            fprintf(fp, "%d %d %d ", Color::to_8bit(srgb.x), Color::to_8bit(srgb.y), Color::to_8bit(srgb.z));
        }
    }
    fclose(fp);
}


void Scene::save_png(std::string img_file/* ="image.png" */) {
    throw std::runtime_error("not implemented, maybe use libpng or opencv");
}