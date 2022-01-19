#include <iostream>
#include <yaml-cpp/yaml.h>


int main() {
    YAML::Node x = YAML::LoadFile("test.yaml");
    YAML::Node y = YAML::Load("beta: 1e5 + 50");
    y["alpha"] = "beta";
    std::cout << x.as<std::string>() << std::endl;
    // std::cout << y.as<std::string>() << std::endl;
    printf("%d, %f", y.IsMap(), y["beta"].as<double>());
    return 0;
}