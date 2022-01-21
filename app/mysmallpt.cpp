#include "scene.h"
#include <cstdio>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::fprintf(stderr, "[*] usage: %s config.yaml [spp]", argv[0]);
        return 0;
    }
    Scene sc(argv[1]);
    if (argc == 3) {
        sc.samples_per_pixel = atoi(argv[2]);
    }
    sc.render();
    sc.save_ppm();
    return 0;
}