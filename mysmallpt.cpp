#include "scene.h"

int main(int argc, char* argv[]) {
    Scene sc = argc == 2 ? Scene(argv[1]) : Scene();
    sc.render();
    sc.save_ppm();
    return 0;
}