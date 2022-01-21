#include "scene.h"
#include <vector>
#include <fstream>
#include <algorithm>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include <stdarg.h>
#include <string>

#define FROM_MASTER 1
#define FROM_WORKER 2

int rank, numprocs;
char processor_name[MPI_MAX_PROCESSOR_NAME];

void log_info(const char* format,...) {
    char time_buf[256];
    time_t curr_time = time(0);
    strftime(time_buf, sizeof(time_buf), "%d-%H:%M:%S", localtime(&curr_time));
    printf("%s [%15.15s:%2d/%2d] ", time_buf, processor_name, rank, numprocs);
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stdout, format, argptr);
    va_end(argptr);
    printf("\n"); 
}

int main(int argc, char* argv[]) {
    int namelen;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processor_name, &namelen);

    Scene sc = Scene(argv[1]);
    int master_id = 0;
    int ymin, ymax;
    
    if (rank == master_id) {
        // master node assign jobs to nodes
        int height_per_node = (sc.height + numprocs - 1) / numprocs;
        for (int id = 1, ymin = height_per_node; id < numprocs && ymin < sc.height; ++id, ymin += height_per_node) {
            int ymax = std::min(ymin + height_per_node, sc.height);
            MPI_Send(&ymin, 1, MPI_INT, id, FROM_MASTER, MPI_COMM_WORLD);
            MPI_Send(&ymax, 1, MPI_INT, id, FROM_MASTER, MPI_COMM_WORLD);
        }

        ymin = 0;
        ymax = height_per_node;
    } else {
        MPI_Recv(&ymin, 1, MPI_INT, master_id, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&ymax, 1, MPI_INT, master_id, FROM_MASTER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    log_info("Get job from master node, ymin: %d, ymax: %d", ymin, ymax);

    std::vector<double> data_raw((ymax - ymin) * sc.width * 3, 0.0);
    sc.render_raw(0, sc.width, ymin, ymax, data_raw);
    log_info("Finish job ymin: %d, ymax: %d", ymin, ymax);
    // collect result
    if (rank == master_id) {
        // fprintf(stderr, "")
        data_raw.resize(sc.height * sc.width * 3);
        for (int id = 1; id < numprocs; ++id) {
            MPI_Recv(&ymin, 1, MPI_INT, id, FROM_WORKER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&ymax, 1, MPI_INT, id, FROM_WORKER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&data_raw[ymin * sc.width * 3], (ymax - ymin) * sc.width * 3, MPI_DOUBLE, id, FROM_WORKER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // save ppm
        FILE* fp = fopen("image.ppm", "w");
        fprintf(fp, "P3\n%d %d\n%d\n", sc.width, sc.height, 255);
        for (int y = sc.height - 1; y >= 0; --y) {
            for (int x = 0; x < sc.width; ++x) {
                int idx = y * sc.width * 3 + x * 3;
                Vec srgb = Color(data_raw[idx], data_raw[idx+1], data_raw[idx+2]).to_sRGB();
                fprintf(fp, "%d %d %d ", Color::to_8bit(srgb.x), Color::to_8bit(srgb.y), Color::to_8bit(srgb.z));
            }
        }
        fclose(fp);

    } else {
        MPI_Send(&ymin, 1, MPI_INT, master_id, FROM_WORKER, MPI_COMM_WORLD);
        MPI_Send(&ymax, 1, MPI_INT, master_id, FROM_WORKER, MPI_COMM_WORLD);
        MPI_Send(&data_raw[0], data_raw.size(), MPI_DOUBLE, master_id, FROM_WORKER, MPI_COMM_WORLD);
    }

    MPI_Finalize();  
    return 0;  
}
