#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#include "heat_eq.h"
#include "grid.h"

// General information about the code.
// Since arrays are row-major in C we will follow this convention.
// row 0, column 0 is in the lower left corner and increases to upwards and to the right respectively.
int main(int argc, char *argv[])
{
	#define REPETITIONS 5

    // Global vairables.
    int weak_scaling = 0; 											// Should the problem size increase with the number of nodes?
    int global_grid_x_size = 1000; 									// Grid size in x-dimension.
    int global_grid_y_size = 1000; 									// Same as above but y.
    int x_nodes = 1;
    int y_nodes = 1;
    int iterations_per_time = 100;									// Number of iterations per timestep.
    double kappa = 0.000001;

	double avg_runtime = 0.0, prev_avg_runtime = 0.0, stddev_runtime = 0.0;
	double start_time, end_time;

	int world_rank;
	
	// Read arguments
    int opt;
    while ((opt = getopt(argc, argv, "x:y:X:Y:st:T:k:r:")) != -1) 
    {
        switch (opt) 
        {
            case 'x':
                global_grid_x_size = atoi(optarg);
                break;
            case 'y':
                global_grid_y_size = atoi(optarg);
                break;
            case 'X':
                x_nodes = atoi(optarg);
                break;
            case 'Y':
                y_nodes = atoi(optarg);
                break;
            case 's':
                weak_scaling = 1;
                break;
            case 't':
                iterations_per_time = atoi(optarg);
                break;
            case 'k':
                kappa = atof(optarg);
                break;
			case 'r':

				break;
        }
    }

	Settings settings =
	{
		weak_scaling,											// Should the problem size increase with the number of nodes?
   		global_grid_x_size, 									// Grid size in x-dimension.
		global_grid_y_size, 									// Same as above but y.
		x_nodes,
		y_nodes,
		iterations_per_time,									// Number of iterations per timestep.
		kappa
	};

	// Print information about the run.
    #ifdef DEBUG
    if (world_rank == 0) 
    {
        printf("\nSolving the 2D heat equation with a grid of size, x:%d, y:%d\n", global_grid_x_size, global_grid_y_size);
        printf("Number of OpenMP threads per node: %d\n", omp_get_num_threads());
        printf("Number of nodes processing nodes, X: %d, Y: %d\n", x_nodes, y_nodes);

        if (weak_scaling)
        {
            printf("Weak scaling is enabled\n");
        }
        else 
        {
            printf("Strong scaling (default) is enabled\n");
        }

        printf("\n");
    }
    #endif

    MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	for (int r = 0; r < REPETITIONS; r++)
	{
		
		start_time = MPI_Wtime();
		solve_heat_equation(settings);
		end_time = MPI_Wtime();	

		if (world_rank == 0)
		{
			printf("run %d: %f s\n", r, end_time - start_time);
			
			prev_avg_runtime = avg_runtime;
			avg_runtime = avg_runtime + ( (end_time - start_time) - avg_runtime ) / (r + 1);
			stddev_runtime = stddev_runtime + ( (end_time - start_time) - avg_runtime) * ( (end_time - start_time) - prev_avg_runtime);
		}
	}

	if (world_rank == 0)
	{
		stddev_runtime = sqrt(stddev_runtime / (REPETITIONS - 1));
		printf("duration\t= %f±%f\n", avg_runtime, stddev_runtime);
	}

    MPI_Finalize();

	return 0;
}