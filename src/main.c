#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "grid.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y)) // Credit: https://stackoverflow.com/questions/3437404/min-and-max-in-c

// General information about the code.
// Since arrays are row-major in C we will follow this convention.
// row 0, column 0 is in the lower left corner and increases to upwards and to the right respectively.

int main(int argc, char *argv[])
{
	#define DEBUG 1
	#define TIME 1000.0
	#define ITERATIONS_PER_TIME 100
	#define KAPPA 0.000001
	#define SIZE 1.0
	#define NDIMS 2
	#define PERIODIC 0
	#define REORDER 1

	// Global vairables.
	int world_size;
	int weak_scaling = 0; // Should the problem size increase with the number of nodes?
	int global_grid_x_size = 1000; // Grid size in x-dimension.
	int global_grid_y_size = 1000; // Same as above but y.
	int x_nodes = 1;
	int y_nodes = 1;
	int number_of_time_steps = (int) TIME * ITERATIONS_PER_TIME; // Default value.
	MPI_Comm cart_comm;
	double** global_grid = NULL;

	// Local variables
	int world_rank;
	int world_rank_2d;
	int world_coord_2d[NDIMS];
	int local_grid_x_size = 0;
	int local_grid_y_size = 0;
	int borders = 0;
	double** local_grid = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Read arguments
	int opt;
	while ((opt = getopt(argc, argv, "x:y:X:Y:s")) != -1) 
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
		}
	}

	// Print information about the run.
	if (world_rank == 0) 
	{
		printf("Solving the 2D heat equation with a grid of size x:%d, y:%d\n", global_grid_x_size, global_grid_y_size);
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

	// Check for configuration errors.
	// Check if the solution will converge according to the expression stated in the pdf.
	double delta_t = 1.0 / (double) ITERATIONS_PER_TIME;
	double h_x_sqrd = pow(SIZE / (double) global_grid_x_size, 2);
	double h_y_sqrd = pow(SIZE / (double) global_grid_y_size, 2);

	if (delta_t > (MIN(h_x_sqrd, h_y_sqrd) / (4 * KAPPA))) 
	{
		printf("ERROR: The time step length is not fine enough to ensure convergence!\n");
		exit(EXIT_FAILURE);
	}

	// Check that the calculation grid is evenly divisible amongst the nodes.
	if (global_grid_x_size % x_nodes != 0 || global_grid_y_size % y_nodes != 0)
	{
		printf("ERROR: The calculation grid cannot be evenly distibuted amongst the computation nodes!\n");
		exit(EXIT_FAILURE);
	}

	// Preprocessing and setting up communication.
	// Setting up the cartesian topology.
	int dimensions[NDIMS] = {y_nodes, x_nodes};
	int periodic[NDIMS] = {PERIODIC, PERIODIC};

	MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dimensions, periodic, REORDER, &cart_comm);
	MPI_Cart_rank(cart_comm, world_coord_2d, &world_rank_2d);

	#ifdef DEBUG
	MPI_Cart_coords(cart_comm, world_rank, NDIMS, world_coord_2d);
	printf("I am %d: (%d, %d); originally rank %d\n", world_rank_2d, world_coord_2d[0], world_coord_2d[1], world_rank);
	#endif

	// Grid generation/distribution
	if (world_rank == 0) 
	{
		if (!weak_scaling) // Strong scaling enabled (default).
		{
			// This means that we should generate a grid and "scatter" it.
			// Allocate memory.
			// TODO: We should distribute inital values in a sparse way instead!
			global_grid = create_two_dimensional_grid(global_grid_y_size, global_grid_x_size);

			// Distribute initial values


			// Split up and distribute the grid.
			// TODO
		}
		else
		{
			// TODO
		}
	}	

	// Perform the numerical solving of the equation.
	local_grid = global_grid;
	local_grid_x_size = global_grid_x_size;
	local_grid_y_size = global_grid_y_size;

	// Loop over time.
	for (int t = 0; t < number_of_time_steps; t++) 
	{
		// Create matrix for next time step.
		double** new_grid = create_two_dimensional_grid(local_grid_y_size, local_grid_x_size);

		// Set boundary values again for the new grid.
		set_boundary_values(new_grid, local_grid_y_size, local_grid_x_size, 1 | 2 | 4 | 8);

		// Peform the numerical solving using Taylor expansion.
		// We should never have a cell on the edge as the center in the iteration.
		// All cells on the edge are either boundary values or ghost cells.
		double kappa_delta_t = KAPPA * delta_t;

		for (int r = 1; r < (local_grid_y_size - 1); r++) 
		{
			for (int c = 1; c < (local_grid_x_size - 1); c++)
			{
				double x_derivative = (local_grid[r][c - 1] + local_grid[r][c + 1] - 2 * local_grid[r][c]) / h_x_sqrd;
				double y_derivative = (local_grid[r - 1][c] + local_grid[r + 1][c] - 2 * local_grid[r][c]) / h_y_sqrd;

				new_grid[r][c] = local_grid[r][c] + kappa_delta_t * (x_derivative + y_derivative);
			}
		}

		free_two_dimensional_grid(local_grid, local_grid_y_size);
		local_grid = new_grid;
	}

	print_grid(local_grid, local_grid_y_size, local_grid_x_size);
	print_grid_to_file("result.csv", local_grid, local_grid_y_size, local_grid_x_size);

	// TODO: Enable this when multinode communication is implemented. Currently frees a already freed array.
	//free_two_dimensional_grid(global_grid, global_grid_y_size);

	MPI_Finalize();

	return 0;
}