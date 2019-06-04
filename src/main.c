#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "grid.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y)) // Credit: https://stackoverflow.com/questions/3437404/min-and-max-in-c

// General information about the code.
// Since arrays are row-major in C we will follow this convention.
// row 0, column 0 is in the lower left corner and increases to upwards and to the right respectively.

int main(int argc, char *argv[])
{
	#define DEBUG
	#define TIME 100.0
	#define ITERATIONS_PER_TIME 20
	#define KAPPA 0.000001
	#define SIZE 1.0
	#define NDIMS 2
	#define PERIODIC 0
	#define REORDER 1

	// Global vairables.
	int world_size;
	int weak_scaling = 0; 											// Should the problem size increase with the number of nodes?
	int global_grid_x_size = 1000; 									// Grid size in x-dimension.
	int global_grid_y_size = 1000; 									// Same as above but y.
	int x_nodes = 1;
	int y_nodes = 1;
	int number_of_time_steps = (int) TIME * ITERATIONS_PER_TIME; 	// Default value.
	MPI_Comm cart_comm;
	double *global_grid;
	
	// Local variables
	int world_rank;
	int world_rank_2d;
	int world_coord_2d[NDIMS];
	int local_grid_x_size = 0;
	int local_grid_y_size = 0;
	int ghost_borders[4]; 												// 0 = N, 1 = E, 2 = S, 3 = W. 
	int bordering_ranks[4]; 											// 0 = N, 1 = E, 2 = S, 3 = W. MPI_PROC_NULL if it's an edge border.
	double *local_grid;

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
	#ifdef DEBUG
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
	#endif

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
	MPI_Cart_coords(cart_comm, world_rank, NDIMS, world_coord_2d);
	MPI_Cart_rank(cart_comm, world_coord_2d, &world_rank_2d);

	// Calculate if a border is a ghost-cell border or if it's border values using the cartesian topology.
	int source[4] = {world_rank_2d, world_rank_2d, world_rank_2d, world_rank_2d}; // MPI_Cart_shift changes the source input so we need a buffer.

	// Get the neighbours for the current rank.
	MPI_Cart_shift(cart_comm, 0, +1, &(source[N]), &(bordering_ranks[N]));
	MPI_Cart_shift(cart_comm, 1, +1, &(source[E]), &(bordering_ranks[E]));
	MPI_Cart_shift(cart_comm, 0, -1, &(source[S]), &(bordering_ranks[S]));
	MPI_Cart_shift(cart_comm, 1, -1, &(source[W]), &(bordering_ranks[W]));

	// Store which borders should be ghost borders.
	for (int i = 0; i < dir_count; i++) 
	{
		ghost_borders[i] = (bordering_ranks[i] != MPI_PROC_NULL);
	}

	#ifdef DEBUG
	printf("I am %d: (%d, %d); originally 2d_rank %d.\n", world_rank_2d, world_coord_2d[0], world_coord_2d[1], world_rank);
	#endif

	// Grid generation/distribution for world_rank 0.
	if (world_rank == 0) 
	{
		if (!weak_scaling) // Strong scaling enabled (default).
		{
			// This means that we should generate a grid and "scatter" it.
			// Allocate memory.
			// TODO: We should distribute inital values in a sparse way instead!
		}
		else // Weak scaling.
		{
			// TODO
		}
	}	

	// Perform the numerical solving of the equation.
	// We should increase the size of the grid if there are ghost cells on the edges.
	// The bitshifting adds one to the dimension if there is a ghost cell on that dimension.
	local_grid_x_size = (global_grid_x_size / x_nodes) + ghost_borders[E] + ghost_borders[W];
	local_grid_y_size = (global_grid_y_size / y_nodes) + ghost_borders[N] + ghost_borders[S];
	local_grid = calloc(local_grid_x_size * local_grid_y_size, sizeof(double));

	#ifdef DEBUG
	printf("Local grid size for rank %d: x: %d, y: %d\n", world_rank_2d, local_grid_x_size, local_grid_y_size);
	#endif

	// Prepare for the halo exchange before entering the loop.
	// The MPI_Neighbor_alltoallv sends/receives data in the following order, S, N, W, E (since we have (Y, X)-cartesian topology).  
	// See https://stackoverflow.com/questions/50608184/what-is-the-correct-order-of-send-and-receive-in-mpi-neighbor-alltoallw
	// For more detail.
	double* sendbuf = malloc((local_grid_x_size * 2 + local_grid_y_size * 2) * sizeof(double));
	double* recvbuf = malloc((local_grid_x_size * 2 + local_grid_y_size * 2) * sizeof(double));

	// Calculate how many elements will be sent/received to each direction.
	int sendcounts[dir_count];
	sendcounts[0] = local_grid_x_size; // S
	sendcounts[1] = local_grid_x_size; // N
	sendcounts[2] = local_grid_y_size; // W
	sendcounts[3] = local_grid_y_size; // E

	// Calculate offsets in the send- and receive-buffer respectively.
	int sdispls[dir_count] = {0, local_grid_x_size, local_grid_x_size * 2, local_grid_x_size * 2 + local_grid_y_size};

	// Loop over time.
	for (int t = 0; t < number_of_time_steps; t++) 
	{
		if (world_rank_2d == 0 && t % 100 == 0) 
		{
			printf("Currently running iteration %d.\n", t);
		}

		// Create matrix for next time step.
		double* new_grid = calloc(local_grid_x_size * local_grid_y_size, sizeof(double));

		// Peform the numerical solving using Taylor expansion.
		// We should never have a cell on the edge as the center in the iteration.
		// All cells on the edge are either boundary values or ghost cells.
		double kappa_delta_t = KAPPA * delta_t;

		for (int r = 1; r < (local_grid_y_size - 1); r++) 
		{
			for (int c = 1; c < (local_grid_x_size - 1); c++)
			{
				double x_derivative = (local_grid[r * local_grid_x_size + c - 1] + local_grid[r * local_grid_x_size + c + 1]
					- 2 * local_grid[r * local_grid_x_size + c]) / h_x_sqrd;

				double y_derivative = (local_grid[(r - 1) * local_grid_x_size + c] + local_grid[(r + 1) * local_grid_x_size + c]
					- 2 * local_grid[r * local_grid_x_size + c]) / h_y_sqrd;

				new_grid[r * local_grid_x_size + c] = local_grid[r * local_grid_x_size + c] + kappa_delta_t * (x_derivative + y_derivative);
			}
		}

		free(local_grid);
		local_grid = new_grid;

		// Send/receive ghost cells.
		// First fill sendbuffer.
		from_grid_to_ghost_array(local_grid, local_grid_y_size, local_grid_x_size, sendbuf);

		// Send/receive
		MPI_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf, sendcounts, sdispls, MPI_DOUBLE, cart_comm);

		// Transfer from receivebuffer to grid.
		from_ghost_array_to_grid(recvbuf, local_grid, local_grid_y_size, local_grid_x_size, ghost_borders);

		// Set boundary values again for the new grid.
		set_boundary_values(local_grid, local_grid_y_size, local_grid_x_size, ghost_borders);
	}

	// Save the data with collective I/O.
	int new_grid_size_x, new_grid_size_y = 0;
	double* new_grid = remove_ghost_cells(local_grid, local_grid_y_size, local_grid_x_size, ghost_borders, &new_grid_size_y, &new_grid_size_x);
	write_output(new_grid, "result.bin", world_rank_2d, world_coord_2d, new_grid_size_y, new_grid_size_x, new_grid_size_y * y_nodes, new_grid_size_x * x_nodes);

	free (new_grid);
	free(local_grid);
	free(sendbuf);
	free(recvbuf);

	MPI_Finalize();

	return 0;
}