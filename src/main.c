#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

void print_usage();
double** create_two_dimensional_array(int);
void free_two_dimensional_array();

int main(int argc, char *argv[])
{
	#define TIME 10.0
	#define KAPPA 0.00001
	#define SIZE 1.0

	// Global vairables.
	int world_size;
	int scaling_grid = 0; // Should the problem size increase with the number of nodes?
	int global_grid_edge_size = 1000; // Default value.
	int number_of_time_steps = 10000; // Default value.
	double** global_grid = NULL;

	// Local variables
	int world_rank;
	int local_grid_size = 0;
	double** local_grid = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Read arguments
	int opt;
	while ((opt = getopt(argc, argv, "N:s")) != -1) 
	{
		switch (opt) 
		{
			case 'N':
				global_grid_edge_size = atoi(optarg);
				if (world_rank == 0) printf("Solving the 2D heat equation with a grid of size %d\n", global_grid_edge_size);
				break;
			case 's':
				scaling_grid = 1;
				break;
		}
	}

	// Since arrays are row-major in C we will follow this convention.
	// row 0, column 0 is in the lower left corner and increases to upwards and to the right respectively.
	// Grid generation/distribution
	if (world_rank == 0) 
	{
		if (!scaling_grid) 
		{
			// This means that we should generate a grid and "scatter" it.
			// Allocate memory.
			global_grid = create_two_dimensional_array(global_grid_edge_size);

			// Set boundary values.
			for (int i = 0; i < global_grid_edge_size; i++)
			{
				global_grid[0][i] = 1.0; 						 // Bottom row.
				global_grid[global_grid_edge_size - 1][i] = 1.0; // Top row.
				global_grid[i][0] = 1.0;						 // Leftmost column.
				global_grid[i][global_grid_edge_size - 1] = 1.0; // Rightmost column.
			}

			// Split up and distribute the grid.
			// TODO
		}
		else
		{
			// TODO
		}
	}

	local_grid = global_grid;
	local_grid_size = global_grid_edge_size;

	// Loop over time.
	for (int t = 0; t < number_of_time_steps; t++) 
	{
		// Create matrix for next time step.
		double** new_grid = create_two_dimensional_array(local_grid_size);

		// Set boundary values again for the new grid.
		for (int i = 0; i < local_grid_size; i++)
		{
			new_grid[0][i] = 1.0; 					// Bottom row.
			new_grid[local_grid_size - 1][i] = 1.0; // Top row.
			new_grid[i][0] = 1.0;					// Leftmost column.
			new_grid[i][local_grid_size - 1] = 1.0; // Rightmost column.
		}

		// Peform the numerical solving using Taylor expansion.
		// We should never have a cell on the edge as the center in the iteration.
		// All cells on the edge are either boundary values or ghost cells.
		for (int r = 1; r < (local_grid_size - 1); r++) 
		{
			for (int c = 1; c < (local_grid_size - 1); c++)
			{
				double x_derivative = (local_grid[r][c - 1] + local_grid[r][c + 1] - 2 * local_grid[r][c]) / pow(SIZE / global_grid_edge_size, 2);
				double y_derivative = (local_grid[r - 1][c] + local_grid[r + 1][c] - 2 * local_grid[r][c]) / pow(SIZE / global_grid_edge_size, 2);

				new_grid[r][c] = local_grid[r][c] + KAPPA * (TIME / number_of_time_steps) * (x_derivative + y_derivative);
			}
		}

		free_two_dimensional_array(local_grid, local_grid_size);
		local_grid = new_grid;

		print_grid(local_grid, local_grid_size);
	}

	print_grid(local_grid, local_grid_size);

	free_two_dimensional_array(global_grid, global_grid_edge_size);

	MPI_Finalize();

	return 0;
}

void print_usage(char *program)
{
	fprintf(stderr, "Usage: %s [-r size of cell grid edge, default 1000]\n [-s enables scaling functionality]\n", program);
}

double** create_two_dimensional_array(int array_size_one_dim)
{
	double** array = malloc(array_size_one_dim * sizeof(int*));
	for (int i = 0; i < array_size_one_dim; i++) 
	{
		array[i] = calloc(array_size_one_dim, sizeof(double));
	}

	return array;
}

void free_two_dimensional_array(double** array, int one_dim_size)
{
	if (array != NULL) 
	{
		// Delete allocated memory.
		for (int i = 0; i < one_dim_size; i++) {
			free(array[i]);
		}
		free(array);
	}
}

void print_grid(double** grid, int one_dim_size) 
{
	for (int r = 0; r < one_dim_size; r++)
	{
		for (int c = 0; c < one_dim_size; c++) 
		{
			printf("%f ", grid[r][c]);
		}
		printf("\n");
	}
	printf("\n");
}
