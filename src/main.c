#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y)) // Credit: https://stackoverflow.com/questions/3437404/min-and-max-in-c

void print_usage();
double** create_two_dimensional_array(int, int);
void free_two_dimensional_array(double**, int);
void print_grid(double**, int, int);
void set_boundary_values(double**, int, int, int);
void print_to_file(char*, double**, int, int);

int main(int argc, char *argv[])
{
	#define TIME 1000.0
	#define ITERATIONS_PER_TIME 100
	#define KAPPA 0.000001
	#define SIZE 1.0

	// Global vairables.
	int world_size;
	int weak_scaling = 0; // Should the problem size increase with the number of nodes?
	int global_grid_x_size = 1000; // Grid size in x-dimension.
	int global_grid_y_size = 1000; // Same as above but y.
	int number_of_time_steps = (int) TIME * ITERATIONS_PER_TIME; // Default value.
	double** global_grid = NULL;

	// Local variables
	int world_rank;
	int local_grid_x_size = 0;
	int local_grid_y_size = 0;
	double** local_grid = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Read arguments
	int opt;
	while ((opt = getopt(argc, argv, "x:y:s")) != -1) 
	{
		switch (opt) 
		{
			case 'x':
				global_grid_x_size = atoi(optarg);
				break;
			case 'y':
				global_grid_y_size = atoi(optarg);
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

	// Check if the solution will converge according to the expression stated in the pdf.
	double delta_t = 1.0 / (double) ITERATIONS_PER_TIME;
	double h_x_sqrd = pow(SIZE / (double) global_grid_x_size, 2);
	double h_y_sqrd = pow(SIZE / (double) global_grid_y_size, 2);

	if (delta_t > (MIN(h_x_sqrd, h_y_sqrd) / (4 * KAPPA))) 
	{
		printf("ERROR: The time step length is not fine enough to ensure convergence!");
		return EXIT_FAILURE;
	}

	// Since arrays are row-major in C we will follow this convention.
	// row 0, column 0 is in the lower left corner and increases to upwards and to the right respectively.
	// Grid generation/distribution
	if (world_rank == 0) 
	{
		if (!weak_scaling) 
		{
			// This means that we should generate a grid and "scatter" it.
			// Allocate memory.
			global_grid = create_two_dimensional_array(global_grid_y_size, global_grid_x_size);

			// Set boundary values.
			set_boundary_values(global_grid, global_grid_y_size, global_grid_x_size, 1 | 2 | 4 | 8);

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
		double** new_grid = create_two_dimensional_array(local_grid_y_size, local_grid_x_size);

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

		free_two_dimensional_array(local_grid, local_grid_y_size);
		local_grid = new_grid;
	}

	print_grid(local_grid, local_grid_y_size, local_grid_x_size);
	print_to_file("result.csv", local_grid, local_grid_y_size, local_grid_x_size);

	// TODO: Enable this when multinode communication is implemented. Currently frees a already freed array.
	//free_two_dimensional_array(global_grid, global_grid_y_size);

	MPI_Finalize();

	return 0;
}

// --------------------- functions ---------------------

void print_usage(char *program)
{
	fprintf(stderr, "Usage: %s [-r size of cell grid edge, default 1000]\n [-s enables scaling functionality]\n", program);
}

double** create_two_dimensional_array(int array_size_y, int array_size_x)
{
	double** array = malloc(array_size_y * sizeof(int*));
	for (int i = 0; i < array_size_y; i++) 
	{
		array[i] = calloc(array_size_x, sizeof(double));
	}

	return array;
}

void free_two_dimensional_array(double** array, int array_size_y)
{
	if (array != NULL)
	{
		// Delete allocated memory.
		for (int i = 0; i < array_size_y; i++) {
			free(array[i]);
		}
		free(array);
	}
}

void print_grid(double** grid, int array_size_y, int array_size_x) 
{
	for (int r = 0; r < array_size_y; r++)
	{
		for (int c = 0; c < array_size_x; c++) 
		{
			printf("%f ", grid[r][c]);
		}
		printf("\n");
	}
	printf("\n");
}

/* 
	Sets the boundaries of the array to 1.0.
	The boundary int tells what edges should be set to one.
	boundaries | 1 is the top row.
	boundaries | 2 is the rightmost column.
	boundaries | 4 is the bottom row.
	boundaries | 8 is the leftmost column.
*/
void set_boundary_values(double** array, int array_size_y, int array_size_x, int boundaries) 
{
	// Set boundary values again for the new grid.
	
	if (boundaries & 1) // Top row.
	{
		for (int i = 0; i < array_size_x; i++)
		{	
			array[array_size_y - 1][i] = 1.0;
		}
	}
	if (boundaries & 2) // Rightmost column.
	{
		for (int i = 0; i < array_size_y; i++) 
		{
			array[i][array_size_x - 1] = 1.0;
		}
	}
	if (boundaries & 4) // Bottom row.
	{
		for (int i = 0; i < array_size_x; i++)
		{
			array[0][i] = 1.0;
		}
	}
	if (boundaries & 8) // Leftmost column.
	{
		for (int i = 0; i < array_size_y; i++) 
		{
			array[i][0] = 1.0;
		}
	}
}

void print_to_file(char* file_name, double** array, int array_size_y, int array_size_x) {
	FILE* file = fopen(file_name, "w");

	if (file == NULL) 
	{
		printf("ERROR: Could not open file for writing. Filename: %s\n", file_name);
		exit(EXIT_FAILURE);
	}

	for (int r = 0; r < array_size_y; r++)
	{
		for (int c = 0; c < array_size_x; c++) 
		{
			if (c != (array_size_x - 1)) // If not last element on row.
			{
				fprintf(file, "%f, ", array[r][c]);
			}
			else 
			{
				// Last cell on row should not have a following comma and should print 
				fprintf(file, "%f\n", array[r][c]);
			}
		}
	}

	fclose(file);
}