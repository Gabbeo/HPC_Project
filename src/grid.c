#include "grid.h"

void print_usage(char *program)
{
	fprintf(stderr, "Usage: %s [-r size of cell grid edge, default 1000]\n [-s enables scaling functionality]\n", program);
}

void print_grid(double* grid, int array_size_y, int array_size_x) 
{
	for (int r = 0; r < array_size_y; r++)
	{
		for (int c = 0; c < array_size_x; c++) 
		{
			printf("%f ", grid[r * array_size_x + c]);
		}
		printf("\n");
	}
	printf("\n");
}

/* 
	Sets the boundaries of the array to 1.0.
	The boundary int array tells what edges are part of the boundary.
*/
void set_boundary_values(double* array, int array_size_y, int array_size_x, int* ghost_borders) 
{
	// Set boundary values again for the new grid.
	
	if (!ghost_borders[N]) // Top row.
	{
		for (int i = 0; i < array_size_x; i++)
		{	
			array[(array_size_y - 1) * array_size_x + i] = 1.0;
		}
	}
	if (!ghost_borders[E]) // Rightmost column.
	{
		for (int i = 0; i < array_size_y; i++) 
		{
			array[i * array_size_x + array_size_x - 1] = 1.0;
		}
	}
	if (!ghost_borders[S]) // Bottom row.
	{
		for (int i = 0; i < array_size_x; i++)
		{
			array[i] = 1.0;
		}
	}
	if (!ghost_borders[W]) // Leftmost column.
	{
		for (int i = 0; i < array_size_y; i++) 
		{
			array[i * array_size_x] = 1.0;
		}
	}
}

void print_grid_to_file(char* file_name, double* array, int array_size_y, int array_size_x) {
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
				fprintf(file, "%f, ", array[r * array_size_x + c]);
			}
			else 
			{
				// Last cell on row should not have a following comma and should print 
				fprintf(file, "%f\n", array[r * array_size_x + c]);
			}
		}
	}
	fclose(file);
}

// The MPI_Neighbor_alltoallv sends/receives data in the following order, S, N, W, E (since we have (Y, X)-cartesian topology).  
// See https://stackoverflow.com/questions/50608184/what-is-the-correct-order-of-send-and-receive-in-mpi-neighbor-alltoallw
// For more detail.
void from_grid_to_ghost_array(double* grid, int grid_size_y, int grid_size_x, double* ghost_array) 
{
	// S
	for (int i = 0; i < grid_size_x; i++)
	{	
		ghost_array[i] = grid[1 * grid_size_x + i];
	}

	// N
	for (int i = 0; i < grid_size_x; i++)
	{
		ghost_array[grid_size_x + i] = grid[(grid_size_y - 2) * grid_size_x + i];
	}

	// W
	for (int i = 0; i < grid_size_y; i++) 
	{
		ghost_array[grid_size_x * 2 + i] = grid[i * grid_size_x + 1];
	}

	// E
	for (int i = 0; i < grid_size_y; i++) 
	{
		ghost_array[grid_size_x * 2 + grid_size_y + i] = grid[i * grid_size_x + grid_size_x - 2];
	}
}

// The MPI_Neighbor_alltoallv sends/receives data in the following order, S, N, W, E (since we have (Y, X)-cartesian topology).  
// See https://stackoverflow.com/questions/50608184/what-is-the-correct-order-of-send-and-receive-in-mpi-neighbor-alltoallw
// For more detail.
void from_ghost_array_to_grid(double* ghost_array, double* grid, int grid_size_y, int grid_size_x, int* borders)
{
	// S
	if (borders[S]) 
	{
		for (int i = 0; i < grid_size_x; i++)
		{	
			grid[i] = ghost_array[i];
		}
	}

	// N
	if (borders[N]) 
	{
		for (int i = 0; i < grid_size_x; i++)
		{
			grid[(grid_size_y - 1) * grid_size_x + i]= ghost_array[grid_size_x + i];
		}
	}

	// W
	if (borders[W])
	{
		for (int i = 0; i < grid_size_y; i++) 
		{
			grid[i * grid_size_x] = ghost_array[grid_size_x * 2 + i];
		}
	}

	// E
	if (borders[E])
	{
		for (int i = 0; i < grid_size_y; i++) 
		{
			grid[i * grid_size_x + grid_size_x - 1] = ghost_array[grid_size_x * 2 + grid_size_y + i];
		}
	}
}

double* remove_ghost_cells(double* old_grid, int old_grid_size_y, int old_grid_size_x, int* ghost_borders, int* new_grid_size_y, int* new_grid_size_x) 
{
	*new_grid_size_y = old_grid_size_y - (ghost_borders[N] + ghost_borders[S]);
	*new_grid_size_x = old_grid_size_x - (ghost_borders[E] + ghost_borders[W]);

	double* new_grid = malloc((*new_grid_size_x) * (*new_grid_size_y) * sizeof(double));
	
	for (int r = ghost_borders[S]; r < ((*new_grid_size_y) + ghost_borders[S]); r++)
	{
		int byte_length = (*new_grid_size_x) * sizeof(double);
		intptr_t from_ptr = (intptr_t) old_grid + (r * old_grid_size_x + ghost_borders[W]) * sizeof(double);
		intptr_t to_ptr = (intptr_t) new_grid + ((r - ghost_borders[S]) * (*new_grid_size_x)) * sizeof(double);

		memcpy((double*) to_ptr, (double*) from_ptr, byte_length);
	}

	return new_grid;
}

//Writes the output data into an output file in binary format using Collective MPI I/O
//The output file format is a <int><int> header with dimensions of the output matrix
//followed by doubles in row major order consisting of the data in the matrix.
void write_output(double* local_array, char* outfile_name, int proc_rank, int *cart_comm_coords,
	int local_array_size_y, int local_array_size_x, int global_array_size_y, int global_array_size_x)
{
	int dimensions[2] = {global_array_size_y, global_array_size_x};
	int local_dimensions[2] = {local_array_size_y, local_array_size_x};
	int starts_coords[2] = {local_array_size_y * cart_comm_coords[0], local_array_size_x * cart_comm_coords[1]};

	MPI_File outfile_handle;
	MPI_File_open(MPI_COMM_WORLD, outfile_name, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &outfile_handle);
	//Rank 0 writes the file header with the matrix dimensions
	if (proc_rank == 0)
	{
		MPI_File_write_at(outfile_handle, 0, dimensions, 2, MPI_INT, MPI_STATUS_IGNORE);
	}

	MPI_Datatype matrix_block;
	MPI_Type_create_subarray(2, dimensions, local_dimensions, starts_coords, MPI_ORDER_C, MPI_DOUBLE, &matrix_block);
	MPI_Type_commit(&matrix_block);
	
	MPI_Offset header_displacement = 2 * sizeof(int);
	MPI_File_set_view(outfile_handle, header_displacement, MPI_DOUBLE, matrix_block, "native", MPI_INFO_NULL);
	MPI_File_write_all(outfile_handle, local_array, local_array_size_x * local_array_size_y, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&outfile_handle); 
}

void stencil_2d_heat_eq(double* old_grid, double* new_grid, int grid_size_x, int r, int c, double h_x_sqrd, double h_y_sqrd, double kappa_delta_t)
{
	double x_derivative = (old_grid[r * grid_size_x + c - 1] + old_grid[r * grid_size_x + c + 1]
				- 2 * old_grid[r * grid_size_x + c]) / h_x_sqrd;

	double y_derivative = (old_grid[(r - 1) * grid_size_x + c] + old_grid[(r + 1) * grid_size_x + c]
		- 2 * old_grid[r * grid_size_x + c]) / h_y_sqrd;

	new_grid[r * grid_size_x + c] = old_grid[r * grid_size_x + c] + kappa_delta_t * (x_derivative + y_derivative);
}

