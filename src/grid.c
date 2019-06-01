#include "grid.h"

void print_usage(char *program)
{
	fprintf(stderr, "Usage: %s [-r size of cell grid edge, default 1000]\n [-s enables scaling functionality]\n", program);
}

double** create_two_dimensional_grid(int array_size_y, int array_size_x)
{
	double** array = malloc(array_size_y * sizeof(int*));
	for (int i = 0; i < array_size_y; i++) 
	{
		array[i] = calloc(array_size_x, sizeof(double));
	}

	return array;
}

void free_two_dimensional_grid(double** array, int array_size_y)
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

void print_grid_to_file(char* file_name, double** array, int array_size_y, int array_size_x) {
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