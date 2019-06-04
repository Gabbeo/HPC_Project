#pragma once

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

enum directions {N, E, S, W, dir_count};

void print_usage();
double** create_two_dimensional_grid(int array_size_y, int array_size_x);
void print_grid(double* grid, int array_size_y, int array_size_x);
void set_boundary_values(double* array, int array_size_y, int array_size_x, int* ghost_borders);
void print_grid_to_file(char* file_name, double* array, int array_size_y, int array_size_x);
void from_grid_to_ghost_array(double* grid, int grid_size_y, int grid_size_x, double* ghost_array);
void from_ghost_array_to_grid(double* ghost_array, double* grid, int grid_size_y, int grid_size_x, int* borders);
double* remove_ghost_cells(double* grid, int old_grid_size_y, int old_grid_size_x, int* ghost_borders, int* new_grid_size_y, int* new_grid_size_x);
void write_output(double* local_array, char* outfile_name, int proc_rank, int *cart_comm_coords,
	int local_array_size_y, int local_array_size_x, int global_array_size_y, int global_array_size_x);