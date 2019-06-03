#pragma once

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

enum directions {N, E, S, W, dir_count};

void print_usage();
double** create_two_dimensional_grid(int, int);
void free_two_dimensional_grid(double**, int);
void print_grid(double**, int, int);
void set_boundary_values(double**, int, int, int*);
void print_grid_to_file(char*, double**, int, int);
void from_grid_to_ghost_array(double**, int, int, double*);
void from_ghost_array_to_grid(double*, double**, int, int, int*);