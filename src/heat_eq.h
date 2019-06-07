#pragma once

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "grid.h"

#define MIN(x, y) (((x) < (y)) ? (x) : (y)) // Credit: https://stackoverflow.com/questions/3437404/min-and-max-in-c

typedef struct
{
    int weak_scaling; 											// Should the problem size increase with the number of nodes?
    int global_grid_x_size; 									// Grid size in x-dimension.
    int global_grid_y_size; 									// Same as above but y.
    int x_nodes;
    int y_nodes;
    int iterations_per_time;									// Number of iterations per timestep.
    double kappa;
} Settings;
void solve_heat_equation(Settings settings);