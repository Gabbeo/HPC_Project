## HPC Project - solving the 2d heat equation using MPI

## Prerequisites

The program needs the module PrgEnv-gnu to be loaded for compilation.

## Compiling

The code can be compiled by using
    
    make

The resulting binary is created in bin/heat_equation.out

## Running 

To run the binary file you need to specify the number of OpenMP threads by running

    export OMP_NUM_THREADS=<number_threads_per_node>
    
And then run the program with

    aprun -n <number_of_nodes> -d <number_of_cores_per_node> ./bin/heat_equation.out -x <grid_size_x> -y <grid_size_y> -X <nodes_x> -Y <nodes_y> -t <number_of_iterations> -k <diffusion_constant>

An additional command '-s' enables weak scaling which means that the grid size scales with the number of nodes.
