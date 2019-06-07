#include "heat_eq.h"

void solve_heat_equation(Settings settings)
{
    #define DEBUG
    #define SIZE 100.0
    #define NDIMS 2
    #define PERIODIC 0
    #define REORDER 1

    // Global vairables.
    int world_size;
    int weak_scaling = settings.weak_scaling; 											// Should the problem size increase with the number of nodes?
    int global_grid_x_size = settings.global_grid_x_size; 									// Grid size in x-dimension.
    int global_grid_y_size = settings.global_grid_y_size; 									// Same as above but y.
    int x_nodes = settings.x_nodes;
    int y_nodes = settings.y_nodes;
    int iterations_per_time = settings.iterations_per_time;									// Number of iterations per timestep.
    double kappa = settings.kappa;
    MPI_Comm cart_comm;

    // Local variables
    int world_rank;
    int world_rank_2d;
    int world_coord_2d[NDIMS];
    int local_grid_x_size = 0;
    int local_grid_y_size = 0;
    int ghost_borders[4]; 												// 0 = N, 1 = E, 2 = S, 3 = W. 
    int bordering_ranks[4];                                             // 0 = N, 1 = E, 2 = S, 3 = W. MPI_PROC_NULL if it's an edge border.

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (weak_scaling == 1) 
    {
        global_grid_x_size *= x_nodes;
        global_grid_y_size *= y_nodes;
    }

    printf("x: %d, y: %d", global_grid_x_size, global_grid_y_size);

    double *local_grid;
    double h_x_sqrd = pow(SIZE / (double) global_grid_x_size, 2);
    double h_y_sqrd = pow(SIZE / (double) global_grid_y_size, 2);
    double delta_t = h_x_sqrd * h_y_sqrd / (2 * kappa * (h_x_sqrd + h_y_sqrd)) ;

    // Check for configuration errors.
    // Check if the solution will converge according to the expression stated in the pdf.
    if (world_rank == 0)
    {
        if (delta_t > (MIN(h_x_sqrd, h_y_sqrd) / (4 * kappa))) 
        {
            fprintf(stderr, "ERROR: The time step length is not fine enough to ensure convergence!\n");
            exit(EXIT_FAILURE);
        }

        // Check that the number of total processes matches the division of nodes.
        if (x_nodes <= 0 || y_nodes <= 0)
        {
            fprintf(stderr, "ERROR: Number of nodes must be greater than 0!\n");
            exit(EXIT_FAILURE);
        }

        if (global_grid_x_size <= 0 || global_grid_y_size <= 0)
        {
            fprintf(stderr, "The size of the grid must be greater than zero in both directions!\n");
            exit(EXIT_FAILURE);
        } 

        if (x_nodes * y_nodes != world_size)
        {
            fprintf(stderr, "The number of nodes in the node-grid must be equal to the nodes assigned in mpirun!\n");
            exit(EXIT_FAILURE);
        }

        // Check that the calculation grid is evenly divisible amongst the nodes.
        if (global_grid_x_size % x_nodes != 0 || global_grid_y_size % y_nodes != 0)
        {
            fprintf(stderr, "ERROR: The calculation grid cannot be evenly distibuted amongst the computation nodes!\n");
            exit(EXIT_FAILURE);
        }
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

    #ifdef DEBUG3
    printf("I am %d: (%d, %d); originally 2d_rank %d.\n", world_rank_2d, world_coord_2d[0], world_coord_2d[1], world_rank);
    #endif

    // Perform the numerical solving of the equation.
    // We should increase the size of the grid if there are ghost cells on the edges.
    // The bitshifting adds one to the dimension if there is a ghost cell on that dimension.
    local_grid_x_size = (global_grid_x_size / x_nodes) + ghost_borders[E] + ghost_borders[W];
    local_grid_y_size = (global_grid_y_size / y_nodes) + ghost_borders[N] + ghost_borders[S];
    local_grid = calloc(local_grid_x_size * local_grid_y_size, sizeof(double));

    #ifdef DEBUG3
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

    // Calculate this outside the loop.
    double kappa_delta_t = kappa * delta_t;

    // Loop over time.
    for (int t = 0; t < iterations_per_time; t++) 
    {
        #ifdef DEBUG
        if (world_rank_2d == 0 && t % 100 == 0) 
        {
            printf("Currently running iteration %d.\n", t);
        }
        #endif

        // Create matrix for next time step.
        double* new_grid = calloc(local_grid_x_size * local_grid_y_size, sizeof(double));

        // Peform the numerical solving using Taylor expansion.
        // We should never have a cell on the edge as the center in the iteration.
        // All cells on the edge are either boundary values or ghost cells.

        // First prepare and send ghost cells with non-blocking group collectives.
        // First fill sendbuffer.
        from_grid_to_ghost_array(local_grid, local_grid_y_size, local_grid_x_size, sendbuf);
        // Send/receive
        MPI_Request request;
        MPI_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf, sendcounts, sdispls, MPI_DOUBLE, cart_comm, &request);

        // First we do the inner cells that do not depend on the ghost cells.
        #pragma omp parallel for
        for (int r = 2; r < (local_grid_y_size - 2); r++) 
        {
            for (int c = 2; c < (local_grid_x_size - 2); c++)
            {
                stencil_2d_heat_eq(local_grid, new_grid, local_grid_x_size, r, c, h_x_sqrd, h_y_sqrd, kappa_delta_t);
            }
        }

        // Set boundary values again for the new grid. Should be done before out edges are calculated.
        set_boundary_values(local_grid, local_grid_y_size, local_grid_x_size, ghost_borders);

        // Wait for communcation to finish.
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        // Transfer from receivebuffer to grid.
        from_ghost_array_to_grid(recvbuf, local_grid, local_grid_y_size, local_grid_x_size, ghost_borders);
        
        // Perform computations on the edge.
        #pragma omp parallel
        {
            #pragma omp for
            for (int r = 1; r < (local_grid_y_size - 1); r++) // Calculate E & W edges.
            {
                int c = 1;
                stencil_2d_heat_eq(local_grid, new_grid, local_grid_x_size, r, c, h_x_sqrd, h_y_sqrd, kappa_delta_t);

                c = local_grid_x_size - 2; // < (local_grid_x_size - 1)
                stencil_2d_heat_eq(local_grid, new_grid, local_grid_x_size, r, c, h_x_sqrd, h_y_sqrd, kappa_delta_t);
            }

            #pragma omp for
            for (int c = 1; c < (local_grid_x_size - 1); c++)
            {
                int r = 1;
                stencil_2d_heat_eq(local_grid, new_grid, local_grid_x_size, r, c, h_x_sqrd, h_y_sqrd, kappa_delta_t);

                r = local_grid_y_size - 2; // < (local_grid_y_size - 1)
                stencil_2d_heat_eq(local_grid, new_grid, local_grid_x_size, r, c, h_x_sqrd, h_y_sqrd, kappa_delta_t);
            }
        }
        free(local_grid); // Delete old grid.
        local_grid = new_grid;
    }

    // Set boundary values one last time for niceness!dd
    set_boundary_values(local_grid, local_grid_y_size, local_grid_x_size, ghost_borders);

    // Save the data with collective I/O.
    int new_grid_size_x, new_grid_size_y = 0;
    double* new_grid = remove_ghost_cells(local_grid, local_grid_y_size, local_grid_x_size, ghost_borders, &new_grid_size_y, &new_grid_size_x);
    write_output(new_grid, "result.bin", world_rank_2d, world_coord_2d, new_grid_size_y, new_grid_size_x, new_grid_size_y * y_nodes, new_grid_size_x * x_nodes);

    free(new_grid);
    free(local_grid);
    free(sendbuf);
    free(recvbuf);
}