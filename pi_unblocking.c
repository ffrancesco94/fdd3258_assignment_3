
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int local_count = 0;
    int rank, num_ranks, i, provided;
    double x, y, z, pi;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    double start_time, stop_time, elapsed_time;
    start_time = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    local_count = 0;
    srand(SEED*rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < NUM_ITER/num_ranks; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            local_count++;
        }
    }
    
    if (rank == 0)
    {
	    int counts[num_ranks-1];
	    MPI_Request requests[num_ranks - 1];
	    int global_count = 0;

	    for (i = 1; i < num_ranks; i++)
	    {
		    MPI_Irecv(&counts[i - 1], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[i - 1]);
	    }
	    MPI_Waitall(num_ranks-1, requests, MPI_STATUSES_IGNORE);
	    
	    global_count += local_count;
	    for (i = 0; i < num_ranks; i++)
	    {
		    global_count += counts[i];
	    }

	    pi = ((double)global_count / (double)NUM_ITER) * 4.0;
    }
    else
    {
	    MPI_Send(&local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } 
    
    stop_time = MPI_Wtime();
    elapsed_time = stop_time - start_time;
    if (rank == 0)
    {
    	    printf("pi: %f\n", pi);
	    printf("Execution_time: %f\n", elapsed_time);
    }

    MPI_Finalize();
    
    return 0;
}

