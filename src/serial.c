#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "strassen.h"
#include <mpi.h>

int main(int argc, char **argv)
{
    /*Initialization of the MPI environment*/
    int p, P, tag, rc;
    double start, end;
    int N = 0;
    tag = 0;
    MPI_Status status;
    srand(time(NULL)); //Init of random seed for this execution

    rc = MPI_Init(NULL, NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &p);

    if (argc == 2) {
        N = atoi(argv[1]);

        // Initialization of the matrices
        char **A = (char **) malloc(N * sizeof(char *));
        char **B = (char **) malloc(N * sizeof(char *));
        for(int i = 0; i < N; i++) {
            A[i] = (char *) malloc(N * sizeof(char));
            B[i] = (char *) malloc(N * sizeof(char));
            for(int j = 0; j < N; j++) {
                B[i][j] = 1;
                A[i][j] = 1;
            }
        }
        start = MPI_Wtime();
        //printf("Init done\n");

        char **C = Strassen(A, B, N);

        char target = N;
        bool correct = true;
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                //printf("%d ", (int) C[i][j]);
                if (C[i][j] != target) {
                    correct = false;
                }
            }
            //printf("\n");
        }

        //printf("%s\n", correct ? "Success" : "Fail");

        free_matrix(A, N);
        free_matrix(B, N);
        free_matrix(C, N);
    }
    if (p == 0) {
        end = MPI_Wtime();
        printf("Strassen N=%d\nTime: %fs\n", N, fabs(end-start));
    }
    rc = MPI_Finalize();
    exit(0);
}
