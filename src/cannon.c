#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>

int main(int argc, char **argv)
{
    if (argc == 2) {
        int N = atoi(argv[1]);

        int p, P, tag, rc;
        tag = 0;
        MPI_Status status;

        rc = MPI_Init(&argc, &argv);
        rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
        rc = MPI_Comm_rank(MPI_COMM_WORLD, &p);

        if(pow(sqrt(P), 2) != P) {
            if (p==0) {
                printf("You must provide a square number of processors\n");
            }
            rc = MPI_Finalize();
            exit(0);
        }

        P = sqrt(P);

        int ip = p / P;
        int jp = p % P;

        int Np = N / P + (N % P == 0) ? 0 : 1;

        char *A = (char *) malloc(Np * Np * sizeof(char));
        char *B = (char *) malloc(Np * Np * sizeof(char));
        char *C = (char *) calloc(Np * Np, sizeof(char));
        for (int i = 0; i < Np * Np; i++) {
            A[i] = 1;
            B[i] = 1;
        }

        for (int c = 0; c < P; c++) {

            for (int i = 0; i < Np * Np; i++) {
                for (int k = 0; k < Np * Np; k++) {
                    for (int j = 0; j < Np * Np; j++) {
                        C[i * Np + j] += A[i * Np + k] * B[k * Np + j];
                    }
                }
            }

            if ((ip + jp) % 2 == 0) {
                //printf("Node(%d, %d): Send A to (%d, %d)\n", ip, jp, ip, (jp + 1) % P);
                rc = MPI_Send(&A, Np*Np, MPI_CHAR, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                //printf("Node(%d, %d): Recv A frm (%d, %d)\n", ip, jp, ip, (jp + P - 1) % P);
                rc = MPI_Recv(&A, Np*Np, MPI_CHAR, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                //printf("Node(%d, %d): Send B to (%d, %d)\n", ip, jp, (ip + 1) % P, jp);
                rc = MPI_Send(&B, Np*Np, MPI_CHAR, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                //printf("Node(%d, %d): Recv B frm (%d, %d)\n", ip, jp, (ip + P - 1) % P, jp);
                rc = MPI_Recv(&B, Np*Np, MPI_CHAR, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
            } else {
                //printf("Node(%d, %d): Recv A frm (%d, %d)\n", ip, jp, ip, (jp + P - 1) % P);
                rc = MPI_Recv(&A, Np*Np, MPI_CHAR, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                //printf("Node(%d, %d): Send A to (%d, %d)\n", ip, jp, ip, (jp + 1) % P);
                rc = MPI_Send(&A, Np*Np, MPI_CHAR, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                //printf("Node(%d, %d): Recv B frm (%d, %d)\n", ip, jp, (ip + P - 1) % P, jp);
                rc = MPI_Recv(&B, Np*Np, MPI_CHAR, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
                //printf("Node(%d, %d): Send B to (%d, %d)\n", ip, jp, (ip + 1) % P, jp);
                rc = MPI_Send(&B, Np*Np, MPI_CHAR, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
            }

        }

        int correct = true;
        for (int i = 0; i < Np * Np; i++) {
            if(C[i] != N) {
                correct = false;
                break;
            }
        }
        printf("Node(%d,%d): %s\n", ip, jp, correct ? "Success" : "Fail");

        rc = MPI_Finalize();
        return 0;
    } else {
        exit(0);
    }
}
