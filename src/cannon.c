#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include "load_matrix.h"

int main(int argc, char **argv)
{
    //int N = 211;


    int p, P, tag, rc;
    tag = 0;
    MPI_Status status;

    rc = MPI_Init(NULL, NULL);
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

    // int Np = N / P + (N % P == 0 ? 0 : 1);
    // int Npr = ip == P - 1 ? N - (N / Np)*Np : Np;
    // int Mpr = jp == P - 1 ? N - (N / Np)*Np : Np;
    //
    // int ip_init = (ip + jp) % P;
    // int jp_init = (ip + jp) % P;
    // int Npr_init = ip_init == P - 1 ? N - (N / Np)*Np : Np;
    // int Mpr_init = jp_init == P - 1 ? N - (N / Np)*Np : Np;
    //
    // double *A = (double *) calloc(Np * Np, sizeof(double));
    // double *A_tmp = (double *) calloc(Np * Np, sizeof(double));
    // double *B = (double *) calloc(Np * Np, sizeof(double));
    // double *B_tmp = (double *) calloc(Np * Np, sizeof(double));
    // double *C = (double *) calloc(Npr * Mpr, sizeof(double));
    //
    // for (int i = 0; i < Npr_init; i++) {
    //     for (int j = 0; j < Mpr; j++) {
    //         B[i * Np + j] = 1;
    //     }
    // }
    // for (int i = 0; i < Npr; i++) {
    //     for (int j = 0; j < Mpr_init; j++) {
    //         A[i * Np + j] = 1;
    //     }
    // }

    int Np, Npr, Mpr;
    double *A = load_A_subpart(argv[1], &Np, &Np, ip, (ip + jp) % P, P);
    double *A_tmp = (double *) calloc(Np * Np, sizeof(double));
    double *B = load_B_subpart(argv[1], &Np, &Np,(ip + jp) % P, jp, P);
    double *B_tmp = (double *) calloc(Np * Np, sizeof(double));
    double *C_res = load_C_subpart(argv[1], &Npr, &Mpr, ip, jp, P);
    double *C = (double *) calloc(Npr * Mpr, sizeof(double));


    if(p==24) {
        for (int i = 0; i < Np*Np; i++) {
            printf("%lf \n", A[i]);
        }
        for (int i = 0; i < Np*Np; i++) {
            printf("%lf \n", B[i]);
        }
        for (int i = 0; i < Npr*Mpr; i++) {
            printf("%lf \n", C_res[i]);
        }
    }

    for (int c = 0; c < P; c++) {
        // if(p==0) {
        //     for (int i = 0; i < Np*Np; i++) {
        //         printf("%lf, \n", A[i]);
        //     }
        // }

        for (int i = 0; i < Npr; i++) {
            for (int k = 0; k < Np; k++) {
                for (int j = 0; j < Mpr; j++) {
                    C[i * Mpr + j] += A[i * Np + k] * B[k * Np + j];
                }
            }
        }

        if ( P % 2 == 0) { // With a grid of even size a simple chestboard exchange pattern is enough
            if ((ip + jp) % 2 == 0) {
                rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
            } else {
                rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
                rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
            }
        } else { // In the odd size case we need to deal with borders
            if ((ip + jp) % 2 == 0) {
                if (jp != P - 1) {
                    // printf("(%d,%d) send A to (%d,%d)\n", ip, jp, ip, (jp + 1) % P);
                    rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                }
                if (ip != P - 1) {
                    rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                }
                if (jp != 0) {
                    // printf("(%d,%d) recv A fr (%d,%d)\n", ip, jp, ip, (jp + P - 1) % P);
                    rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                }
                if (ip != 0) {
                    rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
                }
            } else {
                if (jp != 0) {
                    // printf("(%d,%d) recv A fr (%d,%d)\n", ip, jp, ip, (jp + P - 1) % P);
                    rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                }
                if (ip != 0) {
                    rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
                }
                if (jp != P - 1) {
                    // printf("(%d,%d) send A to (%d,%d)\n", ip, jp, ip, (jp + 1) % P);
                    rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                }
                if (ip != P - 1) {
                    rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                }
            }

            if (jp == P - 1) {
                // printf("(%d,%d) send A to (%d,%d)\n", ip, jp, ip, 0);
                rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P, tag, MPI_COMM_WORLD);
            }
            if (ip == P - 1) {
                rc = MPI_Send(B, Np*Np, MPI_DOUBLE, jp, tag, MPI_COMM_WORLD);
            }
            if (jp == 0) {
                // printf("(%d,%d) recv A fr (%d,%d)\n", ip, jp, ip, P - 1);
                rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + P - 1, tag, MPI_COMM_WORLD, &status);
            }
            if (ip == 0) {
                rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, (P - 1) * P + jp, tag, MPI_COMM_WORLD, &status);
            }
        }
        double *tmp = A;
        A = A_tmp;
        A_tmp = tmp;
        tmp = B;
        B = B_tmp;
        B_tmp = tmp;

    }

    // int correct = true;
    // for (int i = 0; i < Npr; i++) {
    //     for (int j = 0; j < Mpr; j++) {
    //         if(C[i * Mpr + j] != N) {
    //             correct = false;
    //             break;
    //         }
    //     }
    // }
    int correct = true;
    for (int i = 0; i < Npr * Mpr; i++) {
        if (C[i] != C_res[i]) {
            correct = false;
            break;
        }
    }
    // if(p==0) {
    //     for (int i = 0; i < Npr; i++) {
    //         for (int j = 0; j < Mpr; j++) {
    //             printf("%lf, \n", C[i * Mpr + j]);
    //         }
    //     }
    // }
    // for (int i = 0; i < Npr * Mpr; i++) {
    //     if(C[i] != N) {
    //         correct = false;
    //         break;
    //     }
    // }
    printf("Node(%d,%d): %s\n", ip, jp, correct ? "Success" : "Fail");
    free(A);
    free(A_tmp);
    free(B);
    free(B_tmp);
    free(C);
    free(C_res);
    rc = MPI_Finalize();
    exit(0);
}
