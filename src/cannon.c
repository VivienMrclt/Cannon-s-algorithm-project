#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "load_matrix.h"
#include "utils.h"
#include <time.h>

int main(int argc, char **argv)
{
    /*Initialization of the MPI environment*/
    int p, P, tag, rc;
    tag = 0;
    MPI_Status status;
    srand(time(NULL)); //Init of random seed for this execution

    rc = MPI_Init(NULL, NULL);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &p);

    /* Processing of the parameters*/
    bool bd_gt = false; // Do we broadcast the matrices and gather the results
    // It is also change the way we time
    bool timing = false; // Do we time the computing
    double start, end; // Start and end time
    bool verify = false; // Do we verify the result of the operation
    bool verbose = false;// Print the data transfers
    bool testfile = false; // Do we test from a file
    int N = 0;
    char *filename = "";

    parse_param(argc, argv, p, P, &bd_gt, &timing, &verify, &verbose, &testfile, &N, &filename, &rc);

    // Compute the size of the grid PxP
    P = sqrt(P);


    int ip = p / P;
    int jp = p % P;
    double *A_all;
    double *B_all;
    double *C_all;
    double *A;
    double *A_tmp;
    double *B;
    double *B_tmp;
    double *C_res;
    double *C;
    int Np, Npr, Mpr;

    if (bd_gt || (!testfile && verify)) {
        if (timing && bd_gt) {
            // Barrier for timing
            MPI_Barrier(MPI_COMM_WORLD);
            if (p == 0) {
                start = MPI_Wtime();
            }
        }

        //In the of broadcast/gather we need to send the matrices
        if (p==0) {
            if (verify) {
                A_all = (double *) calloc(N * N, sizeof(double));
                B_all = (double *) calloc(N * N, sizeof(double));
                C_all = (double *) calloc(N * N, sizeof(double));
            }

            // We send the correct A and B matrices to all the other processes
            for (int dp = 1; dp < P*P; dp++) {
                // Get the dest process's id
                int idp = dp / P;
                int jdp = dp % P;
                int Ndp;

                // Load the correct A and send it
                double *dA;
                dA = load_A_subpart_rand(N, N, &Ndp, &Ndp, idp, (idp + jdp) % P, P);
                if (verify) {
                    save_subpart(A_all, dA, N, N, idp, (idp + jdp) % P, P, true);
                }

                rc = MPI_Send(dA, Ndp*Ndp, MPI_DOUBLE, dp, tag, MPI_COMM_WORLD);
                verb_sentA(verbose, ip, jp, idp, jdp);
                free(dA);

                // Load the correct B and send it
                double *dB;
                dB = load_B_subpart_rand(N, N, &Ndp, &Ndp, (idp + jdp) % P, jdp, P);
                if (verify) {
                    save_subpart(B_all, dB, N, N, (idp + jdp) % P, jdp, P, true);
                }

                rc = MPI_Send(dB, Ndp*Ndp, MPI_DOUBLE, dp, tag, MPI_COMM_WORLD);
                verb_sentB(verbose, ip, jp, idp, jdp);
                free(dB);
            }

            A = load_A_subpart_rand(N, N, &Np, &Np, 0, 0, P);
            B = load_A_subpart_rand(N, N, &Np, &Np, 0, 0, P);
            if (verify) {
                save_subpart(A_all, A, N, N, ip, jp, P, true);
                save_subpart(B_all, B, N, N, ip, jp, P, true);
            }
            verb_m_loaded(verbose, ip, jp);

        } else {
            Np = N / P + (N % P == 0 ? 0 : 1);

            A = (double *) calloc(Np * Np, sizeof(double));
            B = (double *) calloc(Np * Np, sizeof(double));

            rc = MPI_Recv(A, Np*Np, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
            rc = MPI_Recv(B, Np*Np, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        }
        if (timing && !bd_gt) {
            // Barrier for timing
            MPI_Barrier(MPI_COMM_WORLD);
            if (p == 0) {
                start = MPI_Wtime();
            }
        }
    } else {
        if (testfile) {
            A = load_A_subpart(filename, &Np, &Np, ip, (ip + jp) % P, P);
            B = load_B_subpart(filename, &Np, &Np,(ip + jp) % P, jp, P);
            C_res = load_C_subpart(filename, &Npr, &Mpr, ip, jp, P);
        } else {
            A = load_A_subpart_rand(N, N, &Np, &Np, ip, (ip + jp) % P, P);
            B = load_A_subpart_rand(N, N, &Np, &Np,(ip + jp) % P, jp, P);
        }
        verb_m_loaded(verbose, ip, jp);
        if (timing) {
            // Barrier for timing
            MPI_Barrier(MPI_COMM_WORLD);
            if (p == 0) {
                start = MPI_Wtime();
            }
        }
    }
    // if(p==2) {
    //     printf("A:\n");
    //     for (int i = 0; i < Np*Np; i++) {
    //         printf("%lf, \n", A[i]);
    //     }
    //     printf("B:\n");
    //     for (int i = 0; i < Np*Np; i++) {
    //         printf("%lf, \n", B[i]);
    //     }
    //     printf("C_res:\n");
    //     for (int i = 0; i < Npr*Npr; i++) {
    //         printf("%lf, \n", C_res[i]);
    //     }
    // }


    A_tmp = (double *) calloc(Np * Np, sizeof(double));
    B_tmp = (double *) calloc(Np * Np, sizeof(double));
    Npr = (ip == P - 1 && N % Np != 0) ? N % Np : Np;
    Mpr = (jp == P - 1 && N % Np != 0) ? N % Np : Np;
    C = (double *) calloc(Npr * Mpr, sizeof(double));

    // if(p==1) {
    //     for (int i = 0; i < Np*Np; i++) {
    //         printf("%lf \n", A[i]);
    //     }
    //     for (int i = 0; i < Np*Np; i++) {
    //         printf("%lf \n", B[i]);
    //     }
    //     // for (int i = 0; i < Npr*Mpr; i++) {
    //     //     printf("%lf \n", C_res[i]);
    //     // }
    // }

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
                verb_sentA(verbose, ip, jp, ip, (jp + 1) % P);
                rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                verb_sentB(verbose, ip, jp, (ip + 1) % P, jp);
                rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
            } else {
                rc = MPI_Recv(A_tmp, Np*Np, MPI_DOUBLE, ip * P + (jp + P - 1) % P, tag, MPI_COMM_WORLD, &status);
                rc = MPI_Recv(B_tmp, Np*Np, MPI_DOUBLE, ((ip + P - 1) % P) * P + jp, tag, MPI_COMM_WORLD, &status);
                rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                verb_sentA(verbose, ip, jp, ip, (jp + 1) % P);
                rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                verb_sentB(verbose, ip, jp, (ip + 1) % P, jp);
            }
        } else { // In the odd size case we need to deal with borders
            if ((ip + jp) % 2 == 0) {
                if (jp != P - 1) {
                    // printf("(%d,%d) send A to (%d,%d)\n", ip, jp, ip, (jp + 1) % P);
                    rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P + (jp + 1) % P, tag, MPI_COMM_WORLD);
                    verb_sentA(verbose, ip, jp, ip, (jp + 1) % P);
                }
                if (ip != P - 1) {
                    rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                    verb_sentB(verbose, ip, jp, (ip + 1) % P, jp);
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
                    verb_sentA(verbose, ip, jp, ip, (jp + 1) % P);
                }
                if (ip != P - 1) {
                    rc = MPI_Send(B, Np*Np, MPI_DOUBLE, ((ip + 1) % P) * P + jp, tag, MPI_COMM_WORLD);
                    verb_sentB(verbose, ip, jp, (ip + 1) % P, jp);
                }
            }

            if (jp == P - 1) {
                // printf("(%d,%d) send A to (%d,%d)\n", ip, jp, ip, 0);
                rc = MPI_Send(A, Np*Np, MPI_DOUBLE, ip * P, tag, MPI_COMM_WORLD);
                verb_sentA(verbose, ip, jp, ip, (jp + 1) % P);
            }
            if (ip == P - 1) {
                rc = MPI_Send(B, Np*Np, MPI_DOUBLE, jp, tag, MPI_COMM_WORLD);
                verb_sentB(verbose, ip, jp, (ip + 1) % P, jp);
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

    if (timing && !bd_gt) {
        // Barrier for timing
        MPI_Barrier(MPI_COMM_WORLD);
        if (p == 0) {
            end = MPI_Wtime();
            printf("Cannon N=%d, PxP=%d\nTime: %fs\n", N, P*P, fabs(end-start));
        }
    }

    // Collecting the results if option -bg but also when we need to verify for a random matrix
    if(bd_gt || (!testfile && verify)) {
        if (p==0) {
            for (int dp = 1; dp < P*P; dp++) {
                // Get the dest process's id
                int idp = dp / P;
                int jdp = dp % P;
                int Ndp = N / P + (N % P == 0 ? 0 : 1);
                int Ndpr = (idp == P - 1 && N % Ndp != 0) ? N % Ndp : Ndp;
                int Mdpr = (jdp == P - 1 && N % Ndp != 0) ? N % Ndp : Ndp;

                double *dC = (double *) calloc(Ndpr * Mdpr, sizeof(double));
                rc = MPI_Recv(dC, Ndpr * Mdpr, MPI_DOUBLE, dp, tag, MPI_COMM_WORLD, &status);

                if (verify) {
                    save_subpart(C_all, dC, N, N, idp, jdp, P, false);
                }
                free(dC);
            }
            if (verify) {
                save_subpart(C_all, C, N, N, ip, jp, P, false);
            }


        } else {
            int Np = N / P + (N % P == 0 ? 0 : 1);
            int Npr = (ip == P - 1 && N % Np != 0) ? N % Np : Np;
            int Mpr = (jp == P - 1 && N % Np != 0) ? N % Np : Np;

            rc = MPI_Send(C, Npr*Mpr, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            verb_sentC(verbose, ip, jp, 0, 0);
        }

        if (timing && bd_gt) {
            // Barrier for timing
            MPI_Barrier(MPI_COMM_WORLD);
            if (p == 0) {
                end = MPI_Wtime();
                printf("Cannon N=%d, PxP=%d\nTime: %fs\n", N, P*P, fabs(end-start));
            }
        }
    }


    if (verify) {
        int correct = true;
        if (testfile) {
            for (int i = 0; i < Npr * Mpr; i++) {
                if (C[i] != C_res[i]) {
                    correct = false;
                    break;
                }
            }
            printf("Node(%d,%d): %s\n", ip, jp, correct ? "Success" : "Fail");
        } else {
            if (p==0) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        double c = 0;
                        for (int k = 0; k < N; k++) {
                            c += A_all[i * N + k] * B_all[k * N + j];
                        }
                        if(abs(C_all[i * N + j] - c) > 1e-15) {
                            correct = false;
                            break;
                        }
                    }
                }
                printf("Verification: %s\n", correct ? "Success" : "Fail");
            }


        }
    }

    if (p == 0 && verify && !testfile) {
        free(A_all);
        free(B_all);
        free(C_all);
    }
    free(A);
    free(A_tmp);
    free(B);
    free(B_tmp);
    free(C);
    if (testfile) {
        free(C_res);
    }
    rc = MPI_Finalize();
    exit(0);
}
