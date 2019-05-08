#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "utils.h"
#include "mkl.h"

double mkl_fRand()
{
    return (double)rand() / RAND_MAX;
}

void mkl_load_shape(char * file_name, int *N, int *K, int *M) {
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", M);
    }
}

double *mkl_load_A(char * file_name, int *N, int *K) {
    int M;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", &M);
        double * A =  mkl_malloc(*N * *K * sizeof(double), 64);
        for (int i = 0; i < *N * *K; i++) {
            fscanf(file, "%lf", &A[i]);
        }
        fclose(file);
        return A;
    }
    return NULL;
}

double *mkl_load_B(char * file_name, int *K, int *M) {
    int N;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", &N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", M);
        double * B =  mkl_malloc(*K * *M * sizeof(double), 64);
        for (int i = 0; i < N * *K; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i <*K * *M; i++) {
            fscanf(file, "%lf", &B[i]);
        }
        fclose(file);
        return B;
    }
    return NULL;
}

double *mkl_load_C(char * file_name, int *N, int *M) {
    int K;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", &K);
        fscanf(file, "%d", M);
        double * C =  mkl_malloc(*N * *M * sizeof(double), 64);
        for (int i = 0; i < *N * K; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i <K * *M; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i <*N * *M; i++) {
            fscanf(file, "%lf", &C[i]);
        }
        fclose(file);
        return C;
    }
    return NULL;
}

double *mkl_load_A_subpart(char * file_name, int *Np, int *Kp, int *N, int *K, int ip, int jp, int P) {
    int M;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", &M);
        *Np = size_bloc(ip, *N, P);
        *Kp = size_bloc(jp, *K, P);

        double * A =  mkl_malloc(size_bloc_alloc(*N, P) * size_bloc_alloc(*K, P) * sizeof(double), 64);
        // Go threw the lines before the bloc
        for (int i = 0; i < size_blocs_before(ip, *N, P) * *K; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i < *Np; i++) {
            // Go threw the columns before the bloc on this line
            for (int l = 0; l < size_blocs_before(jp, *K, P); l++) {
                fscanf(file, "%lf", &jump);
            }
            // Scanning the elements of the bloc on this line
            for (int j = 0; j < *Kp; j++) {
                //printf("(%d, %d) -> %d\n", i, j, i * *Kp + j);
                fscanf(file, "%lf", &A[i * *Kp + j]);
            }
            // Go threw the columns after the bloc on this line
            for (int l = size_blocs_before(jp + 1, *K, P); l < *K; l++) {
                fscanf(file, "%lf", &jump);
            }
        }
        fclose(file);
        return A;
    }
    return NULL;
}

double *mkl_load_B_subpart(char * file_name, int *Kp, int *Mp, int *K, int *M, int ip, int jp, int P) {
    int N;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", &N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", M);
        *Kp = size_bloc(ip, *K, P);
        *Mp = size_bloc(jp, *M, P);

        double * B =  mkl_malloc(size_bloc_alloc(*K, P) * size_bloc_alloc(*M, P) * sizeof(double), 64);

        // Go threw A
        for (int i = 0; i < N * *K; i++) {
            fscanf(file, "%lf", &jump);
        }
        // Go threw the lines before the bloc
        for (int i = 0; i < size_blocs_before(ip, *K, P) * *M; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i < *Kp; i++) {
            // Go threw the columns before the bloc on this line
            for (int l = 0; l < size_blocs_before(jp, *M, P); l++) {
                fscanf(file, "%lf", &jump);
            }
            // Scanning the elements of the bloc on this line
            for (int j = 0; j < *Mp; j++) {
                fscanf(file, "%lf", &B[i * *Mp + j]);
            }
            // Go threw the columns after the bloc on this line
            for (int l = size_blocs_before(jp + 1, *M, P); l < *M; l++) {
                fscanf(file, "%lf", &jump);
            }
        }
        fclose(file);
        return B;
    }
    return NULL;
}

double *mkl_load_C_subpart(char * file_name, int *Np, int *Mp, int *N, int *M, int ip, int jp, int P) {
    int K;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", &K);
        fscanf(file, "%d", M);
        *Np = size_bloc(ip, *N, P);
        *Mp = size_bloc(jp, *M, P);

        double * C =  mkl_malloc(*Np * *Mp * sizeof(double), 64);

        // Go threw A
        for (int i = 0; i < *N * K; i++) {
            fscanf(file, "%lf", &jump);
        }
        // Go threw B
        for (int i = 0; i < K * *M; i++) {
            fscanf(file, "%lf", &jump);
        }
        // Go threw the lines before the bloc
        for (int i = 0; i < size_blocs_before(ip, *N, P) * *M; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i < *Np; i++) {
            // Go threw the columns before the bloc on this line
            for (int l = 0; l < size_blocs_before(jp, *M, P); l++) {
                fscanf(file, "%lf", &jump);
            }
            // Scanning the elements of the bloc on this line
            for (int j = 0; j < *Mp; j++) {
                fscanf(file, "%lf", &C[i * *Mp + j]);
            }
            // Go threw the columns after the bloc on this line
            for (int l = size_blocs_before(jp + 1, *M, P); l < *M; l++) {
                fscanf(file, "%lf", &jump);
            }
        }
        fclose(file);
        return C;
    }
    return NULL;
}


double *mkl_load_A_subpart_1(int N, int K, int *Np, int *Kp, int ip, int jp, int P) {
    *Np = size_bloc(ip, N, P);
    *Kp = size_bloc(jp, K, P);

    double *A = (double *)  mkl_malloc(size_bloc_alloc(N, P) * size_bloc_alloc(K, P) * sizeof(double), 64);

    for (int i = 0; i < *Np * *Kp; i++) {
        A[i] = 1;
    }
    return A;
}

double *mkl_load_B_subpart_1(int K, int M, int *Kp, int *Mp, int ip, int jp, int P) {
    *Kp = size_bloc(ip, K, P);
    *Mp = size_bloc(jp, M, P);

    double *B = (double *)  mkl_malloc(size_bloc_alloc(K, P) * size_bloc_alloc(M, P) * sizeof(double), 64);

    for (int i = 0; i < *Kp * *Mp; i++) {
        B[i] = 1;
    }
    return B;
}

double *mkl_load_A_subpart_rand(int N, int K, int *Np, int *Kp, int ip, int jp, int P) {
    *Np = size_bloc(ip, N, P);
    *Kp = size_bloc(jp, K, P);

    double *A = (double *)  mkl_malloc(size_bloc_alloc(N, P) * size_bloc_alloc(K, P) * sizeof(double), 64);

    for (int i = 0; i < *Np * *Kp; i++) {
        A[i] = mkl_fRand();
    }
    return A;
}

double *mkl_load_B_subpart_rand(int K, int M, int *Kp, int *Mp, int ip, int jp, int P) {
    *Kp = size_bloc(ip, K, P);
    *Mp = size_bloc(jp, M, P);

    double *B = (double *)  mkl_malloc(size_bloc_alloc(K, P) * size_bloc_alloc(M, P) * sizeof(double), 64);

    for (int i = 0; i < *Kp * *Mp; i++) {
        B[i] = mkl_fRand();
    }
    return B;
}

// int main(int argc, char **argv)
// {
//     int N, K, M;
//     double *A = load_A("test_files/test_lm_1.txt", &N, &K);
//     double *B = load_B("test_files/test_lm_1.txt", &K, &M);
//     double *C = load_C("test_files/test_lm_1.txt", &N, &M);
//     for (int i = 0; i < N * K; i++) {
//         printf("%lf\n", A[i]);
//     }
//     for (int i = 0; i < K * M; i++) {
//         printf("%lf\n", B[i]);
//     }
//     for (int i = 0; i < N * M; i++) {
//         printf("%lf\n", C[i]);
//     }
//     free(A);
//     free(B);
//     free(C);
//
//     int N00, M00;
//     int N01, M01;
//     int N10, M10;
//     int N11, M11;
//     double *A00 = load_A_subpart("test_files/test_lm_2.txt", &N00, &M00, 0, 0, 2);
//     double *A01 = load_A_subpart("test_files/test_lm_2.txt", &N01, &M01, 0, 1, 2);
//     double *A10 = load_A_subpart("test_files/test_lm_2.txt", &N10, &M10, 1, 0, 2);
//     double *A11 = load_A_subpart("test_files/test_lm_2.txt", &N11, &M11, 1, 1, 2);
//     for (int i = 0; i < N00 * M00; i++) {
//         printf("%lf\n", A00[i]);
//     }
//     for (int i = 0; i < N01 * M01; i++) {
//         printf("%lf\n", A01[i]);
//     }
//     for (int i = 0; i < N10 * M10; i++) {
//         printf("%lf\n", A10[i]);
//     }
//     for (int i = 0; i < N11 * M11; i++) {
//         printf("%lf\n", A11[i]);
//     }
//     double *B00 = load_B_subpart("test_files/test_lm_2.txt", &N00, &M00, 0, 0, 2);
//     double *B01 = load_B_subpart("test_files/test_lm_2.txt", &N01, &M01, 0, 1, 2);
//     double *B10 = load_B_subpart("test_files/test_lm_2.txt", &N10, &M10, 1, 0, 2);
//     double *B11 = load_B_subpart("test_files/test_lm_2.txt", &N11, &M11, 1, 1, 2);
//     for (int i = 0; i < N00 * M00; i++) {
//         printf("%lf\n", B00[i]);
//     }
//     for (int i = 0; i < N01 * M01; i++) {
//         printf("%lf\n", B01[i]);
//     }
//     for (int i = 0; i < N10 * M10; i++) {
//         printf("%lf\n", B10[i]);
//     }
//     for (int i = 0; i < N11 * M11; i++) {
//         printf("%lf\n", B11[i]);
//     }
//     double *C00 = load_C_subpart("test_files/test_lm_2.txt", &N00, &M00, 0, 0, 2);
//     double *C01 = load_C_subpart("test_files/test_lm_2.txt", &N01, &M01, 0, 1, 2);
//     double *C10 = load_C_subpart("test_files/test_lm_2.txt", &N10, &M10, 1, 0, 2);
//     double *C11 = load_C_subpart("test_files/test_lm_2.txt", &N11, &M11, 1, 1, 2);
//     for (int i = 0; i < N00 * M00; i++) {
//         printf("%lf\n", C00[i]);
//     }
//     for (int i = 0; i < N01 * M01; i++) {
//         printf("%lf\n", C01[i]);
//     }
//     for (int i = 0; i < N10 * M10; i++) {
//         printf("%lf\n", C10[i]);
//     }
//     for (int i = 0; i < N11 * M11; i++) {
//         printf("%lf\n", C11[i]);
//     }
//     free(A00);
//     free(A01);
//     free(A10);
//     free(A11);
//     free(B00);
//     free(B01);
//     free(B10);
//     free(B11);
//     free(C00);
//     free(C01);
//     free(C10);
//     free(C11);
// }
