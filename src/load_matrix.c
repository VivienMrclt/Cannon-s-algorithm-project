#include <stdio.h>
#include <stdlib.h>

double *load_A(char * file_name, int *N, int *K) {
    int M;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", &M);
        double * A = calloc(*N * *K, sizeof(double));
        for (int i = 0; i < *N * *K; i++) {
            fscanf(file, "%lf", &A[i]);
        }
        fclose(file);
        return A;
    }
    return NULL;
}

double *load_B(char * file_name, int *K, int *M) {
    int N;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", &N);
        fscanf(file, "%d", K);
        fscanf(file, "%d", M);
        double * B = calloc(*K * *M, sizeof(double));
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

double *load_C(char * file_name, int *N, int *M) {
    int K;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", N);
        fscanf(file, "%d", &K);
        fscanf(file, "%d", M);
        double * C = calloc(*N * *M, sizeof(double));
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

double *load_A_subpart(char * file_name, int *Np, int *Kp, int ip, int jp, int P) {
    int N, K, M;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", &N);
        fscanf(file, "%d", &K);
        fscanf(file, "%d", &M);
        *Np = N / P + (N % P == 0 ? 0 : 1);
        *Kp = K / P + (K % P == 0 ? 0 : 1);

        int Npr = (ip == P - 1 && N % *Np != 0) ? N % *Np : *Np;
        int Kpr = (jp == P - 1 && K % *Kp != 0) ? K % *Kp : *Kp;

        double * A = calloc(*Np * *Kp, sizeof(double));
        // Go threw the lines before the bloc
        for (int i = 0; i < ip * *Np * K; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i < Npr; i++) {
            // Go threw the columns before the bloc on this line
            for (int l = 0; l < jp * *Kp; l++) {
                fscanf(file, "%lf", &jump);
            }
            // Scanning the elements of the bloc on this line
            for (int j = 0; j < Kpr; j++) {
                //printf("(%d, %d) -> %d\n", i, j, i * *Kp + j);
                fscanf(file, "%lf", &A[i * *Kp + j]);
            }
            // Go threw the columns after the bloc on this line
            for (int l = (jp + 1) * *Kp; l < K; l++) {
                fscanf(file, "%lf", &jump);
            }
        }
        fclose(file);
        return A;
    }
    return NULL;
}

double *load_B_subpart(char * file_name, int *Kp, int *Mp, int ip, int jp, int P) {
    int N, K, M;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", &N);
        fscanf(file, "%d", &K);
        fscanf(file, "%d", &M);
        *Kp = K / P + (K % P == 0 ? 0 : 1);
        *Mp = M / P + (M % P == 0 ? 0 : 1);

        int Kpr = (ip == P - 1 && K % *Kp != 0) ? K % *Kp : *Kp;
        int Mpr = (jp == P - 1 && M % *Mp != 0) ? M % *Mp : *Mp;

        double * B = calloc(*Kp * *Mp, sizeof(double));
        // Go threw A
        for (int i = 0; i < N * K; i++) {
            fscanf(file, "%lf", &jump);
        }
        // Go threw the lines before the bloc
        for (int i = 0; i < ip * *Kp * M; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i < Kpr; i++) {
            // Go threw the columns before the bloc on this line
            for (int l = 0; l < jp * *Mp; l++) {
                fscanf(file, "%lf", &jump);
            }
            // Scanning the elements of the bloc on this line
            for (int j = 0; j < Mpr; j++) {
                fscanf(file, "%lf", &B[i * *Mp + j]);
            }
            // Go threw the columns after the bloc on this line
            for (int l = (jp + 1) * *Mp; l < M; l++) {
                fscanf(file, "%lf", &jump);
            }
        }
        fclose(file);
        return B;
    }
    return NULL;
}

double *load_C_subpart(char * file_name, int *Npr, int *Mpr, int ip, int jp, int P) {
    int N, K, M;
    double jump;
    FILE *file = fopen(file_name, "r");
    if (file) {
        fscanf(file, "%d", &N);
        fscanf(file, "%d", &K);
        fscanf(file, "%d", &M);
        int Np = N / P + (N % P == 0 ? 0 : 1);
        int Mp = M / P + (M % P == 0 ? 0 : 1);

        *Npr = (ip == P - 1 && N % Np != 0) ? N % Np : Np;
        *Mpr = (jp == P - 1 && M % Mp != 0) ? M % Mp : Mp;

        double * C = calloc(*Npr * *Mpr, sizeof(double));

        // Go threw A
        for (int i = 0; i < N * K; i++) {
            fscanf(file, "%lf", &jump);
        }
        // Go threw B
        for (int i = 0; i < K * M; i++) {
            fscanf(file, "%lf", &jump);
        }
        // Go threw the lines before the bloc
        for (int i = 0; i < ip * Np * M; i++) {
            fscanf(file, "%lf", &jump);
        }
        for (int i = 0; i < *Npr; i++) {
            // Go threw the columns before the bloc on this line
            for (int l = 0; l < jp * Mp; l++) {
                fscanf(file, "%lf", &jump);
            }
            // Scanning the elements of the bloc on this line
            for (int j = 0; j < *Mpr; j++) {
                fscanf(file, "%lf", &C[i * *Mpr + j]);
            }
            // Go threw the columns after the bloc on this line
            for (int l = (jp + 1) * Mp; l < M; l++) {
                fscanf(file, "%lf", &jump);
            }
        }
        fclose(file);
        return C;
    }
    return NULL;
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
