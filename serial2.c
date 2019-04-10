#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

int max(int a, int b) {
    return a > b ? a : b;
}

int min(int a, int b) {
    return a < b ? a : b;
}

void free_matrix(char **A, int N) {
    if (N > 0) {
        for(int i = 0; i < N; i++) {
            free(A[i]);
        }
        free(A);
    }
}

void print_matrix(char **A, int N, int M) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            printf("%d ", (int) A[i][j]);
        }
        printf("\n");
    }
}

void print_matrix_from(char **A, int ii, int jj, int N, int M) {
    for(int i = ii; i < N; i++) {
        for(int j = jj; j < M; j++) {
            printf("%d ", (int) A[i][j]);
        }
        printf("\n");
    }
}

void add(char **C, int l, char **A, int ia, int ja, char **B, int ib, int jb) {
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < l; j++) {
            C[i][j] = A[(ia-1)*l + i][(ja-1)*l + j] + B[(ib-1)*l + i][(jb-1)*l + j];
        }
    }
}

void sub(char **C, int l, char **A, int ia, int ja, char **B, int ib, int jb) {
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < l; j++) {
            C[i][j] = A[(ia-1)*l + i][(ja-1)*l + j] - B[(ib-1)*l + i][(jb-1)*l + j];
        }
    }
}

void copy(char **C, int l, char **A, int ia, int ja) {
    for(int i = 0; i < l; i++) {
        for(int j = 0; j < l; j++) {
            C[i][j] = A[(ia-1)*l + i][(ja-1)*l + j];
        }
    }
}

char ** Strassen_recc(char **A, char **B, int d) {
    if (d <= 5) {
        int N = 1 << (d+1);
        char **C = (char **) malloc( N * sizeof(char *));
        for(int i = 0; i < N; i++) {
            C[i] = (char *) calloc(N, sizeof(char));
            for(int k = 0; k < N; k++) {
                for(int j = 0; j < N; j++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return C;
    } else {
        int l = 1 << d;

        char **tmp1 = (char **) malloc(l * sizeof(char *));
        char **tmp2 = (char **) malloc(l * sizeof(char *));
        for(int i = 0; i < l; i++) {
            tmp1[i] = (char *) malloc(l * sizeof(char));
            tmp2[i] = (char *) malloc(l * sizeof(char));
        }

        add(tmp1, l, A, 1, 1, A, 2, 2);
        add(tmp2, l, B, 1, 1, B, 2, 2);
        char **m1 = Strassen_recc(tmp1, tmp2, d-1);

        add(tmp1, l, A, 2, 1, A, 2, 2);
        copy(tmp2, l, B, 1, 1);
        char **m2 = Strassen_recc(tmp1, tmp2, d-1);

        copy(tmp1, l, A, 1, 1);
        sub(tmp2, l, B, 1, 2, B, 2, 2);
        char **m3 = Strassen_recc(tmp1, tmp2, d-1);

        copy(tmp1, l, A, 2, 2);
        sub(tmp2, l, B, 2, 1, B, 1, 1);
        char **m4 = Strassen_recc(tmp1, tmp2, d-1);

        add(tmp1, l, A, 1, 1, A, 1, 2);
        copy(tmp2, l, B, 2, 2);
        char **m5 = Strassen_recc(tmp1, tmp2, d-1);

        sub(tmp1, l, A, 2, 1, A, 1, 1);
        add(tmp2, l, B, 1, 1, B, 1, 2);
        char **m6 = Strassen_recc(tmp1, tmp2, d-1);

        sub(tmp1, l, A, 1, 2, A, 2, 2);
        add(tmp2, l, B, 2, 1, B, 2, 2);
        char **m7 = Strassen_recc(tmp1, tmp2, d-1);

        char **C = (char **) malloc(l * 2 * sizeof(char *));
        for (int i = 0; i < l * 2; i++) {
            C[i] = (char *) malloc(l * 2 * sizeof(char));

            if (i < l) {
                for (int j = 0; j < l * 2; j++) {
                    if (j < l) {
                        //C11
                        C[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
                    } else {
                        //C12
                        C[i][j] = m3[i][j-l] + m5[i][j-l];
                    }
                }
            } else {
                for (int j = 0; j < l * 2; j++) {
                    if (j < l) {
                        //C21
                        C[i][j] = m2[i-l][j] + m4[i-l][j];
                    } else {
                        //C22
                        C[i][j] = m1[i-l][j-l] - m2[i-l][j-l] + m3[i-l][j-l] + m6[i-l][j-l];
                    }
                }
            }
        }

        free_matrix(tmp1, l);
        free_matrix(tmp2, l);
        free_matrix(m1, l);
        free_matrix(m2, l);
        free_matrix(m3, l);
        free_matrix(m4, l);
        free_matrix(m5, l);
        free_matrix(m6, l);
        free_matrix(m7, l);

        return C;
    }
}

char ** Strassen(char **A, char **B, int N) {
    int d = floor(log(N)/log(2));
    bool realocated = false;
    if (N > 1 << d) {
        // If the matrices don't have a size power of 2
        // we need to reallocate them
        int Np = 1 << (d+1);
        char ** tmp = (char **) malloc(Np * sizeof(char *));
        for(int i = 0; i < Np; i++) {
            tmp[i] = (char *) calloc(Np, sizeof(char));
            if (i < N) {
                for(int j = 0; j < N; j++) {
                    tmp[i][j] = A[i][j];
                }
            }
        }
        A = tmp;
        tmp = (char **) malloc(Np * sizeof(char *));
        for(int i = 0; i < Np; i++) {
            tmp[i] = (char *) calloc(Np,  sizeof(char));
            if (i < N) {
                for(int j = 0; j < N; j++) {
                    tmp[i][j] = B[i][j];
                }
            }
        }
        B = tmp;
        realocated = true;
        printf("Realocation done\n");
    } else {
        d = d - 1;
    }

    char ** C = Strassen_recc(A, B, d);

    if (realocated) {
        // In the case that we redefined A and B we need to free them
        printf("Realocated\n");
        free_matrix(A, 1 << (d+1));
        free_matrix(B, 1 << (d+1));
        char **tmp = (char **) malloc(N * sizeof(char *));
        for(int i = 0; i < N; i++) {
            tmp[i] = (char *) calloc(N,  sizeof(char));
            for(int j = 0; j < N; j++) {
                tmp[i][j] = C[i][j];
            }
        }
        free_matrix(C, 1 << (d+1));
        C = tmp;
    }

    return C;
}


int main(int argc, char **argv)
{
    if (argc == 2) {
        int N = atoi(argv[1]);

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
        printf("Init done\n");

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

        printf("%s\n", correct ? "Success" : "Fail");

        free_matrix(A, N);
        free_matrix(B, N);
        free_matrix(C, N);
    }
}
