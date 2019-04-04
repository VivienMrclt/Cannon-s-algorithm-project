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

void free_matrix(double **A, int N) {
    if (N > 0) {
        for(int i = 0; i < N; i++) {
            free(A[i]);
        }
        free(A);
    }
}

void print_matrix(double **A, int N, int M) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < M; j++) {
            printf("%d ", (int) A[i][j]);
        }
        printf("\n");
    }
}

void print_matrix_from(double **A, int ii, int jj, int N, int M) {
    for(int i = ii; i < N; i++) {
        for(int j = jj; j < M; j++) {
            printf("%d ", (int) A[i][j]);
        }
        printf("\n");
    }
}

void add(double **A, double **B, double **C, int l, int ia, int ja, int Na, int Ma, int ib, int jb, int Nb, int Mb, int ic, int jc, int *Nc, int *Mc) {

    if ((Na <= ia && Nb <= ib) || (Ma <= ja && Mb <= jb)) {
        *Nc = 0;
        *Mc = 0;
        return;
    }
    // The code is redondant but it avoids unnecessary comparisons
    // This algorithm doesn't allocate the matrix where it is null

    // Computing the 0 zone limits
    int Nmin = min(l, max(min(Na - ia, Nb - ib), 0));
    int Mmin = min(l, max(min(Ma - ja, Mb - jb), 0));
    int Nmax = min(l, max(max(Na - ia, Nb - ib), 0));
    int Mmax = min(l, max(max(Ma - ja, Mb - jb), 0));


    *Nc = ic + Nmax;
    *Mc = jc + Mmax;

    for (int i = 0; i < Nmin; i++) {
        for(int j = 0; j < Mmin; j++) {
            C[ic + i][jc + j] = A[ia + i][ja + j] + B[ib + i][jb + j];
        }
        if(Ma - ja > Mb - jb) {
            for(int j = Mmin; j < Mmax; j++) {
                C[ic + i][jc + j] = A[ia + i][ja + j];
            }
        } else {
            for(int j = Mmin; j < Mmax; j++) {
                C[ic + i][jc + j] = B[ib + i][jb + j];
            }
        }
    }

    if(Na - ia > Nb - ib) {
        for (int i=Nmin; i < Nmax; i++) {
            for(int j = 0; j < Mmin; j++) {
                C[ic + i][jc + j] = A[ia + i][ja + j];
            }
            if(Na - ia > Nb - ib) {
                for(int j = Mmin; j < Mmax; j++) {
                    C[ic + i][jc + j] = A[ia + i][ja + j];
                }
            }
        }
    } else {
        for (int i=Nmin; i < Nmax; i++) {
            for(int j = 0; j < Mmin; j++) {
                C[ic + i][jc + j] = B[ib + i][jb + j];
            }
            if(Na - ia < Nb - ib) {
                for(int j = Mmin; j < Mmax; j++) {
                    C[ic + i][jc + j] = B[ib + i][jb + j];
                }
            }
        }
    }
}

void sub(double **A, double **B, double **C, int l, int ia, int ja, int Na, int Ma, int ib, int jb, int Nb, int Mb, int ic, int jc, int *Nc, int *Mc) {

    if ((Na <= ia && Nb <= ib) || (Ma <= ja && Mb <= jb)) {
        *Nc = 0;
        *Mc = 0;
        return;
    }
    // The code is redondant but it avoids unnecessary comparisons
    // This algorithm doesn't allocate the matrix where it is null

    // Computing the 0 zone limits
    int Nmin = min(l, max(min(Na - ia, Nb - ib), 0));
    int Mmin = min(l, max(min(Ma - ja, Mb - jb), 0));
    int Nmax = min(l, max(max(Na - ia, Nb - ib), 0));
    int Mmax = min(l, max(max(Ma - ja, Mb - jb), 0));

    *Nc = ic + Nmax;
    *Mc = jc + Mmax;

    for (int i = 0; i < Nmin; i++) {
        for(int j = 0; j < Mmin; j++) {
            C[ic + i][jc + j] = A[ia + i][ja + j] - B[ib + i][jb + j];
        }
        if(Ma - ja > Mb - jb) {
            for(int j = Mmin; j < Mmax; j++) {
                C[ic + i][jc + j] = A[ia + i][ja + j];
            }
        } else {
            for(int j = Mmin; j < Mmax; j++) {
                C[ic + i][jc + j] = -B[ib + i][jb + j];
            }
        }
    }

    if(Na - ia > Nb - ib) {
        for (int i=Nmin; i < Nmax; i++) {
            for(int j = 0; j < Mmin; j++) {
                C[ic + i][jc + j] = A[ia + i][ja + j];
            }
            if(Na - ia > Nb - ib) {
                for(int j = Mmin; j < Mmax; j++) {
                    C[ic + i][jc + j] = A[ia + i][ja + j];
                }
            }
        }
    } else {
        for (int i=Nmin; i < Nmax; i++) {
            for(int j = 0; j < Mmin; j++) {
                C[ic + i][jc + j] = -B[ib + i][jb + j];
            }
            if(Na - ia < Nb - ib) {
                for(int j = Mmin; j < Mmax; j++) {
                    C[ic + i][jc + j] = -B[ib + i][jb + j];
                }
            }
        }
    }
}



double ** Strassen(double **A, double **B, int d, int ia, int ja, int Na, int Ma, int ib, int jb, int Nb, int Mb, int *Nc, int *Mc) {
    if (Na <= ia || Nb <= ib || Ma <= ja || Mb <= jb) {
        *Nc = 0;
        *Mc = 0;
        return (double **) NULL;
    }

    if (d == -1) {
        double res = A[ia][ja] * B[ib][jb];
        if (res == 0) {
            *Nc = 0;
            *Mc = 0;
            return (double **) NULL;
        } else {
            double **C = (double **) malloc(sizeof(double *));
            C[0] = (double *) malloc(sizeof(double));
            *Nc=1;
            *Mc=1;
            C[0][0] = res;
            return C;
        }
    }
    int l = 1 << d;

    int N1, M1, N2, M2;

    double **tmp1 = (double **) malloc(l * sizeof(double *));
    double **tmp2 = (double **) malloc(l * sizeof(double *));
    for(int i = 0; i < l; i++) {
        tmp1[i] = (double *) malloc(l * sizeof(double));
        tmp2[i] = (double *) malloc(l * sizeof(double));
    }

    add(A, A, tmp1, l, ia, ja, Na, Ma, ia + l, ja + l, Na, Ma, 0, 0, &N1, &M1);
    add(B, B, tmp2, l, ib, jb, Nb, Mb, ib + l, jb + l, Nb, Mb, 0, 0, &N2, &M2);
    int Nm1, Mm1;
    double **m1 = Strassen(tmp1, tmp2, d-1, 0, 0, N1, M1, 0, 0, N2, M2, &Nm1, &Mm1);

    add(A, A, tmp1, l, ia + l, ja, Na, Ma, ia + l, ja + l, Na, Ma, 0, 0, &N1, &M1);
    int Nm2, Mm2;
    double **m2 = Strassen(tmp1, B, d-1, 0, 0, N1, M1, ib, jb, Nb, Mb, &Nm2, &Mm2);

    sub(B, B, tmp1, l, ib, jb + l, Nb, Mb, ib + l, jb + l, Nb, Mb, 0, 0, &N1, &M1);
    int Nm3, Mm3;
    double **m3 = Strassen(A, tmp1, d-1, ia, ja, Na, Ma, 0, 0, N1, M1, &Nm3, &Mm3);

    sub(B, B, tmp1, l, ib + l, jb, Nb, Mb, ib, jb, Nb, Mb, 0, 0, &N1, &M1);
    int Nm4, Mm4;
    double **m4 = Strassen(A, tmp1, d-1, ia + l, ja + l, Na, Ma, 0, 0, N1, M1, &Nm4, &Mm4);

    add(A, A, tmp1, l, ia, ja, Na, Ma, ia, ja + l, Na, Ma, 0, 0, &N1, &M1);
    int Nm5, Mm5;
    double **m5 = Strassen(tmp1, B, d-1, 0, 0, N1, M1, ib + l, jb + l, Nb, Mb, &Nm5, &Mm5);

    sub(A, A, tmp1, l, ia + l, ja, Na, Ma, ia, ja, Na, Ma, 0, 0, &N1, &M1);
    add(B, B, tmp2, l, ib, jb, Nb, Mb, ib, jb + l, Nb, Mb, 0, 0, &N2, &M2);
    int Nm6, Mm6;
    double **m6= Strassen(tmp1, tmp2, d-1, 0, 0, N1, M1, 0, 0, N2, M2, &Nm6, &Mm6);

    sub(A, A, tmp1, l, ia, ja + l, Na, Ma, ia + l, ja + l, Na, Ma, 0, 0, &N1, &M1);
    add(B, B, tmp2, l, ib + l, jb, Nb, Mb, ib + l, jb + l, Nb, Mb, 0, 0, &N2, &M2);
    int Nm7, Mm7;
    double **m7= Strassen(tmp1, tmp2, d-1, 0, 0, N1, M1, 0, 0, N2, M2, &Nm7, &Mm7);

    free_matrix(tmp1, l);
    free_matrix(tmp2, l);

    *Nc = max(max(max(Nm2, Nm4), max(Nm1, Nm3)), Nm6);
    if (*Nc == 0) {
        *Nc = max(Nm5, Nm7);
    } else {
        *Nc += l;
    }
    *Mc = max(max(max(Mm3, Mm5), max(Mm1, Mm2)), Mm6);
    if (*Mc == 0) {
        *Mc = max(Mm4, Mm7);
    } else {
        *Mc += l;
    }


    double **C = (double **) malloc(*Nc * sizeof(double *));
    for(int i = 0; i < *Nc; i++) {
        C[i] = (double *) malloc(*Mc * sizeof(double));
    }


    //C11
    add(m1, m4, C, l, 0, 0, Nm1, Mm1, 0, 0, Nm4, Mm4, 0, 0, &N1, &M1);
    sub(C, m5, C, l, 0, 0, N1, M1, 0, 0, Nm5, Mm5, 0, 0, &N2, &M2);
    add(C, m7, C, l, 0, 0, N2, M2, 0, 0, Nm7, Mm7, 0, 0, &N1, &M1);

    //C12
    add(m3, m5, C, l, 0, 0, Nm3, Mm3, 0, 0, Nm5, Mm5, 0, l, &N1, &M1);

    //C21
    add(m2, m4, C, l, 0, 0, Nm2, Mm2, 0, 0, Nm4, Mm4, l, 0, &N1, &M1);

    //C22
    sub(m1, m2, C, l, 0, 0, Nm1, Mm1, 0, 0, Nm2, Mm2, l, l, &N1, &M1);
    add(C, m3, C, l, l, l, N1, M1, 0, 0, Nm3, Mm3, l, l, &N2, &M2);
    add(C, m6, C, l, l, l, N2, M2, 0, 0, Nm6, Mm6, l, l, &N1, &M1);


    free_matrix(m1, Nm1);
    free_matrix(m2, Nm2);
    free_matrix(m3, Nm3);
    free_matrix(m4, Nm4);
    free_matrix(m5, Nm5);
    free_matrix(m6, Nm6);
    free_matrix(m7, Nm7);

    return C;
}

int main(int argc, char **argv)
{
    const int N = 256;

    // Initialization of the matrices
    double **A = (double **) malloc(N * sizeof(double *));
    double **B = (double **) malloc(N * sizeof(double *));
    double **C = (double **) malloc(N * sizeof(double *));
    for(int i = 0; i < N; i++) {
        A[i] = (double *) malloc(N * sizeof(double));
        B[i] = (double *) malloc(N * sizeof(double));
        C[i] = (double *) malloc(N * sizeof(double));
        for(int j = 0; j < N; j++) {
            A[i][j] = 1;
            B[i][j] = 1;
        }
    }


    int d = floor(log(N)/log(2));
    if (N == 1 << d) {
        d = d - 1;
    }
    // int Nc, Mc;
    int Nd, Md;
    //sub(A, B, C, N, 0, 0, N, N, 0, 0, N, N, 0, 0, &Nc, &Mc);
    double **D = Strassen(A, B, d, 0, 0, N, N, 0, 0, N, N, &Nd, &Md);

    double target = N;
    bool correct = true;
    for(int i = 0; i < Nd; i++) {
        for(int j = 0; j < Md; j++) {
            printf("%d ", (int) D[i][j]);
            if (D[i][j] != target) {
                correct = false;
            }
        }
        printf("\n");
    }
    if(Nd != N && Md !=N) {
        correct = false;
    }
    printf("Nd:%d, Md:%d\n", Nd, Md);

    printf("%s\n", correct ? "Success" : "Fail");
    free_matrix(A, N);
    free_matrix(B, N);
    free_matrix(C, N);

}
