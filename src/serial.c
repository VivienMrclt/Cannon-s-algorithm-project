#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "strassen.h"


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
