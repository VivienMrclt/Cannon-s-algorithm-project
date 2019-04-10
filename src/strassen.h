#ifndef __STRASSEN__
#define __STRASSEN__

char ** Strassen(char **A, char **B, int N);

void free_matrix(char **A, int N);

void print_matrix(char **A, int N, int M);

#endif
