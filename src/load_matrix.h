#ifndef __LOAD_MATRIX__
#define __LOAD_MATRIX__

void load_shape(char * file_name, int *N, int *K, int *M);

/*
Loads the matrix A from the file file_name. The shape of the matrix is stored
in N, K
A is returned as matrix stored line-wised
*/
double *load_A(char * file_name, int *N, int *K);

/*
Loads the matrix B from the file file_name. The shape of the matrix is stored
in K, M
B is returned as matrix stored line-wised
*/
double *load_B(char * file_name, int *K, int *M);

/*
Loads the matrix C from the file file_name. The shape of the matrix is stored
in N, M
C is returned as matrix stored line-wised
*/
double *load_C(char * file_name, int *N, int *M);

/*
Loads the bloc of matrix A from the file file_name. (ip, jp) is the index of the bloc on a PxP grid. The shape of the matrix is stored in Np, Kp.
The bloc is returned as matrix stored line-wised and is padded with zeros.
*/
double *load_A_subpart(char * file_name, int *Np, int *Kp, int *N, int *K, int ip, int jp, int P);


/*
Loads the bloc of matrix B from the file file_name. (ip, jp) is the index of the bloc on a PxP grid. The shape of the matrix is stored in Kp, Mp.
The bloc is returned as matrix stored line-wised and is padded with zeros.
*/
double *load_B_subpart(char * file_name, int *Kp, int *Mp, int *K, int *M, int ip, int jp, int P);

/*
Loads the bloc of matrix C from the file file_name. (ip, jp) is the index of the bloc on a PxP grid. The shape of the matrix is stored in Np, Mp.
The bloc is returned as matrix stored line-wised and is NOT padded with zeros.
*/
double *load_C_subpart(char * file_name, int *Np, int *Mp, int *N, int *M, int ip, int jp, int P);

double *load_A_subpart_1(int N, int K, int *Np, int *Kp, int ip, int jp, int P);

double *load_B_subpart_1(int K, int M, int *Kp, int *Mp, int ip, int jp, int P);

double *load_A_subpart_rand(int N, int K, int *Np, int *Kp, int ip, int jp, int P);

double *load_B_subpart_rand(int K, int M, int *Kp, int *Mp, int ip, int jp, int P);

double *load_rand(int N, int M);

#endif
