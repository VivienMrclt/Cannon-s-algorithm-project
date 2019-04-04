#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>

int size_sublist(int N, int P, int p) {
    """
    Computes the optimal distribution on the criteria that every process should have the same number of values plus/minus 1.
    """
    if (p<0 || p >= P) {
        return 0;
    }
    int remaining = N;
    int share = 0;
    for (int i = 0; i <= p; i++) {
        share = floor(remaining / (P-i));
        remaining -= share;
    }
    return share ;
}

int min (int a, int b) {
    return a < b ? a : b;
}

void merging_btm(double *merged, double *l1, double *l2, int limit, int limit1, int limit2) {
    """
    Saves the begining of the merge list between l1 and l2 in merged. It stops after having merged limit elements.
    """
    int i1 = 0;
    int i2 = 0;
    for (int i = 0; i < limit; i++) {
        if(i1 < limit1 && i2 < limit2) {
            if (l1[i1] < l2[i2]) {
                merged[i] = l1[i1];
                i1++;
            } else {
                merged[i] = l2[i2];
                i2++;
            }
        } else if (i1 == limit1) {
            merged[i] = l2[i2];
            i2++;
        } else if (i2 == limit2) {
            merged[i] = l1[i1];
            i1++;
        }
    }
}

void merging_top(double *merged, double *l1, double *l2, int limit, int limit1, int limit2) {
    """
    Saves the end of the merge list between l1 and l2 in merged. It stops after having merged limit elements.
    """
    int i1 = limit1 - 1;
    int i2 = limit2 - 1;
    for (int i = limit -1; i >= 0; i--) {
        if(i1 >= 0 && i2 >= 0) {
            if (l1[i1] > l2[i2]) {
                merged[i] = l1[i1];
                i1--;
            } else {
                merged[i] = l2[i2];
                i2--;
            }
        } else if (i1 < 0) {
            merged[i] = l2[i2];
            i2--;
        } else if (i2 < 0) {
            merged[i] = l1[i1];
            i1--;
        }
    }
}

void swap(double *List, int a, int b) {
    double tmp = List[a];
    List[a] = List[b];
    List[b] = tmp;
}

int partition( double *List, int lo, int hi) {
    double pivot = List[hi];
    int i = lo;
    for (int j = lo; j < hi; j++) {
        if (List[j] < pivot) {
            swap(List, i, j);
            i++;
        }
    }
    swap(List, i, hi);
    return i;
}

void quicksort(double *List, int lo, int hi) {
    if (lo < hi) {
        int p = partition (List, lo, hi);
        quicksort(List, lo, p-1);
        quicksort(List, p+1, hi);
    }
}

int main(int argc, char **argv)
{
    int p, P, tag, rc;
    tag = 0;
    MPI_Status status;

    /* Find problem size N from command line */
    if (argc < 2) {
        fprintf(stderr, "No size N given\n");
        exit(1);
    }
    int N = atoi(argv[1]);

    // Init the MPI environment
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &p);

    /*Computing the optimal distribution on the criteria that every process should have the same number of values plus/minus 1.
    */
    int rcounts[P];
    int displs[P];
    int remaining = N;
    int share = 0;
    for (int i = 0; i < P; i++) {
        share = floor(remaining / (P-i));
        rcounts[i] = share;
        displs[i] = N - remaining;
        remaining -= share;
    }

    int Iprev = rcounts[p-1];
    int I = rcounts[p];
    int Inext = rcounts[p+1];

    // Generation of the list
    srand(p+1);
    double list[I];
    double *finalList = NULL;
    if (p == 0) {
        finalList = malloc(sizeof(double) * N);
    }

    for (int i = 0; i < I; i++) {
        list[i] = ((double) rand())/(RAND_MAX);
    }

    // Barrier for timing
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Sorting the list of the process
    quicksort(list, 0, I-1);


    // Even-Odd transposition algorithm
    bool evenprocess = p % 2 == 0;
    bool evenphase = true;

    double list_tmp[I];
    double list_recv[I];

    for (int step = 0; step <= P; step++) {
        if ((evenphase && evenprocess) || (!evenphase && !evenprocess)) {
            if (p != P - 1) {
                rc = MPI_Sendrecv(list, min(I, Inext), MPI_DOUBLE, p+1, tag, list_recv, min(I, Inext), MPI_DOUBLE, p+1, tag, MPI_COMM_WORLD, &status);
                merging_btm(list_tmp, list, list_recv, I, I, min(I, Inext));
                for (int i = 0; i < I; i++) {
                    list[i] = list_tmp[i];
                }
            }
        } else {
            if (p != 0) {
                rc = MPI_Sendrecv(list, min(I, Iprev), MPI_DOUBLE, p-1, tag, list_recv, min(I, Iprev), MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status);
                merging_top(list_tmp, list, list_recv, I, I, min(I, Iprev));
                for (int i = 0; i < I; i++) {
                    list[i] = list_tmp[i];
                }
            }
        }
        evenphase = !evenphase;
    }

    // End of timing
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    if (p == 0) {
        printf("Time: %fs\n",fabs(end-start));
    }

    rc = MPI_Gatherv(list, I, MPI_DOUBLE, finalList, rcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    rc = MPI_Finalize();
    return 0;
}
