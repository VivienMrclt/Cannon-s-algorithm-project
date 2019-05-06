#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>



void verb_m_loaded(bool verbose, int i, int j) {
    if (verbose) {
        printf("Node(%d,%d): Matrice loaded\n", i, j);
    }
}

void verb_sentA(bool verbose, int i1, int j1, int i2, int j2) {
    if (verbose) {
        printf("Node(%d,%d): Sent A to (%d,%d)\n", i1, j1, i2, j2);
    }
}

void verb_sentB(bool verbose, int i1, int j1, int i2, int j2) {
    if (verbose) {
        printf("Node(%d,%d): Sent B to (%d,%d)\n", i1, j1, i2, j2);
    }
}

void verb_sentC(bool verbose, int i1, int j1, int i2, int j2) {
    if (verbose) {
        printf("Node(%d,%d): Sent C to (%d,%d)\n", i1, j1, i2, j2);
    }
}

void save_subpart(double *C, double *Cp, int N, int M, int ip, int jp, int P, bool padded) {
    int Np = N / P + (N % P == 0 ? 0 : 1);
    int Mp = M / P + (M % P == 0 ? 0 : 1);

    int Npr = (ip == P - 1 && N % Np != 0) ? N % Np : Np;
    int Mpr = (jp == P - 1 && M % Mp != 0) ? M % Mp : Mp;

    for (int i = 0; i < Npr; i++) {
        for (int j = 0; j < Mpr; j++) {
            C[(ip * Np + i) * M + jp * Mp + j] = Cp[i * (padded ? Mp : Mpr) + j];
        }
    }

}

void parse_param (int argc, char **argv, int p, int P, bool *bd_gt, bool *timing, bool *verify, bool *verbose, bool *testfile, int *N, char **filename, int *rc) {
    char *help =
"Use:\n\
   mpirun -np #process ./cannon -s matrix_size [-bg]| -t filename [-T|-V|-v]\n\
\n\
Options:\n\
   -s, --shape               Specify the shape of the matrices A and B. Example: -s 5 means multiplication of two matrices of shape 5x5\n\
   -t, --test                Specify what file to use for A and B. It should be followed by the path to the file. For more information on the test files format refer to README.txt\n\
   -h, --help                Display help\n\
   -bg, --broadcast-gather   Broadcast the matrices A and B from the process 0 to the others and gather the result at the end. This option can only be used with -s (not -t)\n\
   -T, --time                Print the computing time. Depending on the flag -bg the time will include or not the broadcast and gather steps\n\
   -v, --verbose             Prints the data transfers\n\
   -V, --verify              Verify the result. If we use -s the result is computed with a serial algo. and compared to the parallel result. If we use -t the result is compared to the result of the test file\n";
    if (argc >= 3) {
        int i = 1;
        bool correct = false;
        //Parsing the arguments
        while (i < argc) {
            if ((strcmp(argv[i],"-s") == 0) | (strcmp(argv[i],"--shape") == 0)) {
                if (i + 1 < argc && !correct) {
                    *N = atoi(argv[i+1]);
                    i++;
                    correct = true;
                } else {
                    correct = false;
                }
            } else if ((strcmp(argv[i],"-t") == 0) | (strcmp(argv[i],"--test") == 0)) {
                if (i + 1 < argc && !correct && !*bd_gt) {
                    *testfile = true;
                    *filename = argv[i+1];
                    i++;
                    correct = true;
                } else {
                    correct = false;
                }
            } else if ((strcmp(argv[i],"-bg") == 0)| (strcmp(argv[i],"--broadcast-gather") == 0)) {
                if (!*testfile) {
                    *bd_gt = true;
                } else {
                    correct = false;
                }
            } else if ((strcmp(argv[i],"-v") == 0) | (strcmp(argv[i],"--verbose") == 0)) {
                *verbose = true;
            } else if ((strcmp(argv[i],"-T") == 0) | (strcmp(argv[i],"--time") == 0)) {
                *timing = true;
            } else if ((strcmp(argv[i],"-V") == 0) | (strcmp(argv[i],"--verify") == 0)) {
                 *verify = true;
            } else {
                if(p==0) {
                    printf("%s", help);
                }
                *rc = MPI_Finalize();
                exit(0);
            }
            i++;
        }
        if (!correct) {
            if(p==0) {
                printf("%s", help);
            }
            *rc = MPI_Finalize();
            exit(0);
        }
    } else {
        if(p==0) {
            printf("%s", help);
        }
        *rc = MPI_Finalize();
        exit(0);
    }

    /*Verify that we can have a square processor grid*/
    if(pow(sqrt(P), 2) != P) {
        if (p==0) {
            printf("You must provide a square number of processors\n");
        }
        *rc = MPI_Finalize();
        exit(0);
    }
}

int size_bloc(int i, int N, int P) {
    printf("SB: %d\n", N / P + (N % P > i ? 1 : 0));
    return N / P + (N % P > i ? 1 : 0);
    // int Np = N / P + (N % P == 0 ? 0 : 1);
    // return (i == P - 1 && N % Np != 0) ? N % Np : Np;
}

int size_blocs_before(int i, int N, int P) {
    printf("SBB: %d\n", N / P * i + (N % P < i ? N % P : i));
    return N / P * i + (N % P < i ? N % P : i);
    // int Np = N / P + (N % P == 0 ? 0 : 1);
    // return i * Np;
}
