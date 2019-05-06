#ifndef __UTILS__
#define __UTILS__


/* Functions useful for verbose messages*/
void verb_m_loaded(bool verbose, int i, int j);

void verb_sentA(bool verbose, int i1, int j1, int i2, int j2);

void verb_sentB(bool verbose, int i1, int j1, int i2, int j2);

void verb_sentC(bool verbose, int i1, int j1, int i2, int j2);

/* Manipulation of matrices*/
void save_subpart(double *C, double *Cp, int N, int M, int ip, int jp, int P, bool padded);

/* Parsing of parameters*/
void parse_param (int argc, char **argv, int p, int P, bool *bd_gt, bool *timing, bool *verify, bool *verbose, bool *testfile, int *N, char **filename, int *rc);

int size_bloc(int i, int N, int P);

int size_blocs_before(int i, int N, int P);


#endif
