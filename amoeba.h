#ifndef AMOEBA_H
#define AMOEBA_H

/*
  amoeba_move should create a new state which is equal to sum of (1+a)*(V[0], V[1]... V[n-1])/n - a*X.
 */

typedef void *amoeba_move_t(void **V, int n, void *X, double a, void *params);
typedef void amoeba_free_t(void *);
typedef double amoeba_E_t(void *, void *);

double amoeba(void **v, amoeba_move_t *, amoeba_free_t *, amoeba_E_t *, int n, void *params);

#endif

