#include "amoeba.h"

struct amoeba_v {
  double E;
  void *v;
};

#include <stdlib.h>
#include <stdio.h>

static int amoeba_VEcmp(const void *a, const void *b) {
  if (((struct amoeba_v *)a)->E > ((struct amoeba_v *)b)->E)
    return 1;
  else if (((struct amoeba_v *)a)->E < ((struct amoeba_v *)b)->E)
    return -1;
  return 0;
}

#define AMOEBA_DEBUG 0

double amoeba(void **v, void *(*expcon)(void **, int n, void *, double, void *), void (*freeV)(void *), double (*calcE)(void *, void *), int n, void *params) {
  struct amoeba_v *V=malloc(sizeof(*V)*n);
  int nshrink=0;
  int i;
  for (i=0; i<n; i++)
    V[i].E=calcE(v[i], params);

  while (nshrink<32) {
    for (i=0; i<n; i++)
      V[i].v=v[i];
    qsort(V, n, sizeof(*V), amoeba_VEcmp);
    for (i=0; i<n; i++)
      v[i]=V[i].v;

    void *refl=expcon(v, n-1, v[n-1], 1., params); // new 1
    double reflE=calcE(refl, params);

    if (reflE >= V[0].E && reflE < V[n-2].E) {
      // reflection
#if AMOEBA_DEBUG
      printf("r");
#endif
      freeV(v[n-1]);
      v[n-1]=refl;                         // own 1
      V[n-1].E=reflE;
      continue;
    }

    if (reflE < V[0].E) {
      // expansion
#if AMOEBA_DEBUG
      printf("e");
#endif
      nshrink--;
      void *expa=expcon(v, n-1, v[n-1], 2., params); // new 2
      double expaE=calcE(expa, params);
      if (expaE < reflE) {
	freeV(v[n-1]);
	freeV(refl);                       // free 1
	v[n-1]=expa;                        // own 2
	V[n-1].E=expaE;
	continue;
      } else {
	freeV(v[n-1]);
	freeV(expa);                        // free 2
	v[n-1]=refl;                       // own 1
	V[n-1].E=reflE;
	continue;
      }
    }

    freeV(refl);                           // free 1

    if (reflE >= V[n-2].E) {
      // contraction
      void *cont=expcon(v, n-1, v[n-1], 0.5, params); // new 2
      double contE=calcE(cont, params);
      if (contE < V[n-1].E) {
#if AMOEBA_DEBUG
	printf("c");
#endif
	nshrink++;
	freeV(v[n-1]);
	v[n-1]=cont;                          // own 2
	V[n-1].E=contE;
	continue;
      } else {
	freeV(cont);                          // free 2
      }

    }

    // shrink
#if AMOEBA_DEBUG
    printf("s");
#endif
    nshrink+=n-1;
    for (i=1; i<n; i++) {
      void *shrink=expcon(v, 1, v[i], -0.5, params);
      double shrinkE=calcE(shrink, params);
      freeV(v[i]);
      v[i]=shrink;
      V[i].E=shrinkE;
    }
  }

  double retval=V[0].E;
#if AMOEBA_DEBUG
  printf("\nE=%g %g\n",V[0].E,V[n-1].E);
#endif
  free(V);
  return retval;
}

