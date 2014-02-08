#ifndef HIPP8_H
#define HIPP8_H

struct star {
  double ra, de, mag;
  double r[3];
  int streak;
};

extern struct star star[];
extern int nstar;

#endif

