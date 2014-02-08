/*

Copyright (c) 2014, Svein Berge
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include "amoeba.h"
#include "hipp8.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define PI 3.1415926535897932384626433832795

struct streak {
  double x;
  double y;
  double angle;
  double len;
  double intensity;
  double r[3];
  double v[3];
  int use;
  int star;
};

#if 0
struct star *loadstar(FILE *mag6, int *nstar) {
#define MAXSTAR 50000
  struct star *star=malloc(sizeof(*star)*MAXSTAR);
  char buf[1024];

  int i;
  for (i=0; i<35; i++)
    fgets(buf, 1024, mag6);

#define MAXMAG 6.

  *nstar=0;

  for (i=0; i<MAXSTAR; i++) {
    float d,m,s;
    fscanf(mag6, "%lg\t", &star[*nstar].ra);
    if (feof(mag6))
      break;
    fscanf(mag6, "%lg\t", &star[*nstar].de);
    star[*nstar].ra *= PI/180.;
    star[*nstar].de *= PI/180.;

    star[*nstar].r[0]=cos(star[*nstar].ra)*cos(star[*nstar].de);
    star[*nstar].r[1]=-sin(star[*nstar].ra)*cos(star[*nstar].de);
    star[*nstar].r[2]=sin(star[*nstar].de);

    int dummy;
    fscanf(mag6, "%d\t", &dummy);

    char sign;
    fscanf(mag6, "%g %g %g\t", &d,&m,&s);
    fscanf(mag6, "%c%g %g %g\t", &sign, &d, &m, &s);

    fscanf(mag6, "%lg", &star[*nstar].mag);

    if (star[*nstar].mag<MAXMAG)
      (*nstar)++;

    fgets(buf, 1024, mag6);
  }
  fclose(mag6);

  return star;
}
#endif

int floodfill(char *fill, unsigned short *pic,int w,int h,int x,int y,int threshold,double *cx, double *cy, double *sqx, double *sqy, double *xy, double *weight) {
  if (pic[x+w*y]<=threshold || fill[x+w*y])
    return 0;
  int cnt=1;
  fill[x+w*y]=1;
  double pw=((double)pic[x+w*y]-threshold);
  *cx+=(double)x*pw;
  *cy+=(double)y*pw;
  *sqx+=(double)x*(double)x*pw;
  *sqy+=(double)y*(double)y*pw;
  *xy+=(double)x*(double)y*pw;
  *weight+=pw;
  if (x<w-1)
    cnt+=floodfill(fill, pic, w, h, x+1, y, threshold, cx, cy, sqx, sqy, xy, weight);
  if (x>0)
    cnt+=floodfill(fill, pic, w, h, x-1, y, threshold, cx, cy, sqx, sqy, xy, weight);
  if (y>0)
    cnt+=floodfill(fill, pic, w, h, x, y-1, threshold, cx, cy, sqx, sqy, xy, weight);
  if (y<h-1)
    cnt+=floodfill(fill, pic, w, h, x, y+1, threshold, cx, cy, sqx, sqy, xy, weight);
  return cnt;
}


unsigned short *loadpic(FILE *picfile, int *_w, int *_h) {
  int w, h, d;
  char l1[80];
  fscanf(picfile, "%10s\n",l1);
  if (strcmp(l1,"P5") && strcmp(l1,"P6")) {
    fprintf(stderr, "Wrong file format\n");
    exit(1);
  }
  int c;
  for (;;) {
    char dummy[1024];
    c=fgetc(picfile);
    if (c=='#')
      fgets(dummy, 1024, picfile);
    else {
      ungetc(c, picfile);
      break;
    }
  }
  fscanf(picfile, "%d %d",&w, &h);
  *_w = w;
  *_h = h;
  fscanf(picfile, "%d",&d);
  unsigned short *pic = malloc(w*h*3*2);
  if (!pic) {
    fprintf(stderr, "Out of memory\n");
    return NULL;
  }
  if (d==255 && !strcmp(l1,"P5")) {
    fread(pic+w*h/2, w*sizeof(char), h, picfile);
    int i;
    for (i=0; i<w*h; i++) {
      pic[i]=256*((unsigned char *)pic)[i+w*h];
    }
  } else if (d==255 && !strcmp(l1,"P6")) {
    fread(pic, 3*w, h, picfile);
    int i;
    for (i=0; i<w*h; i++) {
      pic[i]=256*(((int)(((unsigned char *)pic)[3*i]) + (int)(((unsigned char *)pic)[3*i+1]) + (int)(((unsigned char *)pic)[3*i+2])))/3;
    }
  } else if (d==65535 && !strcmp(l1,"P5")) {
    fread(pic, w*sizeof(short), h, picfile);
  } else if (d==65535 && !strcmp(l1,"P6")) {
    fread(pic, w*sizeof(short)*3, h, picfile);
    int i;
    for (i=0; i<w*h; i++) {
      pic[i]=((pic[3*i] + pic[3*i+1] + pic[3*i+2])/3);
    }
  } else {
    fprintf(stderr, "Wrong file format\n");
    return NULL;
  }
  return pic;
}

struct streak *detect(unsigned short *pic, char *fill, int w, int h, double foclen, double pixpitch, double thresh, int *nstreak) {
  int hist[65536];
  int i,j,k,l;
  int discard=0, olddiscard;
  int threshold=0;
  int background=0;
  int iter=0;

  int *trig=malloc(sizeof(int)*w*h);

  do {

    olddiscard=discard;

    for (i=0; i<65536; i++) {
      hist[i]=0;
    }

    for (i=0; i<h; i++) {
      for (j=0; j<w; j++) {
	if (!fill[i*w+j])
	  hist[pic[i*w+j]]++;
      }
    }
    int sum=0;
    for (i=0; i<65536; i++) {
      sum+=hist[i];
      hist[i]=sum;
    }

  /*
#define THRESHOLD 0.9994
  //#define THRESHOLD 0.9999
#define BACKGROUND 0.99
  */

    int tscount=(double)hist[65535]*(1-thresh);
    int bgcount=(double)hist[65535]*(1-thresh*5); // 20
    
    for (i=0; i<65536; i++) {
      if (hist[i]<tscount)
	threshold=i;
      if (hist[i]<bgcount)
	background=i;
    }

#define BOX 64

    /*
      ...*....*
      .........
      *...*....
      ..*...*..
      ....*****
     */


    discard=0;
    memset(fill,0,w*h);

    int t=0;
    for (i=0; i<w; i++) {
      if (pic[i] > background)
	t++;
      trig[i]=t;
    }

    /*
      000111112
      .........
      *...*....
      ..*...*..
      ....*****
     */

    t=0;
    for (i=1; i<h; i++) {
      if (pic[w*i] > background)
	t++;
      trig[w*i]=t;
    }
    
    /*
      000111112
      0........
      1...*....
      1.*...*..
      1...*****
     */

    for (i=1; i<h; i++)
      for (j=1; j<w; j++) {
	trig[w*i+j]=trig[w*i+j-w]+
	  trig[w*i+j-1]-
	  trig[w*i+j-1-w];
	if (pic[w*i+j] > background)
	  trig[w*i+j]++;
      }

    /*
      000111112
      000111112
      111233334
      112344556
      11235689B
     */

    for (i=0; i<h-BOX; i++)
      for (j=0; j<w-BOX; j++) {
	int t;
	t=trig[(i+BOX)*w+j+BOX]-
	  trig[(i)*w+j+BOX]-
	  trig[(i+BOX)*w+j]+
	  trig[(i)*w+j];
	if (t>BOX*BOX/16)
	  trig[i*w+j]=1;
	else
	  trig[i*w+j]=0;
      }
    
    /*
      000000012
      000000012
      000010034
      112344556
      11235689B
     */

    for (i=0; i<h; i++)
      for (j=0; j<w; j++) {
	if (i>=h-BOX || j>=w-BOX)
	  trig[i*w+j]=0;
      }

    /*
      000000000
      000000000
      000010000
      000000000
      000000000
     */

    int cnt;
    for (i=0; i<h; i++) {
      cnt=0;
      for (j=0; j<w; j++) {
	if (trig[i*w+j])
	  cnt=BOX-1;
	if (cnt) {
	  trig[i*w+j]=1;
	  cnt--;
	}
      }
    }

    /*
      000000000
      000000000
      000011000
      000000000
      000000000
     */

    for (i=0; i<w; i++) {
      cnt=0;
      for (j=0; j<h; j++) {
	if (trig[j*w+i])
	  cnt=BOX-1;
	if (cnt) {
	  fill[j*w+i]=1;
	  //	  pic[j*w+i]=50000;
	  discard++;
	  cnt--;
	}
      }

    /*
      000000000
      000000000
      000011000
      000011000
      000000000
     */

    }

    iter++;
  } while (discard!=olddiscard && iter<3);
  free(trig);

#define MAXSTREAK 10000
  struct streak *orig = malloc(sizeof(*orig)*MAXSTREAK);
  *nstreak=0;

  for (i=1; i<h-1; i++) {
    for (j=2; j<w-3; j++) {
      if (pic[j+i*w]>threshold && !fill[j+i*w]) {
	double cx=0., cy=0., sqx=0., sqy=0., weight=0., xy=0.;
	int n=floodfill(fill, pic, w, h, j, i, background, &cx, &cy, &sqx, &sqy, &xy, &weight);
	if (n<2)
	  continue;
	cx/=weight;
	cy/=weight;
	sqx/=weight;
	sqy/=weight;
	sqx-=cx*cx;
	sqy-=cy*cy;
	xy/=weight;
	xy-=cx*cy;
	double angle = 0.5*atan2(xy*2,sqx-sqy);
	double v1 = cos(angle)*cos(angle)*sqx+sin(2*angle)*xy+sin(angle)*sin(angle)*sqy;
	double v2 = sin(angle)*sin(angle)*sqx-sin(2*angle)*xy+cos(angle)*cos(angle)*sqy;
	double len = sqrt((v1 - v2)*3.)*2;

#if 0
	//#define IMVIEW
#ifdef IMVIEW
	printf("%g %g 130 0 130 130 %g %g 0 3 CHR 134 133 141\n",
	       round(cy),
	       round(cx),
	       round(cx),
	       round(cy)
	       );
	printf("%g %g 130 0 130 130 %g %g 0 3 CHR 134 133 141\n",
	       round(round(cy)+3*len*sin(angle)),
	       round(round(cx)+3*len*cos(angle)),
	       round(round(cx)+3*len*cos(angle)),
	       round(round(cy)+3*len*sin(angle))
	       );
#else
	printf("%g %g %g %g %g\n",cx, cy, angle, len, weight);
#endif
#endif

	if (*nstreak<MAXSTREAK) {
	  orig[*nstreak].x = (cx - w/2)*pixpitch;
	  orig[*nstreak].y = (h/2 - cy)*pixpitch;
	  orig[*nstreak].angle = angle;
	  orig[*nstreak].len = len*pixpitch;
	  orig[*nstreak].intensity = weight;
	  orig[*nstreak].use = 1;
	  orig[*nstreak].star = -1;
	  double l=sqrt(foclen*foclen+orig[*nstreak].x*orig[*nstreak].x+orig[*nstreak].y*orig[*nstreak].y);
	  orig[*nstreak].r[0]=foclen/l;
	  orig[*nstreak].r[1]=orig[*nstreak].x/l;
	  orig[*nstreak].r[2]=orig[*nstreak].y/l;

	  double duv_x = cos(angle) * len*pixpitch;
	  double duv_y = sin(angle) * len*pixpitch;

	  orig[*nstreak].v[0]=-duv_x*orig[*nstreak].x*foclen/l/l/l - duv_y*orig[*nstreak].y*foclen/l/l/l;
	  orig[*nstreak].v[1]=duv_x/l - duv_x*orig[*nstreak].x*orig[*nstreak].x/l/l/l - duv_y*orig[*nstreak].x*orig[*nstreak].y/l/l/l;
	  orig[*nstreak].v[2]=duv_y/l - duv_y*orig[*nstreak].y*orig[*nstreak].y/l/l/l - duv_x*orig[*nstreak].x*orig[*nstreak].y/l/l/l;

	  (*nstreak)++;
	}
      }
    }
  }
  return orig;
}


void rotate(double M[3][3], double v[3], double r[3]) {
  r[0]=M[0][0]*v[0] + M[1][0]*v[1] + M[2][0]*v[2];
  r[1]=M[0][1]*v[0] + M[1][1]*v[1] + M[2][1]*v[2];
  r[2]=M[0][2]*v[0] + M[1][2]*v[1] + M[2][2]*v[2];
}

void rotateT(double M[3][3], double v[3], double r[3]) {
  r[0]=M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2];
  r[1]=M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2];
  r[2]=M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2];
}

// trans(A) * B

void mulT(double A[3][3], double B[3][3], double P[3][3]) {
  int i,j,k;
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      P[i][j]=0.;

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      for (k=0; k<3; k++)
	P[i][j] += A[k][i]*B[k][j];
}


int magcmp(const void *a, const void *b) {
  if (((struct star *)a)->mag > ((struct star *)b)->mag)
    return 1;
  else if (((struct star *)a)->mag < ((struct star *)b)->mag)
    return -1;
  return 0;
}

int intcmp(const void *a, const void *b) {
  if (((struct streak *)a)->intensity > ((struct streak *)b)->intensity)
    return -1;
  else if (((struct streak *)a)->intensity < ((struct streak *)b)->intensity)
    return 1;
  return 0;
}

/*
  refstreak: Identify suitable reference streaks
 */


void refstreak (struct streak *streak, int nstreak, int *r1, int *r2) {
  int i,j;

  double *dneighbour = malloc(sizeof(double) * nstreak);
  for (i=0; i<nstreak; i++) {
    dneighbour[i]=0.;
  }


  for (i=0; i<nstreak; i++) {
    double mind = INFINITY;
    if (streak[i].use) {
      for (j=0; j<nstreak; j++) {
	if (streak[j].use && i!=j) {
	  double dx = streak[i].r[0] - streak[j].r[0];
	  double dy = streak[i].r[1] - streak[j].r[1];
	  double dz = streak[i].r[2] - streak[j].r[2];
	  double d = dx*dx + dy*dy + dz*dz;
	  if (d < mind)
	    mind=d;
	}
      }
    }
    dneighbour[i] = sqrt(mind);
  }

  double bscore = 0.;

  for (i=0; i<nstreak; i++) {
    if (streak[i].use) {
      for (j=i+1; j<nstreak; j++) {
	if (streak[j].use) {
	  double dx = streak[i].r[0] - streak[j].r[0];
	  double dy = streak[i].r[1] - streak[j].r[1];
	  double dz = streak[i].r[2] - streak[j].r[2];
	  // far apart
	  // high intensity
	  // far from nearest neighbours

	  /*
	  double score = sqrt(dx*dx + dy*dy + dz*dz) *
	    dneighbour[i] * dneighbour[j] * streak[i].intensity *
	    streak[j].intensity;
	  */
	  //	  double score = sqrt(dx*dx + dy*dy + dz*dz);
	  double score = streak[i].intensity *
	    streak[j].intensity;
	  
	  if (score > bscore) {
	    bscore=score;
	    *r1=i;
	    *r2=j;
	  }
	}
      }
    }
  }
}

static inline void xproduct(double v1[3], double v2[3], double p[3]) {
  p[0] = v1[1] * v2[2] - v1[2] * v2[1];
  p[1] = v1[2] * v2[0] - v1[0] * v2[2];
  p[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

static inline void vdiff(double v1[3], double v2[3], double p[3]) {
  p[0] = v1[0] - v2[0];
  p[1] = v1[1] - v2[1];
  p[2] = v1[2] - v2[2];
}

static inline void vsum(double v1[3], double v2[3], double p[3]) {
  p[0] = v1[0] + v2[0];
  p[1] = v1[1] + v2[1];
  p[2] = v1[2] + v2[2];
}

static inline void vscale(double v1[3], double a) {
  v1[0]*=a;
  v1[1]*=a;
  v1[2]*=a;
}

static inline double dproduct(double v1[3], double v2[3]) {
  return
    v1[0] * v2[0] +
    v1[1] * v2[1] +
    v1[2] * v2[2];
}

static inline double sqlen(double v[3]) {
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

static inline void vnorm(double v1[3]) {
  double a = sqrt(sqlen(v1));
  vscale(v1, 1./a);
}

/* Innehlder feil? */
void rotmatrix(double M[3][3], double rot[3]) {
  double c[3]={cos(rot[0]), cos(rot[1]),cos(rot[2])};
  double s[3]={sin(rot[0]), sin(rot[1]),sin(rot[2])};
  M[0][0]=c[2]*c[1];
  M[0][1]=c[2]*s[1]*s[0]-s[2]*c[0];
  M[0][2]=c[2]*s[1]*c[0]+s[2]*s[0];
  M[1][0]=s[2]*c[1];
  M[1][1]=s[2]*s[1]*s[0]+c[2]*c[0];
  M[1][2]=s[2]*s[1]*c[0]-c[2]*s[0];
  M[2][0]=-s[2];
  M[2][1]=c[1]*s[0];
  M[2][2]=c[2]*c[0];
}

void matrixrot(double M[3][3], double rot[3]) {
  double c[3];
  double s[3];
  s[2]=-M[2][0];
  c[1]=M[1][0]/s[2];
  s[0]=M[2][1]/c[1];
  c[2]=M[0][0]/c[1];
  c[0]=M[2][2]/c[2];
  s[1]=(M[0][1]+s[2]*c[0])/c[2]/s[0];
  rot[0]=atan2(s[0],c[0]);
  rot[1]=atan2(s[1],c[1]);
  rot[2]=atan2(s[2],c[2]);
}  

void printmatrix(double M[3][3]) {
  printf("%g %g %g\n",M[0][0],M[0][1],M[0][2]);
  printf("%g %g %g\n",M[1][0],M[1][1],M[1][2]);
  printf("%g %g %g\n",M[2][0],M[2][1],M[2][2]);
}

static inline void proj(double a[3], double b[3], double c[3]) {
  memcpy(c, b, sizeof(double)*3);
  vscale(c, dproduct(a,b)/sqlen(b));
}

static inline void gramschmidt(double V[3][3], double E[3][3]) {
  double tmp[3];
  E[0][0]=V[0][0];
  E[0][1]=V[0][1];
  E[0][2]=V[0][2];

  E[1][0]=V[1][0];
  E[1][1]=V[1][1];
  E[1][2]=V[1][2];
  proj(V[1],E[0],tmp);
  E[1][0]-=tmp[0];
  E[1][1]-=tmp[1];
  E[1][2]-=tmp[2];

  E[2][0]=V[2][0];
  E[2][1]=V[2][1];
  E[2][2]=V[2][2];
  proj(V[2],E[0],tmp);
  E[2][0]-=tmp[0];
  E[2][1]-=tmp[1];
  E[2][2]-=tmp[2];
  proj(V[2],E[1],tmp);
  E[2][0]-=tmp[0];
  E[2][1]-=tmp[1];
  E[2][2]-=tmp[2];
  
  vnorm(E[0]);
  vnorm(E[1]);
  vnorm(E[2]);
}


int match(struct streak *streak, int nstreak, struct star *star, int nstar, int r1, int r2, double tol, double foclen, double pixpitch, int w, int h, double rot[3][3]) {
  double ctol = cos(tol);

  double M[3][3];
  vdiff(streak[r1].r, streak[r2].r, M[0]);
  double dr1r2 = sqrt(sqlen(M[0]));
  vsum(streak[r1].r, streak[r2].r, M[1]);
  xproduct(M[0], M[1], M[2]);
  vnorm(M[0]);
  vnorm(M[1]);
  vnorm(M[2]);
  int ret = 0;

  double (*streakr)[3] = malloc(sizeof(*streakr) * nstreak);
  double (*starr)[3] = malloc(sizeof(*starr) * nstar);
  int *starid = malloc(sizeof(int) * nstar);
  
  int i, j, swap;
  for (i=0; i<nstreak; i++) {
    rotateT(M, streak[i].r, streakr[i]);
  }

  /*
  printf("%g %g %g\n",streakr[r1][0], streakr[r1][1], streakr[r1][2]);
  printf("%g %g %g\n",streakr[r2][0], streakr[r2][1], streakr[r2][2]);
  */


  for (i=1; i<nstar && i<400; i++) 
    for (j=0; j<i; j++) {
      // Assume r1 is i
      // Assume r2 is j
      double N[3][3];
      vdiff(star[i].r, star[j].r, N[0]);
      double dij = sqrt(sqlen(N[0]));
      if (fabs(dij-dr1r2) > 2*tol) // Too close or too far apart.
	continue;
      vsum(star[i].r, star[j].r, N[1]);
      xproduct(N[0], N[1], N[2]);
      vnorm(N[0]);
      vnorm(N[1]);
      vnorm(N[2]);
      for (swap=0; swap<2; swap++) {
	if (swap==1) {
	  N[0][0]=-N[0][0];
	  N[0][1]=-N[0][1];
	  N[0][2]=-N[0][2];
	  N[2][0]=-N[2][0];
	  N[2][1]=-N[2][1];
	  N[2][2]=-N[2][2];
	}

	int k,l;

	for (k=0; k<nstar; k++)
	  star[k].streak = -1;
	for (k=0; k<nstreak; k++)
	  streak[k].star = -1;

	int nmiss=0;
#define MAXMISS 10
	double p=1.;

	mulT(N,M,rot);

	int imgstars=0;
	for (l=0; l<nstar; l++) {
	  double s[3];
	  //	  rotate(M, starr[imgstars], s);
	  rotate(rot, star[l].r, s);
	  s[1] *= foclen / s[0] / pixpitch;
	  s[2] *= foclen / s[0] / pixpitch;
	  s[1] += w/2;
          s[2] += h/2;
	  if (s[1]>=0 && s[1]<w && s[2]>=0 && s[2]<h) {
	    rotateT(N, star[l].r, starr[imgstars]);
	    starid[imgstars]=l;
	    imgstars++;
	  }
	}

	double totint=0.;
	int bs;
	int nmatch=0;

	for (k=0; k<nstreak; k++)
	  if (streak[k].use) {
	    double bp = 0;
	    for (l=0; l<imgstars && l<10*k+nstreak/10; l++)
	      if (star[starid[l]].streak==-1) {
		double dp=dproduct(starr[l], streakr[k]);
		if (dp >= bp) {
		  bp = dp;
		  bs = l;
		}
	      }
	    if (bp < ctol) {
	      nmiss++;
	      if (nmiss>MAXMISS)
		break;
	    } else {
	      nmatch++;
	      p *= bp;
	      star[starid[bs]].streak = k;
	      streak[k].star = starid[bs];
	      totint += pow(10, -star[bs].mag);
	    }
	  }
	if (nmatch>nstreak/3) {
	  //	  printf("k=%d nmiss=%d acos=%g mag=%g\n",k, nmiss, 180./PI*acos(p), -log(totint)/log(10));
	  ret = 1;
	  goto exitfunc;
	}
      }
    }

 exitfunc:
  free(streakr);
  free(starr);
  free(starid);

  return ret;
}

void printstreak(struct streak *s) {
  //  printf("x=%g y=%g len=%g int=%g\n", s->x, s->y, s->len, s->intensity);
  printf("x=%g y=%g z=%g len=%g int=%g\n", s->r[0], s->r[1], s->r[2], s->len, s->intensity);
}

void paramcorr(double *x, double *y, double *param, int nparam) {
  int i=2;
  
  if (nparam==0)
    return;

  double rx=*x + param[0];
  double ry=*y + param[1];

  /*
  X = xm + xo;
  Y = ym + yo;
  R2 = X*X + Y*Y;
  c = r + X * rx + Y * ry + R2 * r3 + R2*X * r3x + R2*Y * r3y + R2*R2 * r5 ...
  xc = xm + X*c;
  yc = ym + Y*c;
  */

  double rn=1.;
  double c=0.;
  
  for (;;) {
    c+=rn * param[i];
    i++;
    if (i>=nparam)
      break;
    c+=rn*rx * param[i];
    i++;
    if (i>=nparam)
      break;
    c+=rn*ry * param[i];
    i++;
    if (i>=nparam)
      break;
    rn*=rx*rx+ry*ry;
  }

  *x+=c*rx;
  *y+=c*ry;
}


int rematch(struct streak *streak, int nstreak, struct star *star, int nstar, int r1, int r2, double sqtol, double foclen, double pixpitch, int w, int h, double rot[3][3], double *param, int nparam) {
  int k,l;

  for (l=0; l<nstar; l++) {
    star[l].streak=-1;
  }
  for (l=0; l<nstreak; l++) {
    streak[l].star=-1;
  }

  int imgstars=0;

  int *starid = malloc(sizeof(int)*nstar);
  double *sx = malloc(sizeof(double)*nstar);
  double *sy = malloc(sizeof(double)*nstar);
  for (l=0; l<nstar; l++) {
    double r[3];
    rotate(rot, star[l].r,r);
    double x = r[1]/r[0]*foclen/pixpitch;
    double y = r[2]/r[0]*foclen/pixpitch;
    if (x >= -w/2 && x < w/2 &&
	y >= -h/2 && y < h/2) {
      starid[imgstars] = l;
      sx[imgstars]=x*pixpitch/foclen;
      sy[imgstars]=y*pixpitch/foclen;
      imgstars++;
    }
  }

  int nmatch=0;

  for (k=0; k<nstreak; k++)
    if (streak[k].use) {
      double px = streak[k].x/foclen;
      double py = streak[k].y/foclen;
      paramcorr(&px, &py, param, nparam);
      double bp = INFINITY;
      int bl;
      for (l=0; l<imgstars && l<10*k+nstreak/10; l++)
	if (star[starid[l]].streak==-1) {	
	  double dp = (sx[l]-px)*(sx[l]-px)+(sy[l]-py)*(sy[l]-py);
	  if (dp < bp) {
	    bp = dp;
	    bl = l;
	  }
	}
      if (bp < sqtol) {
	star[starid[bl]].streak = k;
	streak[k].star = starid[bl];
	nmatch++;
      }
    }

  free(starid);
  free(sx);
  free(sy);
  
  return nmatch;
}

struct nmv {
  double rot[3][3];
  double *param;
};

struct params {
  struct streak *streak;
  int nstreak;
  double foclen;
  //  double factor;
  struct star *star;
  int nstar;
  int nparam;
};

void freenmv(struct nmv *X) {
  free(X->param);
  free(X);
}

double calcE(struct nmv *X, struct params *p) {
  double E = 0;
  int k;
  int nuse=0;
  for (k=0; k<p->nstreak; k++)
    if (p->streak[k].star!=-1) {
      double r[3];
      rotate(X->rot, p->star[p->streak[k].star].r,r);
      double sx = r[1]/r[0];
      double sy = r[2]/r[0];
      
      double px = p->streak[k].x/p->foclen;
      double py = p->streak[k].y/p->foclen;

      paramcorr(&px, &py, X->param, p->nparam);

      //      X->E += (sx-rx)*(sx-rx)+(sy-ry)*(sy-ry);
      //      X->E += (sx*sx*ry*ry - 2*sx*sy*rx*ry + sy*sy*rx*rx);
      E +=  (px-sx)*(px-sx)+(py-sy)*(py-sy);
      nuse++;
    }
  E /= nuse;
  return E;
}

struct nmv *expcon(struct nmv **V, int n, struct nmv *X, double alpha, struct params *p) {
  int j,k,l;
  double tmp[3][3];
  struct nmv *test=malloc(sizeof(*test));
  memcpy(test, X, sizeof(*test));
  test->param=malloc(sizeof(double)*p->nparam);

  // 6 = adjustp+1 = nparam+3

  for (j=0; j<3; j++)
    for (k=0; k<3; k++)
      test->rot[j][k]=0.;
  for (j=0; j<p->nparam; j++)
    test->param[j]=0.;
  
  for (l=0; l<n; l++) {
    for (j=0; j<3; j++)
      for (k=0; k<3; k++)
	test->rot[j][k]+=V[l]->rot[j][k]/((double)n)*(1+alpha);
    for (j=0; j<p->nparam; j++)
      test->param[j]+=V[l]->param[j]/((double)n)*(1+alpha);
  }
  
  for (j=0; j<3; j++) {
    for (k=0; k<3; k++)
      test->rot[j][k]-=X->rot[j][k]*alpha;
    //    vnorm(test->rot[j]);
  }
  for (j=0; j<p->nparam; j++)
    test->param[j]-=X->param[j]*alpha;

  gramschmidt(test->rot,tmp);
  memcpy(test->rot,tmp,sizeof(double)*9);

  return test;
}

double adjust(struct streak *streak, int nstreak, double foclen, double pixpitch, int w, int h, struct star *star, int nstar, double rot[3][3], double *param, int nparam) {
  struct nmv **V;
  V=malloc(sizeof(*V)*(nparam+4));

  int i,j,k,l;

  for (i=0; i<nparam+4; i++) {
    V[i]=malloc(sizeof(struct nmv));
    V[i]->param = malloc(sizeof(double)*nparam);
  }

  struct params params;
  params.streak = streak;
  params.nstreak = nstreak;
  params.foclen = foclen;
  //  params.factor = 1./foclen;
  params.star = star;
  params.nstar = nstar;
  params.nparam = nparam;

  double tmp[3][3];

  double a=1e-6;
  for (i=0; i<nparam+4; i++) {
    memcpy(V[i]->rot,rot,sizeof(double)*9);
    memcpy(V[i]->param,param,sizeof(double)*nparam);
  }

  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) { 
      V[i]->rot[j][(i+1)%3] += a*rot[j][(i+2)%3];
      V[i]->rot[j][(i+2)%3] -= a*rot[j][(i+1)%3];
    }
    double tmp[3][3];
    gramschmidt(V[i]->rot,tmp);
    memcpy(V[i]->rot,tmp,sizeof(double)*9);
  }

  for (i=0; i<nparam; i++)
    V[i+3]->param[i] += a;
  

  double retval = amoeba((void**)V, (amoeba_move_t *)expcon, (amoeba_free_t *)freenmv, (amoeba_E_t *)calcE, nparam+4, &params);

  memcpy(param, V[0]->param, sizeof(double)*nparam);
  memcpy(rot, V[0]->rot, sizeof(double)*9);

  for (i=0; i<nparam+4; i++) {
    freenmv(V[i]);
  }

  free(V);
  return retval;
}

void printcalib(struct streak *streak, int nstreak, struct star *star, int nstar, double foclen, double pixpitch, int w, int h, double rot[3][3], double *param, int nparam) {
  /*
  int i;
  for (i=0; i<nstreak; i++)
    if (streak[i].star!=-1) {
      double r[3];
      rotate(rot, star[streak[i].star].r,r);
      double ex = streak[i].x + offset[0];
      double ey = streak[i].y + offset[1];
      double sx = r[1]/r[0]*foclen-ex;
      double sy = r[2]/r[0]*foclen-ey;
      double el = sqrt(ex*ex+ey*ey);
      double exn = ex/el;
      double eyn = ey/el;
      double sr = sx*exn+sy*eyn;
      double st = sx*eyn-sy*exn;
      double corr = paramcorr(ex/pixpitch/(double)w*2, ey/pixpitch/(double)w*2, param, nparam);
      printf("%g %g %g %g %g\n", ex/pixpitch/(double)w*2, ey/pixpitch/(double)w*2, sr/pixpitch, st/pixpitch, corr);
    }
  */  
}

double curveerr(double *param, double *basisfunc, double *target, int nparam, int nstreak) {
  int i,j;
  double e=0.;
  for (i=0; i<nstreak; i++) {
    double le=0.;
    for (j=0; j<nparam; j++) {
      le += basisfunc[i*nparam+j]*param[j];
    }
    le -= target[i];
    e += le*le;
  }
  return e;
}

int doublecmp(const void *a, const void *b) {
  if (*((double *)a) > *((double *)b))
    return 1;
  else if (*((double *)a) < *((double *)b))
    return -1;
  return 0;
}

void fitcont(double *simplex, int nparam, double *ret, double a) {
  int i,j;
  for (i=0; i<nparam+1; i++)
    ret[i]=0.;
  for (i=0; i<nparam; i++)
    for (j=1; j<nparam+1; j++)
      ret[j] += simplex[i*(nparam+1)+j];
  for (i=1; i<nparam+1; i++)
    ret[i] = ret[i]*(1+a)/nparam-a*simplex[nparam*(nparam+1)+i];
}


#define MARKINT 32768
#define mark(x,r) ((x)=(x)>65535-MARKINT*(r)?65535:(x)+MARKINT*(r))

void cross(double x,double y, unsigned short *pic, int w, int h, double a, int f, int t) {
  int j;
  for (j=f; j<t; j++) {
    if (x+j<w)
      mark(pic[(int)x+((int)y)*w+j], a*(1-(y-floor(y))));
    if (x+j<w && y<h-1)
      mark(pic[(int)x+((int)y)*w+j+w], a*(y-floor(y)));
    if (y+j<h)
      mark(pic[(int)x+((int)y)*w+j*w], a*(1-(x-floor(x))));
    if (y+j<h && x<w-1)
      mark(pic[(int)x+((int)y)*w+j*w+1], a*(x-floor(x)));
    if (x-j>=0)
      mark(pic[(int)x+((int)y)*w-j], a*(1-(y-floor(y))));
    if (x-j>=0 && y<h-1)
      mark(pic[(int)x+((int)y)*w-j+w], a*(y-floor(y)));
    if (y-j>=0)
      mark(pic[(int)x+((int)y)*w-j*w], a*(1-(x-floor(x))));
    if (y-j>=0 && x<w-1)
      mark(pic[(int)x+((int)y)*w-j*w+1], a*(x-floor(x)));
  }
}

void writepic(char *outputfile, unsigned short *pic, unsigned char *fill, int w, int h, double foclen, double pixpitch, struct streak *streak, int nstreak, struct star *star, int nstar, double rot[3][3], double *param, int nparam) {
  // Write picture
  FILE *out;
  if (strcmp(outputfile,"-")==0)
    out=stdout;
  else
    out=fopen(outputfile,"wb");
  if (!out) {
    perror("Could not open output file");
    return;
  }

  int i,j;

  for (i=0; i<h; i++)
    for (j=0; j<w; j++) {
      if (!fill[i*w+j] &&
	  (fill[i*w+j+1] ||
	   fill[i*w+j-1] ||
	   fill[i*w+j+w] ||
	   fill[i*w+j-w]))
	mark(pic[i*w+j],1);
    }

  for (i=0; i<nstreak; i++) {
    double x=streak[i].x/pixpitch+w/2;
    double y=h/2-streak[i].y/pixpitch;
    cross(x,y,pic, w,h,0.5,5,30);
    if (streak[i].star!=-1) {
      double a;
      for (a=0; a<2*PI; a+=1/150.) {
	double xr=x+15*cos(a);
	double yr=y+15*sin(a);
	if (xr>=0 && xr<w && yr>=0 && yr<h)
	  mark(pic[(int)xr+((int)yr)*w], 0.1);
      }

      if (param) {
	double r[3];
	rotate(rot, star[streak[i].star].r,r);
	double xs=r[1]/r[0];
	double ys=r[2]/r[0];

	double xp=streak[i].x/foclen;
	double yp=streak[i].y/foclen;
	paramcorr(&xp, &yp,param,nparam);
	
	double re;
	xp-=xs;
	yp-=ys;

	re=sqrt(xp*xp+yp*yp);
	xp/=re;
	yp/=re;

	re*=foclen/pixpitch;//w/2;
	
	for (a=15; a<15+re*15; a+=0.5) {
	  double xr=x+a*xp;
	  double yr=y-a*yp;
	  if (xr>=0 && xr<w && yr>=0 && yr<h)
	    mark(pic[(int)xr+((int)yr)*w], 0.5);
	}
      }
    }
  }

  if (param) {
    //    cross(w/2-offset[0]/pixpitch,h/2+offset[1]/pixpitch, pic, w,h,1,1,50);

    for (i=0; i<h; i++)
      for (j=0; j<w; j++) {
	double x=(double)(j-w/2)*pixpitch/foclen;///(w/2.);
	double y=(double)(h/2-i)*pixpitch/foclen;///(w/2.);
	double px=x;
	double py=y;
	paramcorr(&px,&py,param,nparam);
	double v=sqrt((px-x)*(px-x)+(py-y)*(py-y))*foclen/pixpitch;//*w/2;
	if (v-floor(v) < 0.1)
	  mark(pic[i*w+j], 10*(0.1-(v-floor(v))));
	if (v-floor(v) > 0.99)
	  mark(pic[i*w+j], 100*(v-floor(v)-0.99));
      }
  }
  
  fprintf(out,"P5\n%d %d\n255\n",w,h);

  for (i=0; i<w*h; i++) {
    char v=(pic[i]/256);
    fwrite(&v, 1, 1, out);
  }
  fclose(out);
}

int main(int argc, char **argv) {

  static struct option loptions[] = {
    {"help",0,0,'h'},
    {"foclen",1,0,'f'},
    {"pitch",1,0,'p'},
    {"width",1,0,'w'},
    {"nparam",1,0,'n'},
    {"threshold",1,0,'t'},
    {"quiet",0,0,'q'},
    {"output",0,0,'o'},
    {"catalogue",0,0,'c'}
  };
  
  double foclen=NAN;
  double pixpitch=NAN;
  double dwidth=NAN;
  double threshold=5e-4;
  int nparam=11;
  int quiet=0;
  char *cataloguename=0;
  char *outputfile=0;

  if (argc==1) {
    fprintf(stderr,"%s -h for options\n",argv[0]);
    exit(1);
  }

  while(1) {
    int opt = getopt_long (argc, argv, "hf:p:w:n:x:y:qo:c:t:", loptions, NULL);
    if (opt==-1) break;
    switch(opt) {
    case 'h':
      printf("Usage: starcalib [OPTION] file.pnm\n");
      printf("\n");
      printf("Calibrates a camera and optics by matching an image of the night sky against a star\n");
      printf("catalogue\n");
      printf("\n");
      printf("Options:\n");
      printf("  -h, --help         Print this help\n");
      printf("  -f, --foclen=MM    Focal length in mm (required)\n");
      printf("  -p, --pitch=UM     Detector pixel pitch in microns (-p or -w required)\n");
      printf("  -w, --width=MM     Effective detector width in mm (-p or -w required)\n");
      printf("  -n, --nparam=n     Number of parameters to use in model, must be >= 3\n");
      printf("                     (default = 11)\n");
      printf("  -t, --threshold=T  Detection threshold, as ratio of star pixels to all pixels\n");
      printf("                     (default = 5e-4)\n");
      //      printf("  -c, --catalogue=F  Name of file containing star catalogue\n");
      printf("  -o, --output=F     Name of output image file\n");
      printf("Output:\n");
      printf("  xo yo r rx ry r3 r3x r3y r5 r5x r5y ...\n");
      printf("\n");
      printf("Corrected pixel positions (xc, yc) can be calculated from measued positions (xm, ym)\n");
      printf("using the calibration parameters. Coordinates should be translated such that (0,0)\n");
      printf("is at the center of the image and scaled so that the increment from one pixel to the\n");
      printf("next is pitch/foclen (see options). Positive y-axis points upwards. Sample code:\n");
      printf("\n");
      printf("  X = xm + xo;\n");
      printf("  Y = ym + yo;\n");
      printf("  R2 = X*X + Y*Y;\n");
      printf("  c = r + X * rx + Y * ry + R2 * r3 + R2*X * r3x + R2*Y * r3y + R2*R2 * r5 ...\n");
      printf("  xc = xm + X*c;\n");
      printf("  yc = ym + Y*c;\n");
      exit(0);
    case 'f':
      foclen = atof(optarg)/1000.;
      break;
    case 'p':
      pixpitch = atof(optarg)/1000000.;
      break;
    case 'w':
      dwidth = atof(optarg)/1000.;
      break;
    case 'n':
      nparam = atoi(optarg);
      break;
    case 'q':
      quiet = 1;
      break;
    case 't':
      threshold = atof(optarg);
      break;
    case 'o':
      outputfile = optarg;
      break;
    case 'c':
      cataloguename = optarg;
      break;
    default:
      exit(1);
    }
  }

  int nstreak;
  //  int nstar;
  int w,h;

  FILE *input;

  if (nparam<3) {
    fprintf(stderr,"nparam must be >=3\n");//,(nparam+1)/3*3+1);
    exit(1);
  }
  if (isnan(foclen)) {
    fprintf(stderr,"-f is required\n");
    exit(1);
  }
  if (isnan(dwidth) && isnan(pixpitch)) {
    fprintf(stderr,"-p or -w is required\n");
    exit(1);
  }
  if (optind == argc || strcmp(argv[optind],"-")==0)
    input=stdin;
  else {
    input=fopen(argv[optind],"rb");
    if (!input) {
      perror("Could not open file");
      exit(1);
    }
  }
  if (!quiet)
    fprintf(stderr,"Loading image...\n");
  unsigned short *pic = loadpic(input, &w, &h);
  if (!pic)
    exit(1);
  if (isnan(pixpitch))
    pixpitch = dwidth / (double)w;
  char *fill=calloc(w,h);
  struct streak *streak = detect(pic, fill, w, h, foclen, pixpitch, threshold, &nstreak);
  if (nstreak<10) {
    fprintf(stderr, "Found only %d streaks - too few to continue\n", nstreak);
    exit(1);
  }
  //  printf("nstreak=%d\n",nstreak);
  /*
  FILE *mag8;
  if (cataloguename)
    mag8 = fopen(cataloguename,"rb");
  else {
    mag8 = fopen("/usr/lib/hipp8.tsv","rb");
    if (!mag8)
      mag8 = fopen("/usr/local/lib/hipp8.tsv","rb");
    if (!mag8)
      mag8 = fopen("hipp8.tsv","rb");
  }

  if (!mag8) {
    fprintf(stderr, "Unable to open star catalogue file\n");
    exit(1);
  }
  */

  //  struct star *star = loadstar(mag8, &nstar);
  qsort(star, nstar, sizeof(*star), magcmp);
  qsort(streak, nstreak, sizeof(*streak), intcmp);

  double rot[3][3];

  int r1=-1;
  int r2=-1;

  //    writepic(outputfile, pic, w, h, foclen, pixpitch, streak, nstreak, star, nstar, rot, 0, 0, 0);
  //    exit(0);
  if (!quiet)
    fprintf(stderr,"Searching star catalogue...\n");
  int i;
  for (i=0; i<4; i++) {
    refstreak(streak, nstreak, &r1, &r2);

    if (match(streak, nstreak, star, nstar, r1, r2, 2.*PI/180, foclen, pixpitch, w, h, rot))
      break;
    streak[r1].use=0;
    streak[r2].use=0;
  }
  if (i==4) {
    fprintf(stderr,"No match in star catalogue\n");
    if (outputfile) {
      if (!quiet)
	fprintf(stderr,"Writing output image...\n");
      writepic(outputfile, pic, (unsigned char *)fill, w, h, foclen, pixpitch, streak, nstreak, star, nstar, rot, 0, 0);
    }
    exit(1);
  }

  double *param=malloc(sizeof(double)*nparam);
  for (i=0; i<nparam; i++)
    param[i]=0.;
  //  printmatrix(rot);
  double sqtol;
  double convcrit=0.001*2./w;
  double olderr=INFINITY;
  int oldmatch=0;
  int nstag=0;
  double mintol=1*2./w;
  mintol*=mintol;
  for (;;) {
    double sqerr=adjust(streak, nstreak, foclen, pixpitch, w, h, star, nstar, rot, param, nparam);
    //    double sqerr=0.01*0.01;
    sqtol = sqerr+mintol;
    int nmatch=rematch(streak, nstreak, star, nstar, r1, r2, sqtol, foclen, pixpitch, w, h, rot, param, nparam);
    if (!quiet)
      fprintf(stderr,"Fitting %d stars, std dev=%g pixels\n",nmatch,sqrt(sqerr)*foclen/pixpitch);//w/2.);
    if (oldmatch == nmatch && fabs(olderr - sqerr) < convcrit)
      nstag++;
    else
      nstag=0;
    if (nstag>3)
      break;
    olderr = sqerr;
    oldmatch = nmatch;
  }
  //  printmatrix(rot);

  if (outputfile) {
    if (!quiet)
      fprintf(stderr,"Writing output image...\n");
    writepic(outputfile, pic, (unsigned char *)fill, w, h, foclen, pixpitch, streak, nstreak, star, nstar, rot, param, nparam);
  }

  //  double rots[3];
  //  matrixrot(rot, rots);
  //  printf("A=%g B=%g C=%g\n",rots[0]*180/PI, rots[1]*180/PI, rots[2]*180/PI);

  for (i=0; i<nparam; i++)
    printf("%g ",param[i]);
  printf("\n");

  free(fill);

  //  printcalib(streak, nstreak, star, nstar, foclen, pixpitch, w, h, rot, param, nparam);

  //  printf("y=%g z=%g\n", asin(v[2][2]), atan2(v[1][2], v[0][2]));
  //  orient(streak, nstreak, star, nstar, v, w, h, pixpitch, foclen);

}
