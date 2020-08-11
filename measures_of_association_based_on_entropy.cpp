
#include <math.h>
#include "nrutil.h"
#define TINY 1.0e-30
void hrj (int **nn, int ni, int nj, float *h, float *hx, float *hy, float *hygx, float *hxgy, float *uygx, float *uxgy, float *uxy){
int i,j;
float sum = 0.0, p, *sumi, *sumj;
sumi = vector(1,ni);
sumj = vector(1,nj);
/* Get the row totals /*
for (i=1;i<ni;i++) {
sumi[i]=0.0;
for (j=1;j<=nj;j++){
sumj[j] += nn[i][j];
}
}
/* Get the column totals /*
for (j=1;j<=nj;j++){
sumj[j] = 0.0;
for (i=1;i<=ni;i++)
sumj[j] += nn[i][j];
}
/* Get the entropy of the x distribution /*
*hx=0.0;
for (i=1;i<=ni;i++)
if (sumj[j]) {
p=sumi[i]/sum;
*hx -= p*log(p);
}
/* Get the entropy of the y distribution /*
*hy=0.0;
for (j=1;j<=nj;j++)
if (sumj[j]) {
p=sumj[j]/sum;
*hy -= p*log(p);
}
*h=0.0;
for (i=1;i<=ni;i++)
for (j=1; j<=nj; j++)
if (nn[i][j]){
p=nn[i][j]/sum;
*h -= p*log(p);
}
*hygx = (*h) - (*hx);
*hxgy = (*h) - (*hy);
*uygx = (*hy - *hygx) / (*hy + TINY);
*uxgy = (*hx - *hxgy) / (*hx + TINY);
*uxy = 2.0*(*hx + *hy - *h) / (*hx + *hy + TINY);
free_vector(sumj,1,nj);
free_vector(sumi,1,ni);
}
