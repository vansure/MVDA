#include <R.h>
#include <stdlib.h>
#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#define EPS 1e-4
void SOM_online(double *data, double *codes, double *nhbrdist,
		double *alphas, double *radii, double *changes,
		Sint *pn, Sint *pp, Sint *pncodes, Sint *prlen)
{
  int n = *pn, p = *pp, ncodes = *pncodes, rlen = *prlen;
  int cd, i, j, k, l, nearest, niter, nind;
  double dm, dist, tmp, alpha, threshold;

  double mean_cd = 0.0;         double sum_cd = 0.0;
  double sum_obj = 0.0;         double mean_obj = 0.0;
  double cov = 0.0;             double sum_cov = 0.0;
  double sum_var_obj = 0.0;     double sum_var_cd = 0.0;
  double var_obj = 0.0;         double var_cd = 0.0;
  double pearson = 0.0;

  RANDIN; /*R is column major*/

  niter = rlen * n;

  for (k = 0; k < niter; k++) {
    i = (int)(n * UNIF);

	sum_obj=0.0;            mean_obj = 0.0;
	sum_var_obj=0.0;        var_obj=0.0;
	/* media dell'oggetto i-esimo */
	for(j=0; j<p; j++){ sum_obj +=data[i+j*n]; }
	mean_obj = sum_obj / p;

	for(j=0; j<p; j++){ /*varianza dell'oggetto i-esimo */
		sum_var_obj += (data[i+j*n]-mean_obj) * (data[i+j*n]-mean_obj);
	}

	var_obj = sum_var_obj / (p-1);

    /* find the nearest code 'near' */
        nind = 0; dm = DOUBLE_XMAX; nearest = -1;
        for (cd = 0; cd < ncodes; cd++) {
          sum_cd=0.0;           mean_cd =0.0;
	  sum_cov = 0.0;        sum_var_cd =0.0;
	  cov =0.0;             var_cd =0.0;
	  pearson=0.0;
	   /* media del codes cd-esimo */
	  for(j = 0; j<p; j++){ sum_cd += codes[cd +j * ncodes];}
	  mean_cd = sum_cd /p;

	  for(j=0; j<p; j++){
            sum_cov += (data[i+j*n]-mean_obj)*(codes[cd+j*ncodes]-mean_cd);
            sum_var_cd += (codes[cd+j*ncodes]-mean_cd) * (codes[cd+j*ncodes]-mean_cd);
	  }
	  cov = sum_cov / (p-1); /*cov of codes and object i*/
	  var_cd  = sum_var_cd / (p-1);  /* var of codes cd */

      pearson = cov / (var_obj * var_cd); 
      /* pearson correlation between object i and codes cd */

      dist = 1-pearson; /* dissimilarity between two objects */

      if (dist <= dm * (1 + EPS)) {
        if (dist < dm * (1 - EPS)) {
            nind = 0;
            nearest = cd;
	}else {
            if(++nind * UNIF < 1.0) nearest = cd;
        }
	dm = dist;
      }
    } //end for each codes cd

    if (nearest < 0)
      error("No nearest neighbour found...");

    /* update all codes within threshold of 'nearest'. Linear decrease
       for both radius and learning parameter. */
    threshold = radii[0] - (radii[0] - radii[1]) * (double)k/(double)niter;
    if (threshold < 1.0) threshold = 0.5;
    alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)k/(double)niter;

    l = (int)(k/n);

    for (cd = 0; cd < ncodes; cd++) {
      if(nhbrdist[cd + ncodes*nearest] > threshold) continue;

      for(j = 0; j < p; j++) {
	tmp = data[i + j*n] - codes[cd + j*ncodes];
	codes[cd + j*ncodes] += tmp * alpha;

	if (cd == nearest) changes[l] += tmp * tmp;
      }
    }
  } //end for iterazioni

  for (l = 0; l < rlen; l++) {
    /* mean difference per variable per object */
    changes[l] = sqrt(changes[l]/p)/n;
  }

  RANDOUT;
}
