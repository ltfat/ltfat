#include "ltfatresample.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

/* This is actual structture used to hold info between consecutive blocks */
struct resample_plan_struct
{
   /* target to source sampking rates ratio */
   const double ratio;
   /* Type of polynomial interpolation */
   const resample_type restype;
   /* First sample in block global index */
   size_t inPos;
   /* Buffer for holding overlaps, length depend on interpoaltion technique */
   SAMPLE* overlap;
   /* Function pointer to the block interpolation functon */
   resample_error (*executer)(const resample_plan,
                              const SAMPLE*, const size_t,
                              SAMPLE*);
   /* Filters */
   EMQFfilters ef;

};

resample_plan
resample_init(const resample_type restype, const double ratio)
{
   SAMPLE* overlap;
   struct resample_plan_struct rp = {ratio, restype, 0, NULL, NULL, NULL };
   resample_plan rpret = malloc(sizeof * rpret);


   if (restype == LINEAR)
   {
      overlap = calloc(1, sizeof(SAMPLE));
      rp.executer = &resample_execute_linear;
   }
   else if (restype == LAGRANGE)
   {
      overlap = calloc(5, sizeof(SAMPLE));
      rp.executer = &resample_execute_lagrange;
   }
   else if (restype == HERMITE)
   {
      overlap = calloc(5, sizeof(SAMPLE));
      rp.executer = &resample_execute_hermite;
   }

   rp.overlap = overlap;
   /* When subsampling, do the antialiasing filtering. */
   if (ratio < 0.95)
   {
      rp.ef = emqffilters_init(ratio / 2.0);
   }


   memcpy(rpret, &rp, sizeof * rpret);
   return rpret;
}

/*
 *
 * */
resample_error
resample_execute(const resample_plan rp,
                 const SAMPLE* in, const size_t Lin,
                 SAMPLE* out)
{
   const SAMPLE* tmpIn = in;
   /* Do filtering if initialized */
   if (rp->ef)
   {

      tmpIn = malloc(Lin * sizeof * tmpIn);
      emqffilters_dofilter(rp->ef, in, Lin, (SAMPLE*) tmpIn);

   }

   /* Execute the computation */
   resample_error status = rp->executer(rp, tmpIn, Lin, out);
   /* Advance the "pointer" in data stream */
   resample_advanceby(rp, Lin);

   if (rp->ef)
   {
      free((SAMPLE*)tmpIn);
   }

   return status;
}

size_t
resample_nextoutlen(const resample_plan rp, size_t Lin)
{
   size_t retval = 0;
   const double outSpos = ceil( (rp->inPos) * rp->ratio );
   const double outEpos = ceil( (rp->inPos + Lin) * rp->ratio );
   retval =  (size_t)( outEpos - outSpos);

   return retval;
}

/*  Be carefull, rp is effectivelly a double pointer. */
void
resample_done(resample_plan *rp)
{
   free(*rp);
   *rp = NULL;
}


void
resample_advanceby(const resample_plan rp, size_t Lin)
{
   rp->inPos += Lin;
}


/* INTERNALS */
resample_error
resample_execute_linear(const resample_plan rp,
                        const SAMPLE* in, const size_t Lin,
                        SAMPLE* out)
{
   double truepos;
   ptrdiff_t lowpos, highpos;
   size_t ii, iiThre;
   /* Just to avoid division further*/
   const double oneOverRatio = 1.0 / rp->ratio;
   /* Starting position in the output stream */
   const double outSpos = ceil( (rp->inPos) * rp->ratio );
   /* How many samples will this routine produce */
   size_t Lout = resample_nextoutlen(rp, Lin);

   /* First handle all samples which need overlap */
   iiThre = floor((rp->inPos + 1.0 ) * rp->ratio - outSpos) + 1;

   for (ii = 0; ii < iiThre; ii++)
   {
      truepos = (ii + outSpos) * oneOverRatio - rp->inPos;
      out[ii] = *(rp->overlap) * (1.0 - truepos) + *in * truepos;
   }

   /* All the rest */
   for (ii = iiThre; ii < Lout; ii++)
   {
      truepos = (ii + outSpos) * oneOverRatio - rp->inPos;
      highpos = ceil(truepos);
      lowpos = highpos - 1;

      out[ii] = in[lowpos - 1] * (highpos - truepos) +
                in[lowpos] * (truepos - lowpos);

   }

   /* Store last input sample as overlap */
   *(rp->overlap) = in[Lin - 1];

}

resample_error
resample_execute_lagrange(const resample_plan rp,
                          const SAMPLE* in, const size_t Lin,
                          SAMPLE* out)
{
   double truepos, x;
   ptrdiff_t highpos;
   size_t ii, zz, *iiThre;
   SAMPLE* buf;

   /* Just to avoid division further*/
   const double oneOverRatio = 1.0 / rp->ratio;
   /* Starting position in the output stream */
   const double outSpos = ceil( (rp->inPos) * rp->ratio );
   /* How many samples will this routine produce */
   size_t Lout = resample_nextoutlen(rp, Lin);

   /* 6 thresholds plus leading zero */
   iiThre = calloc(7, sizeof * iiThre);
   buf = calloc(6, sizeof * buf);

   /* First handle all samples which need overlap */
   for (ii = 1; ii < 7; ii++)
   {
      iiThre[ii] = floor((rp->inPos + ((double) ii) ) * rp->ratio - outSpos) + 1;
   }


   for (zz = 0; zz < 6; zz++)
   {
      for (ii = iiThre[zz]; ii < iiThre[zz + 1]; ii++)
      {
         truepos = (ii + outSpos) * oneOverRatio - rp->inPos;
         x = truepos - zz;

         memcpy(buf, rp->overlap + zz, (5 - zz)*sizeof * buf);
         memcpy(buf + (5 - zz), in, (zz + 1)*sizeof * buf );
         out[ii] = lagrange_interp(x, buf);
      }
   }

   for (ii = iiThre[6]; ii < Lout; ii++)
   {
      truepos = (ii + outSpos) * oneOverRatio - rp->inPos;
      highpos = ceil(truepos);
      x = truepos - (highpos - 1);

      memcpy(buf, &in[highpos - 6], 6 * sizeof * buf);
      out[ii] = lagrange_interp(x, buf);
   }

   /* Copy last 5 samples to overlap .*/
   memcpy(rp->overlap, in + Lin - 5, 5 * sizeof * in);
   free(iiThre);
   free(buf);
}

/* y = [y(-2),y(-1),y(0),y(1),y(2),y(3)] */
/*  */
SAMPLE
lagrange_interp(const double x, const SAMPLE* yin)
{
   SAMPLE ym1py1, twentyfourthym2py2, c0, c1, c2, c3, c4, c5;
   const SAMPLE* y = yin + 2;
   ym1py1 = y[-1] + y[1];

   twentyfourthym2py2 = 1 / 24.0 * (y[-2] + y[2]);
   c0 = y[0];
   c1 = 1 / 20.0 * y[-2] - 1 / 2.0 * y[-1] - 1 / 3.0 * y[0] + y[1] -
        1 / 4.0 * y[2] + 1 / 30.0 * y[3];
   c2 = 2 / 3.0 * ym1py1 - 5 / 4.0 * y[0] - twentyfourthym2py2;
   c3 = 5 / 12.0 * y[0] - 7 / 12.0 * y[1] + 7 / 24.0 * y[2] -
        1 / 24.0 * (y[-2] + y[-1] + y[3]);
   c4 = 1 / 4.0 * y[0] - 1 / 6.0 * ym1py1 + twentyfourthym2py2;
   c5 = 1 / 120.0 * (y[3] - y[-2]) + 1 / 24.0 * (y[-1] - y[2]) +
        1 / 12.0 * (y[1] - y[0]);
   return (SAMPLE) ( ((((c5 * x + c4) * x + c3) * x + c2) * x + c1) * x + c0 );
}

resample_error
resample_execute_hermite(const resample_plan rp,
                         const SAMPLE* in, const size_t Lin,
                         SAMPLE* out)
{


}

SAMPLE
hermite_interp(const double x, const SAMPLE* yin)
{


}



/* This is actual structture used to hold info between consecutive blocks */
struct EMQFfilters_struct
{
   const double fc;
   const SAMPLE* beta0;
   const SAMPLE* gamma0;
   SAMPLE** d0;
   const size_t stages0;
   const SAMPLE* beta1;
   const SAMPLE* gamma1;
   SAMPLE** d1;
   const size_t stages1;
   const SAMPLE alpha1;
};

EMQFfilters
emqffilters_init(const double fc)
{
   double alpha, alpha1;
   size_t ii, stages0, stages1;
   SAMPLE *beta0, *gamma0, *beta1, *gamma1, **d0, **d1;
   /* EMQFcoefs is a global variable, defined in filtcoefs.h
    * generated by a matlab script genfiltcoefs.m */
   const double* beta = EMQFcoefs;

   stages0 = (size_t) ceil(EMQFCOEFLEN / 2.0);
   stages1 = (size_t) floor(EMQFCOEFLEN / 2.0);

   beta0 = malloc(stages0 * sizeof * beta0);
   gamma0 = malloc(stages0 * sizeof * gamma0);
   beta1 = malloc(stages1 * sizeof * beta1);
   gamma1 = malloc(stages1 * sizeof * gamma1);
   d0 = malloc(stages0 * sizeof(SAMPLE*));
   d1 = malloc((stages1 + 1) * sizeof(SAMPLE*));

   alpha = -cos(M_PI * fc);
   alpha1 = (1.0 - sqrt(1.0 - alpha * alpha)) / alpha;

   for (ii = 0; ii < stages0; ii++)
   {
      beta0[ii] = (SAMPLE) ((beta[2 * ii] + alpha1 * alpha1) /
                            (beta[2 * ii] * alpha1 * alpha1 + 1.0) );
      gamma0[ii] = (SAMPLE) (alpha * (1.0 + beta0[ii]));
      d0[ii] = calloc(2, sizeof(SAMPLE));
   }

   for (ii = 0; ii < stages1; ii++)
   {
      beta1[ii] = (SAMPLE) ((beta[2 * ii + 1] + alpha1 * alpha1) /
                            (beta[2 * ii + 1] * alpha1 * alpha1 + 1.0) );
      gamma1[ii] = (SAMPLE) (alpha * (1.0 + beta1[ii]));
      d1[ii] = calloc(2, sizeof(SAMPLE));
   }

   d1[stages1] = calloc(1, sizeof(SAMPLE));

   EMQFfilters rv = malloc(sizeof * rv);
   memcpy((double*)&rv->fc, &fc, sizeof(fc));
   memcpy((double*)&rv->alpha1, &alpha1, sizeof(alpha1));
   rv->beta0 = beta0;
   rv->gamma0 = gamma0;
   rv->beta1 = beta1;
   rv->gamma1 = gamma1;
   rv->d0 = d0;
   rv->d1 = d1;

   memcpy((double*)&rv->stages0, &stages0, sizeof(stages0));
   memcpy((double*)&rv->stages1, &stages1, sizeof(stages1));
   return rv;
}

void
emqffilters_dofilter(EMQFfilters ef, const SAMPLE* in, const size_t Lin,
                     SAMPLE* out)
{
   size_t ii, jj;
   SAMPLE x1, y1;

   for (ii = 0; ii < Lin; ii++)
   {
      /* Branch 0 */
      /* Feedig output of one stage to the input of the next stage */
      x1 = in[ii];
      y1 = x1;
      for (jj = 0; jj < ef->stages0; jj++)
      {
         y1 = x1 * ef->beta0[jj] + ef->d0[jj][0];
         ef->d0[jj][0] = ef->gamma0[jj] * (x1 - y1) +
                         ef->d0[jj][1];
         ef->d0[jj][1] = x1 - y1 * ef->beta0[jj];
         x1 = y1;
      }
      /* Store the partial output */
      out[ii] = y1;
      /* And start over with the second branch */
      x1 = in[ii];

      /* Branch 1 */
      for (jj = 0; jj < ef->stages1; jj++)
      {
         y1 = x1 * ef->beta1[jj] + ef->d1[jj][0];
         ef->d1[jj][0] = ef->gamma1[jj] * (x1 - y1) +
                         ef->d1[jj][1];
         ef->d1[jj][1] = x1 - y1 * ef->beta1[jj];
         x1 = y1;
      }

      /* Final all-pass filter in Branch 1  */
      y1 = x1 * ef->alpha1 + ef->d1[ef->stages1][0];
      ef->d1[ef->stages1][0] = x1 - y1 * ef->alpha1;

      /* Add output of the second branch to output */
      out[ii] += y1;
      out[ii] /= 2.0;
   }
}


void
emqffilters_done(EMQFfilters* ef)
{
   size_t ii;
   free((void*)(*ef)->beta0);
   free((void*)(*ef)->gamma0);
   free((void*)(*ef)->beta1);
   free((void*)(*ef)->gamma1);

   for (ii = 0; ii < (*ef)->stages0; ii++)
   {
      free((*ef)->d0[ii]);
   }
   free((*ef)->d0);

   for (ii = 0; ii < (*ef)->stages1 + 1; ii++)
   {
      free((*ef)->d1[ii]);
   }
   free((*ef)->d1);


   free(*ef);
   ef = NULL;
}
