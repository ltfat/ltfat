#include "ltfatresample.h"
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <mex.h>

struct resample_plan_struct
{
   const double ratio;
   const resample_type restype;
   size_t inPos;
   SAMPLE* overlap;
   resample_error (*executer)(const resample_plan,
                              const SAMPLE*, const size_t,
                              SAMPLE*);
};

resample_plan
resample_init(const resample_type restype, const double ratio)
{
   SAMPLE* overlap;
   struct resample_plan_struct rp = {ratio, restype, 0, NULL, NULL };
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
   /* Execute the computation */
   resample_error status = rp->executer(rp, in, Lin, out);
   /* Advance the "pointer" in data stream */
   resample_advanceby(rp, Lin);
   /* rp->inPos += Lin;*/

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
   SAMPLE prd;
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
   c1 = 1 / 20.0 * y[-2] - 1 / 2.0 * y[-1] - 1 / 3.0 * y[0] + y[1] - 1 / 4.0 * y[2] + 1 / 30.0 * y[3];
   c2 = 2 / 3.0 * ym1py1 - 5 / 4.0 * y[0] - twentyfourthym2py2;
   c3 = 5 / 12.0 * y[0] - 7 / 12.0 * y[1] + 7 / 24.0 * y[2] - 1 / 24.0 * (y[-2] + y[-1] + y[3]);
   c4 = 1 / 4.0 * y[0] - 1 / 6.0 * ym1py1 + twentyfourthym2py2;
   c5 = 1 / 120.0 * (y[3] - y[-2]) + 1 / 24.0 * (y[-1] - y[2]) + 1 / 12.0 * (y[1] - y[0]);
   return (SAMPLE) ( ((((c5 * x + c4) * x + c3) * x + c2) * x + c1) * x + c0 );
}

resample_error
resample_execute_hermite(const resample_plan rp,
                         const SAMPLE* in, const size_t Lin,
                         SAMPLE* out)
{


}
