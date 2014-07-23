/*
 * ltfatresample.h
 * This header contains routines for blockwise resampling.
 */
#ifndef _LTFATRESAMPLE_H
#define _LTFATRESAMPLE_H

/* The following defines SAMPLE */
#include "config.h"
#include "stdlib.h"

/* Resampling algorithm */
typedef enum {
   LINEAR = 0, /* Mere lin. intp. */
   LAGRANGE,   /* 6point Lagrange interpolator */
   HERMITE     /* 6point Hermite interpolator */
} resample_type;

/* Error code to check */
typedef enum {
  RESAMPLE_OK = 0,
  RESAMPLE_INVALID_MEMORY
} resample_error;

typedef struct resample_plan_struct *resample_plan;

resample_plan 
resample_init(const resample_type restype, 
              const double ratio);

resample_error
resample_execute(const resample_plan rp,const SAMPLE* in,const size_t Lin, SAMPLE* out);

size_t
resample_nextoutlen(const resample_plan rp, size_t Lin);

void
resample_advanceby(const resample_plan rp, size_t Lin);



void
resample_done(resample_plan *rp);

/* Internal routines */

resample_error
resample_execute_linear(const resample_plan rp,
                        const SAMPLE* in,const size_t Lin,
                        SAMPLE* out);

resample_error
resample_execute_lagrange(const resample_plan rp,
                          const SAMPLE* in,const size_t Lin,
                          SAMPLE* out);

SAMPLE
lagrange_interp(const double x, const SAMPLE *yin);


resample_error
resample_execute_hermite(const resample_plan rp,
                          const SAMPLE* in,const size_t Lin,
                          SAMPLE* out);







#endif
