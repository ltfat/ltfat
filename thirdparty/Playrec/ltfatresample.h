/*
 * ltfatresample.h
 *
 * This is an attempt for a self-containded arbitrary factor resampling API
 * working with streams of blocks of arbitrary lengths.
 *
 * Copyright (C) 2014 Zdenek Prusa <prusa@users.sourceforge.net>.
 * This file is part of LTFAT http://ltfat.sourceforge.net
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _LTFATRESAMPLE_H
#define _LTFATRESAMPLE_H

/* The following defines SAMPLE */
#include "config.h"
#include "stdlib.h"

/* Here we include a generated file containing prototype filters */
/* We use elliptic minimal Q-factors IIR filters (EMQF) from
 * Multirate Filtering for DSP: MATLAB Applications by Ljiljana Milic, chapter V
 * */
#include "filtcoefs.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/*
 *  POLYNOMIAL INTERPOLATION
 *
 *
 * */

/* Resampling algorithm */
typedef enum
{
   LINEAR = 0, /* Mere lin. intp. */
   LAGRANGE,   /* 6point Lagrange interpolator */
   HERMITE     /* 6point Hermite interpolator */
} resample_type;

/* Error code to check */
typedef enum
{
   RESAMPLE_OK = 0,
   RESAMPLE_INVALID_MEMORY
} resample_error;

typedef struct resample_plan_struct *resample_plan;

resample_plan
resample_init(const resample_type restype,
              const double ratio);

resample_error
resample_execute(const resample_plan rp, const SAMPLE* in, const size_t Lin, SAMPLE* out);

size_t
resample_nextoutlen(const resample_plan rp, size_t Lin);

void
resample_advanceby(const resample_plan rp, size_t Lin);

void
resample_done(resample_plan *rp);

/* Internal routines */

resample_error
resample_execute_linear(const resample_plan rp,
                        const SAMPLE* in, const size_t Lin,
                        SAMPLE* out);

resample_error
resample_execute_lagrange(const resample_plan rp,
                          const SAMPLE* in, const size_t Lin,
                          SAMPLE* out);

SAMPLE
lagrange_interp(const double x, const SAMPLE *yin);


resample_error
resample_execute_hermite(const resample_plan rp,
                         const SAMPLE* in, const size_t Lin,
                         SAMPLE* out);



/* Struct for holding EMQF filters */
typedef struct EMQFfilters_struct *EMQFfilters;

EMQFfilters
emqffilters_init(const double fc);

void
emqffilters_dofilter(EMQFfilters ef, const SAMPLE* in, const size_t Lin,
                     SAMPLE* out);

void
emqffilters_done(EMQFfilters* ef);

#endif
