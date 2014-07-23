#ifndef _CONFIG_H
#define _CONFIG_H

/* The mx class equivalent to SAMPLE */
#define mxSAMPLE    mxSINGLE_CLASS
/* Format to be used for samples with PortAudio = 32bit */
typedef float SAMPLE;

/* See enum resample_type in ltfatresample.h for 
 * possible values. */
#define RESAMPLING_TYPE LAGRANGE

#endif
