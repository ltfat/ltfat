#include "config.h"
#include "ltfat.h"
#include "winmanip.h"

/* Compute canonical dual/tight window. This last parameter
 * indicates the type: 0 = dual, 1 = tight.
 */
void fircanon_r(const double *g, const ltfatInt Lg, const ltfatInt L, const ltfatInt a,
		const ltfatInt M, double *gdual, const ltfatInt Ldual, const ltfatInt symm,
		const ltfatInt wintype)
{
   
   double *tmp_fir, *tmp_iir;
   
   tmp_fir = (double*)ltfat_malloc(Lg*sizeof(double));
   tmp_iir = (double*)ltfat_malloc(L*sizeof(double));

   /* Move center of window from the middle of the vector to the beginning. */
   ifftshift_r(g, Lg, tmp_fir);
   
   /* Extend the FIR window to an IIR window. */
   fir2iir_r(tmp_fir, Lg, L, tmp_iir);
         
   if (wintype==0)
   {
     gabdualreal_long(g, L, 1, a, M, tmp_iir);
   }
   else
   {
     gabtightreal_long(g, L, 1, a, M, tmp_iir);
   }
      
   /* Cut dual IIR window to a FIR window. */
   iir2fir_r(tmp_iir, L, Ldual, symm, tmp_fir);
   
   /* Move center of window to the middle of the vector. */
   fftshift_r(tmp_fir, Ldual, gdual);
   
   ltfat_free(gdualf);
   ltfat_free(gf);
   ltfat_free(tmp_iir);
   ltfat_free(tmp_fir);

}


/* Driver routine to calculate dual of FIR window. This routine
 * Input:
 *          g     : pointer to FIR window
 *          Lg    : Length of g
 *          L     : Length of system for which g and gdual should be
 *                  dual windows.
 *          a     : Length of time step (hop size)
 *          M     : Number of channels.
 *          gdual : pointer to dual window
 *          Ldual : Length of dual window
 *          symm  : Symmetry of input window, see the help for iir2fir
 *
 */ 
void firdual_r(const double *g, const ltfatInt Lg, const ltfatInt L, const ltfatInt a,
		const ltfatInt M, double *gdual, const ltfatInt Ldual, const ltfatInt symm)
{
   /* The final 0 indicates that we want the dual window.*/
   fircanon_r(g, Lg, L, a, M, gdual, Ldual,  symm, 0);
}


/* Driver routine to calculate tight window of FIR window. Same input/output
 * parameters as firdual_r
 */
void firtight_r(const double *g, const ltfatInt Lg, const ltfatInt L, const ltfatInt a,
		const ltfatInt M, double *gdual, const ltfatInt Ldual, const ltfatInt symm)
{
   /* The final 1 indicates that we want the tight window.*/
   fircanon_r(g, Lg, L, a, M, gdual, Ldual,  symm, 1);
}


