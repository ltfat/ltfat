function c=comp_nonsepdgt(f,g,a,M,lt,do_timeinv,alg)
%NONSEPDGT  Compute Non-separable Discrete Gabor transform
%   Usage:  c=nonsepdgt(f,g,a,M,lt);
%
%   Input parameters:
%         f     : Input data.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%         lt    : Lattice type
%         do_timeinv : Do a time invariant phase ?
%         alg   : Choose algorithm
%   Output parameters:
%         c     : $M \times N$ array of coefficients.
%
%   `nonsepdgt(f,g,a,M,lt)` computes the non-separable discrete Gabor
%   transform of the input signal *f* with respect to the window *g*,
%   time-shift *a*, number of channels *M* and lattice type *lt*.
%
%     * $alg=0$ : Choose the fastest algorithm
%
%     * $alg=0$ : Always choose multi-win
%
%     * $alg=1$ : Always choose shear
%
%   This is a computational subroutine, do not call it directly.

%   AUTHOR : Nicki Holighaus and Peter L. Soendergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

% Assert correct input.

L=size(f,1);
W=size(f,2);
b=L/M;
N=L/a;

% ----- algorithm starts here, split into sub-lattices ---------------


if (alg==1) || (alg==0 && lt(2)<=2) 

    c=comp_nonsepdgt_multi(f,g,a,M,lt);
    
else

    s=b*lt(1)/lt(2);
    
    [s0,s1,br] = shearfind(a,b,s,L);

    c=comp_nonsepdgt_shear(f,g,a,M,s,s0,s1,br);
    
end;
