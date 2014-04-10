function c=comp_dgt(f,g,a,M,lt,phasetype,algfir,algns)
%COMP_DGT  Compute a DGT
%   Usage:  c=comp_dgt(f,g,a,M,L,phasetype);
%
%   Input parameters:
%         f     : Input data
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : Length of transform to do.
%         lt    : Lattice type
%         phasetype : Type of phase
%         algtype : Select algorithm
%   Output parameters:
%         c     : M*N*W array of coefficients.
%
%   If phasetype is zero, a freq-invariant transform is computed. If
%   phase-type is one, a time-invariant transform is computed.
%
%   The algorithm chooser do the following:
%
%      * $algfir=0$ : Default value, automatically choose the fastest
%        algorithm.
%       
%      * $algfir=1$ : Choose the algorithm depending on the input.
%
%      * $algns=0$  : Default value, automatically choose the fastest
%        algorithm.
%
%      * $algns=1$  : Always choose multiwindow.
%
%      * $algns=2$  : Always choose shear
      

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: OK
%   REFERENCE: OK

L=size(f,1);

if lt(2)==1
        c=comp_sepdgt(f,g,a,M,phasetype);
else
        
    g=fir2long(g,L);
    
    if  (algns==0 && lt(2)<=2) || (algns==1)
        
        c=comp_nonsepdgt_multi(f,g,a,M,lt);
        
    else
        
        [s0,s1,br] = shearfind(L,a,M,lt);
        
        c=comp_nonsepdgt_shear(f,g,a,M,s0,s1,br);
        
    end;
    
    % FIXME : Calls non-comp function 
    if phasetype==1
        c=phaselock(c,a,'lt',lt);
    end;

end;




