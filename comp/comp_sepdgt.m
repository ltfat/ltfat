function [coef]=comp_sepdgt(f,g,a,M,phasetype)
%COMP_SEPDGT  Separable DGT
%   Usage:  c=comp_sepdgt(f,g,a,M);
%  
%   This is a computational routine. Do not call it directly.
%
%   See help on DGT.

%   AUTHOR : Peter L. Søndergaard.

L=size(f,1);
Lwindow=size(g,1);

if Lwindow<L
   % Do the filter bank algorithm
   coef=comp_dgt_fb(f,g,a,M);
else
   % Do the factorization algorithm
   coef=comp_dgt_long(f,g,a,M);
end;


% FIXME : Calls non-comp function 
if phasetype==1
    coef=phaselock(coef,a);
end;
