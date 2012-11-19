function f=comp_igdgt(c,g,a,M,L,c_t,c_f,c_w,timeinv)
%COMP_IGDGT  Compute IGDGT
%   Usage:  f=comp_igdgt(c,g,a,M,L,c_t,c_f,c_w,timeinv);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : length of transform.
%         c_t     : Centering in time of modulation.
%         c_f     : Centering in frequency of modulation.
%         c_w     : Centering in time of window.
%         timeinv : Should we compute a time invariant Gabor system.
%
%   Output parameters:
%         f     : Signal.
%

% AUTHOR : Peter L. Søndergaard.

b=L/M;
N=L/a;

Lwindow=size(g,1);
W=size(c,3);

% Pre-process if c_t is different from 0.
if (c_t~=0)
  halfmod=repmat(exp(2*pi*i*c_t*((0:M-1)+c_f).'/M),1,N*W);

  % The following is necessary because REPMAT does not work for
  % 3D arrays.
  halfmod=reshape(halfmod,M,N,W);
    
  c=c.*halfmod;
end;

% Eventual phaselocking
if timeinv
  c=phaseunlock(c,a);
end;

f=comp_idgt(c,g,a,[0 1],0,0);

% Postprocess to handle c_f different from 0.
if (c_f~=0)
  halfmod=exp(2*pi*i*c_f*(0:L-1).'/M);
  f=f.*repmat(halfmod,1,W);
end;


