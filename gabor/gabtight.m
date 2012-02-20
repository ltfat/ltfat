function gt=gabtight(p1,p2,p3,p4)
%GABTIGHT  Canonical tight window of Gabor frame
%   Usage:  gt=gabtight(a,M,L);
%           gt=gabtight(g,a,M);
%           gt=gabtight(g,a,M,L);
%
%   Input parameters:
%         g     : Gabor window.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         L     : Length of window. (optional)
%   Output parameters:
%         gt    : Canonical tight window, column vector.
%
%   `gabtight(a,M,L)` computes a nice tight window of length *L* for a
%   lattice with parameters *a*, *M*. The window is not an FIR window,
%   meaning that it will only generate a tight system if the system
%   length is equal to *L*.
%
%   `gabtight(g,a,M)` computes the canonical tight window of the Gabor frame
%   with window *g* and parameters *a*, *M*.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%  
%   If the length of *g* is equal to *M*, then the input window is assumed to
%   be a FIR window. In this case, the canonical dual window also has
%   length of *M*. Otherwise the smallest possible transform length is
%   chosen as the window length.
%
%   `gabtight(g,a,M,L)` returns a window that is tight for a system of
%   length *L*. Unless the input window *g* is a FIR window, the returned
%   tight window will have length *L*.
%
%   If $a>M$ then an orthonormal window of the Gabor Riesz sequence with
%   window *g* and parameters *a* and *M* will be calculated.
%
%   See also:  gabdual, gabwin, fir2long, dgt

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK

%% ------------ decode input parameters ------------
  
if numel(p1)==1
  % First argument is a scalar.

  error(nargchk(3,3,nargin));

  a=p1;
  M=p2;
  L=p3;

  g='gauss';
  
else    
  % First argument assumed to be a vector.
    
  error(nargchk(3,4,nargin));
  
  g=p1;
  a=p2;
  M=p3;
    
  if nargin==3
    L=[];
  else
    L=p4;
  end;

end;

[g,L,info] = gabpars_from_window(g,a,M,L);
  
% -------- Are we in the Riesz sequence of in the frame case

scale=1;
if a>M
  % Handle the Riesz basis (dual lattice) case.
  % Swap a and M, and scale differently.
  scale=sqrt(a/M);
  tmp=a;
  a=M;
  M=tmp;
end;

% -------- Compute ------------- 

if (info.gl<=M)
     
  % FIR case
  N_win = ceil(info.gl/a);
  Lwin_new = N_win*a;
  if Lwin_new ~= info.gl
    g_new = fir2long(g,Lwin_new);
  else
    g_new = g;
  end
  weight = sqrt(sum(reshape(abs(g_new).^2,a,N_win),2));
  
  gt = g_new./repmat(weight,N_win,1);
  gt = gt/sqrt(M);
  if Lwin_new ~= info.gl
    gt = long2fir(gt,info.gl);
  end
  
else
  
  % Long window case

  % Just in case, otherwise the call is harmless. 
  g=fir2long(g,L);
  
  gt=comp_gabtight_long(g,a,M)*scale;
  
end;

% --------- post process result -------

%if info.wasreal
%  gt=real(gt);
%end;
      
if info.wasrow
  gt=gt.';
end;
