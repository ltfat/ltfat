function [V,D]=gabmuleigs(K,c,p3,varargin)
%GABMULEIGS  Eigenpairs of Gabor multiplier
%   Usage:  h=gabmuleigs(K,c,g,a);
%           h=gabmuleigs(K,c,a);
%           h=gabmuleigs(K,c,ga,gs,a);
%
%   Input parameters:
%         K     : Number of eigenvector to compute.
%         c     : symbol of Gabor multiplier
%         g     : analysis/synthesis window
%         ga    : analysis window
%         gs    : synthesis window
%         a     : Length of time shift.
%   Output parameters:
%         V     : Matrix containing eigenvectors.
%         D     : Eigenvalues.
%
%   GABMULEIGS(K,c,g,a) computes the K largest eigenvalues and eigen-
%   vectors of the Gabor multiplier with symbol c and time shift _a.
%   The number of channels is deduced from the size of the symbol c.
%   The window g will be used for both analysis and synthesis.
%
%   GABMULEIGS(K,c,ga,gs,a) will do the same using the window the window ga
%   for analysis and gs for synthesis.
%
%   GABMULEIGS(K,c,a) will do the same using the a tight Gaussian window of
%   for analysis and synthesis.
%
%   If K is empty, then all eigenvalues/pairs will be returned.
%
%   See also: gabmul, dgt, idgt, gabdual, gabtight
%

% Change this to 1 or 2 to see the iterative method in action.
printopts=0;

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if nargout==2
  doV=1;
else
  doV=0;
end;

M=size(c,1);
N=size(c,2);

% Tolerance of iterative method
tol=1e-10;

% Maximum number of iterations.
maxit=200;

% When to use an iterative method.
crossover=200;

istight=1;
if numel(p3)==1
  % Usage: h=gabmuleigs(c,K,a);  
  a=p3;
  L=N*a;
  ga=gabtight(a,M,L);
  gs=ga;
  arglist=varargin;
else 
  if numel(varargin{1})==1
    % Usage: h=gabmuleigs(c,K,g,a);  
    ga=p3;
    gs=p3;
    a=varargin{1};
    L=N*a;
    arglist=varargin(2:end);
  else 
    if numel(varargin{2})==1
      % Usage: h=gabmuleigs(c,K,ga,gs,a);  
      ga=p3;
      gs=varargin{1};
      a =varargin{2};
      L=N*a;
      istight=0;
      arglist=varargin(3:end);
    end;    
  end;
end;

crossover=200;

% Do the computation.
% eigs does not work in Octave, and for small problems a direct calculation is just as fast.
if isoctave || L<crossover
  if istight
    gabbas=tfmat('dgt',ga,a,M);
    bigM=gabbas*diag(c(:))*gabbas';
  else
    gabbasa=tfmat('dgt',ga,a,M);
    gabbass=tfmat('dgt',gs,a,M);
    bigM=gabbass*diag(c(:))*gabbasa';      
  end;
  
  if doV
    [V,D]=eig(bigM);
  else
    D=eig(bigM);
  end;

else
  opts.disp=printopts;
  gfa=comp_wfac(ga,a,M);
  if istight
    gfs=gfa;
  else
    gfs=comp_wfac(gs,a,M);
  end;

  if doV
    [V,D] = eigs(@afun_sep,L,K,'LM',opts,c,gfa,gfs,L,a,M);
  else
    D = eigs(@afun_sep,L,K,'LM',opts,c,gfa,gfs,a,M);
  end;

end;

D=diag(D);


function y=afun_sep(x,c,gfa,gfs,L,a,M)
 % Apply ifft to the coefficients.
  c=ifft(c);

  y=comp_idgt_fac(ifft(c.*fft(comp_dgt_fac(x,gfa,a,M))),gfs,L,a,M);
