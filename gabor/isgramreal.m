function [f,relres,iter]=isgramreal(s,g,a,M,varargin)
%ISGRAM  Spectrogram inversion
  %   Usage:  f=isgram(s,g,a,M);
%           f=isgram(s,g,a,M,Ls);
%           [f,relres,iter]=isgram(...);
%
%   Input parameters:
%         c       : Array of coefficients.
%         g       : Window function.
%         a       : Length of time shift.
%         M       : Number of channels.
%         Ls      : length of signal.
%   Output parameters:
%         f       : Signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   ISGRAMREAL(s,g,a) attempts to invert a spectrogram computed by
%
%C     s = abs(dgtreal(f,g,a,M));
%
%   by an iterative method.
%
%   ISGRAMREAL(c,g,a,Ls) does as above but cuts or extends f to length Ls.
%
%   If the phase of the spectrogram is known, it is much better to use
%   IDGTREAL.
%
%   [f,relres,iter]=ISGRAMREAL(...) additionally return the residuals in a
%   vector relres and the number of iteration steps done.
%
%   Generally, if the spectrogram has not been modified, the iterative
%   algorithm will converge slowly to the correct result. If the
%   spectrogram has been modified, the algorithm is not guaranteed to
%   converge at all.  
%
%   ISGRAMREAL takes the following parameters at the end of the line of input
%   arguments:
%
%-    'zero'     - Choose a starting phase of zero. This is the default
%
%-    'rand'     - Choose a random starting phase.
%
%-    'int'      - Construct a starting phase by integration. Only works
%                  for Gaussian windows. This is not implemented yet.
%
%-    'tol',t    - Stop if relative residual error is less than the specified tolerance.  
%
%-    'maxit',n  - Do at most n iterations.
%
%-    'print'    - Display the progress.
%
%-    'quiet'    - Don't print anything, this is the default.
%
%-    'printstep',p - If 'print' is specified, then print every p'th
%                  iteration. Default value is p=10;
%
%   See also:  dgt, idgt
%
%   Demos: demo_isgram
%
%R  griffin1984sem
  
%   AUTHOR : Peter Soendergaard.
%   REFERENCE: OK

% Check input paramameters.

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if numel(g)==1
  error('g must be a vector (you probably forgot to supply the window function as input parameter.)');
end;

definput.keyvals.Ls=[];
definput.keyvals.tol=1e-6;
definput.keyvals.maxit=100;
definput.keyvals.printstep=10;
definput.flags.print={'print','quiet'};
definput.flags.startphase={'zero','rand','int'};

[flags,kv,Ls]=ltfatarghelper({'Ls','tol','maxit'},definput,varargin);

wasrow=0;

if isnumeric(g)
  if size(g,2)>1
    if size(g,1)>1
      error('g must be a vector');
    else
      % g was a row vector.
      g=g(:);
      
      % If the input window is a row vector, and the dimension of c is
      % equal to two, the output signal will also
      % be a row vector.
      if ndims(s)==2
        wasrow=1;
      end;
    end;
  end;
end;

N=size(s,2);
W=size(s,3);

% use assert_squarelat to check a and the window size.
assert_squarelat(a,M,1,'ISGRAMREAL');

L=N*a;

if flags.do_zero
  % Start with a phase of zero.
  c=s;
end;

if flags.do_rand
  c=s.*exp(2*pi*i*rand(M,N));
end;

if flags.do_int
  error('This has not been implemented yet. Please use ISGRAM, which supports this option.');
  %c=constructphase(s,g,a);
end;

g  = gabwin(g,a,M,L,'ISGRAMREAL');

gd = gabdual(g,a,M);

assert_L(L,size(g,1),L,a,M,'ISGRAMREAL');

% For normalization purposes
norm_s=norm(s,'fro');

relres=zeros(kv.maxit,1);
for iter=1:kv.maxit
  f=comp_idgtreal(c,gd,a,M,L,0);
  c=comp_dgtreal(f,g,a,M,L,0);
  
  relres(iter)=norm(abs(c)-s,'fro')/norm_s;
  
  c=s.*exp(i*angle(c));
  
  if flags.do_print
    if mod(iter,kv.printstep)==0
      fprintf('ISGRAMREAL: Iteration %i, residual = %f.\n',iter,relres(iter));
    end;    
  end;
  
  if relres(iter)<kv.tol
    relres=relres(1:iter);
    break;
  end;

end;

% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,wasrow,[0; W]);




