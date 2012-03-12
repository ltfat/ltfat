function [f,relres,iter]=iframeabs(s,F,varargin)
%IFRAMEABS  Reconstruction from magnitude of coefficients
%   Usage:  f=iframeabs(s,F);
%           f=iframeabs(s,F,Ls);
%           [f,relres,iter]=iframeabs(...);
%
%   Input parameters:
%         c       : Array of coefficients.
%         F       : Frame 
%         Ls      : length of signal.
%   Output parameters:
%         f       : Signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `iframeabs(s,F)` attempts to find a signal which has `s` as the absolute
%   value of its frame coefficients ::
%
%     s = abs(framet(f,F));
%
%   using an iterative method.
%
%   `iframeabs(s,F,Ls)` does as above but cuts or extends *f* to length *Ls*.
%
%   If the phase of the coefficients *s* is known, it is much better to use
%   `iframet`.
%
%   `[f,relres,iter]=iframeabs(...)` additionally return the residuals in a
%   vector *relres* and the number of iteration steps *iter*.
%
%   Generally, if the absolute value of the frame coefficients has not been
%   modified, the iterative algorithm will converge slowly to the correct
%   result. If the coeffficients have been modified, the algorithm is not
%   guaranteed to converge at all.
%
%   `iframeabs` takes the following parameters at the end of the line of input
%   arguments:
%
%     'zero'       Choose a starting phase of zero. This is the default
%
%     'rand'       Choose a random starting phase.
%
%     'tol',t      Stop if relative residual error is less than the
%                  specified tolerance.  
%
%     'maxit',n    Do at most n iterations.
%
%     'print'      Display the progress.
%
%     'quiet'      Don't print anything, this is the default.
%
%     'printstep',p  If 'print' is specified, then print every p'th
%                    iteration. Default value is p=10;
%
%   See also:  dgt, idgt
%
%   References: griffin1984sem
  
%   AUTHOR : Remi Decorsiere and Peter Soendergaard.
%   REFERENCE: OK

% Check input paramameters.

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;
  
definput.keyvals.Ls=[];
definput.keyvals.tol=1e-6;
definput.keyvals.maxit=100;
definput.keyvals.printstep=10;
definput.flags.print={'quiet','print'};
definput.flags.startphase={'zero','rand'};
definput.flags.method={'griflim','bfgs'};

[flags,kv,Ls]=ltfatarghelper({'Ls','tol','maxit'},definput,varargin);

% Determine the proper length of the frame
L=framelengthcoef(s,F);   
W=size(s,2);

% Initialize windows to speed up computation
F=frameaccel(F,L);

if flags.do_zero
  % Start with a phase of zero.
  c=s;
end;

if flags.do_rand
  c=s.*exp(2*pi*i*rand(size(s)));
end;

% For normalization purposes
norm_s=norm(s,'fro');

relres=zeros(kv.maxit,1);
if flags.do_griflim
  
  for iter=1:kv.maxit
    f=iframet(c,F);
    c=framet(f,F);
    
    relres(iter)=norm(abs(c)-s,'fro')/norm_s;
    
    c=s.*exp(i*angle(c));
    
    if flags.do_print
      if mod(iter,kv.printstep)==0
        fprintf('IFRAMEABS: Iteration %i, residual = %f.\n',iter,relres(iter));
      end;    
    end;
    
    if relres(iter)<kv.tol
      relres=relres(1:iter);
      break;
    end;
    
  end;
end;
 
    
% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,0,[0; W]);

