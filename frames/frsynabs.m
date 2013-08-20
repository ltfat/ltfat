function [f,relres,iter]=frsynabs(F,s,varargin)
%FRSYNABS  Reconstruction from magnitude of coefficients
%   Usage:  f=frsynabs(F,s);
%           f=frsynabs(F,s,Ls);
%           [f,relres,iter]=frsynabs(...);
%
%   Input parameters:
%         F       : Frame   
%         s       : Array of coefficients.
%         Ls      : length of signal.
%   Output parameters:
%         f       : Signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `frsynabs(F,s)` attempts to find a signal which has `s` as the absolute
%   value of its frame coefficients ::
%
%     s = abs(frana(F,f));
%
%   using an iterative method.
%
%   `frsynabs(F,s,Ls)` does as above but cuts or extends *f* to length *Ls*.
%
%   If the phase of the coefficients *s* is known, it is much better to use
%   `frsyn`.
%
%   `[f,relres,iter]=frsynabs(...)` additionally returns the residuals in a
%   vector *relres* and the number of iteration steps *iter*. The
%   residuals are computed as:
%
%   .. relres = norm(abs(cn)-s,'fro')/norm(s,'fro') 
%
%   .. math:: relres=\frac{\left\Vert \left|c_{n}\right|-s\right\Vert _{fro}}{\left\Vert s\right\Vert _{fro}},
%
%   where $c_n$ is the Gabor coefficients of the signal in iteration *n*.
%
%   Generally, if the absolute value of the frame coefficients has not been
%   modified, the iterative algorithm will converge slowly to the correct
%   result. If the coeffficients have been modified, the algorithm is not
%   guaranteed to converge at all.
%
%   `frsynabs` takes the following parameters at the end of the line of input
%   arguments:
%
%     'input'      Choose the starting phase as the phase of the input
%                  *s*. This is the default
%
%     'zero'       Choose a starting phase of zero.
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
  
%   AUTHOR : Remi Decorsiere and Peter L. Søndergaard.
%   REFERENCE: OK

% Check input paramameters.

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;
  
definput.keyvals.Ls=[];
definput.keyvals.tol=1e-6;
definput.keyvals.maxit=100;
definput.keyvals.printstep=10;
definput.keyvals.alpha=0.99;
definput.flags.print={'quiet','print'};
definput.flags.startphase={'input','zero','rand'};
definput.flags.method={'griflim','bfgs','fgriflim'};

[flags,kv,Ls]=ltfatarghelper({'Ls','tol','maxit'},definput,varargin);

% Determine the proper length of the frame
L=framelengthcoef(F,size(s,1));   
W=size(s,2);

% Initialize windows to speed up computation
F=frameaccel(F,L);


if flags.do_input
  % Start with the phase given by the input.
  c=s;
end;

if flags.do_zero
  % Start with a phase of zero.
  c=abs(s);
end;

if flags.do_rand
  c=abs(s).*exp(2*pi*i*rand(size(s)));
end;

% For normalization purposes
norm_s=norm(s,'fro');

relres=zeros(kv.maxit,1);
if flags.do_griflim
  
  Fs=frameaccel(framedual(F),L);
  
  for iter=1:kv.maxit
    f=Fs.frsyn(c);
    c=F.frana(f);
    
    relres(iter)=norm(abs(c)-s,'fro')/norm_s;
    
    c=s.*exp(i*angle(c));        


    if flags.do_print
      if mod(iter,kv.printstep)==0
        fprintf('FRSYNABS: Iteration %i, residual = %f.\n',iter,relres(iter));
      end;    
    end;
    
    if relres(iter)<kv.tol
      relres=relres(1:iter);
      break;
    end;
    
  end;
end;

if flags.do_fgriflim
  
  Fs=frameaccel(framedual(F),L);

  told=s;

  for iter=1:kv.maxit
    f=Fs.frsyn(c);
    tnew=F.frana(f);

    relres(iter)=norm(abs(tnew)-s,'fro')/norm_s;

    tnew=s.*exp(i*angle(tnew));
    c=tnew+kv.alpha*(tnew-told);
    

    if flags.do_print
      if mod(iter,kv.printstep)==0
        fprintf('FRSYNABS: Iteration %i, residual = %f.\n',iter,relres(iter));
      end;    
    end;
    
    if relres(iter)<kv.tol
      relres=relres(1:iter);
      break;
    end;
    
    told=tnew;
    
  end;
end;

    
% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,0,[0; W]);

