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
%   ISGRAMREAL(s,g,a,M) attempts to invert a spectrogram computed by
%
%C     s = abs(dgtreal(f,g,a,M)).^2;
%
%   by an iterative method.
%
%   ISGRAMREAL(s,g,a,M,Ls) does as above but cuts or extends f to length Ls.
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
%     'griflim'  - Use the Griffin-Lim iterative method, this is the
%                  default.
%
%-    'bfgs'     - Use the limited-memory Broyden Fletcher Goldfarb
%                  Shanno (BFGS) method.  
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
%   To use the BFGS method, please install the minFunc software from
%   http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
%
%   See also:  dgt, idgt
%
%   Demos: demo_isgram
%
%R  griffin1984sem decorsiere2011 liu1989limited
  
%   AUTHOR : Remi Decorsiere and Peter Soendergaard.
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
  definput.flags.method={'griflim','bfgs'};
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
  
  M2=size(s,1);
  N=size(s,2);
  W=size(s,3);
  
  % use assert_squarelat to check a and the window size.
  assert_squarelat(a,M,1,'ISGRAMREAL');
  
  L=N*a;
  
  sqrt_s=sqrt(s);
  
  if flags.do_zero
    % Start with a phase of zero.
    c=sqrt(s);
  end;
  
  if flags.do_rand
    c=sqrt_s.*exp(2*pi*1i*rand(size(s)));
  end;
  
  if flags.do_int
    s2=zeros(M,N);
    s2(1:M2,:)=s;
    if rem(M,2)==0
      s2(M2+1:M,:)=flipud(s(2:end-1,:));
    else
      s2(M2+1:M,:)=flipud(s(2:end));
    end;
    c=constructphase(s2,g,a);
    c=c(1:M2,:);
  end;
  
  g  = gabwin(g,a,M,L,'ISGRAMREAL');
  
  gd = gabdual(g,a,M);
  
  assert_L(L,size(g,1),L,a,M,'ISGRAMREAL');
  
  % For normalization purposes
  norm_s=norm(s,'fro');
  
  relres=zeros(kv.maxit,1);
  if flags.do_griflim
    for iter=1:kv.maxit
      f=comp_idgtreal(c,gd,a,M,L,0);
      c=comp_dgtreal(f,g,a,M,L,0);
      
      relres(iter)=norm(abs(c).^2-s,'fro')/norm_s;
      
      c=sqrt_s.*exp(1i*angle(c));
      
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
    
  end;
  
  if flags.do_bfgs
    if exist('minFunc')~=2
      error(['To use the BFGS method in ISGRAMREAL, please install the minFunc ' ...
             'software from http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html.']);
    end;
    
    % Setting up the options for minFunc
    opts = struct;
    opts.display = kv.printstep;
    opts.maxiter = kv.maxit;
    opts.usemex = 0;

    % Don't limit the number of function evaluations, just the number of
    % time-steps.
    opts.MaxFunEvals = 1e9;
    
    f0 = comp_idgtreal(c,gd,a,M,L,0);
    [f,fval,exitflag,output]=minFunc(@objfun,f0,opts,g,a,M,s);
    % First entry of output.trace.fval is the objective function
    % evaluated on the initial input. Skip it to be consistent.
    relres = sqrt(output.trace.fval(2:end))/norm_s;
    iter = output.iterations;
  end;
  
  % Cut or extend f to the correct length, if desired.
  if ~isempty(Ls)
    f=postpad(f,Ls);
  else
    Ls=L;
  end;
  
  f=comp_sigreshape_post(f,Ls,wasrow,[0; W]);
  
%  Subfunction to compute the objective function for the BFGS method.
function [f,df]=objfun(x,g,a,M,s);
  L=size(s,2)*a;
  c=comp_dgtreal(x,g,a,M,L,0);
  
  inner=abs(c).^2-s;
  f=norm(inner,'fro')^2;
  
  df=4*real(conj(comp_idgtreal(inner.*c,g,a,M,L,0)));



