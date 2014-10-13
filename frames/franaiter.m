function [c,relres,iter]=franaiter(F,f,varargin)
%FRANAITER  Iterative analysis
%   Usage:  c=franaiter(F,f);
%           [c,relres,iter]=franaiter(F,f,...);
%
%   Input parameters:
%         F       : Frame.
%         f       : Signal.
%         Ls      : Length of signal.
%   Output parameters:
%         c       : Array of coefficients.    
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `c=franaiter(F,f)` computes the frame coefficients *c* of the signal *f*
%   using an iterative method such that perfect reconstruction can be
%   obtained using |frsyn|. `franaiter` always works, even when |frana|
%   cannot generate perfect reconstruction coefficients.
%
%   `[c,relres,iter]=franaiter(...)` additionally returns the relative
%   residuals in a vector *relres* and the number of iteration steps *iter*.
%  
%   **Note:** If it is possible to explicitly calculate the canonical dual
%   frame then this is usually a much faster method than invoking
%   `franaiter`.
%
%   `franaiter` takes the following parameters at the end of the line of
%   input arguments:
%
%     'tol',t      Stop if relative residual error is less than the
%                  specified tolerance. Default is 1e-9 (1e-5 for single precision)
%
%     'maxit',n    Do at most n iterations.
%
%     'pg'        Solve the problem using the Conjugate Gradient
%                  algorithm. This is the default.
%
%     'pcg'        Solve the problem using the Preconditioned Conjugate Gradient
%                  algorithm.
%
%     'print'      Display the progress.
%
%     'quiet'      Don't print anything, this is the default.
%
%   Examples
%   --------
%
%   The following example shows how to rectruct a signal without ever
%   using the dual frame:::
%
%      f=greasy;
%      F=frame('dgtreal','gauss',40,60);
%      [c,relres,iter]=franaiter(F,f,'tol',1e-14);
%      r=frsyn(F,c);
%      norm(f-r)/norm(f)
%      semilogy(relres);
%      title('Conversion rate of the CG algorithm');
%      xlabel('No. of iterations');
%      ylabel('Relative residual');
%
%   See also: frame, frana, frsyn, frsyniter
  
% AUTHORS: Peter L. Søndergaard
    
complainif_notenoughargs(nargin,2,'FRANAITER');
complainif_notvalidframeobj(F,'FRANAITER');

tolchooser.double=1e-9;
tolchooser.single=1e-5;

definput.keyvals.Ls=[];
definput.keyvals.tol=tolchooser.(class(f));
definput.keyvals.maxit=100;
definput.flags.alg={'cg','pcg'};
definput.keyvals.printstep=10;
definput.flags.print={'quiet','print'};

[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,~,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));

F=frameaccel(F,Ls);
L=F.L;

%% -- run the iteration 

A=@(x) F.frsyn(F.frana(x));

% An explicit postpad is needed for the pcg algorithm to not fail
f=postpad(f,L);

if flags.do_pcg
    d=framediag(F,L);
    M=spdiags(d,0,L,L);
    
    [fout,flag,~,iter,relres]=pcg(A,f,kv.tol,kv.maxit,M);
else
    
    [fout,flag,~,iter,relres]=pcg(A,f,kv.tol,kv.maxit);          
end;

c=F.frana(fout);

if nargout>1
    relres=relres/norm(fout(:));
end;


%% --- cleanup -----

permutedsize=[size(c,1),permutedsize(2:end)];

c=assert_sigreshape_post(c,dim,permutedsize,order);

