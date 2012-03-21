function [tc,relres,iter,xrec] = gablasso(x,g,a,M,lambda,varargin)
%GABLASSO  LASSO regression in Gabor domain
%   Usage: [tc,xrec] = gablasso(x,a,M,lambda,C,tol,maxit)
%
%   `gablasso` has been deprecated. Please use |framelasso|_ instead.
%
%   A call to `gablasso(x,g,a,M,lambda)` can be replaced by ::
%
%     F=newframe('dgt',[],g,a,M);
%     tc=framelasso(F,lambda);
%
%   Any additional parameters passed to `gablasso` can be passed to
%   |framelasso|_ in the same manner.
%
%   See also: newframe, framelasso

if nargin<5
  error('%s: Too few input parameters.',upper(mfilename));
end;

% Define initial value for flags and key/value pairs.
definput.keyvals.C=[];
definput.keyvals.tol=1e-2;
definput.keyvals.maxit=100;
definput.keyvals.printstep=10;
definput.flags.print={'print','quiet'};
definput.flags.startphase={'zero','rand','int'};

[flags,kv]=ltfatarghelper({'C','tol','maxit'},definput,varargin);

  
%   AUTHOR : Bruno Torresani.  
%   TESTING: OK

%   XXX Removed Remark: When the frame is an orthonormal basis, the solution
%   is obtained by soft thresholding of the basis coefficients, with
%   threshold lambda.  When the frame is a union of orthonormal bases, the
%   solution is obtained by applying soft thresholding cyclically on the
%   basis coefficients (BCR algorithm)

% Determine transform length, and calculate the window.
[x,g,L] = gabpars_from_windowsignal(x,g,a,M,[],'GABLASSO');

% Initialization of thresholded coefficients
c0 = dgt(x,g,a,M);

L=size(c0,2)*a;
if isempty(kv.C)
  [A_dummy,kv.C] = gabframebounds(g,a,M,L);
end;

% Various parameter initializations
threshold = lambda/kv.C;

tc0 = c0;
relres = 1e16;
iter = 0;

% Main loop
while ((iter < kv.maxit)&&(relres >= kv.tol))
    tc = c0 - dgt(idgt(tc0,g,a),g,a,M);
    tc = tc0 + tc/kv.C;
    tc = thresh(tc,threshold,'soft');
    relres = norm(tc(:)-tc0(:))/norm(tc0(:));
    tc0 = tc;
    iter = iter + 1;
    if flags.do_print
      if mod(iter,kv.printstep)==0        
        fprintf('Iteration %d: relative error = %f\n',iter,relres);
      end;
    end;
end

% Reconstruction
if nargout>3
  xrec = idgt(tc,g,a);
end;
