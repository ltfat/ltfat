function [c,relres,iter,frec,cd] = franabp(F,f,lambda,varargin)
%FRANABP Frame Analysis Basis Pursuit
%   Usage: c = franabp(F,f,lambda)
%          c = franabp(F,f,lambda,C)
%          c = franabp(F,f,lambda,C,tol)
%          c = franabp(F,f,lambda,C,tol,maxit)
%
%   Input parameters:
%       F        : Frame definition
%       f        : Input signal
%       lambda   : Regularization parameter, controls sparsity of the solution
%       C        : Step size of the algorithm.
%       tol      : Reative error tolerance.
%       maxit    : Maximum number of iterations.
%   Output parameters:
%       c        : Coefficients.
%       relres   : Last relative error.
%       iter     : Number of iterations done.  
%       frec     : Reconstructed signal such that frec = frsyn(F,c)
%       cd       : The min ||c||_2 solution using the canonical dual frame.
%  
%   `c = franabp(F,f,lambda)` solves the Basis Pursuit problem
% 
%   .. min ||c||_1 subject to Fc = f
%
%   .. math:: ||c||_1 \\ \text{subject to } Fc = f
%
%   for a general frame *F* using SALSA (Split Augmented Lagrangian 
%   Srinkage algorithm) which is an appication of ADMM (Alternating 
%   Direction Method of Multipliers) to a basis pursuit problem.
%
%   The algorithm acts as follows:
%
%   Initialize d,C>0
%   repeat
%      v <- soft(c+d,lambda/C) - d
%      d <- A*(AA*)^(-1)(f - Av)
%      c <- d + v
%   end
%
%   For a quick execution, the function requires compution of the canonical
%   dual frame. If it is not available, the conjugate cradient algorithm is
%   used for inverting the frame operator in each iteration of the algorithm.
%
%   REMARK: `tol` defines tolerance of `relres` which is norm or a relative
%   difference of coefficients obtained in two consecutive iterations of the
%   alforithm.
%
%   Examples:
%   ---------
%
%   The following example shows how franabp produces a sparse
%   representation of a test signal *greasy* and still maintain perfect
%   reconstruction:::
%
%   f = greasy;
%   % Weight enforcing sparsity 
%   lambda = 0.3; 
%   % Gabor frame with redundancy 8
%   F = frame('dgtreal','gauss',64,512);
%   % Solve the basis pursuit problem
%   [c,~,~,frec,cd] = franabp(F,greasy,lambda);
%   % Plot sparse coefficients
%   figure(1);
%   plotframe(F,c,'dynrange',50);
%
%   % Plot coefficients obtained by applying an analysis operator of a
%   % dual Gabor system to *f*
%   figure(2);
%   plotframe(F,cd,'dynrange',50);
%
%   % Check the reconstruction error (should be close do zero).
%   % frec is obtained by applying the synthesis operator of frame *F*
%   % to sparse coefficients *c*.
%   norm(f-frec)
%
%   % Compare decay of coefficients sorted by absolute values 
%   % (compressibility of coefficients)
%   figure(3);
%   semilogx([sort(abs(c),'descend')/max(abs(c)),...
%   sort(abs(cd),'descend')/max(abs(cd))]);
%   legend({'sparsified coefficients','dual system coefficients'});
%
%   References: se14 bopachupeec11

%   AUTHOR: Zdenek Prusa

% Define initial value for flags and key/value pairs.
definput.keyvals.C=1;
definput.keyvals.tol=1e-4;
definput.keyvals.maxit=100;
definput.keyvals.printstep=10;
definput.flags.print={'print','quiet'};
definput.flags.startpoint={'zeros','frana'};
[flags,kv,C]=ltfatarghelper({'C','tol','maxit'},definput,varargin);

if ~isnumeric(lambda) || lambda<0
    error('%s: ''lambda'' parameter must be a positive scalar.',...
          upper(mfilename));
end

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,~,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],[],upper(mfilename));

if W>1
    error('%s: Input signal can be single channel only.',upper(mfilename));
end

% Do a correct postpad so that we can call F.frana and F.frsyn
% directly.
L = framelength(F,Ls);
f = postpad(f,L);

% It is crucial to have the CANONICAL DUAL frame
have_dual = 0;
try
   % Try to create and accelerate the dual frame
   Fd = frameaccel(framedual(F),L);
   have_dual = 1;
catch
   warning(sprintf(['The canonical dual system is not available for a given ',...
       'frame.\n Using franaiter.'],upper(mfilename)));
    % err = lasterror.message;
    % The dual system cannot be created.
    % We will use franaiter instead
end

% Accelerate the frame
F = frameaccel(F,L);

% Cache the constant part
if have_dual
   cd = Fd.frana(f);
else
   cd = franaiter(F,f,'tol',1e-14); 
end

% Intermediate results
d = zeros(size(cd));
% Initial point
if flags.do_frana
   tc0 = F.frana(f);
elseif flags.do_zeros
   tc0 = zeros(size(cd));
end
c = tc0;

threshold = lambda/C;
relres = 1e16;
iter = 0;

  while ((iter < kv.maxit)&&(relres >= kv.tol))
       v = thresh(c + d,threshold,'soft') - d; 
       if have_dual
           d = cd - Fd.frana(F.frsyn(v));
       else
           d = cd - franaiter(F,F.frsyn(v),'tol',1e-14);
       end
       c = d + v;
       relres = norm(c(:)-tc0(:))/norm(tc0(:));
       tc0 = c;
       iter = iter + 1;
       if flags.do_print
         if mod(iter,kv.printstep)==0        
           fprintf('Iteration %d: relative error = %f\n',iter,relres);
         end;
       end;
  end
   
  
if nargout>3
    % Do a reconstruction with the original frame
    frec = postpad(F.frsyn(c),Ls);
    % Reformat to the original shape
    frec = assert_sigreshape_post(frec,dim,permutedsize,order);
end