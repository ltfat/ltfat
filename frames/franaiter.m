function [c,relres,iter]=franaiter(F,c,varargin)
%FRANAITER  Iterative analysis frame inversion
%   Usage:  f=franaiter(F,c);
%
%   Input parameters:
%         F       : Frame   
%         c       : Array of coefficients.
%         Ls      : length of signal.
%   Output parameters:
%         f       : Signal.
%         relres  : Vector of residuals.
%         iter    : Number of iterations done.
%
%   `c=franaiter(F,f)` iteratively inverts the analysis frame of *F* using a
%   least-squares method.
%
%   `[c,relres,iter]=franaiter(...)` additionally returns the residuals in a
%   vector *relres* and the number of iteration steps *iter*.
%  
%   **Note:** If it is possible to explicitly calculate the canonical dual
%   frame then this is usually a much faster method than invoking
%   `franaiter`.
%
%   `franaiter` takes the following parameters at the end of the line of
%   input arguments:
%
%     'tol',t      Stop if relative residual error is less than the
%                  specified tolerance. Default is 1e-9 
%
%     'maxit',n    Do at most n iterations.
%
%     'pcg'        Solve the problem using the Conjugate Gradient
%                  algorithm. This is the default.
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
%      F=newframe('dgtreal','none','gauss',10,20);
%      [c,relres,iter]=franaiter(F,bat);
%      r=frsyn(F,c);
%      norm(bat-r)/norm(bat)
%
%   See also: newframe, frana, frsyn, frsyniter
  
% AUTHORS: Nathanael Perraudin and Peter L. Søndergaard
    
  if nargin<2
    error('%s: Too few input parameters.',upper(mfilename));
  end;
  
  definput.keyvals.Ls=[];
  definput.keyvals.tol=1e-9;
  definput.keyvals.maxit=100;
  definput.flags.alg={'pcg'};
  definput.keyvals.printstep=10;
  definput.flags.print={'quiet','print'};

  [flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);
  
  % Determine L from the first vector, it must match for all of them.
  L=framelengthcoef(F,size(c,1));
    
  if flags.do_pcg

      A=@(x) frsyn(F,frsynadj(F,x));
                  
      % It is possible to specify the initial guess
      [fout,flag,relres,iter]=pcg(A,c,kv.tol,kv.maxit);
      c=frsynadj(F,fout);
  end;
           
end


    