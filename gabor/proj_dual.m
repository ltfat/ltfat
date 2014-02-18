function [ sol ] = proj_dual( x,~, param )
%PROJ_DUAL projection onto the dual windows space
%   Usage:  sol=proj_proj(x, ~, param)
%           [sol, infos]=proj_b2(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   `proj_dual(x,~,param)` solves:
%
%   .. sol = argmin_{z} ||x - z||_2^2   s.t.  A z=y
%
%   .. math::  sol = \min_z ||x - z||_2^2 \hspace{1cm} s.t. \hspace{1cm} A z= y
%
%   param is a Matlab structure containing the following fields:
%
%   * *param.y* : measurements (default: 0).
%
%   * *param.A* : Matrix (default: Id).
%
%   * *param.AAtinv* : $(A A^*)^(-1)$ Define this parameter to speed up computation.
%
%   * *param.verbose* : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%
%   infos is a Matlab structure containing the following fields:
%
%   * *infos.algo* : Algorithm used
%
%   * *infos.iter* : Number of iteration
%
%   * *infos.time* : Time of execution of the function in sec.
%
%   * *infos.final_eval* : Final evaluation of the function
%
%   * *infos.crit* : Stopping critterion used 
%
%   * *infos.residue* : Final residue  
%
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   See also:  prox_l2 proj_b1
%
%   References: fadili2009monotone

%
% Author: Nathanael Perraudin
% Date: Feb 20, 2013
%

% Start the time counter
t1 = tic;

% Optional input arguments
if ~isfield(param, 'y'), param.y = 0; end
if ~isfield(param, 'A'), param.A = eye(length(x)); end
if ~isfield(param, 'AAtinv'), param.AAtinv=pinv(A*At); end
if ~isfield(param, 'verbose'), param.verbose = 1; end


% Projection  
   
sol = x - param.A'*param.AAtinv*(param.A*x-param.y);
crit = 'TOL_EPS'; iter = 0; u = NaN;
    


% Log after the projection onto the L2-ball
error=norm(param.y-param.A *sol );
if param.verbose >= 1
    fprintf(['  Proj. dual windows: ||y-Ax||_2 = %e,', ...
        ' %s, iter = %i\n'],error , crit, iter);
end

% Infos about algorithm
infos.algo=mfilename;
infos.iter=iter;
infos.final_eval=error;
infos.crit=crit;
infos.residue=u;
infos.time=toc(t1);

end


