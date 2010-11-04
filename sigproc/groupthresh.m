function [xo]=groupthresh(xi,lambda,varargin)
%GROUPTHRESH   group (hard/soft) thresholding
%   Usage:  xo=groupthresh(x,lambda);
%
%   GROUPTHRESH(x,lambda) will perform hard group thresholding on x, with
%   threshold lambda xi is a two-dimensional array, the first dimension
%   labelling groups, and the second one labelling members All coefficients
%   within a given group are shrunk according to the value of the L2 norm of
%   the group in comparison to the threshold lambda
%
%   GROUPTHRESH(x,lambda,'soft') will do the same using soft
%   thresholding.
%
%   GROUPTHRESH takes the following flags at the end of the line of input
%   arguments:
%
%-     'hard'   - Perform hard thresholding. This is the default.
%
%-     'soft'   - Perform soft thresholding.  
%
%-     'full'   - Returns the output as a full matrix. This is the default.
%
%-     'sparse' - Returns the output as a sparse matrix.
%  
%   See also:  gabgrouplasso
%
%   Demos:  demo_audioshrink
%
%R  Kowalski08sparsity Kowalski09mixed

%   AUTHOR : Bruno Torresani.  
%   REFERENCE: OK
 
if nargin<2
  error('Too few input parameters.');k
end;

if (prod(size(lambda))~=1 || ~isnumeric(lambda))
  error('lambda must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.flags.iofun={'hard','soft'};
definput.flags.outclass={'full','sparse'};

[flags,keyvals]=ltfatarghelper({},definput,varargin);

NbGroups = size(xi,1);
NbMembers = size(xi,2);

xo = zeros(size(xi));

for g=1:NbGroups,
    threshold = norm(xi(g,:));
    mask = (1-lambda/threshold);
    mask = mask * (mask>0);
    if flags.do_hard
      mask = (mask>0);
    else      
      mask = mask * (mask>0);
    end;
    xo(g,:) = xi(g,:) * mask;
end
