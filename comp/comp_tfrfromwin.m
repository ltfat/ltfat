function [tfr,gamma] = comp_tfrfromwin(g, atheight)


gl = winwidthatheight(g, 1e-10);
w  = winwidthatheight(g, atheight);

Cg = -pi/4*(w/gl)^2/log(atheight);
gamma = Cg*gl^2;
tfr = @(L) gamma/L;


function width = winwidthatheight(gnum,atheight)
%WINWIDTHATHEIGHT Window width at height
%   Usage: width = winwidthatheight(gnum, height)
%
%   Input parameters:
%         gnum      : Window.
%         atheight  : Relative height.
%   Output parameters:
%         width   : Window width in samples.
%
%   `winwidthatheight(gnum,atheight)` computes width of a window *gnum* at
%   the relative height *atheight*. *gnum* must be a numeric vector as
%   returned from |gabwin|. If *atheight* is an array, width will have the
%   same shape with correcpondng values.
%

% AUTHOR: Zdenek Prusa

if ~isnumeric(gnum) || isempty(gnum) || ~isvector(gnum) || ~isreal(gnum)
    error('%s: gnum must be a numeric vector.', upper(mfilename));
end

if isempty(atheight) || any(atheight) > 1 || any(atheight) < 0
    error('%s: h must be in the interval [0-1].', upper(mfilename));
end

width = zeros(size(atheight));
for ii=1:numel(atheight)
    gl = numel(gnum);
    gmax = max(gnum);
    frac=  1/atheight(ii);
    fracofmax = gmax/frac;

    ind =find(gnum(1:floor(gl/2)+1)==fracofmax,1,'first');
    if isempty(ind)
        %There is no sample exactly half of the height
        ind1 = find(gnum(1:floor(gl/2)+1)>fracofmax,1,'last');
        ind2 = find(gnum(1:floor(gl/2)+1)<fracofmax,1,'first');
        if isempty(ind2)
            width(ii) = gl;
        else
            rest = 1-(fracofmax-gnum(ind2))/(gnum(ind1)-gnum(ind2));
            width(ii) = 2*(ind1+rest-1);
        end
    else
        width(ii) = 2*(ind-1);
    end
end