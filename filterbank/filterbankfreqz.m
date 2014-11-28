function gf = filterbankfreqz(g,a,L,varargin)
%FILTERBANKFREQZ  Filterbank frequency responses
%   Usage: gf = filterbankfreqz(g,a,L)
%
%   `gf = filterbankfreqz(g,a,L)` calculates length *L* frequency responses
%   of filters in *g* and returns them as columns of *gf*.
%
%   If an optional parameters 'plot' is passed to `filterbankfreqz`,
%   the frequency responses will be plotted using |plotfft|. Any
%   optional parameter undestood by |plotfft| can be passed in addition
%   to 'plot'.
%

complainif_notenoughargs(nargin,3,'FILTERBANKFREQZ');
complainif_notposint(L,'L','FILTERBANKFREQZ');

% Wrap g if it is not a cell. The format of g will be checked
% further in filterbankwin.
if ~iscell(g)
    g = {g};
end

% The only place we need a
% It is necessary for cases when filterbank is given by e.g.
% {'dual',g}
g = filterbankwin(g,a,L,'normal');
M = numel(g);

G1=comp_transferfunction(g{1},L);
gf = zeros(L,M,class(G1));
gf(:,1) = G1;
for m=2:M
    gf(:,m) = cast(comp_transferfunction(g{m},L),class(G1));
end

% Search for the 'plot' flag
do_plot = any(strcmp('plot',varargin));

if do_plot
    % First remove the 'plot' flag from the arguments
    varargin(strcmp('plot',varargin)) = [];
    % and pass everything else to plotfft
    plotfft(gf,varargin{:});
end;

if nargout<1
    clear gf;
end

