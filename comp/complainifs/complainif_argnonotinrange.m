function complainif_argnonotinrange(fnargin,limitlo,limithi,callfun)
% This is just a replacement for 
% narginchk and nargchk

if nargin<3
    callfun = mfilename;
end

if fnargin<limitlo || fnargin>limithi
   error('%s: Number of arguments is not in the range [%d,%d].',upper(callfun),limitlo,limithi);
end
