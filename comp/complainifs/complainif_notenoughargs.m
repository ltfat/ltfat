function complainif_notenoughargs(fnargin,limit,callfun)

if nargin<3
    callfun = mfilename;
end

if fnargin<limit
   error('%s: Not enough input arguments.',upper(callfun));
end
