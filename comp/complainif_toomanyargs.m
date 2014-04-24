function complainif_toomanyargs(fnargin,limit,callfun)

if nargin<3
    callfun = mfilename;
end

if fnargin>limit
   error('%s: Too many input arguments.',upper(callfun));
end
