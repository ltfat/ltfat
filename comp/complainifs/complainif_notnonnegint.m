function complainif_notnonnegint(var,varname,callfun)

if nargin<3
    callfun = mfilename;
end

if isempty(var) || ~isscalar(var) || ~isnumeric(var) || var<0 || rem(var,1)~=0 
   error('%s: %s should be non-negative integer.',upper(callfun),varname);
end
