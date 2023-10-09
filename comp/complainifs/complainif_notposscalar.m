function complainif_notposscalar(var,varname,callfun)

if nargin<3
    callfun = mfilename;
end

if  isempty(var) || ~isscalar(var) || ~isnumeric(var) || var<=0 || ...
        ( isinf(var) )
   error('%s: %s should be a positive scalar.',upper(callfun),varname);
end
