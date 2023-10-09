function complainif_notposint(var,varname,callfun)

if nargin<3
    callfun = mfilename;
end

if  isempty(var) || ~isscalar(var) || ~isnumeric(var) || var<=0 || ...
        ( ~isinf(var) && rem(var,1)~=0 )
   error('%s: %s should be a positive integer.',upper(callfun),varname);
end
