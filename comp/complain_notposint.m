function complain_notposint(var,varname)

if isempty(var) || ~isscalar(var) || ~isnumeric(var) || var<=0 || rem(var,1)~=0 
   error('%s: %s should be a positive integer.',upper(mfilename),varname)
end
