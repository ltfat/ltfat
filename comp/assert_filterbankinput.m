function [M,longestfilter]=assert_filterbankinput(g,a);

  f=dbstack;  
  callfun=f(2).name;

  if ~iscell(g)
    error('%s: g must be a cell-array.',upper(callfun));
  end;
  
  if ~isnumeric(a) || ~isscalar(a) || a<=0
    error('%s: a must be a positive scalar.',upper(callfun));
  end;

  M=numel(g);
  
  longestfilter=max(cellfun(@numel,g));