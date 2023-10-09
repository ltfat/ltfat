function complainif_notvalidframeobj(F,callfun)

if nargin<2
    callfun = mfilename;
end

if ~isstruct(F) || ~isfield(F,'frana')
  error('%s: Argument F must be a frame definition structure.',...
        upper(callfun));
end;

