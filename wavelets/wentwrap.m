function E = wentwrap(x,fname,varargin)

if(iscell(x))
    xtmp = x;
    x = [];
    for ii=1:numel(xtmp)
       x = [x,xtmp{ii}(:)]; 
    end
end

if(~isempty(varargin))
   E = feval(sprintf('%s%s','went_',fname),x,varargin{:});
else
   E = feval(sprintf('%s%s','went_',fname),x);
end