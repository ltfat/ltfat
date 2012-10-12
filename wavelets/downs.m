function [y] = downs(x,varargin)

if(length(varargin)>0)
    L = varargin{1};
else
    L=2;
end

if(length(varargin)>1)
   sudLich = varargin{2};
else
   sudLich = 0;
end


sudy = mod(sudLich,2);

if(sudy==0)
    y = x(L:L:end);
else
    y = x(1:L:end);
end