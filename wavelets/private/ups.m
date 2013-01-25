function [y] = ups(x,varargin)

if(~isempty(varargin))
    L = max([varargin{1},1]);
else
    L=2;
end

if(length(varargin)>1)
   sudLich = varargin{2};
else
   sudLich = 0;
end


sudy = sudLich;

%xrow =  x(:)';

if(sudy==0)
  y=zeros(L*length(x)+L-1,1);    
  y(L:L:end)=x;    
elseif(sudy==1)
  y=zeros(L*length(x)-(L-1),1);    
  y(1:L:end)=x;    
elseif(sudy==2)
  y=zeros(L*length(x),1);    
  y(L:L:end)=x;  
elseif(sudy==3)
  y=zeros(L*length(x),1);    
  y(1:L:end)=x;  
end

