function [y] = ups(x,varargin)

if(~isempty(varargin))
    L = varargin{1};
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

% if(sudy==0)
%    Lmat=[zeros(L-1,length(x)); xrow];
%     y = [Lmat(:); zeros(L-1,1)];
% elseif(sudy == 1)
%    Lmat=[ xrow; zeros(L-1,length(x));];
%    ytemp = Lmat(:);
%    y = ytemp(1:end-(L-1));
% elseif(sudy==2)
%     Lmat=[ xrow; zeros(L-1,length(x));];
%    y = Lmat(:);
% else
%          
%             
% end