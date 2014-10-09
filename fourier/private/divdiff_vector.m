function f=divdiff_vector(x,m,y)

% f is the divided difference of the columns in data matrix y
% with knots in x. m is an integer vector of multiplicities in x.
% x must be increasing.
%
% Repeated x-values are allowed.
% The rows in y are function, derivative and higher 
% derivative values, respectively.
% 
% Example of input parameters: 
%  x=[1;5;5;5;7;8;8]; m=[0;0;1;2;0;0;1];
%  y(k,:) contains the value of f^(m(k))(x(k))/(k!) 

[l,l1]=size(y); % size of the columns and rows of y
if (l ~= length(x))
   error('x and y must have the same dimension.')
end

index=find(m>0);
ys=y; %save for derivative values
y(index,:)=y(index-m(index),:); % fill function values

for k=1:l-1
   % for differences when consecutive x-values are different
   xdiff=(x(k+1:l)-x(1:l-k));
   ydiff=y(2:end,:)-y(1:end-1,:);
   index=find(m(k+1:l)<k);
   if isempty(index) ~= 1
       y(index,:)=ydiff(index,:)./repmat(xdiff(index),1,l1);
   end
   % for other values insert the proper derivative
   index=find(m>k-1);
   y(index-k,:)=ys(index-m(index)+k,:);
end

f=y(1,:);