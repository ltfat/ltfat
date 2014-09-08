function f=divdiff_vector(x,y)

% f is the divided difference of the data in y
% repeated x-values are allowed if x is increasing
% then the data values in y are taken as function,
% derivative and higher derivative values, respectively.
% y may be a matrix, divdiff is taken per column

x=x(:); % writes all x-values into a column-vector

[x,i]=sort(x);
y=y(i,:); % sorts the data vector/matrix the same way as the x-values

[l,l1]=size(y); % size of the columns and rows of y
if (l ~= length(x))
   error('x and y must have the same dimension.')
end

%use the spline toolbox to find multiplicities
[m,t]=myknt2mlt(x); % see description for knt2mlt
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