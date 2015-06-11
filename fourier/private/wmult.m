function [m,t] = wmult(t)
% find multiplicities in a column vector t of real numbers
% and return the sorted vector t and multiplicities m
% as follows: if t=[3;5;5;5;8;9;9], then 
%                m=[0;0;1;2;0;0;1]
% This numbering is used to indicate the order of 
% derivatives of a function f, which are needed for
% the computation of divided differences with multiple knots


t = sort(t(:));
index = [0;diff(t)==0];
m = cumsum(index);
jump = find(diff(index)<0);

if isempty(jump) 
    return, 
end

z = zeros(length(t),1);
pt = m(jump);
pt(2:end) = diff(pt);
z(jump+1) = pt;
m = m - cumsum(z);
