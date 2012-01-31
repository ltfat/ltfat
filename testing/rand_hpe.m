function f=rand_hpe(varargin)
%RAND_HPE  Random HPE even function.
%   Usage:  f=rand_hpe(s);
%
%   RAND_HPE(s) generates an array of size s, which is HPE along the
%   first dimension.

if isnumeric(varargin);
  s=varargin;
else
  s=cell2mat(varargin);
end;


if length(s)==1
  error('To avoid confusion, the size must be at least two-dimensional.');
end;

shalf=s;

shalf(1)=floor(shalf(1)/2);

f=(randn(shalf)-.5)+(randn(shalf)-.5)*i;

if rem(s(1),2)==0
  f=[f;flipud(conj(f))];
else
  f=[f; ...
     randn([1 s(2:end)])-.5; ...
     flipud(conj(f))];
end;

