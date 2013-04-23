function f=rand_rhpe(varargin)
%RAND_HPE  Random real HPE even function.
%   Usage:  f=rand_rhpe(s);
%
%   RAND_RHPE(s) generates an array of size s, which is real and HPE along the
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

f=tester_rand(shalf)-.5;

if rem(s(1),2)==0
  f=[f;flipud(f)];
else
  f=[f; ...
     tester_rand([1 s(2:end)])-.5; ...
     flipud(f)];
end;

