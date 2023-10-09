function assert_squarelat(a,M,R,callfun,flag)
%ASSERT_SQUARELAT  Validate lattice and window size.
%   Usage:  assert_squarelat(a,M,R,callfun,flag);
%
%   Input parameters:
%         a       : Length of time shift.
%         M       : Number of modulations.
%         R       : Number of multiwindows.
%         callfun : Name of calling function.
%         flag    : See below.
%         
%  if flag>0 test if system is at least critically sampled.
%
%  This routine deliberately checks the validity of M before a, such
%  that it can be used for DWILT etc., where you just pass a=M.
%
%  
if nargin==4
  flag=0;
end;

if  (prod(size(M))~=1 || ~isnumeric(M))
  error('%s: M must be a scalar',callfun);
end;

if (prod(size(a))~=1 || ~isnumeric(a))
  error('%s: a must be a scalar',callfun);
end;

if rem(M,1)~=0
  error('%s: M must be an integer',callfun);
end;

if rem(a,1)~=0
  error('%s: a must be an integer',callfun);
end;

if flag>0
  if a>M*R
    error('%s: The lattice must not be undersampled',callfun);
  end;
end;


