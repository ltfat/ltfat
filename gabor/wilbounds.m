function [AF,BF]=wilbounds(g,M,varargin)
%WILBOUNDS  Calculate frame bounds of Wilson basis
%   Usage:  [AF,BF]=wilbounds(g,M)
%           [AF,BF]=wilbounds(g,M,L)
%
%   Input parameters:
%           g       : Window function.
%           M       : Number of channels.
%           L       : Length of transform to do (optional)
%   Output parameters:
%           AF,BF   : Frame bounds.
%          
%   `wilbounds(g,M)` calculates the frame bounds of the Wilson frame operator
%   of the Wilson basis with window *g* and *M* channels.
%
%   `[A,B]=wilbounds(g,a,M)` returns the frame bounds *A* and *B* instead of
%   just the ratio.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin| for more details.
%
%   If the length of *g* is equal to $2\cdot M$ then the input window is
%   assumed to be a FIR window. Otherwise the smallest possible transform
%   length is chosen as the window length.
%
%   If the optional parameter *L* is specified, the window is cut or
%   zero-extended to length *L*.
%
%   See also: wilwin, gabframebounds

if nargin<2
   error('Too few input parameters.');    
end;

definput.keyvals.L=[];

[flags,keyvals]=ltfatarghelper({'L'},definput,varargin);
L=keyvals.L;

if size(M,1)>1 || size(M,2)>1
  error('M must be a scalar');
end;

if rem(M,1)~=0
  error('M must be an integer.')
end;

[g,info]=wilwin(g,M,L,'WILBOUNDS');
Lwindow=length(g);

[b,N,L]=assert_L(Lwindow,Lwindow,L,M,2*M,'WILBOUNDS');

g=fir2long(g,L);

a=M;

N=L/a;

if rem(N,2)==1
  error('L/M must be even.');
end;


% Get the factorization of the window.
gf=comp_wfac(g,a,2*M);

% Compute all eigenvalues.
lambdas=comp_gfeigs(gf,L,a,2*M);

% Min and max eigenvalue.
AF=lambdas(1);
BF=lambdas(size(lambdas,1));

% Divide by 2 (only difference to gfeigs).
AF=AF/2;
BF=BF/2;

if nargout<2
  AF=BF/AF;
end;
