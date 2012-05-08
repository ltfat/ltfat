function [c,Ls,g]=dgt(f,g,a,M,s,varargin)
%NONSEPDGT  Non-separable Discrete Gabor transform
%   Usage:  c=nonsepdgt(f,g,a,M,s);
%           c=nonsepdgt(f,g,a,M,s,L);
%           [c,Ls]=nonsepdgt(f,g,a,M,s);
%           [c,Ls]=nonsepdgt(f,g,a,M,s,L);
%
%   Input parameters:
%         f     : Input data.
%         g     : Window function.
%         a     : Length of time shift.
%         s     : Shear of lattice.
%         M     : Number of channels.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : $M \times N$ array of coefficients.
%         Ls    : Length of input signal.
%
%   `nonsepdgt(f,g,a,M,s)` computes the non-separable discrete Gabor
%   transform of the input signal *f* with respect to the window *g*,
%   time-shift *a*, number of channels *M* and shear *s*. The output is a
%   vector/matrix in a rectangular layout.
%
%   The shear *s* is a $1 \times 2$ vector $\[s_1 s_2\]$ denoting an
%   irreducible fraction $s_1/s_2$. This fraction describes the distance
%   in frequency (counted in frequency channels) that each coefficient is
%   offset when moving in time by the time-shift of *a*. Some examples:
%   `s=[0 1]` defines a square lattice, `s=[1 2]` defines the quinqux
%   (almost hexagonal) lattice, `s=[1 3]` describes a lattice with a
%   $1/3$ frequency offset for each time shift and so forth.
%
%   The length of the transform will be the smallest multiple of *a* and *M*
%   that is larger than the signal. *f* will be zero-extended to the length of
%   the transform. If *f* is a matrix, the transformation is applied to each
%   column. The length of the transform done can be obtained by
%   `L=size(c,2)*a;`
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%
%   `nonsepdgt(f,g,a,M,s,L)` computes the Gabor coefficients as above, but
%   does a transform of length *L*. f will be cut or zero-extended to length
%   *L* before the transform is done.
%
%   `[c,Ls]=nonsepdgt(...)` additionally returns the length of the input
%   signal *f*. This is handy for reconstruction::
%
%               [c,Ls]=nonsepdgt(f,g,a,M,s);
%               fr=inonsepdgt(c,gd,a,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is,
%   provided that *gd* is a dual window of *g*.
%
%   `[c,Ls,g]=nonsepdgt(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   The non-separable discrete Gabor transform is defined as follows:
%   Consider a window *g* and a one-dimensional signal *f* of length *L* and
%   define $N=L/a$.  The output from `c=nonsepdgt(f,g,a,M,s)` is then given
%   by: FIXTHIS
%
%   ..              L-1 
%      c(m+1,n+1) = sum f(l+1)*conj(g(l-a*n+1))*exp(-2*pi*i*m*l/M), 
%                   l=0  
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)\overline{g(l-an+1)}e^{-2\pi ilm/M}
%
%   where $m=0,\ldots,M-1$ and $n=0,\ldots,N-1$ and $l-an$ is computed modulo *L*.
%
%   Coefficient layout:
%   -------------------
%
%   The following code generates plots which show the coefficient layout
%   and enumeration of the first 4 lattices in the time-frequecy plane:::
%
%     a=6;
%     M=6;
%     N=6;
%     L=N*a;
%     b=L/M;
%     
%     [x,y]=meshgrid(a*(0:N-1),b*(0:M-1));
%
%     s1=[0 1 1 2];
%     s2=[0 2 3 3];
%
%     for fignum=1:4
%       figure;
%       z=y;
%       if s2(fignum)>0
%         z=z+mod(s1(fignum)*x/s2(fignum),b);
%       end;
%       for ii=1:M*N
%         text(x(ii),z(ii),sprintf('%2.0i',ii));
%         rectangle('Curvature',[1 1], 'Position',[x(ii)-.5,z(ii)-.5,3,3]);
%       end;
%       axis([0 L 0 L]);
%       axis('square');
%       title(sprintf('s=[%i %i]',s1(fignum),s2(fignum)));
%     end;
%
%   See also:  inonsepdgt, nonsepgabwin, nonsepgabdual

% Assert correct input.

if nargin<5
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
[flags,kv]=ltfatarghelper({'L'},definput,varargin);

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

if ~isempty(L)
  if (prod(size(L))~=1 || ~isnumeric(L))
    error('%s: L must be a scalar',callfun);
  end;
  
  if rem(L,1)~=0
    error('%s: L must be an integer',callfun);
  end;
end;

% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,callfun,0);

if isnumeric(g)
  if ~isvector(g)
    error('%s: g must be a vector',upper(callfun));
  end;
  Lwindow=length(g);
else
  Lwindow=0;
end;


if isempty(L)
  % Smallest length transform. FIXME
  Lsmallest=lcm(a,M);

  % Choose a transform length larger than both the length of the
  % signal and the window.
  L=ceil(Ls/Lsmallest)*Lsmallest;
else

  if rem(L,M)~=0
    error('%s: The length of the transform must be divisable by M = %i',...
          callfun,M);
  end;

  if rem(L,a)~=0
    error('%s: The length of the transform must be divisable by a = %i',...
          callfun,a);
  end;

  if L<Lwindow
    error('%s: Window is too long.',callfun);
  end;

end;

b=L/M;
N=L/a;

[g,info] = comp_window(g,a,M,L,s,'NONSEPDGT');

if (info.isfir)  
  if info.istight
    g=g/sqrt(2);
  end;  
end;

f=postpad(f,L);

% If the signal is single precision, make the window single precision as
% well to avoid mismatches.
if isa(f,'single')
  g=single(g);
end;


% FIXME ----- algorithm starts here ---------------

% simple algorithm: split into sublattices
c=zeros(M,N,W);
for ii=1:s(2)
  % Modulate the window
  gw=g.*expwave(L,a*s1
  
end;

