function [c,Ls,g]=nonsepdgt(f,g,a,M,lt,varargin)
%NONSEPDGT  Non-separable Discrete Gabor transform
%   Usage:  c=nonsepdgt(f,g,a,M,lt);
%           c=nonsepdgt(f,g,a,M,lt,L);
%           [c,Ls]=nonsepdgt(f,g,a,M,lt);
%           [c,Ls]=nonsepdgt(f,g,a,M,lt,L);
%
%   Input parameters:
%         f     : Input data.
%         g     : Window function.
%         a     : Length of time shift.
%         lt    : Lattice type
%         M     : Number of channels.
%         L     : Length of transform to do.
%   Output parameters:
%         c     : $M \times N$ array of coefficients.
%         Ls    : Length of input signal.
%
%   `nonsepdgt(f,g,a,M,lt)` computes the non-separable discrete Gabor
%   transform of the input signal *f* with respect to the window *g*,
%   time-shift *a*, number of channels *M* and lattice type *lt*. The output is a
%   vector/matrix in a rectangular layout.
%
%   The lattice type *lt* is a $1 \times 2$ vector $[lt_1,lt_2]$ denoting an
%   irreducible fraction $lt_1/lt_2$. This fraction describes the distance
%   in frequency (counted in frequency channels) that each coefficient is
%   offset when moving in time by the time-shift of *a*. Some examples:
%   `lt=[0 1]` defines a square lattice, `lt=[1 2]` defines the quinqux
%   (almost hexagonal) lattice, `lt=[1 3]` describes a lattice with a
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
%   `nonsepdgt(f,g,a,M,lt,L)` computes the Gabor coefficients as above, but
%   does a transform of length *L*. f will be cut or zero-extended to length
%   *L* before the transform is done.
%
%   `[c,Ls]=nonsepdgt(...)` additionally returns the length of the input
%   signal *f*. This is handy for reconstruction::
%
%               [c,Ls]=nonsepdgt(f,g,a,M,lt);
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
%   define $N=L/a$.  The output from `c=nonsepdgt(f,g,a,M,lt)` is then given
%   by:
%
%   ..              L-1 
%      c(m+1,n+1) = sum f(l+1)*conj(g(l-a*n+1))*exp(-2*pi*i*(m+w(n))*l/M), 
%                   l=0  
%
%   .. math:: c\left(m+1,n+1\right)=\sum_{l=0}^{L-1}f(l+1)\overline{g(l-an+1)}e^{-2\pi il(m+w(n))/M}
%
%   where $m=0,\ldots,M-1$ and $n=0,\ldots,N-1$ and $l-an$ is computed
%   modulo *L*.  The additional offset $w$ is given by
%   $w(n)=\mod(n\cdot lt_1,lt_2)/lt_2$ in the formula above.
%
%   Coefficient layout:
%   -------------------
%
%   The following code generates plots which show the coefficient layout
%   and enumeration of the first 4 lattices in the time-frequecy plane:::
%
%     a=6;
%     M=6;
%     L=36;
%     b=L/M;
%     N=L/a;
%     cw=3;
%     ftz=12;
%     
%     [x,y]=meshgrid(a*(0:N-1),b*(0:M-1));
%
%     lt1=[0 1 1 2];
%     lt2=[1 2 3 3];
%
%     for fignum=1:4
%       figure;
%       z=y;
%       if lt2(fignum)>0
%         z=z+mod(lt1(fignum)*x/lt2(fignum),b);
%       end;
%       for ii=1:M*N
%         text(x(ii)-cw/4,z(ii),sprintf('%2.0i',ii),'Fontsize',ftz);
%         rectangle('Curvature',[1 1], 'Position',[x(ii)-cw/2,z(ii)-cw/2,cw,cw]);
%       end;
%       axis([-cw L -cw L]);
%       axis('square');
%       title(sprintf('s=[%i %i]',lt1(fignum),lt2(fignum)),'Fontsize',ftz);
%     end;
%
%   See also:  nonsepgabdual

%   AUTHOR : Nicki Holighaus and Peter L. Soendergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

% Assert correct input.

if nargin<5
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.flags.nsalg={'multiwin','shear'};
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);

if isempty(L)
  L=nonsepdgtlengthsignal(Ls,a,M,lt);
else
  Lcheck=nonsepdgtlengthsignal(L,a,M,lt);
  if Lcheck~=L
    error('%s: Invalid transform size L',upper(mfilename));
  end;
end;

b=L/M;
N=L/a;

[g,info] = comp_window(g,a,M,L,lt,'NONSEPDGT');

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


% ----- algorithm starts here, split into sub-lattices ---------------

if flags.do_multiwin
    mwin=comp_nonsepwin2multi(g,a,M,lt);
    
    % simple algorithm: split into sublattices
    c=zeros(M,N,W);
    
    for ii=0:lt(2)-1
        c(:,ii+1:lt(2):end,:)=comp_dgt(f,mwin(:,ii+1),lt(2)*a,M,L,0);% .* e;
    end;
    
    % phase factor correction 
    % Explanation:
    %         For c(m,l*lt(2)+n), the missing phase factor is given by 
    %       exp(-2*pi*i*l*lt(2)*a*rem(n*lt(1),lt(2))/lt(2)/M),
    %       independently of m. Furthermore lt(2)/lt(2) cancels.
    %         The Kronecker products in the following setup a matrix 
    %       the size of the coefficient matrix and compute the product.
    
    E = exp(-2*pi*i*a*kron(0:N/lt(2)-1,ones(1,lt(2))).*...
            rem(kron(ones(1,N/lt(2)),0:lt(2)-1)*lt(1),lt(2))/M);
    
    for w=1:W
        c(:,:,w) = c(:,:,w).*repmat(E,M,1);
    end;
    
end;

if flags.do_shear

    b=L/M;
    N=L/a;
    s=b*lt(1)/lt(2);
    
    [s0,s1,X] = shearfind(a,b,s,L);
    
    if s1 ~= 0
        g = pchirp(L,s1).*g;
        f = repmat(pchirp(L,s1),1,W).*f;
    end
    
    if s0 ~= 0
        g = ifft(pchirp(L,-s0).*fft(g));
        f = ifft(repmat(pchirp(L,-s0),1,W).*fft(f));
    end
    
    br = X;
    ar = a*b/X;
    Mr = L/br;
    c = comp_dgt(f,g,ar,Mr,L,0);
    
    ind = [ar 0; 0 br]*[kron((0:L/ar-1),ones(1,L/br));kron(ones(1,L/ar),(0:L/br-1))];
    phs = reshape(mod((-s1*(-s0*ind(2,:)+ind(1,:)).^2-s0*ind(2,:).^2)*(L+1),2*L),L/br,L/ar);
    phs = exp(-1*pi*1i*phs/L);
    
    c = phs.*c;
    
    ind_final = [1 0;-s1 1]*[1 -s0;0 1]*ind;
    ind_final = mod(ind_final,L);

    c2 = zeros(M,N);
    
    % The code line below this comment executes the commented for-loop
    % using Fortran indexing.
    %
    % for jj = 1:size(ind,2)        
    %     c2(floor(ind_final(2,jj)/b)+1, ind_final(1,jj)/a+1) = ...
    %         c(ind(2,jj)/br+1, ind(1,jj)/ar+1);
    % end
    c2(floor(ind_final(2,:)/b)+1+(ind_final(1,:)/a)*M) = c(ind(2,:)/br+1+(ind(1,:)/ar)*Mr);

    c = c2;    
    
end;
