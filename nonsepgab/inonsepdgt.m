function [f,g]=inonsepdgt(coef,g,a,lt,varargin)
%INONSEPDGT  Inverse discrete Gabor transform
%   Usage:  f=inonsepdgt(c,g,a,lt);
%           f=inonsepdgt(c,g,a,lt,Ls);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         lt    : Lattice type
%         Ls    : length of signal.
%   Output parameters:
%         f     : Signal.
%
%   `inonsepdgt(c,g,a,lt)` computes the Gabor expansion of the input
%   coefficients *c* with respect to the window *g*, time shift *a* and
%   lattice type *lt*. The number of channels is deduced from the size of
%   the coefficients *c*.
%
%   `inonsepdgt(c,g,a,lt,Ls)` does as above but cuts or extends *f* to length
%   *Ls*.
%
%   `[f,g]=inonsepdgt(...)` additionally outputs the window used in the
%   transform. This is useful if the window was generated from a description
%   in a string or cell array.
%
%   For perfect reconstruction, the window used must be a dual window of the
%   one used to generate the coefficients.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%
%   If *g* is a row vector, then the output will also be a row vector. If *c* is
%   3-dimensional, then `inonsepdgt` will return a matrix consisting of one column
%   vector for each of the TF-planes in *c*.
%
%   Assume that `f=inonsepdgt(c,g,a,lt,L)` for an array *c* of size $M\times N$.
%   Then the following holds for $k=0,\ldots,L-1$:
% 
%   ..          N-1 M-1          
%     f(l+1)  = sum sum c(m+1,n+1)*exp(2*pi*i*m*l/M)*g(l-a*n+1)
%               n=0 m=0          
%
%   .. math:: f(l+1) = \sum_{n=0}^{N-1}\sum_{m=0}^{M-1}c(m+1,n+1)e^{2\pi iml/M}g(l-an+1)
%
%   `inonsepdgt` takes the following flags at the end of the line of input
%   arguments:
%
%     'freqinv'  Compute an `inonsepdgt` using a frequency-invariant phase. This
%                is the default convention described above.
%
%     'timeinv'  Compute an `inonsepdgt` using a time-invariant phase. This
%                convention is typically used in filter bank algorithms.
%
%   See also:  dgt, gabwin, dwilt, gabtight

%   AUTHOR : Nicki Holighaus and Peter Soendergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: OK

% Check input paramameters.

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if numel(g)==1
  error('g must be a vector (you probably forgot to supply the window function as input parameter.)');
end;

definput.keyvals.Ls=[];
definput.flags.nsalg={'multiwin','shear'};
[flags,kv,Ls]=ltfatarghelper({'Ls'},definput,varargin);

wasrow=0;

M=size(coef,1);
N=size(coef,2);
W=size(coef,3);

% use assert_squarelat to check a and the window size.
assert_squarelat(a,M,1,'INONSEPDGT');

L=N*a;

[g,info] = comp_window(g,a,M,L,lt,'INONSEPDGT');

assert_L(L,size(g,1),L,a,M,'INONSEPDGT');

if flags.do_multiwin
    % ----- algorithm starts here, split into sub-lattices ---------------
    
    mwin=comp_nonsepwin2multi(g,a,M,lt);
    
    % phase factor correction (backwards), for more information see 
    % analysis routine
    
    E = exp(2*pi*i*a*kron(0:N/lt(2)-1,ones(1,lt(2))).*...
            rem(kron(ones(1,N/lt(2)), 0:lt(2)-1)*lt(1),lt(2))/M);
    for w=1:W
        coef(:,:,w) = coef(:,:,w).*repmat(E,M,1);
    end;
    
    % simple algorithm: split into sublattices and add the result from eacg
    % sublattice.
    f=zeros(L,W);
    for ii=0:lt(2)-1
        % Extract sublattice
        sub=coef(:,ii+1:lt(2):end);
        f=f+comp_idgt(sub,mwin(:,ii+1),lt(2)*a,M,L,0);  
    end;

end;


if flags.do_shear
    V = latticetype2matrix(L,a,M,lt);
    
    b = L/M;
    [s0,s1,X] = shearfind(a,b,V(2,1),L);
    
    br = X;
    ar = a*b/X;
    
    ind = [ar 0; 0 br]*[kron((0:L/ar-1),ones(1,L/br));kron(ones(1,L/ar),(0:L/br-1))];
    phs = reshape(mod((s1*(ind(1,:)-s0*ind(2,:)).^2+s0*ind(2,:).^2)*(L+1),2*L),L/br,L/ar);
    phs = exp(-pi*1i*phs/L);
    
    ind_final = [1 0;-s1 1]*[1 -s0;0 1]*ind;
    ind_final = mod(ind_final,L);
    
    c2 = zeros(L/br,L/ar);
    
    %This loop should be made a single step, but
    %   c2(floor(ind_final(2,:)/b)+1,ind_final(1,:)/a+1) = c(ind(2,:)/br+1,ind(1,:)/ar+1);
    %seems to result in rubbish.
    
    for jj = 1:size(ind,2)
        c2(ind(2,jj)/br+1,ind(1,jj)/ar+1) = coef(floor(ind_final(2,jj)/b)+1,ind_final(1,jj)/a+1);
    end
    
    coef = c2;
    coef = phs.*coef;
    
    if s1 ~= 0
        g = pchirp(L,s1).*g;
    end
    
    if s0 ~= 0
        g = ifft(pchirp(L,-s0).*fft(g));    
    end
    
    f = comp_idgt(coef,g,a*b/X,L/X,L,0);
    
    if s0 ~= 0
        f = ifft(pchirp(L,s0).*fft(f));    
    end
    
    if s1 ~= 0
        f = pchirp(L,-s1).*f;
    end        
end;
    
% Cut or extend f to the correct length, if desired.
if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

f=comp_sigreshape_post(f,Ls,wasrow,[0; W]);
