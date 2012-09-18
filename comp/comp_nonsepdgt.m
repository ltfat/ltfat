function c=comp_nonsepdgt(f,g,a,M,lt,do_timeinv,alg)
%NONSEPDGT  Compute Non-separable Discrete Gabor transform
%   Usage:  c=nonsepdgt(f,g,a,M,lt);
%
%   Input parameters:
%         f     : Input data.
%         g     : Window function.
%         a     : Length of time shift.
%         M     : Number of channels.
%         lt    : Lattice type
%         do_timeinv : Do a time invariant phase ?
%         alg   : Choose algorithm
%   Output parameters:
%         c     : $M \times N$ array of coefficients.
%
%   `nonsepdgt(f,g,a,M,lt)` computes the non-separable discrete Gabor
%   transform of the input signal *f* with respect to the window *g*,
%   time-shift *a*, number of channels *M* and lattice type *lt*.
%
%     * $alg=0$ : Choose the fastest algorithm
%
%     * $alg=0$ : Always choose multi-win
%
%     * $alg=1$ : Always choose shear
%
%   This is a computational subroutine, do not call it directly.

%   AUTHOR : Nicki Holighaus and Peter L. Soendergaard
%   TESTING: TEST_NONSEPDGT
%   REFERENCE: REF_NONSEPDGT

% Assert correct input.

L=size(f,1);
W=size(f,2);
b=L/M;
N=L/a;

% ----- algorithm starts here, split into sub-lattices ---------------

if (alg==1) || (alg==0 && lt(2)<=2) 

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
    
else

    b=L/M;
    N=L/a;
    s=b*lt(1)/lt(2);
    
    [s0,s1,X] = shearfind(a,b,s,L);

    c = comp_nonsepdgt_shear(f,g,a,M,s0,s1,X);
    
end;
