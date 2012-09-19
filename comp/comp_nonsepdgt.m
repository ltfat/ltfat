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

c=zeros(M,N,W);

if (alg==1) || (alg==0 && lt(2)<=2) 

    mwin=comp_nonsepwin2multi(g,a,M,lt);
    
    % simple algorithm: split into sublattices
    
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

    s=b*lt(1)/lt(2);
    
    [s0,s1,br] = shearfind(a,b,s,L);
            
    ar = a*b/br;
    Mr = L/br;
    Nr = L/ar;
            
    ind = [ar 0; 0 br]*[kron((0:L/ar-1),ones(1,L/br));kron(ones(1,L/ar), ...
                                                      (0:L/br-1))];
    phs = reshape(mod((s1*(ind(1,:)-s0*ind(2,:)).^2+s0*ind(2,:).^2)*(L+1) ...
                      -2*(s0 ~= 0)*ind(1,:).*ind(2,:),2*L),L/br,L/ar);
    phs = exp(pi*1i*phs/L);
    
    ind_final = [1 0;-s1 1]*[1 -s0;0 1]*ind;
    ind_final = mod(ind_final,L);
    
    if s1 ~= 0
        g = comp_pchirp(L,s1).*g;
        f = repmat(comp_pchirp(L,s1),1,W).*f;
    end
    
    if s0 ~= 0
        g = comp_pchirp(L,-s0).*fft(g)/L;
        f = repmat(comp_pchirp(L,-s0),1,W).*fft(f);
        
        c_rect = comp_dgt(f,g,br,Nr,L,0);
        
        for w=0:W-1
            
            c(floor(ind_final(2,:)/b)+1+(ind_final(1,:)/a)*M+w*M*N) = ...
                c_rect(ind(1,[1:Mr,end:-1:Mr+1])/ar+1+(ind(2,:)/br)*Nr+w*M*N).* ...
                phs(ind(2,:)/br+1+(ind(1,:)/ar)*Mr);
        end;
    else 
        
        c_rect = comp_dgt(f,g,ar,Mr,L,0);
        
        for w=1:W
            c_rect(:,:,w) = phs.*c_rect(:,:,w);
        end;
        
        % The code line below this comment executes the commented for-loop
        % using Fortran indexing.
        %
        % for jj = 1:size(ind,2)        
        %     c(floor(ind_final(2,jj)/b)+1, ind_final(1,jj)/a+1) = ...
        %         c_rect(ind(2,jj)/br+1, ind(1,jj)/ar+1);
        % end
        for w=0:W-1
            c(floor(ind_final(2,:)/b)+1+(ind_final(1,:)/a)*M+w*M*N) = ... 
                c_rect(ind(2,:)/br+1+(ind(1,:)/ar)*Mr+w*M*N);
        end;
    end;
    
end;
