function f=comp_idgt(coef,g,a,lt,phasetype,algns)
%COMP_IDGT  Compute IDGT
%   Usage:  f=comp_idgt(c,g,a,lt,phasetype);
%
%   Input parameters:
%         c     : Array of coefficients.
%         g     : Window function.
%         a     : Length of time shift.
%         lt    : Lattice type
%         phasetype : Type of phase
%   Output parameters:
%         f     : Signal.
%
%   Value of the algorithm chooser
%
%     * $algns=0$ : Choose the fastest algorithm
%
%     * $algns=0$ : Always choose multi-win
%
%     * $algns=1$ : Always choose shear

% AUTHOR : Peter L. Søndergaard.

M=size(coef,1);
N=size(coef,2);
W=size(coef,3);

L=N*a;

Lwindow=size(g,1);

% FIXME : Calls non-comp function 
if phasetype==1
    coef=phaseunlock(coef,a,'lt',lt);
end;

if lt(2)==1

    if L==Lwindow
        % Do full-window algorithm.
        
        % Get the factorization of the window.
        gf = comp_wfac(g,a,M);      

        % FIXME: This line is necessary because the mex and oct interfaces expect
        % a matrix as input.
        coef=reshape(coef,M,prod(size(coef))/M);

        
        % Call the computational subroutine.
        f  = comp_idgt_fac(coef,gf,L,a,M);
        
    else
        
        % Do filter bank algorithm.
        % Call the computational subroutine.

        % FIXME: This line is necessary because the mex and oct interfaces expect
        % a matrix as input.
        coef=reshape(coef,M,prod(size(coef))/M);

        f=comp_idgt_fb(coef,g,L,a,M);
        
    end;
    
else
    
    if (algns==1) || (algns==0 && lt(2)<=2) 
        
        % ----- algorithm starts here, split into sub-lattices ---------------
        
        mwin=comp_nonsepwin2multi(g,a,M,lt,L);
        
        % phase factor correction (backwards), for more information see 
        % analysis routine
        
        E = exp(2*pi*i*a*kron(0:N/lt(2)-1,ones(1,lt(2))).*...
                rem(kron(ones(1,N/lt(2)), 0:lt(2)-1)*lt(1),lt(2))/M);

        coef=bsxfun(@times,coef,E);
        
        % simple algorithm: split into sublattices and add the result from each
        % sublattice.
        f=zeros(L,W);
        for ii=0:lt(2)-1
            % Extract sublattice
            sub=coef(:,ii+1:lt(2):end,:);
            f=f+comp_idgt(sub,mwin(:,ii+1),lt(2)*a,[0 1],0,0);  
        end;

    else

        g=fir2long(g,L);
      
        [s0,s1,br] = shearfind(L,a,M,lt);
        
        b=L/M;
        ar = a*b/br;
        Mr = L/br;
        Nr = L/ar;
        
        ind = [ar 0; 0 br]*[kron((0:L/ar-1),ones(1,L/br));kron(ones(1,L/ar), ...
                                                          (0:L/br-1))];
        phs = reshape(mod((s1*(ind(1,:)-s0*ind(2,:)).^2+s0*ind(2,:).^2)*(L+1) ...
                          -2*(s0 ~= 0)*ind(1,:).*ind(2,:),2*L),L/br,L/ar);    
        phs = exp(-pi*1i*phs/L);

        ind_final = [1 0;-s1 1]*[1 -s0;0 1]*ind;
        ind_final = mod(ind_final,L);
        
        if s1 ~= 0
            g = comp_pchirp(L,s1).*g;
        end
        
        if s0 ~= 0
            
            c_rect = zeros(Nr,Mr,W);
            g = comp_pchirp(L,-s0).*fft(g);
            for w=0:W-1
                c_rect(ind(1,[1:Mr,end:-1:Mr+1])/ar+1+(ind(2,:)/br)*Nr+w*M*N) = ...
                    coef(floor(ind_final(2,:)/b)+1+(ind_final(1,:)/a)*M+w*M* ...
                         N).*phs(ind(2,:)/br+1+(ind(1,:)/ar)*Mr);
            end;
            f = comp_idgt(c_rect,g,br,[0 1],0,0);
            f = ifft(bsxfun(@times,comp_pchirp(L,s0),f));   
            
        else
            
            c_rect = zeros(Mr,Nr,W);
            for w=0:W-1
                c_rect(ind(2,:)/br+1+(ind(1,:)/ar)*Mr+w*M*N) = ... 
                    coef(floor(ind_final(2,:)/b)+1+(ind_final(1,:)/a)*M+w*M*N);       
                c_rect(:,:,w+1) = phs.*c_rect(:,:,w+1);
            end;
            f = comp_idgt(c_rect,g,ar,[0 1],0,0);
            
        end
        
        if s1 ~= 0
            f = bsxfun(@times,comp_pchirp(L,-s1),f);
        end            
        
    end;

end;    