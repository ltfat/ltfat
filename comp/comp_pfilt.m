function h=comp_pfilt(f,g,a,do_time);
%COMP_PFILT  Compute filtering
%
%   If *do_time* is 1, the routine will use the time-side algorithm for
%   FIR filters, otherwise it will always do the multiplication in the
%   frequency domain.

[L,W]=size(f);

if numel(a)==1
    N=L/a;
else
    N=L/a(1)*a(2);    
end;


h=zeros(N,W,assert_classname(f));

l=(0:L-1).'/L;

if isfield(g,'h') && do_time
    % Use a direct algorithm
    g_time=circshift(postpad(g.h,L),g.offset).*...
           exp(2*pi*1i*round(g.centre*L/2)*l);            
    if g.realonly
        g_time=real(g_time);
    end;
    g_time=conj(involute(g_time));
    for n=0:N-1
        h(n+1,:)=sum(bsxfun(@times,f,circshift(g_time,a*n)));
    end;
    
else
    
        
    % Zero-extend and use a full length fft algorithm
    % This case can be further optimized

    G=comp_transferfunction(g,L);

    if numel(a)==1
        for w=1:W
            F=fft(f(:,w));
            h(:,w)=ifft(sum(reshape(F.*G,N,a),2))/a;
        end;
    else
        Llarge=ceil(L/N)*N;
        amod=Llarge/N;
        
        for w=1:W
            
            F=fft(f(:,w));
            h(:,w)=ifft(sum(reshape(middlepad(F.*G,Llarge),N,amod),2))/amod;
            
        end;
    end;
    
end;



