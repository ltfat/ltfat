function H=comp_transferfunction(g,L)
%COMP_TRANSFERFUNCTION  Compute the transfer function

l=(0:L-1).'/L;
if isfield(g,'h')
        
    g_time=circshift(postpad(g.h,L),g.offset).*...
           exp(2*pi*1i*round(g.centre*L/2)*l);
        
    if g.realonly
        g_time=real(g_time);
    end;
        
    H=fft(g_time);
    
else
    if ~isnumeric(g.H)
        g.H=g.H(L);
    end;

    H=circshift(postpad(g.H,L),g.foff(L)).*exp(-2*pi*1i*round(g.delay)*l);
    
    if g.realonly
        H=(H+involute(H))/2;
    end;
        
end;
