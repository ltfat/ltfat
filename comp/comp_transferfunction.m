function H=comp_transferfunction(g,L)
%COMP_TRANSFERFUNCTION  Compute the transfer function

l=(0:L-1).'/L;
if isfield(g,'h')
    
    % This is not safe when already having imp. resp. of length L
    % with zero delay (periodically wrapped).
 
    g_time=circshift(postpad(g.h,L),g.offset);

    
    if isfield(g,'fc')
       g_time = g_time.*exp(2*pi*1i*round(g.fc*L/2)*l);
    end
        
    if isfield(g,'realonly') && g.realonly
        g_time=real(g_time);
    end;
        
    H=fft(g_time);
    
elseif isfield(g,'H')
    if ~isnumeric(g.H)
        g.H=g.H(L);
        g.foff=g.foff(L);
    end;
    
    H=circshift(postpad(g.H,L),g.foff);
    
    if isfield(g,'delay')
       H = H.*exp(-2*pi*1i*round(g.delay)*l);
    end
    
    if isfield(g,'realonly') && g.realonly
        H=(H+involute(H))/2;
    end;
else
    error('%s: Unrecognized filter format. The struct should have either .h or .H field.',upper(mfilename));    
end;
