function H=comp_transferfunction(g,L)
%COMP_TRANSFERFUNCTION  Compute the transfer function
%
%  `comp_transferfunction(g,L)` computes length *L* transfer function 
%  (frequency response) of a single filter *g*. This function can only
%  handle filters in a proper internal format i.e. already processed by
%  |filterbankwin|.

l=(0:L-1).'/L;
if isfield(g,'h')
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
    else
        % If g.H is already numeric, we must handle transforming it to a
        % different length.
        if isfield(g,'L') 
            if g.L~=L
               g.H = fft(middlepad(ifft(circshift(postpad(g.H(:),g.L),g.foff)),L));
               g.foff = 0;
            end
        else
            % We do not know which L was g.H created with, there is no way
            % how to handle this properly.
            error('%s: No way of knowing which L was used for g.H.',...
                  upper(mfilename));
        end
    end;
    
    H=circshift(postpad(g.H(:),L),g.foff);
    
    if isfield(g,'delay')
       H = H.*exp(-2*pi*1i*round(g.delay)*l);
    end
    
    if isfield(g,'realonly') && g.realonly
        H=(H+involute(H))/2;
    end;
else
    error('%s: Unrecognized filter format. The struct should have either .h or .H field.',upper(mfilename));    
end;
