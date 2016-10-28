function gout=freqfilter(winname,bw,varargin)
% 
% Authors: Nicki Holighaus & Zdenek Prusa
% Date: September 15

if ~iscell(winname), winname = {winname}; end

% Define initial value for flags and key/value pairs.
definput.import={'normalize'};
definput.importdefaults={'energy'};
definput.keyvals.delay=0;
definput.keyvals.fc=0;
definput.keyvals.fs=2;
%definput.keyvals.order=4;
definput.keyvals.scal=1;
definput.keyvals.min_win=1;
%definput.keyvals.trunc_at=10^(-5);
definput.keyvals.bwtruncmul = 4;
definput.flags.pedantic = {'pedantic','nopedantic'};
definput.flags.real={'complex','real'};

[flags,kv]=ltfatarghelper({'fc'},definput,varargin);

[bw,kv.fc,kv.delay,kv.scal]=scalardistribute(bw,kv.fc,kv.delay,kv.scal);

% Sanitize
kv.fc=modcent(2*kv.fc/kv.fs,2);

Lw = @(L,bw) min(ceil(bw*kv.bwtruncmul*L/kv.fs),L);
    
fsRestricted = @(L,bw) kv.fs/L*Lw(L,bw);
if flags.pedantic
    fc_offset = @(L,fc) L/2*fc-round(L/2*fc);
else
    fc_offset = @(L,fc) 0;
end

Nfilt = numel(bw);
gout = cell(Nfilt,1);
for ii=1:Nfilt
    g=struct();
    
    if flags.do_1 || flags.do_area 
        g.H=@(L)    fftshift(freqwin(winname,Lw(L,bw(ii)),bw(ii),fsRestricted(L,bw(ii)),'shift',fc_offset(L,kv.fc(ii))))*kv.scal(ii)*L;                
    end;
    
    if  flags.do_2 || flags.do_energy
        g.H=@(L)    fftshift(freqwin(winname,Lw(L,bw(ii)),bw(ii),fsRestricted(L,bw(ii)),'shift',fc_offset(L,kv.fc(ii))))*kv.scal(ii)*sqrt(L);                        
    end;
        
    if flags.do_inf || flags.do_peak
        g.H=@(L)    fftshift(freqwin(winname,Lw(L,bw(ii)),bw(ii),fsRestricted(L,bw(ii)),'shift',fc_offset(L,kv.fc(ii))))*kv.scal(ii);           
    end;
    
    g.foff=@(L) round(L/2*kv.fc(ii)) - floor(Lw(L,bw(ii))/2);                    
                    
    g.realonly=flags.do_real;
    g.delay=kv.delay(ii);
    g.fs=kv.fs;
    gout{ii}=g;
end;

if Nfilt==1
    gout=g;
end;