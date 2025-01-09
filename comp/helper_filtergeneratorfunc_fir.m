function [filterfunc,winbw] = helper_filtergeneratorfunc_fir(wintype,winCell,fs,bwmul,min_win,trunc_at,audscale,do_subprec,do_symmetric,do_warped)
firwinflags=getfield(arg_firwin,'flags','wintype');
freqwinflags=getfield(arg_freqwin,'flags','wintype');
probeLs = fs;
probeLg = 5000;

subprecflag = 'pedantic';
if ~do_subprec, subprecflag = 'nopedantic'; end

switch wintype
    case firwinflags
        % probe prototype window full length
        g_probe = fir2long(firwin(wintype,probeLg),probeLs);
        % peak normalize
        gf_probe = fft(g_probe)/max(abs(fft(g_probe)));
        % compute ERB-type bandwidth of the prototype
        winbw = norm(gf_probe).^2*probeLg/probeLs/4; % in normalized frequency
        %winbw = norm(gf_probe).^2/(probebw*probeLg/2);

        if do_symmetric
            filterfunc = @(fsupp,fc)... 
                         firfilter(winCell{1},fsupp,fc/fs*2,'1');
        else
            fsupp_scale=1/winbw*kv.bwmul;
            filterfunc = @(fsupp,fc,scal)...
                         warpedblfilter(winCell,fsupp_scale,fc,fs,...
                                        @(freq) freqtoaud(freq,flags.audscale),@(aud) audtofreq(aud,flags.audscale),'scal',scal,'inf');
        end     
end


function width = winwidthatheight(gnum,atheight)

width = zeros(size(atheight));
for ii=1:numel(atheight)
    gl = numel(gnum);
    gmax = max(gnum);
    frac=  1/atheight(ii);
    fracofmax = gmax/frac;

    ind =find(gnum(1:floor(gl/2)+1)==fracofmax,1,'first');
    if isempty(ind)
        %There is no sample exactly half of the height
        ind1 = find(gnum(1:floor(gl/2)+1)>fracofmax,1,'last');
        ind2 = find(gnum(1:floor(gl/2)+1)<fracofmax,1,'first');
        rest = 1-(fracofmax-gnum(ind2))/(gnum(ind1)-gnum(ind2));
        width(ii) = 2*(ind1+rest-1);
    else
        width(ii) = 2*(ind-1);
    end
end
