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
