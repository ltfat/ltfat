function test_failed = tets_filterbankconstphase()
firwinflags=getfield(arg_firwin,'flags','wintype');
freqwinflags=getfield(arg_freqwin,'flags','wintype');

[f,fs] = gspi;
win = {'gauss','blackman','hann'};

for winId = 1:numel(win)
    bwmul = 1/4;
    [g,a,fc,L,info]=audfilters(fs,numel(f),'fractional','bwmul',bwmul,'spacing',1/9,'redtar',8,win{winId},'subprec');
    aud_red = sum(a(:,2))/L;

    corig = filterbank(f,g,a);

    c=filterbankconstphase(corig,a,info.fc,info.tfr);

    cproj = filterbank(ifilterbankiter(c,g,a,'pcg','tol',1e-10),g,a);
    Cdb = 20*log10( norm(abs(cell2mat(corig)) - abs(cell2mat(cproj)) )/norm( abs(cell2mat(corig))) );

    %figure(1);plotfilterbankphasediff(corig,c,1e-4,a);

    fprintf('AUDFILTERS win=%8s red=%.2f, C=%.2f dB\n', win{winId}, aud_red, Cdb);

    redmul=2;
    switch win{winId}
    case firwinflags
        redmul=1;
    end

    [g,a,fc,L,info]=cqtfilters(fs,100,fs/2 - 100,48,numel(f),'fractional',win{winId},'redmul',redmul,'Qvar',4,'subprec');
    cqt_red = sum(a(:,2))/L;

    corig = filterbank(f,g,a);
    c=filterbankconstphase(corig,a,info.fc,info.tfr);
   % figure(2);plotfilterbankphasediff(corig,c,1e-4,a);

    cproj = filterbank(ifilterbankiter(c,g,a,'pcg','tol',1e-10),g,a);
    Cdb = 20*log10( norm(abs(cell2mat(corig)) - abs(cell2mat(cproj)) )/norm( abs(cell2mat(corig))) );
    fprintf('CQTFILTERS win=%8s red=%.2f, C=%.2f dB\n', win{winId}, cqt_red, Cdb);
end
