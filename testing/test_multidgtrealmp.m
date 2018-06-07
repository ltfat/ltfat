function test_failed=test_multidgtrealmp()
test_failed = 0;

dicts = {...
    {'hann',128,512}
    {'hann',128,512,'tria',64,256}
    {'hann',128,512,'tria',64,256,'tria',256,1024}
    };

errdbArr = [-10,-20,-40];
maxatArr = [100,200,300];
algs = {{'mp'},{'mp','pedanticsearch'},{'selfprojmp'},{'cyclicmp'}};


f = greasy;
for pconv = {'timeinv','freqinv'}
for algId= 1:numel(algs)
for dIdx = 1:numel(dicts)
    for errdbIdx = 1:numel(errdbArr)
        terrdb = errdbArr(errdbIdx);
        [c,frec,info] = multidgtrealmp(f,dicts{dIdx},'errdb',terrdb,algs{algId}{:},pconv{1});

        errdb = 20*log10(norm(f-frec)/norm(f));

        frec2 = sum(info.synthetize(c),2);
        errdb2 = 20*log10(norm(f-frec2)/norm(f));

        errdb3 = 20*log10(info.relres);
        
        [test_failed,fail]=ltfatdiditfail( abs(terrdb - errdb) > 1 ||...
                                           abs(terrdb - errdb2) > 1 || ...
                                           abs(terrdb - errdb3) > 1, test_failed);
        fprintf(['MULTIDGTREALMP ERR  %s  dict:%3i terr=%.2f, itno=%d, %s\n'],strcat(algs{algId}{:}),dIdx,terrdb,info.iter,fail); 

    end
end
end
end

for dIdx = 1:numel(dicts)
    for maxatIdx = 1:numel(maxatArr)
        tmaxat = maxatArr(maxatIdx);
        [c,frec,info] = multidgtrealmp(f,dicts{dIdx},'maxit',tmaxat);
    
        [test_failed,fail]=ltfatdiditfail( info.iter ~= tmaxat || ...
            sum(cellfun(@(cEl) numel(find(abs(cEl))),c)) ~= info.atoms, test_failed);
        fprintf(['MULTIDGTREALMP MAXIT  dict:%3i %s\n'],dIdx,fail); 

    end
end