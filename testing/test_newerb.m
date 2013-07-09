warpname={'symmetric','warped'};

test_failed=0;

for warpidx=1:2
    warping=warpname{warpidx};
    
    for testcase=1:4
        
        isuniform=0;
        switch testcase
          case 1
            % Fails for L=5699
            L=5700
            [g,a]=erbfilters(16000,L,'fractional',warping);
            L=filterbanklength(L,a);
            isreal=1;
          case 2
            L=4000;
            [g,a]=erbfilters(16000,warping);
            L=filterbanklength(L,a);
            isreal=1;
          case 3
            L=6000;
            [g,a]=erbfilters(16000,'uniform',warping);
            L=filterbanklength(L,a);
            isreal=1;
            isuniform=1;
          case 4
            % Fails for L=5699
            L=5700
            [g,a]=erbfilters(10000,L,'fractional',warping);
            L=filterbanklength(L,a);
            isreal=1;
            
        end;
        
        disp('testcase')
        testcase
        
        % Test it
        if 0
            f=greasy+randn(L,1);
        else        
            f=[1;zeros(L-1,1)];
            if 0
                ff=fft(f);
                ff(1999:2003)=0;
                f=ifft(ff);
            end;
        end;
        
        % Inspect it: Dual windows, frame bounds and the response
        disp('Frame bounds:')
        [A,B]=filterbankrealbounds(g,a,L);
        A
        B
        B/A
        filterbankresponse(g,a,L,'real','plot');
        gd=filterbankrealdual(g,a,L);
        
        if isuniform
            c=ufilterbank(f,g,a);
        else
            c=filterbank(f,g,a);
        end;
        r=ifilterbank(c,gd,a);
        if isreal
            r=2*real(r);
        end;
        
        disp('Reconstruction:')
        res=norm(f-r);
    
        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        s=sprintf(['ERBFILTER DUAL  %s testcase: %3i L:%3i %0.5g %s'],warping,testcase,L,res,fail);    
        disp(s);
        
        
        gt=filterbankrealtight(g,a,L);
        if isuniform
            ct=ufilterbank(f,gt,a);
        else
            ct=filterbank(f,gt,a);
        end;
        rt=ifilterbank(ct,gt,a);
        if isreal
            rt=2*real(rt);
        end;

        res=norm(f-r)

        [test_failed,fail]=ltfatdiditfail(res,test_failed);
        s=sprintf(['ERBFILTER TIGHT %s testcase: %3i L:%3i %0.5g %s'],warping,testcase,L,res,fail);    
        disp(s);
        
    end;
    
end;