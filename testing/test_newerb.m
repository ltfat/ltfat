warpname={'warped','symmetric'};

test_failed=0;

for warpidx=1:2
    warping=warpname{warpidx};
    
    for fracidx=1:2
        if fracidx==1
            fractional={'regsampling'};
            fracname='regsamp';
        else
            fractional={'fractional','L',10000};
            fracname='fractional';
        end;
        
        for uniformidx=1:2
            if uniformidx==1                
                isuniform=0;
                uniform='nonuniform';
            else
                isuniform=1;
                uniform='uniform';
            end;
            
            [g,a]=erbfilters(16000,fractional{:},warping,uniform,'redmul',1);
            L=filterbanklength(5000,a);
            isreal=1;
            
            
            f=randn(L,1);
            % Test it
            if 0
                f=[1;zeros(L-1,1)];
                ff=fft(f);
                ff(1999:2003)=0;
                f=ifft(ff);
            end;
            
            if 0
                % Inspect it: Dual windows, frame bounds and the response
                disp('Frame bounds:')
                [A,B]=filterbankrealbounds(g,a,L);
                A
                B
                B/A
                filterbankresponse(g,a,L,'real','plot');
            end;
            
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
            
            res=norm(f-r);
            
            [test_failed,fail]=ltfatdiditfail(res,test_failed);
            s=sprintf(['ERBFILTER DUAL  %s %s %s L:%3i %0.5g %s'],warping,fracname,uniform,L,res,fail);    
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
            
            res=norm(f-r);
            
            [test_failed,fail]=ltfatdiditfail(res,test_failed);
            s=sprintf(['ERBFILTER TIGHT %s %s %s L:%3i %0.5g %s'],warping,fracname,uniform,L,res,fail);    
            disp(s);
            
        end;
        
    end;
    
end;