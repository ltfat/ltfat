ar=2:10;
Mr=3:10;
Lmodr=1:3;
lt1r=0:9;
lt2r=1:10;

test_failed=0;

for lt2=lt2r
    for lt1=lt1r
        if lt1>=lt2
            continue;
        end;
        if gcd(lt1,lt2)>1
            continue
        end;
        for M=Mr            
            for a=ar
                if a>=M
                    continue;
                end;                
                for Lmod=Lmodr

                    L=nonsepdgtlengthsignal(1,a,M,[lt1,lt2]);
                    lt=[lt1,lt2];
                    g=crand(L,1);
                    
                    gd       = nonsepgabdual(g,a,M,lt);
                    gd_smith = nonsepgabdual(g,a,M,lt,'smith');
                    gd_shear = nonsepgabdual(g,a,M,lt,'shear');

                    res=norm(gd-gd_smith)/norm(g);
                    [test_failed,fail]=ltfatdiditfail(res,test_failed);
                    stext=sprintf(['DUAL SMITH L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                                   '%s'], L,a,M,lt(1),lt(2),res,fail);
                    disp(stext)

                    if numel(fail)>0
                        error('Failed test');
                    end;
                    
                    res=norm(gd-gd_shear)/norm(g);
                    [test_failed,fail]=ltfatdiditfail(res,test_failed);
                    stext=sprintf(['DUAL SHEAR L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                                   '%s'], L,a,M,lt(1),lt(2),res,fail);
                    disp(stext)

                    
                    if numel(fail)>0
                        error('Failed test');
                    end;

                end;
            end;
        end;
    end;
end;
