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
                    
                    L=dgtlength(1,a,M,[lt1,lt2]);
                    lt=[lt1,lt2];
                                       
                    b=L/M;
                    s=b*lt(1)/lt(2);

                    [s0,s1,X] = shearfind(a,b,s,L);
                              
                    f=crand(L,1);                                        
                    g=crand(L,1);
                    
                    gd       = nonsepgabdual(g,a,M,lt);
                    gd_shear = nonsepgabdual(g,a,M,lt,'shear');
                    
                    res=norm(gd-gd_shear)/norm(g);
                    [test_failed,fail]=ltfatdiditfail(res,test_failed);
                    stext=sprintf(['DUAL SHEAR L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                                   '%s'], L,a,M,lt(1),lt(2),res,fail);
                    disp(stext)

                    
                    if numel(fail)>0
                        error('Failed test');
                    end;

                    
                    cc = nonsepdgt(f,g,a,M,lt);
                    
                    cc_shear = nonsepdgt(f,g,a,M,lt,'shear');
                    
                    res = norm(cc(:)-cc_shear(:))/norm(cc(:));
                    [test_failed,fail]=ltfatdiditfail(res,test_failed);
                    stext=sprintf(['DGT  SHEAR L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
                                   '%s'], L,a,M,lt(1),lt(2),res,fail);
                    disp(stext)

                    if numel(fail)>0
                        error('Failed test');
                    end;


                    
                    
                    r=inonsepdgt(cc_shear,gd,a,lt,'shear');  
                    res=norm(f-r,'fro');
                    
                    [test_failed,fail]=ltfatdiditfail(res,test_failed);
                    stext=sprintf(['REC  SHEAR L:%3i a:%3i M:%3i lt1:%2i lt2:%2i %0.5g ' ...
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
