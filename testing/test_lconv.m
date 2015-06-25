function test_failed=test_lconv
Lf=[9, 10  9, 9, 1];
Wf=[1,  1, 3, 1, 1];

Lg=[9, 10  9, 9, 1];
Wg=[1,  1, 1, 3, 1];


ctypes={'default','r','rr'};

test_failed=0;

disp(' ===============  TEST_LCONV ==============');

for jj=1:length(Lf)

    for ii=1:3
        for type = {'real','complex'}
            ctype=ctypes{ii};
            
            
            if strcmp(type{1},'complex')
                f=tester_crand(Lf(jj), Wf(jj));
                g=tester_crand(Lg(jj), Wg(jj));
            else
                f=tester_rand(Lf(jj), Wf(jj));
                g=tester_rand(Lg(jj), Wg(jj));
            end
            
            h1=lconv(f,g,ctype);
            h2cell = {};

            if Wf(jj) == 1
                for wId = 1:Wg(jj)
                    h2cell{end+1}=ref_lconv(f,g(:,wId),ctype);
                end
            else
                for wId = 1:Wf(jj)
                    h2cell{end+1}=ref_lconv(f(:,wId),g,ctype);
                end
            end
             

            h2 = cell2mat(h2cell);
            
            res=norm(h1(:)-h2(:));
            [test_failed,fail]=ltfatdiditfail(res,test_failed);
            s=sprintf('PCONV %3s %6s  %0.5g %s',ctype,type{1},res,fail);
            disp(s);
        end
    end
end

