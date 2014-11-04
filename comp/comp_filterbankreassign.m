function sr = comp_filterbankreassign(s,tgrad,fgrad,a,cfreq)

N = cellfun(@numel,s);
M = numel(s);



oneover2 = 1/2;

sr = cell(M,1);
cfreq2 = zeros(M,1);

for mm = 1:M
    sr{mm} = zeros(N(mm),1);
    cfreq2(mm) = cfreq(mm) - floor(cfreq(mm)*oneover2)*2;
end

%% Compute reassigned frequencies and times
for mm = M:-1:1

    fgradIdx = zeros(N(mm),1);
    tgradIdx = zeros(N(mm),1);
    
    cfreqm = cfreq2(mm);
 
    for jj = 1:N(mm)
        
        fgradmjj = fgrad{mm}(jj) + cfreqm;
        
        oldfgrad = 10;
        
        
        fgradIdx(jj) = 0;
        
        if fgrad{mm}(jj) > 0
            pos = mm;
            for ii = mm:M
                pos = ii;
                tmpfgrad = cfreq2(ii) - fgradmjj;
                
                if tmpfgrad >= 0
                    if abs(tmpfgrad) < abs(oldfgrad)
                        
                        fgradIdx(jj) = pos;
                        
                    else
                        
                        fgradIdx(jj) = pos-1;
                       
                    end
                    
                    break;
                   
                end
                
                oldfgrad = tmpfgrad;
             
            end
            
            if pos == M && tmpfgrad < 0
                    
                for ii = 1:mm
                    
                   pos = ii;
                    
                    tmpfgrad = cfreq2(ii) - fgradmjj + 2;
                    
                    if tmpfgrad >= 0
                        if abs(tmpfgrad) < abs(oldfgrad)
                        
                            fgradIdx(jj) = pos;
                            
                        else
                        
                            fgradIdx(jj) = pos-1;
                       
                        end
                    
                        break;
                    
                    end
                
                    oldfgrad = tmpfgrad;
                   
                end
            end
            
            if fgradIdx(jj) < 1
                
                fgradIdx(jj) = M;
                
            end
        else
            pos = mm;
            for ii = mm:-1:1
                pos = ii;                
                tmpfgrad = cfreq2(ii) - fgradmjj;
                
                if tmpfgrad <= 0
                    if abs(tmpfgrad) < abs(oldfgrad)
                        
                        fgradIdx(jj) = pos;
                        
                    else
                        
                        fgradIdx(jj) = pos+1;
                       
                    end
                    
                    break;
                   
                end
                
             oldfgrad = tmpfgrad;
             
            end
            
            if pos == 1 && tmpfgrad > 0
                    
                for ii = M:-1:mm
                    
                    pos = ii;
                    
                    tmpfgrad = cfreq2(ii) - fgradmjj - 2;
                    
                    if tmpfgrad <= 0
                        if abs(tmpfgrad) < abs(oldfgrad)
                        
                            fgradIdx(jj) = pos;
                            
                        else
                        
                            fgradIdx(jj) = pos+1;
                       
                        end
                    
                        break;
                    
                    end
                
                    oldfgrad = tmpfgrad;
                   
                end
                
            end
                
            if fgradIdx(jj) >= M+1
                
                fgradIdx(jj) = 1;
             
            end
        end
    end
            
    for jj = 1:N(mm)
       tmpIdx = fgradIdx(jj);
       tgradIdx(jj) = mod(round((tgrad{mm}(jj) + a(mm)*(jj-1))./a(tmpIdx)),N(tmpIdx))+1;
    end
           
    for jj=1:N(mm)
       sr{fgradIdx(jj)}(tgradIdx(jj)) = sr{fgradIdx(jj)}(tgradIdx(jj))+s{mm}(jj);
    end
  
end