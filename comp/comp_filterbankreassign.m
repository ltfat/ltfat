function [sr,repos] = comp_filterbankreassign(s,tgrad,fgrad,a,cfreq)

Lc = cellfun(@numel,s);
M = numel(s);

oneover2 = 1/2;

sr = cell(M,1);
cfreq2 = zeros(M,1);

for mm = 1:M
    sr{mm} = zeros(Lc(mm),1);
    cfreq2(mm) = cfreq(mm) - floor(cfreq(mm)*oneover2)*2;
end

if nargout>1
   chan_pos = [0;cumsum(Lc)];
   repos = cell(sum(cellfun(@numel,sr)),1);
end

%% Compute reassigned frequencies and times
for mm = M:-1:1

    tgradIdx = zeros(Lc(mm),1);
    fgradIdx = zeros(Lc(mm),1);

    cfreqm = cfreq2(mm);

    for jj = 1:Lc(mm)

        tgradmjj = tgrad{mm}(jj) + cfreqm;
        oldtgrad = 10;
        tgradIdx(jj) = 0;

        if tgrad{mm}(jj) > 0
            pos = mm;
            for ii = mm:M
                pos = ii;
                tmptgrad = cfreq2(ii) - tgradmjj;

                if tmptgrad >= 0
                    if abs(tmptgrad) < abs(oldtgrad)

                        tgradIdx(jj) = pos;

                    else

                        tgradIdx(jj) = pos-1;

                    end

                    break;

                end

                oldtgrad = tmptgrad;

            end

            if pos == M && tmptgrad < 0

                for ii = 1:mm

                   pos = ii;

                   tmptgrad = cfreq2(ii) - tgradmjj + 2;

                   if tmptgrad >= 0
                        if abs(tmptgrad) < abs(oldtgrad)

                            tgradIdx(jj) = pos;

                        else

                            tgradIdx(jj) = pos-1;

                        end

                        break;

                    end

                    oldtgrad = tmptgrad;

                end
            end

            if tgradIdx(jj) < 1

                tgradIdx(jj) = M;

            end
        else
            pos = mm;
            for ii = mm:-1:1
                pos = ii;
                tmptgrad = cfreq2(ii) - tgradmjj;

                if tmptgrad <= 0
                    if abs(tmptgrad) < abs(oldtgrad)

                        tgradIdx(jj) = pos;

                    else

                        tgradIdx(jj) = pos+1;

                    end

                    break;

                end

                oldtgrad = tmptgrad;

            end

            if pos == 1 && tmptgrad > 0

                for ii = M:-1:mm

                    pos = ii;

                    tmptgrad = cfreq2(ii) - tgradmjj - 2;

                    if tmptgrad <= 0
                        if abs(tmptgrad) < abs(oldtgrad)

                            tgradIdx(jj) = pos;

                        else

                            tgradIdx(jj) = pos+1;

                        end

                        break;

                    end

                    oldtgrad = tmptgrad;

                end

            end

            if tgradIdx(jj) >= M+1

                tgradIdx(jj) = 1;

            end
        end
    end

    for jj = 1:Lc(mm)
       tmpIdx = tgradIdx(jj);
       fgradIdx(jj) = mod(round((fgrad{mm}(jj) + a(mm)*(jj-1))./a(tmpIdx)),Lc(tmpIdx))+1;
    end

    for jj=1:Lc(mm)
       sr{tgradIdx(jj)}(fgradIdx(jj)) = sr{tgradIdx(jj)}(fgradIdx(jj))+s{mm}(jj);
    end

    if nargout>1
       for jj=1:Lc(mm)
          repos{chan_pos(tgradIdx(jj))+fgradIdx(jj)}(end+1,1) = chan_pos(mm)+jj;
       end
    end



end
