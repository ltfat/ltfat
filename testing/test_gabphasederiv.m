function test_failed = test_gabphasederiv()
test_failed = 0;

% Tests whether all the algorithms return the same values
disp('-----------------test_gabphasederiv----------------------');

L = 100;
l = 0:L-1; l = l(:);

f1 = exp(1i*pi*0.1*l);
f2 = pchirp(L,1);
f3 = expchirp(L,0.2,2.2,'fc',-1.2);
f4 = exp(1i*pi*0.5*L*(l./(L)).^2);
f5 = exp(1i*pi*0.1*L*(l./(L)).^2);
f6 = zeros(L,1); f6(floor(L/2):floor(L/2)+5) = 1;

f = {f2,f1,f2};

a = [4,4,4];
M = [16,8,8];
%g = {{'gauss',1},'gauss',{'hann',8}};
g = {{'gauss',1}};

phaseconvCell = {'freqinv','timeinv','symphase','relative'};

dflags = {'t','f','tt','ff','tf','ft'};

for ii =1:numel(g)
    for phaseconvId=1:numel(phaseconvCell)
        phaseconv=phaseconvCell{phaseconvId};
        for dflagId = 1:numel(dflags)
            dflag = dflags{dflagId};
            
            c = dgt(f{ii},g{ii},a(ii),M(ii));
            [~,info]=gabwin(g{ii},a(ii),M(ii),dgtlength(numel(f{ii}),a(ii),M(ii)));
            
            minlvl = max(abs(c(:)))*1e-4;
            algArgs = {
                {'dgt',f{ii},g{ii},a(ii),M(ii)}
                {'phase',angle(c),a(ii)}
                {'cross',f{ii},g{ii},a(ii),M(ii)} 
                {'abs',abs(c),g{ii},a(ii)}   
                };

            algRange = 1:numel(algArgs);
            if ~info.gauss
                algRange = 1:numel(algArgs)-1;
            end
            pderivs = cell(numel(algRange),1);
            
            for algId=algRange
                alg = algArgs{algId};
                pderivs{algId}=gabphasederiv(dflag,alg{:},phaseconv);
                % Mask by big enough coefficients
                pderivs{algId}(abs(c)<minlvl)=0;
            end
            % MSE
            res = (cellfun(@(pEl) norm(pEl(:)-pderivs{1}(:)).^2/numel(pEl(:)),pderivs(2:end)));
            [test_failed,failstring]=ltfatdiditfail( sum(res),...
                           test_failed,1e-8);
             fprintf('GABPHASEDERIV dflag:%2s phaseconv:%8s %d %s\n',dflag,phaseconv,sum(res),failstring);
        end
    end
end
