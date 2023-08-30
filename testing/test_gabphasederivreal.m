function test_failed = test_gabphasederivreal()
test_failed = 0;

% Tests whether all the algorithms return the same values
disp('-----------------test_gabphasederivreal----------------------');

L = 300;
l = 0:L-1; l = l(:);

f1 = zeros(L,1); f6(floor(L/2):floor(L/2)+5) = 1; %Box function
modul = 2*sin(8*pi*l./L);
f2 = sin(2*pi*(20*l./L+modul)); % Frequency-modulated sinusoid
X = mod(mod(mod(mod(l./L*.65,2*L).*l./L,2*L).*l./L,2*L)*(L+1),2*L);
f3 = sin(pi*X); % Real-valued quadratic chirp

% Compare the phase derivatives only for coefficients bigger than
% relminlvl*max(c(:)) and away from the borders..
relminlvldb = 20;
relminlvl = 10^(-relminlvldb/20);

a = [1,1,4,4,  1,  4,  2];
M = [L,L,L,L,L/4,L/2,L/2];


%g = {{'gauss',1},'gauss',{'hann',8}};
%g = {{'hann',8}};
g = {{'gauss',1},{'gauss',4},...
    'gauss',{'gauss',4},{'gauss',1},'gauss'};%,{'hann',80}};

f = {f3,f3,f3,f3,f3,f3,f3,f3,f3};

phaseconvCell = {'relative','freqinv','timeinv','symphase'};

dflags = {'t','f','tt','ff','tf'};

for ii =1:numel(g)
    for phaseconvId=1:numel(phaseconvCell)
        phaseconv=phaseconvCell{phaseconvId};
        for dflagId = 1:numel(dflags)
            dflag = dflags{dflagId};
            
            c = dgtreal(f{ii},g{ii},a(ii),M(ii));
            [~,info]=gabwin(g{ii},a(ii),M(ii),dgtlength(numel(f{ii}),a(ii),M(ii)));
            
            minlvl = max(abs(c(:)))*relminlvl;
            algArgs = {
                {'dgt',f{ii},g{ii},a(ii),M(ii)}     
                {'phase',angle(c),a(ii),M(ii)}
                {'cross',f{ii},g{ii},a(ii),M(ii)} 
                {'abs',abs(c),g{ii},a(ii),M(ii)}   
                };

            algRange = 1:numel(algArgs);
            if ~info.gauss
                algRange = 1:numel(algArgs)-1;
            end
            pderivs = cell(numel(algRange),1);
            
            for algId=algRange
                alg = algArgs{algId};
                pderivs{algId}=gabphasederivreal(dflag,alg{:},phaseconv);
                
                % Make it independent of L
                if any(strcmp(dflag,{'t','f'}))
                    pderivs{algId} = pderivs{algId}./L;
                end
                % Mask by big enough coefficients
                pderivs{algId}(abs(c)<minlvl)=0;
            end
            % MSE
            
            N = L/a(ii);
            nfirst = ceil(N/3); nlast = ceil(2*N/3);
            pderivs = cellfun(@(pEl) pEl(:,nfirst:nlast),pderivs,'UniformOutput',0);
            
            % MSE
            res = (cellfun(@(pEl) norm(pEl(:)-pderivs{1}(:)).^2/numel(pEl(:)),pderivs(2:end)));
            [test_failed,failstring]=ltfatdiditfail( sum(res)/numel(algRange),...
                           test_failed,1e-3);
                       
             algStr = cellfun(@(cEl) cEl{1},algArgs(algRange),'UniformOutput',0);
                    

             fail2string = cellfun(@(aEl,rEl)sprintf('%s=%d, ',aEl,rEl),algStr(2:end),num2cell(res),'UniformOutput',0);
             fail2string = strcat(fail2string{:});
             fprintf('GABPHASEDERIVREAL a=%d,M=%d, dflag:%2s pc:%8s MSE %d, %s,  %s\n',a(ii),M(ii),dflag,phaseconv,sum(res),fail2string,failstring);
             if ~isempty(failstring)
                 clim=2;figure(1);plotdgtreal(pderivs{1},1,M(ii),'lin','clim',[-clim,clim]);
                 figure(2);plotdgtreal(pderivs{2},1,M(ii),'lin','clim',[-clim,clim]);
                 figure(3);plotdgtreal(pderivs{3},1,M(ii),'lin','clim',[-clim,clim]);
                 figure(4);plotdgtreal(pderivs{1}-pderivs{3},1,M(ii),'lin','clim',[-clim,clim]);
                 figure(5);plotdgtreal(pderivs{1}-pderivs{2},1,M(ii),'lin','clim',[-clim,clim]);
                 prd = 2;
             end
        end
    end
end
