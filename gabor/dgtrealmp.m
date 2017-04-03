function [ c, relres, fhat, cres] = dgtrealmp( f, a, M, varargin)
%DGTREALMP Matching Pursuit using 
%   Detailed explanation goes here

definput.keyvals.tol=1e-4;
definput.keyvals.maxit=[];
definput.keyvals.relresstep=[];
definput.keyvals.printstep=[];
definput.flags.print = {'quiet','printerr','printdb'};
definput.keyvals.g1 = [];
definput.keyvals.g2 = [];
definput.keyvals.g3 = [];
definput.keyvals.g4 = [];
definput.keyvals.g5 = [];
definput.keyvals.g6 = [];
definput.keyvals.g7 = [];
definput.keyvals.g8 = [];
definput.keyvals.g9 = [];
definput.keyvals.g10 = [];
definput.flags.method={'mp','comp','locomp'};
[flags,kv]=ltfatarghelper({'maxit','tol'},definput,varargin);

Ls = size(f,1);
L = dgtlength(Ls,a,M);

if isempty(kv.maxit), kv.maxit = floor(0.1*L); end

if isempty(kv.relresstep) || kv.relresstep > kv.maxit
    kv.relresstep = kv.maxit; 
end

maxwins = 10;
g = cell(maxwins,1);
for ii=1:10
    g{ii} = kv.(sprintf('g%i',ii));
end
g = g(~cellfun(@isempty,g));

if isempty(g), g{1} = 'gauss'; end

g = cellfun(@(gEl) normalize(gabwin(gEl,a,M,L),'2'),g(:)',...
             'UniformOutput',0);

N = L/a;
superN = numel(g)*N;
M2 = floor(M/2) + 1;

kerns = cell(numel(g),numel(g));
masks = cell(numel(g),numel(g));

if flags.do_mp || flags.do_locomp

for ii = 1:numel(g)
    for jj = 1:numel(g)
        [kerns{ii,jj}, masks{ii,jj}] = projkernels(g{ii},g{jj},a,M,L,kv.tol);
    end
end

cresfull = cell2mat(cellfun(@(gEl) dgt(f,gEl,a,M),g(:)','UniformOutput',0));

elseif flags.do_comp
    gd = gabdual(cell2mat(g),a,M,L);
    gd = num2cell(gd,1);
    atnorms = zeros(numel(g));
    
    for ii = 1:numel(g)
        for jj = 1:numel(g)
            atnorms(ii,jj) = g{ii}.'*gd{ii};
            
            kerns{ii,jj} = projkernels(g{ii},gd{jj},a,M,L,kv.tol);
            for kk = 1:size(kerns{ii,jj},3)
                kerns{ii,jj}(:,:,kk) = kerns{ii,jj}(:,:,kk)./atnorms(ii,jj);
            end
        end  
    end
    
    cresfull = cell2mat(cellfun(@(gEl,nEl) 1/nEl*dgt(f,gEl,a,M),gd(:)',num2cell(diag(atnorms))','UniformOutput',0));
end

c = zeros(M,superN);

kernwidx = cellfun(@(kEl) fftshift(fftindex(size(kEl,2))),kerns,'UniformOutput',0);
kernhidx = cellfun(@(kEl) fftshift(fftindex(size(kEl,1))),kerns,'UniformOutput',0);
bigkernwidx = cellfun(@(kEl) fftshift(fftindex(min([2*size(kEl,2)-1,N]))),kerns,'UniformOutput',0);
bigkernhidx = cellfun(@(kEl) fftshift(fftindex(min([2*size(kEl,1)-1,M]))),kerns,'UniformOutput',0);
kernno = size(kerns{1,1},3);

s = abs(  cresfull );
norm_f = norm(f);
relres = zeros(ceil(kv.maxit/kv.relresstep),1);


tic
iter = 1;
[maxcols,maxcolspos] = max(s(1:M2,:));

if flags.do_mp

    while (iter <= kv.maxit)
       %% 1) Selection
       [~,n] = max(maxcols); n = n-1;
       m = maxcolspos(n+1); m=m-1;  
       winIdx = floor(n/N) + 1;
       nstart = (winIdx - 1)*N;
       nloc = n - nstart;

       cval = cresfull(m+1,n+1);

       %% 2) Residual update
       for secondwinIdx = 1:numel(g)
            idxn = N*(secondwinIdx-1) + mod(nloc + kernwidx{winIdx,secondwinIdx},N)+1;
            idxm = mod(m + kernhidx{winIdx,secondwinIdx},M)+1;
           
            currkern = kerns{winIdx,secondwinIdx}(:,:,1+mod(n,kernno));
    %        currmask = masks{winIdx,secondwinIdx};
            
            cresfulltmp = cresfull(idxm,idxn);
     %       cresfulltmp(currmask) = cresfulltmp(currmask) - cval*currkern(currmask);
            cresfulltmp = cresfulltmp - cval*currkern;
            cresfull(idxm,idxn) = cresfulltmp;

            s(idxm,idxn) = abs(cresfull(idxm,idxn));
            [maxcolupd,maxcolposupd] = max(s(1:M2,idxn));

            maxcols(idxn) = maxcolupd;
            maxcolspos(idxn) = maxcolposupd;
       end

       %% 3) Coefficient update 
       cvalold = c(m+1,n+1);
       % A single atom can be selected more than once
       c(m+1,n+1) = cvalold + cval;

       %% 4) Rest
       if mod(iter,kv.printstep) == 0
            figure(1);plotdgtreal(cresfull(1:M2,:),a,M,'clim',[-120,10]);shg;
            figure(2);plotdgtreal(c(1:M2,:),a,M,'clim',[-120,10]);shg;
       end

       if mod(iter,kv.relresstep) == 0 || iter == kv.maxit
            [fhatsum,fhat] = reccoefs(c(1:M2,:),g,a,M,L,Ls);
            relrescurr = norm(f-fhatsum)/norm_f;
            relres(ceil(iter/kv.relresstep)) = relrescurr;

            if flags.do_printerr || flags.do_printdb
                if flags.do_printdb
                    relrescurr = 20*log10(relrescurr);
                    fprintf('Iteration %d, RMS: %.2f dB\n',iter, relrescurr);
                else
                    fprintf('Iteration %d, RMS: %.2d %s\n',iter, relrescurr);
                end 
            end
       end

       iter = iter+1;
    end
elseif flags.do_locomp
    suppind = false(size(cresfull));
    
    while (iter <= kv.maxit)
       %% 1) Selection
       [~,n] = max(maxcols); n = n-1;
       m = maxcolspos(n+1); m=m-1;  
       winIdx = floor(n/N) + 1;
       nstart = (winIdx - 1)*N;
       nloc = n - nstart;

       suppind(m+1,n+1) = 1;

       %% 2) Residual update
       for secondwinIdx = 1:numel(g)
            % Index range in the output
            idxn = N*(secondwinIdx-1) + mod(nloc + kernwidx{winIdx,secondwinIdx},N)+1;
            idxm = mod(m + kernhidx{winIdx,secondwinIdx},M)+1;
  
            cresfulltmp = cresfull(idxm,idxn);
            % Mask for active atoms
            suppidtmp = suppind(idxm,idxn);
            %suppidtmp(:) = 0;
            %suppidtmp(floor(end/2)+1) = 1;
            
            % Coefficients of active atoms
            cval = cresfulltmp(suppidtmp);
            cvaln = repmat(idxn',numel(idxm),1);
            cvaln = cvaln(suppidtmp);
            cvalm = repmat(idxm,1,numel(idxn));
            cvalm = cvalm(suppidtmp);
            
            G = eye(numel(cval));
            
            for ii = 1:numel(cval)
                currkern = kerns{winIdx,secondwinIdx}(:,:,1+mod(cvaln(ii)-1,kernno));
                [kernh,kernw] = size(currkern);
                mmid = floor(kernh/2) + 1 - cvalm(ii);
                nmid = floor(kernw/2) + 1 - cvaln(ii);
                
                for jj = [1:ii-1,ii+1:numel(cval)]
                    muse = mmid + cvalm(jj);
                    nuse = nmid + cvaln(jj);
                    
                    if muse >=1  && muse <= kernh && nuse>=1 && nuse <=kernw
                        G(ii,jj) = conj(currkern(muse, nuse));
                    end
                end
            end
            

       end
       
       % Update the solution
       cvalinv = G\cval;

       ctmp = c(idxm,idxn);
       ctmp(suppidtmp) = ctmp(suppidtmp) + cvalinv;
       c(idxm,idxn) = ctmp;
          
       for ii=1:numel(cvalinv)
           idxn = mod(cvaln(ii)-1 + kernwidx{winIdx,winIdx},N)+1;
           idxm = mod(cvalm(ii)-1 + kernhidx{winIdx,winIdx},M)+1;   
            
           currkern = kerns{winIdx,secondwinIdx}(:,:,1+mod(cvaln(ii)-1,kernno));
           cresfull(idxm,idxn) = cresfull(idxm,idxn) - cvalinv(ii)*currkern;
       end
        
       idxn = mod(n + bigkernwidx{winIdx,secondwinIdx},N)+1;
       idxm = mod(m + bigkernhidx{winIdx,secondwinIdx},M)+1;
       s(idxm,idxn) = abs(cresfull(idxm,idxn));
       [maxcolupd,maxcolposupd] = max(s(1:M2,idxn));

       maxcols(idxn) = maxcolupd;
       maxcolspos(idxn) = maxcolposupd;


       %% 4) Rest
       if mod(iter,kv.printstep) == 0
            figure(1);plotdgtreal(cresfull(1:M2,:),a,M,'clim',[-120,10]);shg;
            figure(2);plotdgtreal(c(1:M2,:),a,M,'clim',[-120,10]);shg;
       end

       if mod(iter,kv.relresstep) == 0 || iter == kv.maxit
            [fhatsum,fhat] = reccoefs(c(1:M2,:),g,a,M,L,Ls);
            relrescurr = norm(f-fhatsum)/norm_f;
            relres(ceil(iter/kv.relresstep)) = relrescurr;

            if flags.do_printerr || flags.do_printdb
                if flags.do_printdb
                    relrescurr = 20*log10(relrescurr);
                    fprintf('Iteration %d, RMS: %.2f dB\n',iter, relrescurr);
                else
                    fprintf('Iteration %d, RMS: %.2d %s\n',iter, relrescurr);
                end 
            end
       end

       iter = iter+1;
    end

    
    
end

t = toc;
fprintf('%.2f atoms/s\n',kv.maxit/t);

cres = cresfull(1:M2,:);
c = c(1:M2,:);

%% Helper functions
function [fhat,fhats] = reccoefs(c,g,a,M,L,Ls)
N = L/a;

fhats = zeros(Ls,numel(g));

for ii=1:numel(g)
    fhats(:,ii) = idgtreal(c(:,(ii-1)*N+1:ii*N),g{ii},a,M,Ls);
end

fhat = sum(fhats,2);

function [kerns,mask] = projkernels(g1,g2,a,M,L,tol)
N = L/a;
M2 = floor(M/2) + 1;

projfnc = @(c) dgt(middlepad(g1,L),g2,a,M);
ctmp = zeros(M,N); ctmp(1,1) = 1;
kern = projfnc(ctmp);

ksize=findsmallkernelsize(kern(1:M2,:),tol);
if ksize(1)>M, ksize(1) = M; end
kernh= ksize(1); kernw = ksize(2);

kern = middlepad(kern,kernh);
kern = middlepad(kern,kernw,2);

kernno = lcm(M,a)/a;
kerns = zeros([size(kern),kernno]);
for n=1:kernno
    kerns(:,:,n) = fftshift(phasekernfi(kern,n-1,a,M));
end

thr = tol*max(abs(kern(:)));
mask = abs(kerns(:,:,n)) > thr;

function ksize=findsmallkernelsize(kern,relthr)
[M2,N] = size(kern);
thr = relthr*max(abs(kern(:)));

lastrow = 1;
for n=1:N
    newlastrow = find(abs(kern(:,n))>=thr,1,'last');
    if newlastrow > lastrow
        lastrow = newlastrow;
    end
end

lastcol1 = 1;
% Searching from the other end is not neccesary since the kenel
% is always symmetric in the horizontal direction too
% lastcol2 = 1;
for m=1:M2
    newlastcol = find(abs(kern(m,1:ceil(end/2)))>=thr,1,'last');
    if newlastcol > lastcol1
        lastcol1 = newlastcol;
    end

    %     newlastcol = find(abs(kern(m,end:-1:floor(end/2)))>=thr,1,'last');
    %     if newlastcol > lastcol2
    %         lastcol2 = newlastcol;
    %     end
end

ksize = [2*lastrow-1, 2*lastcol1 - 1];

function kernm = phasekernfi(kern,n,a,M)

kernh = size(kern,1);
l = -2*pi*n*fftindex(kernh)*a/M;

kernm = bsxfun(@times,kern,exp(1i*l));

% M/a periodic in m
function kernm = phasekernti(kern,m,a,M)

[kernh, kernw] = size(kern);
l = 2*pi*m*fftindex(kernw)'*a/M;

kernm = bsxfun(@times,kern,exp(1i*l));