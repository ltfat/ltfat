function test_failed=test_fbwarped_framebounds
test_failed = 0;
Ls = 44100;
fs = 44100;
fmax = fs/2;
bins = 1;
fac = [1,7/8,3/4,5/8,1/2];

eigstol = 1e-4;
eigsmaxit = 100;

pcgmaxit = 150;
pcgtol = 1e-4;

warpfun = cell(4,1);
invfun = cell(4,1);

% ERBlet warping
 warpfun{3} = @freqtoerb;
 invfun{3} = @erbtofreq;
% constant-Q warping
 warpfun{4} = @(x) 10*log(x);
 invfun{4} = @(x) exp(x/10);
% sqrt-warping
 warpfun{2} = @(x) sign(x).*((1+abs(x)).^(1/2)-1);
 invfun{2} = @(x) sign(x).*((1+abs(x)).^2-1);
% Linear warping
 warpfun{1} = @(x) x/100;
 invfun{1} = @(x) 100*x;
 
 fmin = [0,0,0,50];

A = zeros(4,length(fac));
B = ones(4,length(fac));
red = A;

for jj = 1:4
    [g,a,fc,L]=warpedfilters(warpfun{jj},invfun{jj},fs,fmin(jj),fmax,bins,Ls,'bwmul',1.5,'fractional','complex');
    
    gf=filterbankresponse(g,a,Ls); framebound_ratio = max(gf)/min(gf);
    disp(['Painless system frame bound ratio: ', num2str(framebound_ratio)]);
    
    for kk = 1:length(fac)
        %[jj,kk]
        atemp = a;
        idx = [2:length(fc)/2,length(fc)/2+2:length(fc)];
        atemp(idx,2) = ceil(atemp(idx,2).*fac(kk));
        red(jj,kk) = sum(atemp(:,2)./atemp(:,1));
        
        [g,asan]=filterbankwin(g,atemp,L,'normal');
        gtemp=comp_filterbank_pre(g,asan,L,10);
        gtemp{1}.H = gtemp{1}.H.*sqrt(fac(kk));
        gtemp{length(fc)/2+1}.H = gtemp{length(fc)/2+1}.H.*sqrt(fac(kk));

        F = frame('filterbank',gtemp,asan,numel(gtemp));
        [A(jj,kk),B(jj,kk)] = framebounds(F,Ls,'tol',eigstol,'pcgtol',pcgtol,'maxit',eigsmaxit,'pcgmaxit',pcgmaxit);

     end
end

disp('This is a ratio B/A. Rows - warping sunction, Cols - redundancy compared to te minimal painless case')
B./A
