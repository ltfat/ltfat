%TIME_TCONV  Comparison of implementations for tconv
%   Usage: time_tconv;
%
%   Development test for comparison of the different implemenations of 
%   tconv in implementationstconv.m. See code for details.

% FJ : Here are the duration results on my computer under ubuntu 8.04 for
% n=400 and density=0.001
% Matlab version : 7.0.4.352 (R14) Service Pack 2
% Octave version : 3.0.0
% Algo | Matlab   | Octave
%-----------------------------
% Ref  | 0.39663  | 41.211
% 1    | 0.81334  | 15.575
% 2    | 0.041054 | 3.4537


n=400; % matrix size is nxn
density=0.001; % density of the sparse matrices
sparsef=sprand(n,n,density); % sparse version of matrix f
sparseg=sprand(n,n,density); % sparse version of matrix g
fullf=full(sparsef); % full version of the matrix f
fullg=full(sparseg); % full version of the matrix g

% reference implementation
tic;
refh=tconv(fullf, fullg);
duration=toc;
disp(['Algorithm ref: duration ' num2str(duration)])

% comparison of algorithms for sparse matrices

tic;
h=ref_tconv_1(sparsef, sparseg);
duration=toc;
errorMax=max(abs(refh(:)-h(:)));
densityh=nnz(h)/numel(h);
disp(['Algorithm 1  : duration ' num2str(duration) ', maximum of error '...
    num2str(errorMax) ', density ' num2str(densityh)])
clear('h');

tic;
h=ref_tconv_2(sparsef, sparseg);
duration=toc;
errorMax=max(abs(refh(:)-h(:)));
densityh=nnz(h)/numel(h);
disp(['Algorithm 2  : duration ' num2str(duration) ', maximum of error '...
    num2str(errorMax) ', density ' num2str(densityh)])
clear('h');

