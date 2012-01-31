%TIMEPREADADJ  Comparison of implementations for spreadadj
%   Usage: timespreadadj;
%
%
% FJ : Here are the duration results on my computer under ubuntu 8.04 for
% n=300
% Matlab version : 7.0.4.352 (R14) Service Pack 2
% Octave version : 3.0.0
% Algo | Matlab   | Octave
%-----------------------------
% Ref  | 0.16746  | 14.048
% 1    | 0.038398 | 13.173
% 2    | 0.10909  | 14.582
% 3    | 0.073976 | 9.0905
% 4    | 0.012999 | 0.032307   Best speed for full matrix
% 5    | 0.2108   | 2.957
% 6    | 0.004917 | 0.0045241  Best speed for sparse matrix

n=300; % matrix size is nxn
sparseMat=sprand(n,n,0.1); % sparse version of the test matrix
fullMat=full(sparseMat); % full version of the test matrix

routinemax=6;

% reference implementation
tic;
refcadj=spreadadj(fullMat);
duration=toc;
disp(['Algorithm ref: duration ' num2str(duration)])

% comparison of algorithms for full matrices
tic;
cadj=ref_spreadadj_1(fullMat);
duration=toc;
disp(['Algorithm 1  : duration ' num2str(duration)])
clear('cadj');

tic;
cadj=ref_spreadadj_2(fullMat);
duration=toc;
disp(['Algorithm 2  : duration ' num2str(duration)])
clear('cadj');

tic;
cadj=ref_spreadadj_3(fullMat);
duration=toc;
disp(['Algorithm 3  : duration ' num2str(duration)])
clear('cadj');

tic;
cadj=ref_spreadadj_4(fullMat);
duration=toc;
disp(['Algorithm 4  : duration ' num2str(duration)])
clear('cadj');

% comparison of algorithms for full sparses matrices

tic;
cadj=ref_spreadadj_5(sparseMat);
duration=toc;
disp(['Algorithm 5  : duration ' num2str(duration)])
clear('cadj');

tic;
cadj=ref_spreadadj_6(sparseMat);
duration=toc;
disp(['Algorithm 6  : duration ' num2str(duration)])
clear('cadj');


