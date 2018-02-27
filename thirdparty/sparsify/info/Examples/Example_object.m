% Example using operator specified as function
%
% Copyright (c) 2007 Thomas Blumensath
%
% The University of Edinburgh
% Email: thomas.blumensath@ed.ac.uk
% Comments and bug reports welcome
%
% This file is part of sparsity Version 0.1
% Created: April 2007
%
% Part of this toolbox was developed with the support of EPSRC Grant
% D000246/1
%
% Please read COPYRIGHT.m for terms and conditions.

% Problem formulation
m=1024; % Sparse domain dimension
n=256;  % observation domain dimension
k=50;   % 50 Non-zero elements



% Make an indicator vector to select FFT sub-matrix
rp                  = randperm(m/2);
SubSet              = zeros(m/2,1);
SubSet(rp(1:n/2))   = 1;
SubSet              = logical(SubSet);

% Create object

MyOb = MyObjectName(SubSet);

% Generate sparse signal
s=zeros(m,1);
RANP=randperm(m);
s(RANP(1:k))=randn(k,1);
x=MyOb*s;


% Use any of the algorithms provided

[s_est, err_mse, iter_time]  = greed_gp(x,MyOb,m,'stopTol',k);

% See how we have done
 

figure(1)
subplot(4,1,1)
stem(s,'k')
title('original')
subplot(4,1,2)
stem(s_est,'k')
title('estimate')
sig_ms=x'*x/n;
figure(1)
subplot(4,1,3)
stem(20*log10(sig_ms./err_mse),'k')
title('SNR_x')
subplot(4,1,4)
stem(iter_time,'k')
title('time')
