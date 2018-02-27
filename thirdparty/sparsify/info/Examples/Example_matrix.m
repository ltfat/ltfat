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
n=1024; % Sparse domain dimension
m=256;  % observation domain dimension
k=50;   % 50 Non-zero elements



% Make an indicator vector to select FFT su-bmatrix
RANP=randperm(n);
FFT_subset=zeros(n,1);
FFT_subset(RANP(1:m))=1;
% Define functions

P = randn(m,n);
% Normalise
for i=1:n
    P(:,i)=P(:,i)/norm(P(:,i));
end

% Generate sparse signal
s=zeros(n,1);
RANP=randperm(n);
s(RANP(1:k))=randn(k,1);
x=MyOp ( s );


% Use any of the algorithms provided

s_approx = greed_omp(MyOp,x,n,'stopTol',k);

% See how we have done

figure(1)
subplot(2,1,1)
stem(s,'k')
subplot(2,1,2)
stem(s_approx,'k')
