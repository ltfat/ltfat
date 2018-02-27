function TestSparsify(Algo,Option,P_format,k,m,n)
% Function to test the algorithms in sparsify Matlab toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage
%   TestSparsify(Algo,Option,P_format,k,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%       Algo    Algorithm name
%               Choose from:
%                   OMP
%                   OMP_qr
%                   OMP_chol
%                   OMP_cg
%                   OMP_cgp
%                   OMP_pinv
%                   OMP_linsolve
%                   OLS
%                   MP
%                   GP
%                   NOMP
%                   l0_Mterm
%                   l0_reg
%
%       Option  Options vector for each of the algorithms
%
%       P_format    Choose from:
%                       Matrix
%                       Object
%                       Function
%               
%   Optional (default):
%       k       k is sparsity of vector (4)
%       n x m   is size of P (64 128)
%       
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
% Please read COPYRIGHT.m for terms and conditions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Default values and initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    k=4;
end
if nargin < 5
    m=128;
end
if nargin < 6
    n=64;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Set up original sparse vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_SubSet            = zeros(m,1);
rp                  = randperm(m);
% rp = [ 4 18 42 56];
s_SubSet(rp(1:k))   = 1;
s_SubSet            = logical(s_SubSet);
s_original          = zeros(m,1);
s_original(s_SubSet)= randn(k,1);
% s_original(s_SubSet)= [ 1 -1 -0.5 0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        What problem to solve?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rp                  = randperm(m/2);
SubSet              = zeros(m/2,1);
SubSet(rp(1:n/2))   = 1;
% if SubSet(1)
%     SubSet(rp(n/2)) = 1;
% else
%     SubSet(1)=1;
% end
SubSet              = logical(SubSet);

% Pob

    Pob             = MyObjectName(SubSet);

% Pfunc

    Pfunc           = @(z)  MyOp_witharg(z,SubSet);
    Ptfunc          = @(z)  MyOpTranspose_witharg(z,SubSet);

    
% Pmat
if strcmp(P_format,'Matrix')
    Pmat = zeros(n,m);
    for i=1:m
        mask=zeros(m,1);
        mask(i)=1;
        Pmat(:,i)=Pfunc(mask);
    end
end


% Test function and object normalisation


%     Pfuncmat = zeros(n,m);
%     for i=1:m
%         mask=zeros(m,1);
%         mask(i)=1;
%         Pfuncmat(:,i)=Pfunc(mask);
%     end
%     
%     Ptfuncmat = zeros(m,n);
%     for i=1:n
%         mask=zeros(n,1);
%         mask(i)=1;
%         Ptfuncmat(:,i)=Ptfunc(mask);
%     end
%     
%     Pobmat = zeros(n,m);
%     for i=1:m
%         mask=zeros(m,1);
%         mask(i)=1;
%         Pobmat(:,i)=Pob*mask;
%     end
%     
%     Ptobmat = zeros(m,n);
%     for i=1:n
%         mask=zeros(n,1);
%         mask(i)=1;
%         Ptobmat(:,i)=Pob'*mask;
%     end
%     

    
% Create operator and observation

    
if strcmp(P_format,'Matrix')
    P = Pmat;
    x = P*s_original;
        
elseif strcmp(P_format,'Object')
   
   
    P  = Pob;
    x                   = P*s_original;
    
elseif strcmp(P_format,'Function')


    P  = Pfunc;
    Pt = Ptfunc;
    x  = P(s_original); 
    % If we use function handle, we need to include transpose in Option cell.
    lo=length(Option);
    Option{lo+1}='P_trans';
    Option{lo+2}=Pt;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        What algorithm to use?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(Algo,'OMP')
    [s_est, err_mse, iter_time]=greed_omp(x,P,m,Option);
    [s_est, err_mse]=greed_omp(x,P,m,Option);
    [s_est]=greed_omp(x,P,m,Option);
elseif strcmp(Algo,'MP')
    [s_est, err_mse, iter_time]=greed_mp(x,P,m,Option);
    [s_est, err_mse]=greed_mp(x,P,m,Option);
    s_est=greed_mp(x,P,m,Option);
elseif strcmp(Algo,'GP')
    [s_est, err_mse, iter_time]=greed_gp(x,P,m,Option);
    [s_est, err_mse]=greed_gp(x,P,m,Option);
    [s_est]=greed_gp(x,P,m,Option);
elseif strcmp(Algo,'NOMP') 
    [s_est, err_mse, iter_time]=greed_nomp(x,P,m,Option);
    [s_est, err_mse]=greed_nomp(x,P,m,Option);
    [s_est]=greed_nomp(x,P,m,Option);
elseif strcmp(Algo,'OMP_qr')
    [s_est, err_mse, iter_time]=greed_omp_qr(x,P,m,Option);
    [s_est, err_mse]=greed_omp_qr(x,P,m,Option);
    [s_est]=greed_omp_qr(x,P,m,Option);
elseif strcmp(Algo,'OMP_chol')
    [s_est, err_mse, iter_time]=greed_omp_chol(x,P,m,Option);
    [s_est, err_mse]=greed_omp_chol(x,P,m,Option);
    [s_est]=greed_omp_chol(x,P,m,Option);
elseif strcmp(Algo,'OMP_cgp')
    [s_est, err_mse, iter_time]=greed_omp_cgp(x,P,m,Option);
    [s_est, err_mse]=greed_omp_cgp(x,P,m,Option);
    [s_est]=greed_omp_cgp(x,P,m,Option);
elseif strcmp(Algo,'OMP_cg')
    [s_est, err_mse, iter_time]=greed_omp_cg(x,P,m,Option);
    [s_est, err_mse]=greed_omp_cg(x,P,m,Option);
    [s_est]=greed_omp_cg(x,P,m,Option);
elseif strcmp(Algo,'OMP_pinv')
    [s_est, err_mse, iter_time]=greed_omp_pinv(x,P,m,Option);
    [s_est, err_mse]=greed_omp_pinv(x,P,m,Option);
    [s_est]=greed_omp_pinv(x,P,m,Option);
elseif strcmp(Algo,'OMP_linsolve')
    [s_est, err_mse, iter_time]=greed_omp_linsolve(x,P,m,Option);
    [s_est, err_mse]=greed_omp_linsolve(x,P,m,Option);
    [s_est]=greed_omp_linsolve(x,P,m,Option);
elseif strcmp(Algo,'hard_l0_reg')
    lam=0.04;
    [s_est, err_mse, iter_time]=hard_l0_reg(x,P,m,lam,Option);
    [s_est, err_mse]=hard_l0_reg(x,P,m,lam,Option);
    [s_est]=hard_l0_reg(x,P,m,lam,Option);
elseif strcmp(Algo,'OLS')
    [s_est, err_mse, iter_time]=greed_ols(x,P,m,Option);
    [s_est, err_mse]=greed_ols(x,P,m,Option);
    [s_est]=greed_ols(x,P,m,Option);
elseif strcmp(Algo,'hard_l0_Mterm')
    M=4;
    [s_est, err_mse, iter_time]=hard_l0_Mterm(x,P,m,M,Option);
    [s_est, err_mse]=hard_l0_Mterm(x,P,m,M,Option);
    [s_est]=hard_l0_Mterm(x,P,m,M,Option);
else
    error('Unrecognised Algorithm specified.')
end

% s_est([4 18 42 56])

figure(1)
subplot(4,1,1)
stem(s_original,'k')
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
