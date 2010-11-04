%DEMO_GABLASSO  Sparse regression by Lasso method

% Signals
x0 = sin(2*pi*64*linspace(0,1,512));
x=x0';
x=x+randn(size(x))/2;

% DCGT parameters
g = 'tight';
a=2;
M=128;
M2 = floor(M/2);

% Regression parameters
lambda = 0.08;

% LASSO
[tcl,relres,iter,xrecl] = gablasso(x,g,a,M,lambda);

% GLASSO
[tcgl,relres,iter,xrecgl] = gabgrouplasso(x,g,a,M,lambda);

% Displays
figure(1);
subplot(2,2,1);plot(x0); axis tight; grid; title('Original')
subplot(2,2,2);plot(x); axis tight; grid; title('Noisy')
subplot(2,2,3);plot(real(xrecl)); axis tight; grid; title('LASSO')
subplot(2,2,4);plot(real(xrecgl)); axis tight; grid; title('GLASSO')

figure(2);
tmp = abs(dgt(x0,g,a,M));
subplot(2,2,1); imagesc(tmp(1:M2,:)); axis xy; title('Original')
tmp = abs(dgt(x,g,a,M));
subplot(2,2,2); imagesc(tmp(1:M2,:)); axis xy; title('Noisy')
subplot(2,2,3); imagesc(abs(tcl(1:M2,:))); axis xy; title('LASSO')
subplot(2,2,4); imagesc(abs(tcgl(1:M2,:))); axis xy; title('GLASSO')
