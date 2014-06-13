function demo_introframes % RUNASSCRIPT
%DEMO_INTROFRAMES Introduction to frames in finite dimension
%
%   This demonstration explains some basic concepts of frame theory using
%   examples in $\mathbb{R}^L$ ($L=3$,$L=2$). 
%
%   References: 


% Naimark theorem Parseval tight frame is a projection of a ONB from a
% higher dimension
% genral frame is a projection of a Riez basis froma higher dimension

% Tightness

% Equal norm

% Maximally robust (every n x n submatrix of F* is invertible)

% Equiangularity

% Symmetry https://www.math.auckland.ac.nz/~waldron/Preprints/Frame-symmetries/CA-02-060-RD.pdf

% If Φ is a frame then_
% (1) V ΦU is a frame for any invertible matrices U, V .
% (2) If Φ is tight frame/unit-norm tight frame, then aVΦU is tight
% frame/unit-norm tight frame for any unitary matrices U, V
% and a = 0.
%(3) If Φ is equal-norm, then aDΦU is equal-norm for any
%diagonal unitary matrix D, unitary matrix U , and a = 0.
%(4) If Φ is maximally robust, then DΦU is maximally robust
%for any invertible diagonal matrix D and any invertible
%matrix U .
%(5) If Φ is unit-norm tight frame and maximally robust, then
%DΦU is unit-norm tight frame and maximally robust for any
%unitary diagonal matrix D and any unitary matrix U .

F1 = eye(3);
plotfinframe(F1);

function plotfinframe(A,B)
% Funkce vykresli v R^3 system generatoru obsazeny ve sloupcich A,
% pokud je zadano, tak i B, pak jinou barvou vykresli i druhy system
% generatoru.
% h....handle na osy, v nichz probehlo vykresleni

%% Kresleni
[~,n] = size(A);
% puvodni frejm
zer = zeros(n,1);
quiver3(zer,zer,zer,A(1,:)',A(2,:)',A(3,:)','b');
hold on;
for cnt = 1:n
    text(A(1,cnt), A(2,cnt), A(3,cnt),['e_' num2str(cnt)]);
end

xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

% druhy frejm
if nargin>1
    quiver3(zer,zer,zer,B(1,:)',B(2,:)',B(3,:)','r');
    for cnt = 1:n
        text(B(1,cnt), B(2,cnt), B(3,cnt),['f_' num2str(cnt)]);
    end
end
axis equal;
hold off;


function plotCurrent(A,x,xhat)
%% Vykresleni v grafu
[m,n] = size(A);
if m ~= 3
    disp('POZOR: Graf je pro R^3, vstupni matice nevyhovuje, graf nebude vykreslen')
    return
end

figure(1);
clf;
% dva frejmy
h = plotfinframe(A);
zer = zeros(3,1);

% vykresleni puvodniho a rekonstruovaneho vektoru
bJsouTotozne = norm(x-xhat) < 1e-3;

axes(h)
h = plot3(x(1),x(2),x(3),'Xk'); % vykreslení pùvodního vektoru
set(h,'LineWidth',2);
set(h,'MarkerSize',12);
if bJsouTotozne
    text(x(1)+0.05,x(2),x(3),['x=xhat']) %jsou
else
    text(x(1)+0.05,x(2),x(3),['x']) %nejsou
end

% pokud se x a xhat lisi, pak se vykresli jeste to xhat
if ~bJsouTotozne
    h = plot3(xhat(1),xhat(2),xhat(3),'Xg'); % vykreslení rekonstruovaného vektoru
    text(xhat(1)+0.05,xhat(2),xhat(3),['xhat'])
    set(h,'LineWidth',2);
    set(h,'MarkerSize',12);
end
axis auto
axis equal

function [xhat,dualsystem,coefs] = demo_frejmy(A,f)

% Demo funkce, která ukazuje báze a frejmy jednak výpoètem, jednak v grafu
% (pøed prvním sputìním je tøeba naèíst data: load('frejmy-data.mat') )

%% Vypocty
% výpoèet souradnic v primarnim systému generátorù pomoci skalarnich
% soucinu s duálním systemem
[coefs,B] = dualcoefs(A,f);
dualsystem = B;

% rekonstruovany vektor
xhat = A*coefs;


%% Vykresleni v grafu
[m,n] = size(A);
if m ~= 3
    disp('POZOR: Graf je pro R^3, vstupni matice nevyhovuje, graf nebude vykreslen')
    return
end

figure
% dva frejmy
h = plotfinframe(A,B);

% vykresleni puvodniho a rekonstruovaneho vektoru
bJsouTotozne = norm(f-xhat) < 10e-3;

axes(h)
h = plot3(f(1),f(2),f(3),'Xk'); % vykreslení pùvodního vektoru
set(h,'LineWidth',2);
set(h,'MarkerSize',12);
if bJsouTotozne
    text(f(1)+0.05,f(2),f(3),['x=xhat']) %jsou
else
    text(f(1)+0.05,f(2),f(3),['x']) %nejsou
end

% pokud se x a xhat lisi, pak se vykresli jeste to xhat
if ~bJsouTotozne
    h = plot3(xhat(1),xhat(2),xhat(3),'Xg'); % vykreslení rekonstruovaného vektoru
    text(xhat(1)+0.05,xhat(2),xhat(3),['xhat'])
    set(h,'LineWidth',2);
    set(h,'MarkerSize',12);
end
axis auto
axis equal

