function [ fplanes,cplanes,info ] = tfjigsawsep( f,varargin )
%TFJIGSAWSEP separation of tonal, transient and noisy parts of a signal
%    Usage:  fplanes = tfjigsawsep(f);
%            [fplanes, cplanes, info] = tfjigsawsep(f);
%
%    Input parameters:
%            f        : Input data
%            t1       : Significance level of the tonal layer refered to
%                       a white noise reference
%            t2       : Same for the transient layer
%            p        : Length of the shorter side of the supertiles
%    
%    Output parameters:
%           fplanes   : signallength-by-3 array containing the 3 produces
%                       layers, tonal in fplanes(:,1), transient in
%                       fplanes(:,2) and noisy residual in fplanes(:,3).
%           cplanes   : ?-by-? array containing the gabor coefficients for
%                       the single layers
%
%           fplanes = tfjigsawsep(f) applies the separation algorithm
%           on the wav-file sig.wav with the default values t1=t2=0.95 and 
%           p=4.
%
%           [...] = tfjigsawsep(f,t1,t2,p) is the most recommended way
%           of using this function. One has to tweak and play with the
%           parameters a bit in order to get good resluts. t1 is usually
%           best for values within [0.85,0.95], t2 sometimes has to be
%           chosen larger (~ < 1.05), eg. for percussions in music signals. 
%           p has to be a positive integer, it works good for p=1,..,10. It
%           depends very much on the type of signal how to choose p, eg.
%           for speech signals, higher values yield better results.
%
%           [...] = tfjigsawsep(...,'ver2') uses an alternative algorithm
%           where usually the tonal layer is very nice.
%
%           [...] = tfjigsawsep(...,'plot') produces a figure with phase
%           and TF-plots and buttons to automatically play the single
%           layers.
%
%
%           [fplanes,cplanes,info] = tfjigsawsep(f) additionally
%           produces an infostructure.
%
%
%           See also: demo_jigsaw1, demojigsaw2, demojigsaw2synth
%         
%           References:
%             F. Jaillet, B. Torresani. Time-Frequency Jigsaw Puzzle:
%             Adaptive Multiwindow and Multilayered Gabor Expansion.
%             IJWMIP, vol.5, pp. 293-315, 2007.
%
%           


% ltfatarghelper:

% no input (except x) -> default
% one input -> t = t1 = t2
% two inputs -> t1 and t2
% three inputs -> t1,t2 and p (error for 0 > p => min(L/a1,M2/2+1) and non_int)
% string 'ver2' for the alternative algorithm
% string 'plot' for the figure
% string 'wintype' for different window ??
% string 'speech' with special large epsilon

definput.keyvals.t1=0.95;
definput.keyvals.t2=0.95;
definput.keyvals.p=4;
definput.keyvals.maxit=15;
definput.keyvals.fs = [];
definput.flags.algver={'ver1','ver2'};
definput.flags.plot={'noplot','plot'};
[flags,kv]=ltfatarghelper({},definput,varargin);

t1 = kv.t1;
t2 = kv.t2;
p = kv.p;

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);

if W>1
    error('%s: Multichannel inputs are not supported.',upper(mfilename));
end

% window settings
wintype = 'hann';
winsize1 = 4096;
winsize2 = 256;
g1 = {wintype,winsize1};
g2 = {wintype,winsize2};

% length of time shifts and channels
M1 = winsize1;
a1 = winsize1/8;
M2 = winsize2;
a2 = winsize2/8;

%adjust signal length
% L2 = dgtlength(Ls,a2,M2); % don't need that??
% L = dgtlength(L2,a1,M1);
%L1 = dgtlength(Ls,a1,M1);
L = dgtlength(Ls,lcm(a1,a2),lcm(M1,M2));

complainif_notposint(p,'p',mfilename);
if ~( p > 0 && p < min([L/a1,M2/2+1]))
    error('%s: p must be in range ]0,%i]',...
          upper(mfilename),min([L/a1,M2/2+1]));
end

% entropy reference from noise signal
[ref1,ref2] = noisest(L,g1,g2,a1,a2,M1,M2,p);
tau1 = ref1*t1;
tau2 = ref2*t2;

% initialization
R = postpad(f,L);
l1 = zeros(L,1);
l2 = zeros(L,1);

T1 = p;
F1 = (a1/a2)*p;
T2 = (M1/M2)*p;
F2 = p;

k = 1;
kmin = 5;
kmax = kv.maxit;
epsilon = 0;

% main loop
while all(abs(R) < epsilon) ~= 1

    switch flags.algver
       case 'ver1'
        [x1,x2]=jigsaw1(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2);
       case 'ver2'
        [x1,x2]=jigsaw2(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2);
    end
    
    R = R-x1-x2;
    l1 = l1+x1;
    l2 = l2+x2;
    
    if k == kmax
        disp('max. number of iteration reached: condition for residual is not fulfilled!')
        disp('Try to increase t2 slightly or change p')
        break
    end
    
    if k == kmin
        epsilon = sort(abs(R),'ascend');
        epsilon = 5/2*epsilon(round(0.98*numel(epsilon)));
        disp(['The current maximum in the residual is ',num2str(max(abs(R)))])
        disp(['An upper bound is estimated to be ',num2str(epsilon)])
    end
    
    k = k+1;
end

% ltfatarghelper?

fplanes = [l1,l2,R];
fplanes = postpad(fplanes,Ls);

info.g1 = {wintype, winsize1};
info.g2 = {wintype, winsize2};
info.g3 = {wintype, 2048};
info.M1 = M1;
info.M2 = M2;
info.M3 = 2048;
info.a1 = a1;
info.a2 = a2;
info.a3 = 512;
info.supertilesizes_g1 = [F1,T1];
info.supertilesizes_g2 = [F2,T2];

if nargout>1
        % it doesn't work as array since c1,c2,c3 have different sizes
        cplanes = cell(3,1);
        cplanes{1} = dgtreal(l1,g1,a1,M1);
        cplanes{2} = dgtreal(l2,g2,a2,M2);
        cplanes{3} = dgtreal(R,info.g3,info.a3,info.M3);
end


if flags.do_plot
    plottfjigsawsep(fplanes,cplanes,info,'fs',kv.fs);
end

function [ref1,ref2] = noisest (L,g1,g2,a1,a2,M1,M2,p)

T1 = p;
F1 = (a1/a2)*p;
T2 = (M1/M2)*p;
F2 = p;

n = noise(L,'white');
n1 = dgtreal(n,g1,a1,M1);
n2 = dgtreal(n,g2,a2,M2);
r1 = n1(1:F1,1:T1);
r2 = n2(1:F2,1:T2);
ref1 = renyientropy(r1);
ref2 = renyientropy(r2);


function [ r ] = renyientropy( t,alpha )
% Renyi Entropy


switch nargin
    case 1
        alpha = 0.5;
        r = (1/(1-alpha))*log2(sum(sum(abs(t).^(2*alpha))).*norm(t(:),2)^(-2*alpha));
    case 2
        if alpha <= 0 || alpha >= 1
            error('alpha must be chosen within (0,1).');
        else
            r = (1/(1-alpha))*log2(sum(sum(abs(t).^(2*alpha))).*norm(t(:),2)^(-2*alpha));
        end
    otherwise
        error('This function takes at least one input argument.');
end


function [x1,x2] = jigsaw1(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2)

if nargin <7
    error('Not enough input arguments.')
elseif iscolumn(R) ~= 1
    error('R has to be a column vector!')
elseif mod(a1/a2,1)~=0 || mod(M1/M2,1)~=0
    error('a2 has to divide a1 and M2 has to divide M1!')
elseif p <= 0 || mod(p,1)~=0
    error('p has to be a positive integer!')
end

% signal length already must have been adjusted!!
L = length(R);

% supertile sizes
T1 = p;
F1 = (a1/a2)*p;
T2 = (M1/M2)*p;
F2 = p;

% size of coef. matrix
d1 = M1/2+1;
e1 = L/a1;
d2 = M2/2+1;
e2 = L/a2;

% tilling
ft1 = floor(d1/F1);
tt1 = floor(e1/T1);
ft2 = floor(d2/F2);
tt2 = floor(e2/T2);

c1 = dgtreal(R,g1,a1,M1,L);
c2 = dgtreal(R,g2,a2,M2,L);

C1 = mat2cell(c1,[F1*ones(1,ft1),mod(d1,F1)],[T1*ones(1,tt1),mod(e1,T1)]);
C2 = mat2cell(c2,[F2*ones(1,ft2),mod(d2,F2)],[T2*ones(1,tt2),mod(e2,T2)]);

for j=1:ceil(e1/T1)
    for i=1:ceil(d1/F1)

        E1 = renyientropy(C1{i,j});
        E2 = renyientropy(C2{i,j});

        if E1 > tau1 && E2 > tau2
            C1{i,j}(:) = 0;               
            C2{i,j}(:) = 0;
        else
            if (min([E1,E2]) == E1 && E1 < tau1) || (min([E1,E2]) == E2 && E2 > tau2 && E1 < tau1)
                C2{i,j}(:) = 0;
            elseif (min([E1,E2]) == E2 && E2 < tau2) || (min([E1,E2]) == E1 && E1 > tau1 && E2 < tau2)
                C1{i,j}(:) = 0;
            else
                warning('something went wrong with the decision criteria!');
            end
        end
    end
end

c1 = cell2mat(C1);
c2 = cell2mat(C2);

x1 = idgtreal(c1,{'dual',g1},a1,M1,L);
x2 = idgtreal(c2,{'dual',g2},a2,M2,L);


function [x1,x2] = jigsaw2(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2)

if nargin <7
    error('Not enough input arguments.')
elseif iscolumn(R) ~= 1
    error('R has to be a column vector!')
elseif mod(a1/a2,1)~=0 || mod(M1/M2,1)~=0
    error('a2 has to divide a1 and M2 has to divide M1!')
elseif p <= 0 || mod(p,1)~=0
    error('p has to be a positive integer!')
end

% signal length already must have been adjusted!!
L = length(R);

% supertile sizes
T1 = p;
F1 = (a1/a2)*p;
T2 = (M1/M2)*p;
F2 = p;

% size of coef. matrix
d1 = M1/2+1;
e1 = L/a1;
d2 = M2/2+1;
e2 = L/a2;

% tilling
ft1 = floor(d1/F1);
tt1 = floor(e1/T1);
ft2 = floor(d2/F2);
tt2 = floor(e2/T2);

c1 = dgtreal(R,g1,a1,M1,L);
c2 = dgtreal(R,g2,a2,M2,L);

C1 = mat2cell(c1,[F1*ones(1,ft1),mod(d1,F1)],[T1*ones(1,tt1),mod(e1,T1)]);
C2 = mat2cell(c2,[F2*ones(1,ft2),mod(d2,F2)],[T2*ones(1,tt2),mod(e2,T2)]);
    
for j=1:ceil(e1/T1)
    for i=1:ceil(d1/F1)

        E1 = renyi(C1{i,j});
        E2 = renyi(C2{i,j});

        if min([E1,E2]) == E2 || E1 > tau1
            C1{i,j} = C1{i,j}-C1{i,j};
        end
    end
end

c1 = cell2mat(C1(:,:));
x1 = idgtreal(c1,{'dual',g1},a1,M1,L);
RR = R-x1;
    

c1 = dgtreal(RR,g1,a1,M1,L);
c2 = dgtreal(RR,g2,a2,M2,L);
C1 = mat2cell(c1,[F1*ones(1,ft1),mod(d1,F1)],[T1*ones(1,tt1),mod(e1,T1)]);
C2 = mat2cell(c2,[F2*ones(1,ft2),mod(d2,F2)],[T2*ones(1,tt2),mod(e2,T2)]);

for j=1:ceil(e1/T1)
    for i=1:ceil(d1/F1)

        E1 = renyi(C1{i,j});
        E2 = renyi(C2{i,j});

        if min([E1,E2]) == E1 || E2 > tau2
            C2{i,j} = C2{i,j}-C2{i,j};
        end
    end
end

c2 = cell2mat(C2(:,:));
x2 = idgtreal(c2,{'dual',g2},a2,M2,L);







