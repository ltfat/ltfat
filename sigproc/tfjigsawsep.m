function [ fplanes,cplanes,info ] = tfjigsawsep( f, varargin )
%TFJIGSAWSEP Time frequency jigsaw puzzle tonal-transient separation
%   Usage:  fplanes = tfjigsawsep(f);
%           fplanes = tfjigsawsep(f,t1,t2);
%           fplanes = tfjigsawsep(f,t1,t2,p);
%           [fplanes, cplanes, info] = tfjigsawsep(...);
%
%   Input parameters:
%            f        : Input signal
%            t1       : Significance level of the tonal layer refered to
%                       a white noise reference
%            t2       : Same for the transient layer
%            p        : Length of the shorter side of the supertiles
%    
%   Output parameters:
%           fplanes   : signallength-by-3 array containing the 3 produced
%                       layers, tonal in fplanes(:,1), transient in
%                       fplanes(:,2) and the noisy residual in fplanes(:,3).
%           cplanes   : 3-by-1 cellarray containing the Gabor coefficients
%                       for the single layers
%
%   `tfjigsawsep(f)` applies the separation algorithm on the input signal *f*
%   and returns the tonal, the ransient and the residual parts.
%   The default default values of the parameters is *t1=t2=0.95* and *p=4*. 
%   The parameters of the 3 Gabor systems used are the following:
%   
%       "Tonal" system:     g1={'hann',4096}; a1 = 512; M1 = 4096;
%   
%       "Transient" system: g2={'hann',256};  a2 = 32;  M2 = 256;
%
%       "Residual" system:  g3={'hann',2048}; a3 = 512; M3 = 2048;
%   
%   `tfjigsawsep(f,t1,t2)` works as before, but allows changing threshold
%   parameters *t1* and *t2*. Good values are in range [0.85,0.95]. 
%   *t2* sometimes has to be chosen larger (~ 1.05), eg. for 
%   percussions in music signals.
%
%   `tfjigsawsep(f,t1,t2,p)` additionally allows changing the size of
%   the supertiles *p*. It must be a positive integer, good values are
%   in the range of p=1,..,10, but it depends very much on the type of signal.
%   E.g. for speech signals, higher values yield better results.
%
%   `[fplanes, cplanes, info] = tfjigsawsep(...)` additionally returns
%   a 3 element cell array *cplanes* containing |dgtreal| coefficients
%   of the respective separated layers and a structure *info*, with the
%   parameters used in the algorithm. The relationship between *fplanes*
%   and *cplanes* is the following:

%   Additional parameters:
%   ----------------------
%
%   The function accepts the following flags and key-value pairs:
%
%       'ver2' Use the second version of the algorithm.
%
%       'plot' Plot the separated waveforms and the coefficients.
%
%       'maxit',maxit  Maximum number of iterations. The default value is 15.   
%
%
%   Algorithm:
%   ----------
%
%   The algorithm is based on [1]. It transforms a signal into a two-windowed
%   Gabor expansion such that one wide window shall lead to a high frequency
%   resolution (tonal layer is represented well) and a narrow one to a high
%   time resolution (transient layer is repr. well). The resulting Gabor
%   coefficients in the time-frequency plane are grided adaptively into rectangular
%   'supertiles'. An entropy criterion chooses those tiles, where
%   the tonal part of the signal is represented better and is below a given
%   threshold. The rest is set to zero. The leftover Gabor coefficients are
%   transformed back and substracted from the original signal. Then the
%   same is applied again to choose those tiles, where the transient
%   part is represented better. After that, one gets the first
%   approximation of the two layers. By applying this procedure iteratively
%   on the residual, tonal and transient layers emerge.
%
%   Examples:
%   ---------
%   
%   The following example shows the basic usage of the function:::
%     
%     % Load the glockenspiel test signal and add some noise
%     [f,fs] = gspi; f = f + 0.001*randn(size(f));
%     % Setup the parameters
%     p = 2;
%     t1 = 0.92;
%     t2 = 0.96;
%     [fplanes,cplanes,info] = tfjigsawsep(f,t1,t2,p,'plot','ver2');
%     
%
%   See also: dgtreal plottfjigsawsep
%         
%   References: jato07

%AUTHOR: Daniel Haider, 2017 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Remarks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A few improvements in the algorithm can still be done such as         %
%   designing the window settings, the tf-tilling and the residual        %
%   condition a bit more flexible.                                        %
%   It would also be useful to provide parameter settings for specific    %
%   types of signals (ie. speech, music,...) .                            %
%                                                                         %
%   Version 1 of the algorithm (default) works particularly well for      %
%   speech signals, but also for depicting the transient layer            %
%   (eg. percussive elements) nicely in musical signals.                  %
%   Version 2 works particularly well for depicting a nice tonal layer.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.t1=[];
definput.keyvals.t2=[];
definput.keyvals.p=4;
definput.keyvals.maxit=15;
definput.keyvals.fs = [];
definput.flags.algver={'ver1','ver2'};
definput.flags.plot={'noplot','plot'};
[flags,kv]=ltfatarghelper({'t1','t2','p'},definput,varargin);

t1 = kv.t1;
t2 = kv.t2;
p = kv.p;

if xor(isempty(t1),isempty(t2))
    error('%s: Both t1 and t2 must be defined.',upper(mfilename));
else
    if isempty(t1), t1 = 0.95; end
    if isempty(t2), t2 = 0.95; end
end

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

% length of time shifts and number of frequency channels
M1 = winsize1;
a1 = winsize1/8;
M2 = winsize2;
a2 = winsize2/8;

% adjust signal length
% to the smallest multiple of a1 and M1 AND a2 and M2
% L2 = dgtlength(Ls,a2,M2);
% L = dgtlength(L2,a1,M1);
%L1 = dgtlength(Ls,a1,M1);
L = dgtlength(Ls,lcm(a1,a2),lcm(M1,M2));

complainif_notposint(p,'p',mfilename);
if ~( p > 0 && p < min([L/a1,M2/2+1]))
    error('%s: p must be in range ]0,%i]',...
          upper(mfilename),min([L/a1,M2/2+1]));
end

% supertile sizes are adapted to p in order to respresent the same 
% part of the signal
T1 = p;
F1 = (a1/a2)*p;
T2 = (M1/M2)*p;
F2 = p;

% entropy reference from noise signal
% thresholds tau1,tau2 for tonal resp. transient layer are chosen to have
% a certain significance with respect to the estimated reference
[ref1,ref2] = noisest(L,g1,g2,a1,a2,M1,M2,p);
tau1 = ref1*t1;
tau2 = ref2*t2;

% initialization of residual, layers, min and max number of iterations
% and epsilon, the upper limit for the residual condition, which is
% computed at the kmin-th iteration
R = postpad(f,L);
l1 = zeros(L,1);
l2 = zeros(L,1);
k = 1;
kmin = 5;
kmax = kv.maxit;
epsilon = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% runs until all values in R are below epsilon to exclude single peaks 
while all(abs(R) < epsilon) ~= 1

    switch flags.algver
       case 'ver1'
        [x1,x2]=jigsaw1(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2);
       case 'ver2'
        [x1,x2]=jigsaw2(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2);
    end
    
    % residual parts of the signal
    R = R-x1-x2;
    l1 = l1+x1;
    l2 = l2+x2;
    
    if k == kmax
        disp('max. number of iteration reached: condition for residual is not fulfilled!')
        disp('Try to increase t2 slightly or change p')
        break
    end
    
    if k == kmin
        % epsilon is computed as upper limit, corresponding to
        % an empirical p-quantile
        epsilon = sort(abs(R),'ascend');
        epsilon = 5/2*epsilon(round(0.995*numel(epsilon)));
        disp(['The current maximum in the residual is ',num2str(max(abs(R)))])
        disp(['An upper bound is estimated to be ',num2str(epsilon)])
    end
    
    k = k+1;
end

% fplanes as Ls-by-3 array containg the layers
fplanes = [l1,l2,R];
fplanes = postpad(fplanes,Ls);

% info structure
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
info.noiseentropy_tonal = ref1;
info.threshold_tonal = tau1;
info.noiseentropy_transient = ref2;
info.threshold_transient = tau2;

% cplanes as 3-by-1 cell array containing
% the gabor coefficients corresp. to the layers
if nargout > 1
    cplanes = cell(3,1);
    cplanes{1} = dgtreal(l1,g1,a1,M1);
    cplanes{2} = dgtreal(l2,g2,a2,M2);
    cplanes{3} = dgtreal(R,info.g3,info.a3,info.M3);
end

% option for plots
if flags.do_plot
    plottfjigsawsep(fplanes,cplanes,info,'fs',kv.fs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% compiling functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ r ] = renyientropy( t,alpha )

% computes the Renyi entropy for an array
% yields high values for peaky data and
% low values for almost constant 
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

function [ref1,ref2] = noisest (L,g1,g2,a1,a2,M1,M2,p)

% computes the entropy for a white noise signal within one supertile
% as estimation reference for the tresholds tau1,tau2
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


function [x1,x2] = jigsaw1(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2)

%%%%%%%%%%%%%%% Version 1 of the jigsaw puzzle algorithm %%%%%%%%%%%%%%%%%%

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

% size of gabor coefficient matrices
d1 = M1/2+1;
e1 = L/a1;
d2 = M2/2+1;
e2 = L/a2;

% tilling of the tf-plane
ft1 = floor(d1/F1);
tt1 = floor(e1/T1);
ft2 = floor(d2/F2);
tt2 = floor(e2/T2);

% gabor transformation for real valued signals
c1 = dgtreal(R,g1,a1,M1,L);
c2 = dgtreal(R,g2,a2,M2,L);

% block matrix formation corresp. to the supertile sizes (with rest)
C1 = mat2cell(c1,[F1*ones(1,ft1),mod(d1,F1)],[T1*ones(1,tt1),mod(e1,T1)]);
C2 = mat2cell(c2,[F2*ones(1,ft2),mod(d2,F2)],[T2*ones(1,tt2),mod(e2,T2)]);

% decision criteria
%   case 1: both entropies are above their thresholds -> set them to zero
%   case 2: g1 has a better repr. and the entropy is below its
%           threshold -> set the tile corr. to g2 to zero
%   case 3: vice versa
%   case 4: g1 has a better repr. and its entropy is above its
%           threshold tau1 but despite of that, the entropy corr. to
%           g2 is below its threshold tau2 -> set the tile corr. to
%           g1 to zero
%   case 5: vice versa
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

% coefficients back into matrix form
c1 = cell2mat(C1);
c2 = cell2mat(C2);

% synthesis with the canonical dual windows
x1 = idgtreal(c1,{'dual',g1},a1,M1,L);
x2 = idgtreal(c2,{'dual',g2},a2,M2,L);


function [x1,x2] = jigsaw2(R,g1,g2,a1,a2,M1,M2,p,tau1,tau2)

%%%%%%%%%%%%%%% Version 2 of the jigsaw puzzle algorithm %%%%%%%%%%%%%%%%%%

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

% size of coefficient matrices
d1 = M1/2+1;
e1 = L/a1;
d2 = M2/2+1;
e2 = L/a2;

% tilling of the tf-plane
ft1 = floor(d1/F1);
tt1 = floor(e1/T1);
ft2 = floor(d2/F2);
tt2 = floor(e2/T2);

% gabor tranformation for real valued signals
c1 = dgtreal(R,g1,a1,M1,L);
c2 = dgtreal(R,g2,a2,M2,L);

% block matrix formation corresp. to the supertile sizes (with remainer)
C1 = mat2cell(c1,[F1*ones(1,ft1),mod(d1,F1)],[T1*ones(1,tt1),mod(e1,T1)]);
C2 = mat2cell(c2,[F2*ones(1,ft2),mod(d2,F2)],[T2*ones(1,tt2),mod(e2,T2)]);

% decision criteria:
% keep only those tiles where there is a better tonal representation of the
% signal
for j=1:ceil(e1/T1)
    for i=1:ceil(d1/F1)
        E1 = renyientropy(C1{i,j});
        E2 = renyientropy(C2{i,j});
        if min([E1,E2]) == E2 || E1 > tau1
            C1{i,j} = C1{i,j}-C1{i,j};
        end
    end
end

% reminding coefficients back into matrix form
% synthesis with the canonical dual window of g1
% RR is a first instance of the residual
c1 = cell2mat(C1(:,:));
x1 = idgtreal(c1,{'dual',g1},a1,M1,L);
RR = R-x1;

% the same procedure is now applied with respect to the narrow window g2
c1 = dgtreal(RR,g1,a1,M1,L);
c2 = dgtreal(RR,g2,a2,M2,L);
C1 = mat2cell(c1,[F1*ones(1,ft1),mod(d1,F1)],[T1*ones(1,tt1),mod(e1,T1)]);
C2 = mat2cell(c2,[F2*ones(1,ft2),mod(d2,F2)],[T2*ones(1,tt2),mod(e2,T2)]);

for j=1:ceil(e1/T1)
    for i=1:ceil(d1/F1)

        E1 = renyientropy(C1{i,j});
        E2 = renyientropy(C2{i,j});

        if min([E1,E2]) == E1 || E2 > tau2
            C2{i,j} = C2{i,j}-C2{i,j};
        end
    end
end

c2 = cell2mat(C2(:,:));
x2 = idgtreal(c2,{'dual',g2},a2,M2,L);
