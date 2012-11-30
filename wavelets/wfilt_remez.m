function [h,g,a]=wfilt_remez(L,K,B)

%REMEZFLT Generates orthonormal wavelet filters based on the Remez
%	      exchange algorithm.
%
%	      [H,G,RH,RG] = REMEZFLT(L,K,B) returns a set of wavelet
%	      filters. You can control regularity, frequency selectivity,
%	      and length of the filters. It works performing a factorization
%	      based on the complex cepstrum of the polynomial returned by
%	      REMEZWAV.
%
%	      L is the length of the filters. K is degree of flatness at
%	      z=-1. B is the normalized transition bandwidth.
%
%	      See also: REMEZWAV, FC_CCEPS.
%
%	      References: O. Rioul and P. Duhamel, "A Remez Exchange Algorithm
%			  for Orthonormal Wavelets", IEEE Trans. Circuits and
%			  Systems - II: Analog and Digital Signal Processing,
%			  41(8), August 1994

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es



poly=remezwav(L,K,B);
rh=fc_cceps(poly);

[g{1},g{2},h{1},h{2}]=rh2rg(rh);
a= [2;2];

function [p,r]=remezwav(L,K,B)

%REMEZWAV    P=REMEZWAV(L,K,B) gives impulse response of maximally
%	     frequency selective P(z), product filter of paraunitary
%	     filter bank solution H(z) of length L satisfying K flatness
%	     constraints (wavelet filter), with normalized transition
%	     bandwidth B (optional argument if K==L/2).
% 
%	     [P,R]=REMEZWAV(L,K,B) also gives the roots of P(z) which can
%	     be used to determine H(z).
%
%	     See also: REMEZFLT, FC_CCEPS.
%
%	     References: O. Rioul and P. Duhamel, "A Remez Exchange Algorithm
%			 for Orthonormal Wavelets", IEEE Trans. Circuits and
%			 Systems - II: Analog and Digital Signal Processing,
%			 41(8), August 1994
%                                                                          
%       Author: Olivier Rioul, Nov. 1, 1992 (taken from the
%		above reference)
%  Modified by: Jose Martin Garcia
%       e-mail: Uvi_Wave@tsc.uvigo.es
%--------------------------------------------------------


computeroots=(nargout>1);

%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(L,2), error('L must be even'); end
if rem(L/2-K,2), K=K+1; end
N=L/2-K;
%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2  %%%%%%%%%%%%%%%%%%%%%%%%%%
% Daubechies solution
% PK(z)=z^(-2K-1))+AK(z^2)
if K==0, AK=0;
else
   binom=pascal(2*K,1);
   AK=binom(2*K,1:K)./(2*K-1:-2:1);
   AK=[AK AK(K:-1:1)];
   AK=AK/sum(AK);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2' %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daubechies factor
% PK(z)=((1+z^(-1))/2)^2*K QK(z)
if computeroots & K>0
   QK=binom(2*K,1:K);
   QK=QK.*abs(QK);
   QK=cumsum(QK);
   QK=QK./abs(binom(2*K-1,1:K));
   QK=[QK QK(K-1:-1:1)];
   QK=QK/sum(QK)*2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output Daubechies solution PK(z)
if K==L/2
   p=zeros(1,2*L-1);
   p(1:2:2*L-1)=AK; p(L)=1;
   if computeroots
      r=[roots(QK); -ones(L,1)];
   end
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Daubechies polinomial
% PK(x)=1+x*DK(x^2)
if K==0, DK=0;
else
   binom=pascal(K,1);
   binom=binom(K,:);
   DK=binom./(1:2:2*K-1);
   DK=fliplr(DK)/sum(DK);
end

wp=(1/2-B)*pi;  % cut-off frequency
gridens=16*(N+1);  % grid density
found=0;  % boolean for Remez loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial estimate of yk
a=min(4,K)/10;
yk=linspace(0,1-a,N+1);
yk=(yk.^2).*(3+a-(2+a)*yk);
yk=1-(1-yk)*(1-cos(wp)^2);
ykold=yk;

iter=0;
while 1  % REMEZ LOOP
iter=iter+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute delta
Wyk=sqrt(yk).*((1-yk).^K);
Dyk=(1-sqrt(yk).*polyval(DK,yk))./Wyk;
for k=1:N+1
   dy=yk-yk(k); dy(k)=[];
   dy=dy(1:N/2).*dy(N:-1:N/2+1);
   Lk(k)=prod(dy);
end
invW(1:2:N+1)=2./Wyk(1:2:N+1);
delta=sum(Dyk./Lk)/sum(invW./Lk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute R(y) on fine grid
Ryk=Dyk-delta.*invW; Ryk(N+1)=[];
Lk=(yk(1:N)-yk(N+1))./Lk(1:N);
y=linspace(cos(wp)^2,1-K*1e-7,gridens);
yy=ones(N,1)*y-yk(1:N)'*ones(1,gridens);
% yy contain y-yk on each line
ind=find(yy==0);  % avoid division by 0
if ~isempty(ind)
   yy(ind)=1e-30*ones(size(ind));
end
yy=1./yy;
Ry=((Ryk.*Lk)*yy)./(Lk*yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP IV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find next yk
Ey=1-delta-sqrt(y).*(polyval(DK,y)+((1-y).^K).*Ry);
k=find(abs(diff(sign(diff(Ey))))==2)+1;
% N extrema
if length(k)>N
% may happen if L and K are large 
   k=k(1:N);
end
yk=[yk(1) y(k)];
% N+1 extrema including wp
if K==0, yk=[yk 1]; end
% extrema at y==1 added
if all(yk==ykold), break; end
ykold=yk;

end  % REMEZ LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  STEP A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute impulse response
w=(0:2*N-2)*pi/(2*N-1);
y=cos(w).^2;
yy=ones(N,1)*y-yk(1:N)'*ones(1,2*N-1);
ind=find(yy==0);
if ~isempty(ind)
   yy(ind)=1e-30*ones(size(ind));
end
yy=1./yy;
Ry=((Ryk.*Lk)*yy)./(Lk*yy);
Ry(2:2:2*N-2)=-Ry(2:2:2*N-2);
r=Ry*cos(w'*(2*(0:N-1)+1));
% partial real IDFT done
r=r/(2*N-1);
r=[r r(N-1:-1:1)];
p1=[r 0]+[0 r];
pp=p1;  % save p1 for later use
for k=1:2*K
   p1=[p1 0]-[0 p1];
end
if rem(K,2), p1=-p1; end
p1=p1/2^(2*K+1);
p1(N+1:N+2*K)=p1(N+1:N+2*K)+AK;
% add Daubechies response:
p(1:2:2*L-1)=p1; p(L)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP A' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute roots
if computeroots
   Q(1:2:2*length(pp)-1)=pp;
   for k=1:2*K
     Q=[Q 0]-[0 Q];
   end
   if rem(K,2), Q=-Q; end
   Q=Q/2;
   if K>0  % add Daubechies factor QK
      Q(2*N+1:L-1)=Q(2*N+1:L-1)+QK;
   else
      Q(L)=1;
   end
   r=[roots(Q); -ones(2*K,1)];
end






