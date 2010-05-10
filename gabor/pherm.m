function g=pherm(L,order,tfr)
%PHERM  Periodized Hermite function.
%   Usage: g=pherm(L,order);
%          g=pherm(L,order,tfr);
% 
%   Input parameters:
%        L : Length of vector.
%    order : Order of Hermite function.
%      tfr : ratio between time and frequency support.
%   Output parameters:
%        g : The periodized Gaussian(s).
%
%   PHERM(L,order,tfr) computes samples of a periodized Hermite function
%   of order _order. order is counted from 0, so the zeroth order Hermite
%   function is the Gaussian.
%
%   The returned functions are eigenvectors of the DFT. The first four 
%   Hermite functions are orthonormal, but in general they are not.    
%
%   The parameter tfr determines the ratio between the effective
%   support of g and the effective support of the DFT of g. If tfr>1 then
%   g has a wider support than the DFT of g.
%
%   PHERM(L,order) does the same setting tfr=1.
%
%   If _order is a vector, PHERM will return a matrix, where each column
%   is a Hermite function with the corresponding order.
%
%   If tfr is a vector, PHERM will return a matrix, where each column
%   is a Hermite function with the corresponding tfr.
%
%   If both _order and tfr are vectors, they must have the same length,
%   and the values will be paired.
%
%   See also:  hermbasis, pgauss, psech

% AUTHORs: Thomasz Hrycak and Peter Soendergaard.
% 

error(nargchk(1,3,nargin));
  
if nargin==2
  tfr=1;
end;

if size(L,1)>1 || size(L,2)>1
  error('L must be a scalar');
end;

if rem(L,1)~=0
  error('L must be an integer.')
end;

% Parse tfr and order.
if sum(1-(size(tfr)==1))>1
  error('tfr must be a scalar or vector');
end;

if sum(1-(size(order)==1))>1
  error('"order" must be a scalar or vector');
end;

Ltfr=length(tfr);
Lorder=length(order);

if Ltfr>1 && Lorder>1 && Ltfr~=Lorder
  error('order and tfr must be vectors of same length');
end;

% Figure out how many vectors to compute: W
W=max(Ltfr,Lorder);
tfr=tfr(:);
order=order(:);

% Repeat tfr and order so they both have length W
if Ltfr==1
  tfr=repmat(tfr,1,W);
end;

if Lorder==1
  order=repmat(order,1,W);
end;

% Calculate W windows.
g=zeros(L,W);
for w=1:W
  
  thisorder=order(w);
  thisw=tfr(w);
    
  % These numbers have been computed numerically.
  if thisorder<=6
    safe=4;
  else if thisorder<=18
      safe=5;
    else if thisorder<=31
	safe=6;
      else if thisorder<=46
	  safe=7;
	else
	  % Anything else, use a high number.
	  safe=12;
	end;
      end;
    end;
  end;
  
  % Outside the interval [-safe,safe] then H(thisorder) is numerically zero.
  nk=ceil(safe/sqrt(L/sqrt(thisw)));
  
  sqrtl=sqrt(L);
  
  lr=(0:L-1).';
  for k=-nk:nk
    xval=(lr/sqrtl-k*sqrtl)/sqrt(thisw);
    g(:,w)=g(:,w)+comp_hermite(thisorder, sqrt(2*pi)*xval);
  end;
  
  % Normalize it.
  g(:,w)=g(:,w)/norm(g(:,w));
  
end;
  