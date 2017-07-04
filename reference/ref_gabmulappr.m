function sym=ref_gabmulappr(T,p2,p3,p4,p5);
%GABMULMAT  Best Approximation by a Gabor multiplier.
%   Usage:  sym=gabmulappr(T,a,M);
%           sym=gabmulappr(T,g,a,M);
%           sym=gabmulappr(T,ga,gs,a,M);
%
%   Input parameters:
%         T     : matrix to be approximated
%         g     : analysis/synthesis window
%         ga    : analysis window
%         gs    : synthesis window
%         a     : Length of time shift.
%         M     : Number of channels.
%
%   Output parameters:
%         sym   : symbol
%
%   GABMULAPPR(T,g,a,M) will calculate the best approximation of the given
%   matrix T in the frobenius norm by a Gabor multiplier determined by the
%   symbol sym over the rectangular time-frequency lattice determined by a
%   and M.  The window g will be used for both analysis and synthesis.
%   IMPORTANT: The chosen Gabor system has to be a frame sequence!
%
%   GABMULAPPR(T,a,M) will do the same using an optimally concentrated,
%   tight Gaussian as window function.
%
%   GABMULAPPR(T,gs,ga,a) will do the same using the window ga for analysis
%   and gs for synthesis.
%
%   SEE ALSO: GABMUL, GABMULINV
% 
%   In this algorithm first the 'lower symbol' is calcualted, then the
%   so-called 'upper symbol'.
% 
%   The lower symbol is the inner product < T , P_\lambda > where P_\lambda 
%   are the projections gamma_\lambda \otimes g_\lambda . These operators
%   form a Bessel sequence and the lower symbol, the lower symbol is the 
%   analysis sequence of T using this Bessel sequence.
%
%   The upper symbol is the inner product < T , Q_\lambda > where 
%   Q_\lambda are the dual projections operator. Therefore the upper 
%   symbol is the analysis with the dual sequence (if the P have formed a 
%   frame). Because the 
%
%   References :
%   * P. Balazs, Basic Definition and Properties of Bessel Multipliers, 
%     Journal of Mathematical Analysis and Applications 325(1):571--585, 
%     January 2007.
%   * P.~Balazs, Hilbert-Schmidt operators and frames - classification, best
%     approximation by multipliers and algorithms,  International Journal 
%     of Wavelets, Multiresolution and Information Processing, accepted, 
%     to appear.
%   * H. G. Feichtinger, M. Hampjes, G. Kracher, "Approximation of matrices 
%     by Gabor multipliers" , IEEE Signal Procesing Letters Vol. 11, 
%     Issue 11, pp 883-- 886 (2004)
%
%   Author: P. Balazs (XXL)

complainif_argnonotinrange(nargin,3,5,mfilename);

L=size(T,1);

if size(T,2)~=L
  error('T must be square.');
end;

if nargin==3
  % Usage: sym=gabmulappr(T,a,M);
  a=p2;
  M=p3;
  ga=gabtight(a,M,L);
  gs=ga;
end;

if nargin==4
  % Usage: sym=gabmulappr(T,g,a,M);
  ga=p2;
  gs=p2;
  a=p3;
  M=p4;
end;
  
if nargin==5
  % Usage: sym=gabmulappr(T,ga,gm,a,M);
  ga=p2;
  gs=p3;
  a=p4;
  M=p5;
end;

N=L/a;
b=L/M;


%       gg = GBa(:,ii+jj*M);  % Element of analysis frame
%       hh = GBs(:,ii+jj*M);  % Element of synthesis frame
%       HS inner product : < T , g \tensor h> = < T h , g > =
part1=reshape(dgt(T',ga,a,M),M*N,L);
part2=reshape(dgt(part1',gs,a,M),M*N,M*N).';
lowsym = reshape(diag(part2),M,N);

GBa = frsynmatrix(frame('dgt',ga,a,M),length(ga));
% Gabor frame synthesis matrix
if ga ~= gs
   GBs = frsynmatrix(frame('dgt',gs,a,M),length(gs));
else
  GBs = GBa;
end

if ga == gs
  % HS Gram matrix
  Gram = abs(GBa'*GBa).^2;
else
  Gram = (GBs'*GBs) .* conj(GBa'*GBa); 
end

% The Gram matrix is square and Toeplitz.
%iGram=inv(Gram);
sym = reshape(Gram\lowsym(:),M,N);
% sym = involute(sym); % strange but true ?
 



