function ggd = pgrpdelay(g,L)
%PGRPDELAY Group delay of a filter with periodic boundaries
%   Usage:  ggd = pgrpdelay(g,L);
%
%   `pgrpdelay(g,L)` computes group delay of filter *g* as a negative
%   derivative of the phase frequency response of filter *g* assuming
%   periodic (cyclic) boundaries i.e. the delay may be a negative number.
%   The derivative is calculated using the second order centered difference
%   approximation. 
%   The resulting group delay is in samples.
%
%   Example:
%   --------
%
%   The following example shows a group delay of causal, moving average
%   6tap FIR filter and it's magnitude frequency response for comparison.
%   The dips in the group delay correspond to places where modulus of the 
%   frequency response falls to zero.:::
%
%      g = struct(struct('h',ones(6,1),'offset',0)); 
%      L = 512;
%      figure(1);
%      subplot(2,1,1);
%      plot(-L/2+1:L/2,fftshift(pgrpdelay(g,512)));
%      axis tight;ylim([0,4]);
%
%      subplot(2,1,2);
%      magresp(g,L,'nf');
%      


g = comp_fourierwindow(g,L,upper(mfilename));

H = comp_transferfunction(g,L);
Harg = angle(H);

% Forward approximation
tgrad_1 = Harg-circshift(Harg,-1);
tgrad_1 = tgrad_1 - 2*pi*round(tgrad_1/(2*pi));
% Backward approximation
tgrad_2 = circshift(Harg,1)-Harg;
tgrad_2 = tgrad_2 - 2*pi*round(tgrad_2/(2*pi));
% Average
ggd = (tgrad_1+tgrad_2)/2;
 
ggd = ggd/(2*pi)*L;

