%DEMO_NEXTFASTFFT  Next fast FFT number
%
%   This demo shows the behaviour of the |nextfastfft| function.
%
%   .. figure::
%
%      Benchmark of the FFT routine
%
%      The figure shows the sizes returned by the |nextfastfft| function
%      compared to using `nextpow2`. As can be seen, the |nextfastfft|
%      approach gives FFT sizes that are much closer to the input size.
%
%   .. figure::
%
%      Efficiency of the table
%
%      The figure show the highest output/input ratio for varying input
%      sizes. As can be seen, the efficiency is better for larger input
%      values, where the output size is at most a few percent larger than
%      the input size.
%
%   See also: nextfastfft

% Range to use for testing.
% It is important for this script that
% range_max = 2^nextpow2(range_max) = nextfastfft(range_max), so it must
% be a power of 2.
range_min=100;
range_max=1024;
r=range_min:range_max;

% r2 contains the next higher sizes using nextpow2
r2=2.^nextpow2(r);

% r3 contains the next higher sizes using nextfastfft
r3=nextfastfft(r);

figure(1);
plot(r,r,r,r2,r,r3);
xlabel('Input size.');
ylabel('FFT size.');
legend('Same size','nextpow2','nextfastfft','Location','SouthEast');

%% Efficiency analysis of the table
[dummy,table]=nextfastfft(1);

eff=table(2:end)./(table(1:end-1)+1);

figure(2);
semilogx(table(2:end),eff);
xlabel('Input size.');
ylabel('Output/input ratio.');
mean(eff)
