function c=crestfactor(insig)
%CRESTFACTOR  Crest factor of input signal in dB
%   Usage:  c=crestfactor(insig);
%
%   `crestfactor(insig)` computes the crest factor of the input signal
%   *insig*. The output is measured in dB.
%
%   See also: rms, gaindb

c=20*log10(norm(insig,Inf)/rms(insig));

