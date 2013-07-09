function H=comp_warpedfreqresponse(wintype,fc,bw,fs,L,freqtoscale,scaletofreq,varargin)
%COMP_WARPEDFREQRESPONSE  Transfer function of warped filter
%   Usage: H=comp_warpedfreqresponse(wintype,fc,bw,fs,L,freqtoscale);
%          H=comp_warpedfreqresponse(wintype,fc,bw,fs,L,freqtoscale,normtype);
%
%   Input parameters:
%      wintype     : Type of window (from firwin)
%      fc          : Centre frequency, in scale units.
%      bw          : Bandwith, in scale units.
%      fs          : Sampling frequency in Hz.
%      L           : Transform length (in samples).
%      freqtoscale : Function to convert Hz into scale units.
%      scaletofreq : Function to convert scale units into Hz.
%      normtype    : Normalization flag to pass to |normalize|.

definput.import={'normalize'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Compute the values in Aud of the channel frequencies of an FFT of
% length L.
bins_lo   = freqtoscale(modcent(fs*(0:L-1)/L,fs)).';

% This one is necessary to represent the highest frequency filters, which
% overlap into the negative frequencies.
nyquest2  = 2*freqtoscale(fs/2);
bins_hi   = nyquest2+bins_lo;

% firwin makes a window of width 1 centered around 0 on the scale, so we rescale the
% bins in order to pass the correct width to firwin and subtract fc
bins_lo=(bins_lo-fc)/bw;
bins_hi=(bins_hi-fc)/bw;

% pos_lo is the same as foff
pos_lo=floor(scaletofreq(fc-.5*bw)/fs*L);
pos_hi=ceil(scaletofreq(fc+.5*bw)/fs*L);

if pos_hi>L/2
    % Filter is high pass and spilling into the negative frequencies
    pos_hi=ceil(scaletofreq(fc+.5*bw-nyquest2)/fs*L);    
end;

win_lo=firwin(wintype,bins_lo);
win_hi=firwin(wintype,bins_hi);

H=win_lo+win_hi;
   
H=normalize(H,flags.norm);

% Testing
H=circshift(H,-pos_lo);
H=H(1:modcent(pos_hi-pos_lo,L));
