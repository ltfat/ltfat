function H=comp_warpedfreqresponse(wintype,fc,bw,fs,L,freqtoscale,varargin)
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
%      normtype    : Normalization flag to pass to |normalize|.

definput.import={'normalize'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Compute the values in Aud of the channel frequencies of an FFT of
% length L.
bins_lo   = freqtoscale(modcent(fs*(0:L-1)/L,fs)).';

% This one is necessary to represent the highest frequency filters, which
% overlap into the negative frequencies.
bins_hi   = 2*freqtoscale(fs/2)+bins_lo;

% firwin makes a window of width 1 centered around 0 on the scale, so we rescale the
% bins in order to pass the correct width to firwin and subtract fc
bins_lo=(bins_lo-fc)/bw;
bins_hi=(bins_hi-fc)/bw;

win_lo=firwin(wintype,bins_lo);
win_hi=firwin(wintype,bins_hi);

H=win_lo+win_hi;
   
H=normalize(H,flags.norm);
