function outsig=expchirp(L,fstart,fend,varargin)
%EXPCHIRP  Exponential chirp
%   Usage: outsig=expchirp(L,fstart,fend)
% 
%   `expchirp(L,fstart,fend)` computes an exponential chirp of length *L*
%   starting at normalized frequency *fstart* and ending at frequency *fend*.
%
%   `expchirp` takes the following parameters at the end of the line of input
%   arguments:
%
%     'fs',fs    Use a sampling frequency of *fs* Hz. If this option is
%                specified, *fstart* and *fend* will be measured in Hz.
%
%     'phi',phi  Starting phase of the chirp. Default value is 0.
%
%   See also: pchirp
    
% AUTHORS:  Piotr Majdak, Peter L. Søndergaard.

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.phi=0;
definput.keyvals.fs=[];

[flags,kv]=ltfatarghelper({},definput,varargin);

if ~isempty(kv.fs)
  fstart=fstart/kv.fs*2;
  fend  =  fend/kv.fs*2;
end;

w1=pi*fstart*L;
w2=pi*fend*L;

A=w1/log(w2/w1);
tau=1/log(w2/w1);

t=((0:L-1)/L).';
outsig=exp(i*A*(exp(t/tau)-1)+kv.phi);

