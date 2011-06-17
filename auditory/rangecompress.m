function [outsig, sigweight] = rangecompress(insig,varargin)
%RANGECOMPRESS   Compress the dynamic range of a signal 
%   Usage: [outsig, sigweight] = rangecompress(insig,mu);
%   
%   [outsig, sigweight]=RANGECOMPRESS(insig,mu) mu-law rangecompresss the input
%   signal insig using mu-law rangecompressing with parameters mu.
%
%R  jano90

% AUTHOR: Bruno Torresani and Peter L. Soendergaard

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.flags.method={'mulaw','alaw'}
definput.keyvals.mu=255;
[flags,kv]=ltfatarghelper({'L'},definput,varargin);

if flags.do_mulaw
  tmp = log(1+kv.mu);
  
  sigweight = max(abs(insig(:)));
  insig = insig/sigweight;
  
  outsig = sign(insig) .* log(1+kv.mu*abs(insig))/tmp;

end;

if flags.do_alaw
  error('Not implemented yet.');
end;
