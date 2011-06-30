function outsig = rangeexpand(insig,varargin);
%RANGEEXPAND  Expand the dynamic range of a signal
%   Usage:  sig = rangeexpand(insig,mu,sigweight);
%
%   RANGEEXPAND(insig,mu,sigweight) inverts a previously
%   applied mu-law companding to the signal insig. The parameters
%   mu and sigweight must match those from the call to MULAWENCODE
%
%   RANGEEXPAND takes the following optional arguments:
%
%-     'mulaw'  - Do mu-law compression, this is the default.
%
%-     'alaw'   - Do A-law compression.
%
%-     'mu',mu  - mu-law parameter. Default value is 255.
%  
%R  jano90

% AUTHOR: Bruno Torresani and Peter Soendergaard

definput.flags.method={'mulaw','alaw'}
definput.keyvals.mu=255;
[flags,kv]=ltfatarghelper({'L'},definput,varargin);

if flags.do_mulaw

  cst = (1+kv.mu);
  outsig = cst.^(abs(insig));
  outsig = sign(insig) .* (sig-1);
  outsig = sig * sigweight/kv.mu;

end;

if flags.do_alaw
  error('Not implemented yet.');
end;

