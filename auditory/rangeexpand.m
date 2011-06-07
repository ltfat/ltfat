function sig = rangeexpand(codedsig,mu,sigweight);
%RANGEEXPAND  Expand the dynamic range of a signal
%   Usage:  sig = rangeexpand(codedsig,mu,sigweight);
%
%   RANGEEXPAND(codedsig,mu,sigweight) inverts a previously
%   applied mu-law companding to the signal codedsig. The parameters
%   mu and sigweight must match those from the call to MULAWENCODE
%  
%R  jano90

% AUTHOR: Bruno Torresani and Peter Soendergaard

definput.flags.method={'mulaw','alaw'}
definput.keyvals.mu=255;
[flags,kv]=ltfatarghelper({'L'},definput,varargin);

if flags.do_mulaw

  cst = (1+kv.mu);
  sig = cst.^(abs(codedsig));
  sig = sign(codedsig) .* (sig-1);
  sig = sig * sigweight/kv.mu;

end;