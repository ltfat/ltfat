function sig = mulawdecode(codedsig,mu,sigweight);
%MULAWDECODE  Inverse of Mu-Law companding
%   Usage:  sig = mulawdecode(codedsig,mu,sigweight);
%
%   MULAWDECODE(codedsig,mu,sigweight) inverts a previously
%   applied mu-law companding to the signal codedsig. The parameters
%   mu and sigweight must match those from the call to MULAWENCODE
%  
%R  jano90

% AUTHOR: Bruno Torresani and Peter Soendergaard

error(nargchk(3,3,nargin));

cst = (1+mu);
sig = cst.^(abs(codedsig));
sig = sign(codedsig) .* (sig-1);
sig = sig * sigweight/mu;
