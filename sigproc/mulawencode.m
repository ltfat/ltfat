function [outsig, sigweight] = mulawencode(insig,mu)
%MULAWENCODE   Mu-Law compand signal 
%   Usage: [outsig, sigweight] = mulawencode(insig,mu);
%   
%   [outsig, sigweight]=MULAWENCODE(insig,mu) mu-law compands the input
%   signal insig using mu-law companding with parameters mu.
%
%R  jano90

% AUTHOR: Bruno Torresani

error(nargchk(2,2,nargin));

tmp = log(1+mu);

sigweight = max(abs(insig(:)));
insig = insig/sigweight;

outsig = sign(insig) .* log(1+mu*abs(insig))/tmp;

