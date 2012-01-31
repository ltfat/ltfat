function out=ref_bincoeff(u,v)
%REF_BINCOEFF  Binomial coefficients, possibly rational
%
%  Compted by lambda functions.
%
%  See formula 1.2 in Unsers paper: "Fractional splines and wavelets"

out=gamma(u+1)./(gamma(v+1).*gamma(u-v+1));




