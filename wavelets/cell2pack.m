function [cvec,Lc] = cell2pack(ccell)
%CELL2PACK Changes wavelet coefficients storing format
%   Usage:  [cvec,Lc] = cell2pack(ccell);
%
%   Input parameters:
%         ccell    : Coefficients stored in J+1 cell-array.
%
%   Output parameters:
%         cvec     : Coefficients in packed format.
%         lengths  : J+1 element vector containing coefficients lengths.
%
%   *cvec* is collumn vector or matrix with *W* collumns for multi-channel inputs containing
%   coefficients in the packed format. Coefficients are stored as folows:
%   cvec(1:lengths(1),w) - approximation coefficients at level *J* of the channel *w*,
%   cvec(1+sum(lengths(1:j-1)):sum(lengths(1:j),w) for *j=2,...,J+1* - detail coefficients at
%   level *J-2+j*.
%

JJ = length(ccell);

[cLen, W] = size(ccell{end});


Lc = zeros(JJ,1);
for jj=1:JJ
   Lc(jj) =  length(ccell{jj});
end

cvec = zeros(sum(Lc),W);

lenSum = 0;
for jj=1:JJ
    for w=1:W
        cvec(1+lenSum:Lc(jj)+lenSum,w) = ccell{jj}(:,w);
    end
    lenSum = lenSum + Lc(jj);
end