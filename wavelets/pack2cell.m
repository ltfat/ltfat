function ccell = pack2cell(cvec,Lc)
%PACK2CELL Changes wavelet coefficients storing format
%   Usage:  ccell = pack2cell(cvec,lengths);
%
%   Input parameters:
%         cvec     : Coefficients in packed format.
%         lengths  : J+1 element vector containing coefficients lengths.
%
%   Output parameters:
%         ccell    : Coefficients stored in J+1 cell-array.
%
%   *cvec* is collumn vector or matrix with *W* collumns for multi-channel inputs containing
%   coefficients in the packed format. Coefficients are stored as folows:
%   cvec(1:lengths(1),w) - approximation coefficients at level *J* of the channel *w*,
%   cvec(1+sum(lengths(1:j-1)):sum(lengths(1:j),w) for *j=2,...,J+1* - detail coefficients at
%   level *J-2+j*.
%

JJ = length(Lc);

[cLen, W] = size(cvec);

ccell = cell(JJ,1);


for jj=1:JJ
   ccell{jj} = zeros(Lc(jj),W);
end

lenSum = 0;
for jj=1:JJ
    for w=1:W
        ccell{jj}(:,w) = cvec(1+lenSum:Lc(jj)+lenSum,w);
    end
    lenSum = lenSum + Lc(jj);
end