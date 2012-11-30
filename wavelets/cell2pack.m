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

[JJ bands] = size(ccell);
JJtotal = (JJ-1)*bands + 1;

[cLen, W] = size(ccell{end});


Lc = zeros(JJtotal,1);
Lc(1) =  length(ccell{1});
jjIdx = 2;
for jj=2:JJ
    for bb=1:bands
       Lc(jjIdx) =  length(ccell{jj,bb}(:,1));
       jjIdx = jjIdx +1;
    end
end


cvec = zeros(sum(Lc),W);

for w=1:W
    cvec(1:Lc(1),w) = ccell{1}(:,w);
end


for w=1:W
 lenSumIdx = 1;
 lenSum = Lc(lenSumIdx);

 for jj=2:JJ
    for bb=1:bands 
        cvec(1+lenSum:Lc(lenSumIdx+1)+lenSum,w) = ccell{jj,bb}(:,w);
        lenSum = lenSum+Lc(lenSumIdx+1);
        lenSumIdx=lenSumIdx+1;
    end
  end
end