function [cvec,Lc] = wavcell2pack(ccell)
%WAVCELL2PACK Changes wavelet coefficients storing format
%   Usage:  [cvec,Lc] = wavcell2pack(ccell);
%
%   Input parameters:
%         ccell    : Coefficients stored in a collumn cell-array.
%
%   Output parameters:
%         cvec     : Coefficients in packed format.
%         Lc       : Vector containing coefficients lengths.
%
%   *cvec* is collumn vector or matrix with *W* collumns for multi-channel inputs containing
%   coefficients in the packed format. Coefficients are stored as folows:
%   cvec(1:Lc(1),w) - approximation coefficients at level *J* of the channel *w*,
%   cvec(1+sum(Lc(1:j-1)):sum(Lc(1:j),w) for *j>1*.
%

JJtotal = length(ccell);
[cLen, W] = size(ccell{end});


Lc = zeros(JJtotal,1);
for jj=1:JJtotal
   Lc(jj) =  length(ccell{jj});
end


cvec = zeros(sum(Lc),W);

 for w=1:W
   lenSumIdx = 1;
   lenSum = 0;
   for jj=1:JJtotal
      cvec(1+lenSum:Lc(lenSumIdx)+lenSum,w) = ccell{jj}(:,w);
      lenSum = lenSum+Lc(lenSumIdx);
      lenSumIdx=lenSumIdx+1;
   end
 end

