function ccell = wavpack2cell(cvec,Lc)
%WAVPACK2CELL Changes wavelet coefficients storing format
%   Usage:  
%          ccell = wavpack2cell(cvec,Lc);
%
%   Input parameters:
%         cvec     : Coefficients in packed format.
%         Lc       : Vector containing coefficients lengths.
%
%   Output parameters:
%         ccell    : Coefficients stored in a collumn cell-array.
%
%   *cvec* is collumn vector or matrix with *W* collumns for multi-channel inputs containing
%   coefficients in the packed format. Coefficients are stored as folows:
%   cvec(1:Lc(1),w) - approximation coefficients at level *J* of the channel *w*,
%   cvec(1+sum(Lc(1:j-1)):sum(Lc(1:j),w) for *j>1*.
%

JJ = length(Lc);

[cLen, W] = size(cvec);

% ALLOCATING OUTPUT
ccell = cell(JJ,1);

for jj=1:JJ
     ccell{jj} = zeros(Lc(jj),W);
end

% DO THE COPY

for w=1:W
    lenSumIdx = 1;
    lenSum = 0;
    for jj=1:JJ
       ccell{jj}(:,w) = cvec(1+lenSum:Lc(lenSumIdx)+lenSum,w);
       lenSum = lenSum+Lc(lenSumIdx);
       lenSumIdx=lenSumIdx+1;
    end
end












