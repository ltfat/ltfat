function ccell = pack2cell(cvec,Lc)
%PACK2CELL Changes wavelet coefficients storing format
%   Usage:  
%          ccell = pack2cell(cvec,Lc);
%          ccell = pack2cell(cvec,Lc,bands);
%
%   Input parameters:
%         cvec     : Coefficients in packed format.
%         Lc       : J+1 element vector containing coefficients lengths.
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


% check if (length(Lc)-1)/bands is integer
JJ = length(Lc);

[cLen, W] = size(cvec);

% ALLOCATING OUTPUT
ccell = cell(JJ,1);

for jj=1:JJ
     ccell{jj} = zeros(Lc(jj),W);
end

% ccell{1} = zeros(Lc(1),W);
% jjIdx = 2;
% for jj=2:JJ
%   for bb=1:bands
%      ccell{jj,bb} = zeros(Lc(jjIdx),W);
%      jjIdx = jjIdx +1;
%   end
% end

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


% for w=1:W
%     ccell{1}(:,w) = cvec(1:Lc(1),w);
% end
% 
% for w=1:W
%    lenSumIdx = 1;
%    lenSum = Lc(lenSumIdx);
%     for jj=2:JJ
%       for bb=1:bands
%         ccell{jj,bb}(:,w) = cvec(1+lenSum:Lc(lenSumIdx+1)+lenSum,w);
%         lenSum = lenSum+Lc(lenSumIdx+1);
%         lenSumIdx=lenSumIdx+1;
%       end
%     end
% end










