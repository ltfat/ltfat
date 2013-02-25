function [cvec,Lc] = wavcell2pack(ccell,varargin)
%WAVCELL2PACK Changes wavelet coefficients storing format
%   Usage:  [cvec,Lc] = wavcell2pack(ccell);
%           [cvec,Lc] = wavcell2pack(ccell,dim);
%
%   Input parameters:
%         ccell    : Coefficients stored in a collumn cell-array.
%         dim      : Identifies dimension along which the data were transformed. 
%
%   Output parameters:
%         cvec     : Coefficients in packed format.
%         Lc       : Vector containing coefficients lengths.
%
%   This function transformes coefficients stored as elements of cell-array *ccell*
%   collumns/matrices to a single vector/matrix.
%   *cvec* is vector or matrix containing coefficients in the packed format.
%   For dim==1, the coefficients are stored as collumns:
%   cvec(1:Lc(1),w) - approximation coefficients at level *J* of the channel *w*,
%   cvec(1+sum(Lc(1:j-1)):sum(Lc(1:j),w) for *j>1*. For dim==2 as rows in
%   the similar manner.
%

if(nargin<1)
    error('%s: Too few input parameters.',upper(mfilename));
end

definput.keyvals.dim = 1;
[flags,kv,dim]=ltfatarghelper({'dim'},definput,varargin);
if(dim>2)
    error('%s: Multidimensional data is not accepted.',upper(mfilename));
end


JJ = numel(ccell);
Lc = zeros(JJ,1);
for jj=1:JJ
   Lc(jj) =  size(ccell{jj},1);
end

W = size(ccell{end},2);
% ALLOCATING OUTPUT
cvec = zeros(sum(Lc),W);
% DO THE COPY
LcStart = 1 + cumsum([0;Lc(1:end-1)]); 
LcEnd = cumsum(Lc); 
for jj=1:JJ
  cvec(LcStart(jj):LcEnd(jj),:) = ccell{jj};
end

% Reshape back to rows
if(dim==2)
    cvec = cvec.';
end



