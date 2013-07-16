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
%   This function transforms coefficients stored as elements of the cell-array *ccell*
%   columns/matrices to a single vector/matrix.
%   *cvec* is vector or matrix containing coefficients in the packed format.
%   For $dim=1$, the coefficients are stored as columns::
%
%     cvec(1:Lc(1),w)
%
%   is the approximation coefficients at level *J* of the channel *w*, ::
%
%     cvec(1+sum(Lc(1:j-1)):sum(Lc(1:j),w)
%
%   for *j>1*. For *dim=2* as rows in the similar manner.
%

if(nargin<1)
    error('%s: Too few input parameters.',upper(mfilename));
end

definput.keyvals.dim = 1;
[flags,kv,dim]=ltfatarghelper({'dim'},definput,varargin);
if(dim>2)
    error('%s: Multidimensional data is not accepted.',upper(mfilename));
end

% Actual computation
Lc = cellfun(@(x) size(x,1), ccell);
cvec = cell2mat(ccell);

% Reshape back to rows
if(dim==2)
    cvec = cvec.';
end



