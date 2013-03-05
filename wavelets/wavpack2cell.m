function [ccell,dim] = wavpack2cell(cvec,Lc,varargin)
%WAVPACK2CELL Changes wavelet coefficients storing format
%   Usage:  
%          ccell = wavpack2cell(cvec,Lc);
%          ccell = wavpack2cell(cvec,Lc,dim);
%
%   Input parameters:
%         cvec     : Coefficients in packed format.
%         Lc       : Vector containing coefficients lengths.
%         dim      : Identifies dimension along which the data were transformed. 
%
%   Output parameters:
%         ccell    : Coefficients stored in a cell-array. Each element is
%         collum or matrix.
%         dim      : Return used dim. Usefull as an input of the
%         complementary function |wavcell2pack|_.
%
%   This function transformes coefficients in vector/matrix *cvec* to the
%   cell array *ccell* of collumn vectors/matrices.
%   *cvec* is vector or matrix containing coefficients in the packed format.
%   For dim==1, the coefficients are expected to be stored as collumns:
%   cvec(1:Lc(1),w) - approximation coefficients at level *J* of the channel *w*,
%   cvec(1+sum(Lc(1:j-1)):sum(Lc(1:j),w) for *j>1*, for dim==2 as rows in
%   the similar manner. 
%   The resulting *ccell* has `length(Lc)` elements.


if(nargin<2)
    error('%s: Too few input parameters.',upper(mfilename));
end

definput.keyvals.dim = [];
[flags,kv,dim]=ltfatarghelper({'dim'},definput,varargin);

%If dim is not specified use first non-singleton dimension.
if(isempty(dim))
    dim=find(size(cvec)>1,1);
end

if(dim>2)
    error('%s: Multidimensional data is not accepted.',upper(mfilename));
end

if(dim==2)
    cvec = cvec.';
end

if(sum(Lc)~=size(cvec,1))
    error('%s: Sum of elements of Lc is not equal to vector length along dimension %d. Possibly wrong dim?',upper(mfilename),dim);
end

% Actual computaion
ccell = mat2cell(cvec,Lc);


% JJ = length(Lc);
% % ALLOCATING OUTPUT
% ccell = cell(JJ,1);
% % DO THE COPY
% LcEnd = cumsum(Lc); 
% LcStart = 1 + cumsum([0;Lc(1:end-1)]); 
% for jj=1:JJ
%   ccell{jj} = cvec(LcStart(jj):LcEnd(jj),:);
% end















