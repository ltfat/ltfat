function [ccell,dim] = wavpack2cell(cvec,Lc,varargin)
%WAVPACK2CELL Changes wavelet coefficients storing format
%   Usage:  
%          ccell = wavpack2cell(cvec,Lc);
%          ccell = wavpack2cell(cvec,Lc,dim);
%
%   Input parameters:
%         cvec     : Coefficients in packed format.
%         Lc       : Vector containing coefficients lengths.
%         dim      : Dimension along which the data were transformed. 
%
%   Output parameters:
%         ccell    : Coefficients stored in a cell-array. Each element is
%                    a column vector or a matrix.
%         dim      : Return used dim. Usefull as an input of the
%                    complementary function |wavcell2pack|.
%
%   `ccell = wavpack2cell(cvec,Lc)` copies coefficients from a single column
%   vector or columns of a matrix *cvec* of size `[sum(Lc), W]` to the cell
%   array *ccell* of length `length(Lc)`. Size of *j*-th element of *ccell*
%   is `[Lc(j), W]` and it is obtained by::
% 
%      ccell{j}=cvec(1+sum(Lc(1:j-1)):sum(Lc(1:j),:);
%
%   `ccell = wavpack2cell(cvec,Lc,dim)` allows specifying along which
%   dimension the coefficients are stored in *cvec*. *dim==1* (default)
%   considers columns (as above) and *dim==2* rows to be coefficients 
%   belonging to separate channels. Other values are not supported. For 
%   *dim=2*, *cvec* size is `[W, sum(Lc)]`, Size of *j*-th element of *ccell*
%   is `[Lc(j), W]` and it is obtained by::
% 
%      ccell{j}=cvec(:,1+sum(Lc(1:j-1)):sum(Lc(1:j)).';
%
%   See also: wavcell2pack, fwt, wfbt, wpfbt

% AUTHOR: Zdenek Prusa


if(nargin<2)
    error('%s: Too few input parameters.',upper(mfilename));
end

if(~isnumeric(cvec))
    error('%s: *cvec* is not a numeric array.',upper(mfilename));
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

% Actual computation
ccell = mat2cell(cvec,Lc);
