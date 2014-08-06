function c = vect2cell(x,idx)
%VECT2CELL Vector to cell
%
%   Works exactly like mat2cell(x,idx,size(x,2))
%   but it is faster.

if sum(idx)~=size(x,1)
    error('%s: Sizes do not comply.',upper(mfilename));
end

idxEnd = cumsum(idx(:));
idxStart = [1;1+idxEnd(1:end-1)];
c = arrayfun(@(idS,idE) x(idS:idE,:),idxStart,idxEnd,'UniformOutput',0);
