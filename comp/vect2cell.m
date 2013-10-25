function c = vect2cell(x,idx)

idxEnd = cumsum(idx(:));
idxStart = [1;1+idxEnd(1:end-1)];
c = arrayfun(@(idS,idE) x(idS:idE,:),idxStart,idxEnd,'UniformOutput',0);

