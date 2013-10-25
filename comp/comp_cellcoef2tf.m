function coef = comp_cellcoef2tf(coef,maxLen)
%COMP_CELLCOEF2TF Cell to a tf-layout
%   Usage: coef = comp_cellcoef2tf(coef,dim)
%


coefLenMax = max(cellfun(@(cEl)size(cEl,1),coef));
coefTmp = zeros(coefLenMax,numel(coef),size(coef{1},2),class(coef{1}));

if nargin>1
   coefLenMax = min([coefLenMax,maxLen]);
end

for ii=1:numel(coef)
   if size(coef{ii},1) == 1
      coefTmp(:,ii) = coef{ii};
      continue;
   end
   coefTmp(:,ii) = interp1(coef{ii},linspace(1,size(coef{ii},1),coefLenMax),'nearest');
end
coef = coefTmp';