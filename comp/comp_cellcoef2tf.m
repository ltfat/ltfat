function coef = comp_cellcoef2tf(coef,maxLen)
%COMP_CELLCOEF2TF Cell to a tf-layout
%   Usage: coef = comp_cellcoef2tf(coef,maxLen)
%




coefLenMax = max(cellfun(@(cEl)size(cEl,1),coef));

if nargin>1
   coefLenMax = min([coefLenMax,maxLen]);
end

coefTmp = zeros(coefLenMax,numel(coef),size(coef{1},2),class(coef{1}));

for ii=1:numel(coef)
   if size(coef{ii},1) == 1
      coefTmp(:,ii) = coef{ii};
      continue;
   end
   if ~isoctave
       coefTmp(:,ii) = interp1(coef{ii},linspace(1,size(coef{ii},1),...
                       coefLenMax),'nearest');
   else
       coefRe = interp1(real(coef{ii}),linspace(1,size(coef{ii},1),...
                       coefLenMax),'nearest');
                   
        coefIm = interp1(imag(coef{ii}),linspace(1,size(coef{ii},1),...
                       coefLenMax),'nearest');  
                   
        coefTmp(:,ii) = coefRe + 1i*coefIm;
   end
end
coef = coefTmp';
