function f = comp_iufwt(c,g,J,a)
%COMP_IFWT_ALL Compute Inverse DWT
%   Usage:  f = comp_ifwt_all(c,g,J,Ls,type,ext);
%
%   Input parameters:
%         c     : L*M*W array of coefficients, M=J*(filtNo-1)+1.
%         g     : Synthesis wavelet filters-Cell-array of length *filtNo*.
%         J     : Number of filterbank iterations.
%         a     : Upsampling factors - array of length *filtNo*.
%
%   Output parameters:
%         f     : Reconstructed data - L*W array.
%
% 

% see comp_ufwt for explanantion
assert(a(1)==a(2),'First two elements of a are not equal. Such wavelet filterbank is not suported.');

% Output length
L = size(c,1);

% For holding the impulse responses.
filtNo = length(g);
gDel = cellfun(@(gEl) gEl.d,g(:));
%Change format to a matrix
gMat = cell2mat(cellfun(@(gEl) gEl.h(:),g(:)','UniformOutput',0));
%Divide each column (filter) by a
gMat = bsxfun(@rdivide,gMat,sqrt(a(:)'));

% Read top-level appr. coefficients.
ca = squeeze(c(:,1,:));
cRunPtr = 2;
for jj=1:J
   % Upsampling the filters.
   filtUps = a(1)^(J-jj); 
   %gMatUps = comp_ups(gMat,filtUps,1);
   % Zero index position of the upsampled filetrs.
   skip = filtUps.*gDel - filtUps; 
   % Run the filterbank
   ca=comp_iatrousfilterbank_td([reshape(ca,size(ca,1),1,size(ca,2)),c(:,cRunPtr:cRunPtr+filtNo-2,:)],gMat,filtUps,skip); 
   % Bookkeeping
   cRunPtr = cRunPtr + filtNo -1;
end
% Copy to the output.
f = ca;


% for ch=1:W
%   tempca = c(:,1,ch);
%   runPtr = 2;
%   for jj=1:J
%      filtUps = a(1)^(J-jj); 
%      skip = floor(filtUps*g{1}.d) - filtUps; 
%      tempca = comp_upconv({tempca},Ls,{tmpg{1}},upFac,skip,'per',filtUps); 
%      for ff=2:filtNo
%          skip = floor(filtUps*g{ff}.d) - filtUps; 
%          tempca = tempca + comp_upconv({c(:,runPtr,ch)},Ls,{tmpg{ff}},upFac,skip,'per',filtUps);
%          runPtr = runPtr + 1;
%      end
%   end
%   f(:,ch) = tempca;
% end

    
    