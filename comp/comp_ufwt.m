function c = comp_ufwt(f,h,J,a)
%COMP_UFWT Compute undecimated FWT
%   Usage:  c=comp_ufwt(f,h,J);
%
%   Input parameters:
%         f     : Input data - L*W array.
%         h     : Analysis Wavelet filters - cell-array of length *filtNo*.
%         J     : Number of filterbank iterations.
%         a     : Subsampling factors - array of length *filtNo*.
%
%   Output parameters:
%         c     : L*M*W array of coefficients, where M=J*(filtNo-1)+1.
%


% This could be removed with some effort. The question is, are there such
% wavelet filters? If your filterbank has different subsampling factors after first two filters, please send a feature request.
assert(a(1)==a(2),'First two elements of a are not equal. Such wavelet filterbank is not suported.');

% For holding the impulse responses.
filtNo = length(h);
hDel = cellfun(@(hEl) hEl.d,h(:));
%Change format to a matrix
hMat = cell2mat(cellfun(@(hEl) hEl.h(:),h(:)','UniformOutput',0));
%Divide each column (filter) by a element of a
hMat = bsxfun(@rdivide,hMat,sqrt(a(:)'));

% Allocate output
[L, W] = size(f);
M = J*(filtNo-1)+1;
c = zeros(L,M,W);

ca = f;
runPtr = size(c,2) - (filtNo-2);
for jj=1:J
    % Upsampling the filters.
    hMatUps = comp_ups(hMat,a(1)^(jj-1),1);
    % Zero index position of the upsampled filters.
    skip = ceil(a(1)^(jj-1).*(hDel - 1));
    % Run filterbank.
    ca=comp_ufilterbank_td(ca,hMatUps,1,skip,'per');
    % Bookkeeping
    c(:,runPtr:runPtr+filtNo-2,:)=ca(:,2:end,:);
    ca = squeeze(ca(:,1,:));
    runPtr = runPtr - (filtNo - 1);
end
% Saving final approximation coefficients.
c(:,1,:) = ca;


% sub = 1;
% for ch=1:W
%   tempca = f(:,ch);  
%   runPtr = 0;
%   for jj=1:J
%      for ff=filtNo:-1:2
%         skip = ceil(a(1)^(jj-1)*(h{ff}.d - 1));
%         c(:,end-runPtr,ch) = comp_convsub(tempca,inLen,{tmph{ff}},sub,skip,'per',a(1)^(jj-1)); 
%         runPtr = runPtr + 1;
%      end
%      skip = ceil(a(1)^(jj-1)*(h{1}.d - 1));
%      tempca = comp_convsub(tempca,inLen,{tmph{1}},sub,skip,'per',a(1)^(jj-1));
%   end
%   c(:,1,ch) = tempca;
% end 


