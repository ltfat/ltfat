function c = comp_ufwt(f,h,J,a)
%COMP_UFWT Compute undecimated FWT
%   Usage:  c=comp_ufwt(f,h,J);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet filters.
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c     : Coefficients stored in a matrix.
%

% This could be removed with some effort. The question is, are there such
% wavelet filters? If your filterbank has different subsampling factors after first two filters, please send a feature request.
assert(a(1)==a(2),'First two elements of a are not equal. Such wavelet filterbank is not suported.');

filtNo = length(h);
[inLen, W] = size(f);
% Allocate output
c = zeros(inLen,J*(filtNo-1)+1,W);

% For holding the impulse responses.
tmph = cell(filtNo,1);
for ff=1:filtNo
  tmph{ff} =  h{ff}.h/sqrt(a(ff)); 
end

sub = 1;
for ch=1:W
  tempca = f(:,ch);  
  runPtr = 0;
  for jj=1:J
     for ff=filtNo:-1:2
        skip = ceil(a(1)^(jj-1)*(h{ff}.d - 1));
        c(:,end-runPtr,ch) = comp_convsub(tempca,inLen,{tmph{ff}},sub,skip,'per',a(1)^(jj-1)); 
        runPtr = runPtr + 1;
     end
     skip = ceil(a(1)^(jj-1)*(h{1}.d - 1));
     tempca = comp_convsub(tempca,inLen,{tmph{1}},sub,skip,'per',a(1)^(jj-1));
  end
  c(:,1,ch) = tempca;
end 


