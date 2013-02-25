function c = comp_fwt(f,h,J,a,Lc,ext)
%COMP_FWT_ALL Compute DWT
%   Usage:  c=comp_fwt_all(f,h,J,type,ext);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet filters.
%         J     : Number of filterbank iterations.
%         ext   : 'per','zero','even','odd' Type of the forward transform boundary handling.
%
%   Output parameters:
%         c     : Coefficients stored in J*length(h)+1 cell-array.
%

% This could be removed with some effort. The question is, are there such
% wavelet filters? If your filterbank has different subsampling factors following the first two filters, please send a feature request.
assert(a(1)==a(2),'First two elements of a are not equal. Such wavelet filterbank is not suported.');


filtNo = length(h);
LcStart = 1 + cumsum([0;Lc(1:end-1)]); 
LcEnd = cumsum(Lc); 

[inLen, W] = size(f);
c = zeros(sum(Lc),W);
%c = cell(length(Lc),1);

% For holding the impulse responses.
tmph = cell(filtNo,1);
% For holding the 
skip = zeros(filtNo,1);
for ff=1:filtNo
  tmph{ff} =  h{ff}.h; 
  if(strcmp(ext,'per'))
     % Initial shift of the filter to compensate for it's delay.
     % "Zero" delay transform is produced
     skip(ff) = h{ff}.d-1;
  else
     % No compensation for the filter delay (filters are made causal with respect to the output sequences).
     % This creates relative shift between levels of coefficients.
     % Initial shift determines type of subsampling. 
     % This is even subsampling. e.g. subs. [1,2,3,4,5,6] by a factor 3 becomes [3,6]
     % The number of output coefficients depends on it.
     skip(ff) = a(ff)-1;
     % For odd subsampling skip(ff) = 0; but it requires slight touches
     % elsewhere
  end
end

for ch=1:W
  tempca = f(:,ch);
  runPtr = 0;
  for jj=1:J
     for ff=filtNo:-1:2
        c(LcStart(end-runPtr):LcEnd(end-runPtr),ch) = comp_convsub(tempca,Lc(end-runPtr),{tmph{ff}},a(ff),skip(ff),ext,0);
        runPtr = runPtr + 1;
     end
     tempca = comp_convsub(tempca,Lc(end-runPtr+1),{tmph{1}},a(1),skip(1),ext,0);
  end
  c(LcStart(1):LcEnd(1),ch) = tempca;
end

       



