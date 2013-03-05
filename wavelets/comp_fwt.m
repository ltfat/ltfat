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
assert(a(1)==a(2),'First two elements of *a* are not equal. Such wavelet filterbank is not suported.');

filtNo = length(h);

% Impulse responses.
filtLen = length(h{1}.h);
hMat = zeros(filtLen,filtNo);
skip = zeros(filtNo,1);
for ff=1:filtNo
  hMat(:,ff) =  h{ff}.h(:);
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
     % elsewhere.
  end
end

ccell = cell(numel(Lc),1);
%Non-uniform filterbank
runPtr = numel(Lc)-filtNo+2;
ctmp = f;
for jj=1:J
    % Run filterbank
    ctmp = comp_filterbank_td(ctmp,hMat,a,skip,ext);
    % Bookkeeping
    [ccell{runPtr:runPtr+filtNo-2}] = ctmp{2:end};
    ctmp = ctmp{1};
    runPtr = runPtr - (filtNo - 1);
end
% Save final approximation coefficients
ccell{1} = ctmp;
c = wavcell2pack(ccell);



       



