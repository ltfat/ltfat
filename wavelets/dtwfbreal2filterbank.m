function [g,a,info] = dtwfbreal2filterbank( dualwt, varargin)



complainif_notenoughargs(nargin,1,'DTWFB2FILTERBANK');


dtw = dtwfbinit({'strict',dualwt},varargin{:});


wtPath = 1:numel(dtw.nodes);
wtPath(noOfNodeOutputs(1:numel(dtw.nodes),dtw)==0)=[];
rangeLoc = rangeInLocalOutputs(wtPath,dtw);
rangeOut = rangeInOutputs(wtPath,dtw);
[g,a] = nodesMultid(wtPath,rangeLoc,rangeOut,dtw);

dtw.nodes = dtw.dualnodes;
[g2,a2] = nodesMultid(wtPath,rangeLoc,rangeOut,dtw);

% Align filter offsets so they can be summed
for ii = 1:numel(g)
   % Sanity checks
   assert(g{ii}.offset<=0,sprintf('%s: Invalid wavelet filters.',upper(mfilename)));
   assert(g2{ii}.offset<=0,sprintf('%s: Invalid wavelet filters.',upper(mfilename)));
    
   offdiff = g{ii}.offset-g2{ii}.offset;
   if offdiff>0
       g{ii}.offset = g{ii}.offset - offdiff;
       g{ii}.h = [zeros(offdiff,1);g{ii}.h(:)];
   elseif offdiff<0
       g2{ii}.offset = g2{ii}.offset + offdiff;
       g2{ii}.h = [zeros(-offdiff,1);g2{ii}.h(:)];
   end
   
   lendiff = numel(g{ii}.h) - numel(g2{ii}.h);
   if lendiff~=0
       maxLen = max(cellfun(@(gEl) numel(gEl.h),g));
       g{ii}.h = postpad(g{ii}.h,maxLen);
       g2{ii}.h = postpad(g2{ii}.h,maxLen);
   end
end

% Return the filterbanks only after they have been properly aligned
info.g1 = g;
info.g2 = g2;

gpos = cellfun(@(gEl,g2El) setfield(gEl,'h',(gEl.h+1i*g2El.h)),g,g2,'UniformOutput',0);

%gpos([1,end]) = cellfun(@(gEl) setfield(gEl,'h',conj(gEl.h)),gpos([1,end]),'UniformOutput',0);

%gneg = cellfun(@(gEl,g2El) setfield(gEl,'h',(gEl.h-1i*g2El.h)),g,g2,'UniformOutput',0);
gneg = [];
g = [gpos;gneg(end:-1:1)];
%a = [a;a(end:-1:1)];

