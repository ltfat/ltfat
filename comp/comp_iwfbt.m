function f=comp_iwfbt(c,wtreePath,outLens,rangeLoc,rangeOut,Ls,ext)
%COMP_IWFBT Compute Inverse Wavelet Filter-Bank Tree
%   Usage:  f=comp_iwfbt(c,wt,Ls,type,ext)
%
%   Input parameters:
%         c     : Coefficients stored in the cell array.
%         wt    : Structure defining the filterbank tree (see |wfbtinit|_) 
%         type  : 'dec','undec' Type of the wavelet transform.
%         ext   : 'per','zpd','sym','symw','asym','asymw','ppd','sp0' Type of the forward transform boundary handling.
%
%   Output parameters:
%         f     : Reconstructed data.
%

% Do non-expansve transform if ext='per'
doPer = strcmp(ext,'per');

ca = {};
for jj=1:length(wtreePath)
    %Node filters to a cell array
    gCell = cellfun(@(gEl) gEl.h(:),wtreePath{jj}.filts(:),'UniformOutput',0);
    %Node filters subs. factors
    a = wtreePath{jj}.a;
    %Node filters initial skips
    if(doPer)
       skip = cellfun(@(gEl) gEl.d-1,wtreePath{jj}.filts);
    else
       skip = cellfun(@(gEl) numel(gEl),gCell) - 1 - (a - 1);
    end
    filtNo = numel(gCell);
    
    %Prepare input arrays
    catmp = cell(filtNo,1);
    %Read from output
    catmp(rangeLoc{jj}) = c(rangeOut{jj});
    diffRange = setdiff(1:filtNo,rangeLoc{jj});
    %Read from intermediate outputs
    catmp(diffRange(end:-1:1)) = ca(1:numel(diffRange));
    
    %Run filterbank
    catmp = comp_ifilterbank_td(catmp,gCell,a,outLens(jj),skip,ext);
    %Save intermediate output
    ca = [ca(numel(diffRange)+1:end);catmp];
end
f = catmp;
% number of channels
% chans = size(c{end},2);  
% 
%     for ch=1:chans
%        tempca = {};
%        for jj=1:length(wtreePath)
%            tmpfilt = wtreePath{jj}.filts;
%            tmpa = wtreePath{jj}.a;
%            tmpOutLen = outLens(jj);
%            tmpNoOfFilters = length(tmpfilt);
%            
%            % find all unconnected outputs of the node
%            tmpOutRange = rangeLoc{jj}(end:-1:1);
%            tmpCoefOutRange = rangeOut{jj}(end:-1:1);
%            
%            tempOut = zeros(tmpOutLen,1);
%            % first, do the filtering that uses c 
%            for ff=1:length(tmpCoefOutRange)
%                ffTmpFilt = tmpfilt{tmpOutRange(ff)}.h;
%                ffTmpFiltD = tmpfilt{tmpOutRange(ff)}.d;
%                ffTmpFlen = length(ffTmpFilt);
%                ffTmpa = tmpa(tmpOutRange(ff));
%                if(doPer)
%                   tmpSkip = ffTmpFiltD-1;
%                else
%                   tmpSkip = ffTmpFlen-2;
%                end
%                tempOut = tempOut + comp_upconv({c{tmpCoefOutRange(ff)}(:,ch)},...
%                    tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,doPer,0);
%            end
%            
%            % second, do other filtering
%            tmpOtherRange = 1:tmpNoOfFilters;tmpOtherRange(tmpOutRange)=0;tmpOtherRange=tmpOtherRange(tmpOtherRange~=0);
%            tmpOtherRange = tmpOtherRange(end:-1:1);
%            tmpOtherOutNo = length(tmpOtherRange);
%            for ff=1:tmpOtherOutNo
%                ffTmpFilt = tmpfilt{tmpOtherRange(ff)}.h;
%                ffTmpFiltD = tmpfilt{tmpOtherRange(ff)}.d;
%                ffTmpFlen = length(ffTmpFilt);
%                ffTmpa = tmpa(tmpOtherRange(ff));
%                if(doPer)
%                   tmpSkip = ffTmpFiltD-1;
%                else
%                   tmpSkip = ffTmpFlen-2;
%                end
%                tempOut = tempOut + comp_upconv({tempca{ff}},...
%                    tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,doPer,0);
%            end
%         
%            % store the intermediate results
%            tempca = {tempca{tmpOtherOutNo+1:end},tempOut};
%        end
%        f(:,ch) = tempca{end};
%     end
