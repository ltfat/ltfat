function c=comp_wfbt(f,wtNodes,rangeLoc,rangeOut,ext)
%COMP_WFBT Compute output of the Wavelet Filterbank Tree
%   Usage:  c=comp_wfbt(f,wt,type,ext);
%
%   Input parameters:
%         f     : Input data.
%         wt    : Structure defining the filterbank tree (see |wfbtinit|_)
%         ext   : 'per','zpd','sym','symw','asym','asymw','ppd','sp0' Type of the forward transform boundary handling.
%
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%

% Do non-expansve transform if ext='per'
doPer = strcmp(ext,'per');
% Pre-allocated output
c = cell(sum(cellfun(@(rEl) numel(rEl),rangeOut)),1);

 ca = {f};
 % Go over all nodes in breadth-first order
 for jj=1:numel(wtNodes)
    %Node filters to a cell array
    hCell = cellfun(@(hEl) hEl.h(:),wtNodes{jj}.filts(:),'UniformOutput',0);
    %Node filters subs. factors
    a = wtNodes{jj}.a;
    %Node filters initial skips
    if(doPer)
       skip = cellfun(@(hEl) hEl.d-1,wtNodes{jj}.filts);
    else
       skip = a-1;
    end

    % Run filterbank
    catmp=comp_filterbank_td(ca{1},hCell,a,skip,ext);
    % Goes directy to the output
    c(rangeOut{jj}) = catmp(rangeLoc{jj});
    % Is saved for next iterations
    ca = [ca(2:end);catmp(setdiff(1:numel(hCell),rangeLoc{jj}))];
 end        
%tmpfilt = wtNodes{jj}.filts;
%tmpa = wtNodes{jj}.a;
%tmpInLen = inLens(jj);
%         tmpOutRange = rangeLoc{jj};
%         tmpNoOfFilters = length(tmpfilt);
%         tmpCoefOutRange = rangeOut{jj};
% 
% 
%         % first, do the filtering that goes directly to c 
%         for ff=1:length(tmpCoefOutRange)
%             ffTmpFilt = tmpfilt{tmpOutRange(ff)}.h;
%             ffTmpFiltD = tmpfilt{tmpOutRange(ff)}.d;
%             ffTmpFlen = length(ffTmpFilt);
%             ffTmpa = tmpa(tmpOutRange(ff));
%             if(doPer)
%                tmpOutLen = ceil(tmpInLen/ffTmpa);
%                tmpSkip = ffTmpFiltD-1;
%             else
%                tmpOutLen = floor((tmpInLen+(ffTmpFlen-1))/ffTmpa); 
%                tmpSkip = 1;
%            end
%            c{tmpCoefOutRange(ff)}(:,ch) =...
%                comp_convsub(ca{1},tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,ext,0);
%         end
% 
%         % store the other outputs
%         tmpOtherRange = 1:tmpNoOfFilters;tmpOtherRange(tmpOutRange)=0;tmpOtherRange=tmpOtherRange(tmpOtherRange~=0);
%         tmpOtherOutNo = length(tmpOtherRange);
%         tmpOut = cell(1,tmpOtherOutNo);
%         for ff=1:tmpOtherOutNo
%             ffTmpFilt = tmpfilt{tmpOtherRange(ff)}.h;
%             ffTmpFiltD = tmpfilt{tmpOtherRange(ff)}.d;
%             ffTmpFlen = length(ffTmpFilt);
%             ffTmpa = tmpa(tmpOtherRange(ff));
%             if(doPer)
%                tmpOutLen = ceil(tmpInLen/ffTmpa);
%                tmpSkip = ffTmpFiltD-1;
%             else
%                tmpOutLen = floor((tmpInLen+(ffTmpFlen-1))/ffTmpa); 
%                tmpSkip = 1;
%             end
%             tmpOut{ff} =...
%                comp_convsub(ca{1},tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,ext,0);
%         end
% 
%         ca = {ca{2:end},tmpOut{:}};
%    end




