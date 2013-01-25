function c=comp_wfbt(f,wt,type,ext)
%COMP_WFBT Compute output of the Wavelet Filterbank Tree
%   Usage:  c=comp_wfbt(f,wt,type,ext);
%
%   Input parameters:
%         f     : Input data.
%         wt    : Structure defining the filterbank tree (see |wfbtinit|_)
%         type  : 'dec','undec' Type of the wavelet transform.
%         ext   : 'per','zpd','sym','symw','asym','asymw','ppd','sp0' Type of the forward transform boundary handling.
%
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%

% Do non-expansve transform if ext='per'
if(strcmp(ext,'per'))
    doNoExt = 1;
    ext='perdec';
else
    doNoExt = 0;
end

[inLen, chans] = size(f);
% Nodes in a BF order
treePath = nodesBForder(wt);
nodesInLengths = zeros(length(treePath),1);

% for storing lengths of the output coefficient vectors
noOfCoeff = noOfOutputs(wt);
c = cell(noOfCoeff,1);

% The decimated case
if(strcmp(type,'dec'))
    % initialize input data lengths of all nodes to be processed
    for ii=1:length(treePath)
       nodesInLengths(ii) =  nodeInLen(treePath(ii),inLen,doNoExt,wt);
    end
    
    for ch=1:chans
       tempca = {f(:,ch)};
       for jj=1:length(treePath)
           tmpfilt = wt.nodes{treePath(jj)}.h;
           if(isempty(tmpfilt))
               error('%s: Analysis filters not defined.',mfilename);
           end
           tmpa = wt.nodes{treePath(jj)}.a;
           tmpInLen = nodesInLengths(jj);
           tmpCoefOutRange = []; 
           tmpOutRange = rangeInLocalOutputs(treePath(jj),wt);
           tmpNoOfFilters = length(tmpfilt);
           if(~isempty(tmpOutRange))
               tmpCoefOutRange = rangeInOutputs(treePath(jj),wt);
           end
           
           % first, do the filtering that goes directly to c 
           for ff=1:length(tmpCoefOutRange)
               ffTmpFilt = tmpfilt{tmpOutRange(ff)};
               ffTmpFlen = length(ffTmpFilt);
               ffTmpa = tmpa(tmpOutRange(ff));
               if(doNoExt)
                  tmpOutLen = ceil(tmpInLen/ffTmpa);
                  tmpSkip = floor(ffTmpFlen/2);
               else
                  tmpOutLen = floor((tmpInLen+(ffTmpFlen-1))/ffTmpa); 
                  tmpSkip = 1;
              end
              c{tmpCoefOutRange(ff)}(:,ch) =...
                  conv_td_sub(tempca{1},tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,ext,0);
           end
           
           % store the other outputs
           tmpOtherRange = 1:tmpNoOfFilters;tmpOtherRange(tmpOutRange)=0;tmpOtherRange=tmpOtherRange(tmpOtherRange~=0);
           tmpOtherOutNo = length(tmpOtherRange);
           tmpOut = cell(1,tmpOtherOutNo);
           for ff=1:tmpOtherOutNo
               ffTmpFilt = tmpfilt{tmpOtherRange(ff)};
               ffTmpFlen = length(ffTmpFilt);
               ffTmpa = tmpa(tmpOtherRange(ff));
               if(doNoExt)
                  tmpOutLen = ceil(tmpInLen/ffTmpa);
                  tmpSkip = floor(ffTmpFlen/2);
               else
                  tmpOutLen = floor((tmpInLen+(ffTmpFlen-1))/ffTmpa); 
                  tmpSkip = 1;
               end
               tmpOut{ff} =...
                  conv_td_sub(tempca{1},tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,ext,0);
           end

           tempca = {tempca{2:end},tmpOut{:}};
       end
    end
elseif(strcmp(type,'undec'))
    
end