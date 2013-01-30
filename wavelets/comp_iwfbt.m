function f=comp_iwfbt(c,wt,Ls,type,ext)
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
if(strcmp(ext,'per'))
    doNoExt = 1;
else
    doNoExt = 0;
end

% number of channels
chans = size(c{end},2);  

% Nodes in the reverse BF order
treePath = nodesBForder(wt);
treePath = treePath(end:-1:1);

% 


% The decimated case
if(strcmp(type,'dec'))
    nodesOutLengths = zeros(length(treePath),1);
    for ii=1:length(treePath)
       nodesOutLengths(ii) =  nodeInLen(treePath(ii),Ls,doNoExt,wt);
    end
    for ch=1:chans
       tempca = {};
       for jj=1:length(treePath)
           tmpfilt = wt.nodes{treePath(jj)}.filts;
           if(isempty(tmpfilt))
               error('%s: Synthesis filters not defined.',upper(mfilename));
           end
           tmpa = wt.nodes{treePath(jj)}.a;
           tmpOutLen = nodesOutLengths(jj);
           tmpNoOfFilters = length(tmpfilt);
           
           % find all unconnected outputs of the node
           tmpOutRange = rangeInLocalOutputs(treePath(jj),wt);
           tmpOutRange = tmpOutRange(end:-1:1);
           tmpCoefOutRange = []; 
           if(~isempty(tmpOutRange))
               tmpCoefOutRange = rangeInOutputs(treePath(jj),wt);
               tmpCoefOutRange = tmpCoefOutRange(end:-1:1);
           end
           
           tempOut = zeros(tmpOutLen,1);
           % first, do the filtering that uses c 
           for ff=1:length(tmpCoefOutRange)
               ffTmpFilt = tmpfilt{tmpOutRange(ff)}.h;
               ffTmpFiltD = tmpfilt{tmpOutRange(ff)}.d;
               ffTmpFlen = length(ffTmpFilt);
               ffTmpa = tmpa(tmpOutRange(ff));
               if(doNoExt)
                  tmpSkip = ffTmpFiltD-1;
               else
                  tmpSkip = ffTmpFlen-2;
               end
               tempOut = tempOut + up_conv_td({c{tmpCoefOutRange(ff)}(:,ch)},...
                   tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,doNoExt,0);
           end
           
           % second, do other filtering
           tmpOtherRange = 1:tmpNoOfFilters;tmpOtherRange(tmpOutRange)=0;tmpOtherRange=tmpOtherRange(tmpOtherRange~=0);
           tmpOtherRange = tmpOtherRange(end:-1:1);
           tmpOtherOutNo = length(tmpOtherRange);
           for ff=1:tmpOtherOutNo
               ffTmpFilt = tmpfilt{tmpOtherRange(ff)}.h;
               ffTmpFiltD = tmpfilt{tmpOtherRange(ff)}.d;
               ffTmpFlen = length(ffTmpFilt);
               ffTmpa = tmpa(tmpOtherRange(ff));
               if(doNoExt)
                  tmpSkip = ffTmpFiltD-1;
               else
                  tmpSkip = ffTmpFlen-2;
               end
               tempOut = tempOut + up_conv_td({tempca{ff}},...
                   tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,doNoExt,0);
           end
        
           % store the intermediate results
           tempca = {tempca{tmpOtherOutNo+1:end},tempOut};
       end
       f(:,ch) = tempca{end};
    end
elseif(strcmp(type,'undec'))
   for ii=1:length(treePath)
      nodesOutLengths(ii) =  Ls;
   end
   for ch=1:chans
       tempca = {};
       for jj=1:length(treePath)
           tmpfilt = wt.nodes{treePath(jj)}.filts;
           tmpa = wt.nodes{treePath(jj)}.a;
           for ii = 1:numel(tmpfilt)
              tmpfilt{ii}.h = tmpfilt{ii}.h/sqrt(tmpa(ii));
           end
           tmpOutLen = nodesOutLengths(jj);
           tmpNoOfFilters = length(tmpfilt);
           
           % find all unconnected outputs of the node
           tmpOutRange = rangeInLocalOutputs(treePath(jj),wt);
           tmpOutRange = tmpOutRange(end:-1:1);
           tmpCoefOutRange = []; 
           if(~isempty(tmpOutRange))
               tmpCoefOutRange = rangeInOutputs(treePath(jj),wt);
               tmpCoefOutRange = tmpCoefOutRange(end:-1:1);
           end
           
           tempOut = zeros(tmpOutLen,1);
           % first, do the filtering that uses c 
           for ff=1:length(tmpCoefOutRange)
               ffTmpFilt = tmpfilt{tmpOutRange(ff)}.h;
               ffTmpFiltD = tmpfilt{tmpOutRange(ff)}.d;
               ffTmpFlen = length(ffTmpFilt);
               %ffTmpa = tmpa(tmpOutRange(ff));
               ffTmpUpFac = nodeFiltUps(treePath(jj),wt);
               if(doNoExt)
                  tmpSkip = ffTmpUpFac*(ffTmpFiltD-1);
               else
                  error('Not done yet!');
                  %tmpSkip = ffTmpFlen-2;
               end
               tempOut = tempOut + up_conv_td({c{tmpCoefOutRange(ff)}(:,ch)},...
                   tmpOutLen,{ffTmpFilt},1,tmpSkip,doNoExt,ffTmpUpFac);
           end
           
           % second, do other filtering
           tmpOtherRange = 1:tmpNoOfFilters;tmpOtherRange(tmpOutRange)=0;tmpOtherRange=tmpOtherRange(tmpOtherRange~=0);
           tmpOtherRange = tmpOtherRange(end:-1:1);
           tmpOtherOutNo = length(tmpOtherRange);
           for ff=1:tmpOtherOutNo
               ffTmpFilt = tmpfilt{tmpOtherRange(ff)}.h;
               ffTmpFiltD = tmpfilt{tmpOtherRange(ff)}.d;
               ffTmpFlen = length(ffTmpFilt);
               %ffTmpa = tmpa(tmpOtherRange(ff));
               ffTmpUpFac = nodeFiltUps(treePath(jj),wt);
               if(doNoExt)
                  tmpSkip = ffTmpUpFac*(ffTmpFiltD-1);
               else
                  tmpSkip = ffTmpFlen-2;
               end
               tempOut = tempOut + up_conv_td({tempca{ff}},...
                   tmpOutLen,{ffTmpFilt},1,tmpSkip,doNoExt,ffTmpUpFac);
           end
        
           % store the intermediate results
           tempca = {tempca{tmpOtherOutNo+1:end},tempOut};
       end
       f(:,ch) = tempca{end};
    end
end