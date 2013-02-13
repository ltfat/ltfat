function c=comp_wpfbt(f,wt,type,ext)
%COMP_WPFBT Compute Wavelet Packet Coefficients given Filterbank Tree
%   Usage:  c=comp_wpfbt(f,wt,type,ext);
%
%   Input parameters:
%         f     : Input data.
%         wt    : Structure defining the filterbank tree (see |wfbtinit|_)
%         type  : 'dec','undec' Type of the wavelet transform.
%         ext   : 'per','zpd','sym','symw','asym','asymw','ppd','sp0' Type of the forward transform boundary handling.
%
%   Output parameters:
%         c     : Coefficients stored in cell-array.
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
% Number of nodes to go trough
nodesNo = length(treePath);
nodesInLengths = zeros(nodesNo,1);

% for storing lengths of the output coefficient vectors
noOfCoeff = 0;
for ii=1:nodesNo
    noOfCoeff = noOfCoeff + length(wt.nodes{treePath(ii)}.filts);
end
c = cell(noOfCoeff,1);

% The decimated case
if(strcmp(type,'dec'))
    % initialize input data lengths of all nodes to be processed
    for ii=1:nodesNo
       nodesInLengths(ii) =  nodeInLen(treePath(ii),inLen,doNoExt,wt);
    end
    
    for ch=1:chans
       tempa = f(:,ch);
       cOutRunIdx = 1;
       cInRunIdxs = [1];
       for jj=1:length(treePath)
           tmpfilt = wt.nodes{treePath(jj)}.filts;
           tmpfiltNo = length(tmpfilt);
           tmpa = wt.nodes{treePath(jj)}.a;
           tmpInLen = nodesInLengths(jj);
           
           for ff=1:tmpfiltNo
               ffTmpFilt = tmpfilt{ff}.h;
               ffTmpFiltD = tmpfilt{ff}.d;
               ffTmpFlen = length(ffTmpFilt);
               ffTmpa = tmpa(ff);
               if(doNoExt)
                  tmpOutLen = ceil(tmpInLen/ffTmpa);
                  tmpSkip = ffTmpFiltD-1;
               else
                  tmpOutLen = floor((tmpInLen+(ffTmpFlen-1))/ffTmpa); 
                  tmpSkip = 1;
               end
               c{cOutRunIdx+ff-1}(:,ch) =...
                  conv_td_sub(tempa,tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,ext,0);
           end
           
           cInRunIdxs = [cInRunIdxs(2:end),cOutRunIdx:cOutRunIdx+tmpfiltNo-1 ];
           tempa = c{cInRunIdxs(1)};
           cOutRunIdx = cOutRunIdx + tmpfiltNo;

       end
    end
elseif(strcmp(type,'undec'))
    error('Not done yet!');
    % initialize input data lengths of all nodes to be processed
    if(doNoExt)
       for ii=1:length(treePath)
          nodesInLengths(ii) =  inLen;
       end
    else
        error('NOT supported yet.');
    end
    
   for ch=1:chans
       tempca = {f(:,ch)};
       for jj=1:length(treePath)
           tmpfilt = wt.nodes{treePath(jj)}.filts;
           tmpa = wt.nodes{treePath(jj)}.a;
           for ii = 1:numel(tmpfilt)
              tmpfilt{ii}.h = tmpfilt{ii}.h/sqrt(tmpa(ii));
           end
           tmpInLen = nodesInLengths(jj);
           tmpCoefOutRange = []; 
           tmpOutRange = rangeInLocalOutputs(treePath(jj),wt);
           tmpNoOfFilters = length(tmpfilt);
           if(~isempty(tmpOutRange))
               tmpCoefOutRange = rangeInOutputs(treePath(jj),wt);
           end
           
           % first, do the filtering that goes directly to c 
           for ff=1:length(tmpCoefOutRange)
               ffTmpFilt = tmpfilt{tmpOutRange(ff)}.h;
               ffTmpFiltD = tmpfilt{tmpOutRange(ff)}.d;
               ffTmpFlen = length(ffTmpFilt);
               %ffTmpa = tmpa(tmpOutRange(ff));
               ffTmpUpFac = nodeFiltUps(treePath(jj),wt);
               if(doNoExt)
                  tmpOutLen = tmpInLen;
                  tmpSkip = ceil(ffTmpUpFac*(ffTmpFiltD-1));
               else
                  error('NOT done yet');
                  %tmpOutLen = floor((tmpInLen+(ffTmpFlen-1))/ffTmpa); 
                  %tmpSkip = 1;
              end
              c{tmpCoefOutRange(ff)}(:,ch) =...
                  conv_td_sub(tempca{1},tmpOutLen,{ffTmpFilt},1,tmpSkip,ext,ffTmpUpFac);
           end
           
           % store the other outputs
           tmpOtherRange = 1:tmpNoOfFilters;tmpOtherRange(tmpOutRange)=0;tmpOtherRange=tmpOtherRange(tmpOtherRange~=0);
           tmpOtherOutNo = length(tmpOtherRange);
           tmpOut = cell(1,tmpOtherOutNo);
           for ff=1:tmpOtherOutNo
               ffTmpFilt = tmpfilt{tmpOtherRange(ff)}.h;
               ffTmpFiltD = tmpfilt{tmpOtherRange(ff)}.d;
               ffTmpFlen = length(ffTmpFilt);
               %ffTmpa = tmpa(tmpOtherRange(ff));
               ffTmpUpFac = nodeFiltUps(treePath(jj),wt);
               if(doNoExt)
                  tmpOutLen = tmpInLen;
                  tmpSkip = ceil(ffTmpUpFac*(ffTmpFiltD-1));
               else
                  error('NOT done yet');
               end
               tmpOut{ff} =...
                  conv_td_sub(tempca{1},tmpOutLen,{ffTmpFilt},1,tmpSkip,ext,ffTmpUpFac);
           end

           tempca = {tempca{2:end},tmpOut{:}};
       end
    end 

    
    
end


