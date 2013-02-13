function f=comp_iwpfbt(c,wt,Ls,type,ext)
%COMP_IWFBT Compute Inverse Wavelet Packet Filter-Bank Tree
%   Usage:  f=comp_iwpfbt(c,wt,Ls,type,ext)
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
f = zeros(Ls,chans);

% Nodes in the reverse BF order
treePath = nodesBForder(wt,'rev');
[pOutIdxs,chOutIdxs] = rangeWpBF(wt,'rev');
noOfCoeff = length(c);

% The decimated case
if(strcmp(type,'dec'))
    nodesOutLengths = zeros(length(treePath),1);
    for ii=1:length(treePath)
       nodesOutLengths(ii) =  nodeInLen(treePath(ii),Ls,doNoExt,wt);
    end
    for ch=1:chans
       % exclude root
       for jj=1:length(treePath)-1
           tmpfilt = wt.nodes{treePath(jj)}.filts;
           tmpfiltNo = length(tmpfilt);
           tmpa = wt.nodes{treePath(jj)}.a;

           % HOWTO combine redundant coefficients
           c{pOutIdxs(jj)} = c{pOutIdxs(jj)}/2;
           for ff=1:tmpfiltNo
               ffTmpFilt = tmpfilt{ff}.h;
               ffTmpFiltD = tmpfilt{ff}.d;
               ffTmpFlen = length(ffTmpFilt);
               ffTmpa = tmpa(ff);
               if(doNoExt)
                  tmpSkip = ffTmpFiltD-1;
               else
                  tmpSkip = ffTmpFlen-2;
               end
               c{pOutIdxs(jj)} = c{pOutIdxs(jj)} + 0.5*up_conv_td({c{chOutIdxs{jj}(ff)}(:,ch)},...
                   length(c{pOutIdxs(jj)}),{ffTmpFilt},ffTmpa,tmpSkip,doNoExt,0);
           end
       end
       
       tmpfilt = wt.nodes{treePath(end)}.filts;
       tmpfiltNo = length(tmpfilt);
       tmpa = wt.nodes{treePath(end)}.a;
           %tmpOutLen = nodesOutLengths(jj);
           
       for ff=1:tmpfiltNo
          ffTmpFilt = tmpfilt{ff}.h;
          ffTmpFiltD = tmpfilt{ff}.d;
          ffTmpFlen = length(ffTmpFilt);
          ffTmpa = tmpa(ff);
          if(doNoExt)
            tmpSkip = ffTmpFiltD-1;
          else
            tmpSkip = ffTmpFlen-2;
          end
          f(:,ch) = f(:,ch) + up_conv_td({c{chOutIdxs{end}(ff)}(:,ch)},...
               Ls,{ffTmpFilt},ffTmpa,tmpSkip,doNoExt,0);
       end
       
       
       
    end
elseif(strcmp(type,'undec'))
    error('Not done yet!');
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