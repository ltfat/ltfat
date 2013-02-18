function f=comp_iwfbt(c,wtreePath,outLens,rangeLoc,rangeOut,Ls,type,ext)
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


% The decimated case
if(strcmp(type,'dec'))
    for ch=1:chans
       tempca = {};
       for jj=1:length(wtreePath)
           tmpfilt = wtreePath{jj}.filts;
           tmpa = wtreePath{jj}.a;
           tmpOutLen = outLens(jj);
           tmpNoOfFilters = length(tmpfilt);
           
           % find all unconnected outputs of the node
           tmpOutRange = rangeLoc{jj}(end:-1:1);
           tmpCoefOutRange = rangeOut{jj}(end:-1:1);
           
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
               tempOut = tempOut + comp_upconv({c{tmpCoefOutRange(ff)}(:,ch)},...
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
               tempOut = tempOut + comp_upconv({tempca{ff}},...
                   tmpOutLen,{ffTmpFilt},ffTmpa,tmpSkip,doNoExt,0);
           end
        
           % store the intermediate results
           tempca = {tempca{tmpOtherOutNo+1:end},tempOut};
       end
       f(:,ch) = tempca{end};
    end
elseif(strcmp(type,'undec'))
   for ch=1:chans
       tempca = {};
       for jj=1:length(wtreePath)
           tmpfilt = wt.nodes{treePath(jj)}.filts;
           tmpa = wt.nodes{treePath(jj)}.a;
           for ii = 1:numel(tmpfilt)
              tmpfilt{ii}.h = tmpfilt{ii}.h/sqrt(tmpa(ii));
           end
           tmpOutLen = outLens(jj);
           tmpNoOfFilters = length(tmpfilt);
           
           % find all unconnected outputs of the node
           tmpOutRange = rangeLoc{jj}(end:-1:1);
           tmpCoefOutRange = rangeOut{jj}(end:-1:1);
           
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