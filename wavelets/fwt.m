function [c,Ls] = fwt(f,h,J,varargin)
%FWT   Fast Wavelet Transform 
%   Usage:  c = fdwt(f,h,J);
%           c = fdwt(f,h,J,L);
%           c = fdwt(f,h,J,'type');
%           c = fdwt(f,h,J,L,'type');
%           [c,Ls] = fdwt(...);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet filters.
%         J     : Number of filterbank iterations.
%         L     : Length of transform to do.
%         type  : 'dec,'undec' Type of the wavelet transform.
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%         Ls    : Length of input signal.
%
%   See also:
%
%   Demos:
%
%   References:

%   AUTHOR : Zdenek Prusa.
%   TESTING: TEST_FWT
%   REFERENCE: REF_FWT

if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

definput.keyvals.L=[];
definput.flags.type = {'dec','undec'};
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1 && rem(a,1)~=0)
   error('%s: J must be a positive integer.',upper(callfun)); 
end

if(~iscell(h) || length(h)<2)
   error('%s: h is expected to be a cell array containing two wavelet filters.',upper(callfun)); 
end
flen = length(h{1});

%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);

%% ------ step 2: Determine legal length
if isempty(L)
    L=fwtlength(Ls,J);
else
    Luser=fwtlength(L,J);
    if Luser~=L
        error(['%s: Incorrect transform length L=%i specified. Next valid length ' ...
               'is L=%i. See the help of FWTLENGTH for the requirements.'],...
              upper(mfilename),L,Luser);
    end;

end;

%% ----- step 3 : Check whether the 
eqflen = (2^J-1)*(flen-1) + 1;

if L<eqflen
  error('%s: Equivalent filter is too long.',upper(mfilename));
  % Or possibly make the filter circ-all 
end;

%% ----- step 4: final cleanup ---------------

f=postpad(f,L);

if(flags.do_dec)
 c = comp_fwt(f,h,J);
elseif(flags.do_undec)
 c = comp_fwt_undec(f,h,J);    
end

% % always even lengths
% apprLen = @(j) ceil(ceil(Ls/2^min([j,J]))/2)*2;
% coefLen = @(j) ceil(Ls/2^min([j,J]));
% upfLen = @(j) max([flen, 2^j*flen-(2^j-1)]);


% hpad = cell(J,length(h));
% c = cell(J+1,1);
% 
% if(flags.do_dec)
%     maxeqLen = (2^J-1)*(length(h{1})-1);
%     if(maxeqLen>sigLen)
%         % TODO: circle-add filter when it is longer then the input signal
%         error('Not implemented yet');
%     else
%         for hidx=1:length(h)
%             for j=0:J-1
%                hpad{j+1,hidx}=zeros(apprLen(j),1);
%                hpad{j+1,hidx}(1:flen) = h{hidx}(:);
%                % circshift to be consistent with Matlab Wavelet Toolbox
%                hpad{j+1,hidx} = circshift(hpad{j+1,hidx},-ceil(flen/2)+1);
%            end
%         end
%     end
% 
%       for j=J+1:-1:1
%           c{end-j+1} = zeros(coefLen(j),W);
%       end
% 
% 
%     for w=1:W    
%         tempa = f(:,w);
%         for j=1:J
%             if(mod(length(tempa),2)==1)
%                 tempa = [tempa; tempa(end)];
%             end
%             c{end+1-j}(:,w) = downs(pconv(tempa,hpad{j,2}));
%             tempa = downs(pconv(tempa,hpad{j,1}));
%         end
%         c{1}(:,w) = tempa;
%     end
% elseif(flags.do_undec)
%     maxeqLen = upfLen(J);
%     if(maxeqLen>sigLen)
%         % TODO: circle-add filter when it is longer then the input signal
%         error('Not implemented yet');
%     else
% 
%         for hidx=1:length(h)
%             for j=0:J-1
%                flenTemp = upfLen(j);
%                hpad{j+1,hidx}=zeros(apprLen(0),1);
%                %flenTemp = flenTemp +1;
% %                if mod(flenTemp,2)==1
% %                    range = 2:flenTemp+1;
% %                    flenTemp = flenTemp +1;
% %                else
%                    range = 1:flenTemp;
%             %   end
%                hpad{j+1,hidx}(range) = ups(h{hidx}(:),2^j,1);
% 
%               % circshift to be consistent with Matlab Wavelet Toolbox
%             %  if(mod(flenTemp,2)==1)
% 
%              hpad{j+1,hidx} = circshift(hpad{j+1,hidx},-floor(flenTemp/2));
% 
%         end
%        end
%     end  
%     
%     for w=1:W    
%         tempa = f(:,w);
%         for j=1:J
%             c{end+1-j}(:,w) = pconv(tempa,hpad{j,2});
%             tempa = pconv(tempa,hpad{j,1});
%         end
%         c{1}(:,w) = tempa;
%     end 
%  
%     
%   
% end