function c = fwt(f,h,J,varargin)
%FWT   Fast Wavelet Transform 
%   Usage:  c = fdwt(f,h,J)
%
%   Input parameters:
%         f     : Input data.
%         h     : Wavelet filters.
%         J     : Number of filterbank iterations.
%   Output parameters:
%         c     : Coefficients stored in J+1 cell-array.
%
%
%   See also:
%
%   Demos:
%
%   References:

%   AUTHOR : Zdenek Prusa.
%   TESTING: TEST_FWT
%   REFERENCE: REF_FWT

[sigLen,W] = size(f);
flen = length(h{1});
if(sigLen==1)
    sigLen=W;
    W = 1;
end



definput.flags.type = {'dec','undec'};
[flags,kv]=ltfatarghelper({},definput,varargin);

if(J<1)
   error('%s: J must be a positive integer.',upper(callfun)); 
end

% always even lengths
apprLen = @(j) ceil(ceil(sigLen/2^min([j,J]))/2)*2;
coefLen = @(j) ceil(sigLen/2^min([j,J]));
upfLen = @(j) max([flen, 2^j*flen-(2^j-1)]);


hpad = cell(J,length(h));
c = cell(J+1,1);

if(flags.do_dec)
    maxeqLen = (2^J-1)*(length(h{1})-1);
    if(maxeqLen>sigLen)
        % TODO: circle-add filter when it is longer then the input signal
        error('Not implemented yet');
    else
        for hidx=1:length(h)
            for j=0:J-1
               hpad{j+1,hidx}=zeros(apprLen(j),1);
               hpad{j+1,hidx}(1:flen) = h{hidx}(:);
               % circshift to be consistent with Matlab Wavelet Toolbox
               hpad{j+1,hidx} = circshift(hpad{j+1,hidx},-ceil(flen/2)+1);
           end
        end
    end

      for j=J+1:-1:1
          c{end-j+1} = zeros(coefLen(j),W);
      end


    for w=1:W    
        tempa = f(:,w);
        for j=1:J
            if(mod(length(tempa),2)==1)
                tempa = [tempa; tempa(end)];
            end
            c{end+1-j}(:,w) = downs(pconv(tempa,hpad{j,2}));
            tempa = downs(pconv(tempa,hpad{j,1}));
        end
        c{1}(:,w) = tempa;
    end
elseif(flags.do_undec)
    maxeqLen = upfLen(J);
    if(maxeqLen>sigLen)
        % TODO: circle-add filter when it is longer then the input signal
        error('Not implemented yet');
    else

        for hidx=1:length(h)
            for j=0:J-1
               flenTemp = upfLen(j);
               hpad{j+1,hidx}=zeros(apprLen(0),1);
               %flenTemp = flenTemp +1;
%                if mod(flenTemp,2)==1
%                    range = 2:flenTemp+1;
%                    flenTemp = flenTemp +1;
%                else
                   range = 1:flenTemp;
            %   end
               hpad{j+1,hidx}(range) = ups(h{hidx}(:),2^j,1);

              % circshift to be consistent with Matlab Wavelet Toolbox
            %  if(mod(flenTemp,2)==1)

             hpad{j+1,hidx} = circshift(hpad{j+1,hidx},-floor(flenTemp/2));

        end
       end
    end  
    
    for w=1:W    
        tempa = f(:,w);
        for j=1:J
            c{end+1-j}(:,w) = pconv(tempa,hpad{j,2});
            tempa = pconv(tempa,hpad{j,1});
        end
        c{1}(:,w) = tempa;
    end 
 
    
  
end