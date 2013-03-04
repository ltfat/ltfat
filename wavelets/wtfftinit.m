function [H,a] = wtfftinit(w,octaves,voices,Ls)
%WTFFTINIT
%
%
%

% Possible inputs
% {{'db',10},J}
% 
% {'morl',octaves,scalesPerOctave}
%
% cell array of structs containing FIR filters
%
% cell array of structs containing FREQ filters (pass trough)

assert(octaves>0 && rem(octaves,1)==0 && isscalar(octaves));
assert(voices>=1 && rem(voices,1)==0 && isscalar(voices));



% check
prefixFreq = 'wfreq_';
prefixFunc = 'wfunc_';
prefixFilt = 'wfilt_';
if(~iscell(w))
   w = {w};
end
wname = w{1};

expScales = 2.^((0:voices-1)/voices);
expOctaves = 2.^(1:octaves);
s = expScales.'*expOctaves;
s = s(:);
scales = length(s);
H=zeros(Ls,scales);


do_freq = 1;
do_func = 1;
do_filt = 1;

try
   feval([prefixFreq,wname]);
catch err
    do_freq = 0;
    if(strfind(err.identifier,'UndefinedFunction'))
       try 
          feval([prefixFunc,wname]); 
       catch err2
          if(strfind(err2.identifier,'UndefinedFunction'))
             do_func = 0;
             try 
                feval([prefixFilt,wname]); 
             catch err3
                if(strfind(err3.identifier,'UndefinedFunction'))
                    error('Function not found. Tryed %s, %s, %s.',[prefixFreq,wname],[prefixFunc,wname],[prefixFilt,wname]);
                end
             end
          end
       end
    end
end


if(do_freq)
   for ss=1:scales
      H(:,end+1-ss) = fftshift(feval([prefixFreq,wname],s(ss),Ls));
   end
elseif(do_func)
   for ss=1:scales
      wh = feval([prefixFunc,wname],s(ss));
      whlen = length(wh);
      whtemp = zeros(Ls,1);
      whtemp(1:whlen)= wh(:);
      H(:,end+1-ss) = fft(circshift(whtemp,-floor(whlen/2)));
      H(:,end+1-ss) = H(:,end+1-ss)/max(H(:,end+1-ss));
   end
elseif(do_filt)
      h = fwtinit(w,'syn');
      filtLen = length(h.filts{2}.h);
      [wvals,~,xvals] = wavfun(h);
      wvalsInt = cumsum([zeros(1,size(wvals,2));wvals;zeros(1,size(wvals,2))]);
      
   for ss=1:scales
      samples = (filtLen-1)/s(ss);
      t = (0:1/s(ss):filtLen-1);
      wh = interp1(wvalsInt,0:length(wvalsInt)-1,t); 
      whlen = length(wh);
      whtemp = zeros(Ls,1);
      whtemp(1:whlen)= wh(:);
      H(:,end+1-ss) = fft(circshift(whtemp,-floor(whlen/2)));
      H(:,end+1-ss) = H(:,end+1-ss)/max(H(:,end+1-ss));
   end 
else
    error('Should not ever get here.');
end