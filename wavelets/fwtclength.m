function [Lc,L]=fwtclength(Ls,h,J,varargin) 
%FWTCLENGTH FWT subbands length fom signal
%   Usage: L=fwtlength(Ls,h,J);
%          L=fwtlength(Ls,h,J,...);
%
%   `Lc=fwtclength(Ls,h,J)` returns the lengths of the wavelet coefficient
%   subbands for a signal of length *Ls*. Please see the help on |fwt|_ for
%   an explanation of the parameters *h* and *J*. In addition, the function
%   returns the next legal length of the input signal for the given extension
%   type.
%
%   See also: fwt, fwtlength


h = fwtinit(h,'ana');
definput.import = {'fwt'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Get the next legal length
L = fwtlength(Ls,h,J,flags.ext);

filtNo = length(h.filts);
subbNo = (filtNo-1)*J+1;
Lc = zeros(subbNo,1);
runPtr = 0; 
levelLen = L;
if(flags.do_per)
  % Non-expansive case
  for jj=1:J
     for ff=filtNo:-1:2
        Lc(end-runPtr) = ceil(levelLen/h.a(ff));
        runPtr = runPtr + 1;
     end
     levelLen = ceil(levelLen/h.a(1));
  end
else
  % Expansive case
  filts = h.filts;
  for jj=1:J
    for ff=filtNo:-1:2
       Lc(end-runPtr) = floor((levelLen+(length(filts{ff}.h)-1))/h.a(ff));
       runPtr = runPtr + 1;
    end
    levelLen = floor((levelLen+(length(filts{1}.h)-1))/h.a(1));
 end
end
Lc(1)=levelLen; 