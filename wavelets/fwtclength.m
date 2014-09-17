function [Lc,L]=fwtclength(Ls,w,J,varargin)
%FWTCLENGTH FWT subbands lengths from a signal length
%   Usage: Lc=fwtclength(Ls,w,J);
%          [Lc,L]=fwtclength(...);
%
%   `Lc=fwtclength(Ls,w,J)` returns the lengths of the wavelet coefficient
%   subbands for a signal of length *Ls*. Please see the help on |fwt| for
%   an explanation of the parameters *w* and *J*.
%
%   `[Lc,L]=fwtclength(...)` additianally the function returns the next 
%   legal length of the input signal for the given extension type.
%
%   The function support the same boundary-handling flags as the |fwt|
%   does.
%
%   See also: fwt, fwtlength

% AUTHOR: Zdenek Prusa

complainif_notposint(Ls,'Ls','FWTCLENGTH');
complainif_notposint(J,'J','FWTCLENGTH');

w = fwtinit(w);

definput.import = {'fwtext'};
[flags,kv]=ltfatarghelper({},definput,varargin);

% Get the next legal length
L = fwtlength(Ls,w,J,flags.ext);

filtNo = length(w.g);
subbNo = (filtNo-1)*J+1;
Lc = zeros(subbNo,1);
runPtr = 0;
levelLen = L;
if flags.do_per
  % Non-expansive case
  for jj=1:J
     for ff=filtNo:-1:2
        Lc(end-runPtr) = ceil(levelLen/w.a(ff));
        runPtr = runPtr + 1;
     end
     levelLen = ceil(levelLen/w.a(1));
  end
% elseif flags.do_valid
%   % Valid coef. case
%   filts = w.g;
%   for jj=1:J
%      for ff=filtNo:-1:2
%         Lc(end-runPtr) = floor((levelLen-(length(filts{ff}.h)-1))/w.a(ff));
%         runPtr = runPtr + 1;
%      end
%      levelLen = floor((levelLen-(length(filts{1}.h)-1))/w.a(1));
%   end
else
  % Expansive case
  filts = w.g;
  for jj=1:J
    for ff=filtNo:-1:2
       skip = w.a(ff) - 1;
       Lc(end-runPtr) = ceil((levelLen+(length(filts{ff}.h)-1)-skip)/w.a(ff));
       runPtr = runPtr + 1;
    end
    skip = w.a(1) - 1;
    levelLen = ceil((levelLen+(length(filts{1}.h)-1)-skip)/w.a(1));
 end
end
Lc(1)=levelLen;
