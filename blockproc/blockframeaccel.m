function Fo = blockframeaccel(F, Lb, varargin)
%BLOCKFRAMEACCEL Precompute structures for block processing
%   Usage: F = blockframeaccel(F,Lb);
%
%   `F=blockframeaccel(F,Lb)` have to be called for each frame object prior 
%   entering the main loop where |blockana| and |blocksyn| are called.
%   The function work entirely like |frameaccel| but in addition, it prepares
%   structures for the processing of a consecutive stream of blocks.
%
%      `'sliwin',sliwin`   : Slicing window. `sliwin` have to be a window
%                            of length *2Lb* or a string accepted
%                            by the |firwin| function. It is used only in
%                            the slicing window approach. The default is 
%                            `'hann'`.
%
%      `'transa',transa`   : Transition area length. Number of zeros to be 
%                            used to pad the slicing window to *2Lb* from
%                            both sides. The option is used as padding
%                            only in the slicing window approach using 
%                            window defined as a string.



definput.flags.blockalg = {'naive','sliced','segola'};
definput.keyvals.sliwin = [];
definput.keyvals.transa = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

assert(~(~flags.do_sliced && ~isempty(kv.sliwin) && ~isempty(kv.transa)),...
   sprintf('%s: Definig slicing window without setting the ''silced'' flag.',...
   mfilename));

if flags.do_sliced 
   if isempty(kv.sliwin)
      kv.sliwin = 'hann';
   end

   if isempty(kv.transa)
      kv.transa = 0;
    end
   
   if ~isnumeric(kv.sliwin)
      kv.sliwin = fftshift(firwin(kv.sliwin,2*Lb));
   end

   Fo = frameaccel(F,2*Lb+2*kv.transa);
   Fo.sliwin = kv.sliwin;
   Fo.transa = kv.transa;
elseif flags.do_segola
   % Determine window length without calling frameaccel
   % Fo = frameaccel(F,Lb);
   winLen = framefirlen(F);

   if winLen==-1
      error(['%s: Segment overlap cannot be used with this frame.,'...
             ' It does not have FIR windows.'],upper(mfilename));
   end
    
   switch(F.type) 
      case {'dgt','dgtreal'}
         Fo = frameaccel(F,Lb+winLen-1+F.a);
      case {'filterbank','filterbankreal','ufilterbank','ufilterbankreal'}
         lcma =  filterbanklength(1,F.a(:,1));
         Fo = frameaccel(F,Lb+winLen-1+lcma);
         assert(all(Fo.a(:,2)==1), '%s: Fractional subsampling is not supported',upper(mfilename) );
         Fo.lcma =  lcma;
      case {'dwilt'}
         Fo = frameaccel(F,Lb+winLen-1+2*F.M);
         Fo.a = 2*Fo.M;
      case {'wmdct'}
         Fo = frameaccel(F,Lb+winLen-1+F.M);
         Fo.a = Fo.M;
      otherwise
	 error('%s: Unsupported frame for segola.',upper(mfilename));
   end
   
   % This is important otherwise we would get 0 coefficients for some
   % blocks.
   assert(max(Fo.a(:,1)) <= Lb ,sprintf(['%s: Time step %i is bigger than the',...
      ' block length %i.'],upper(mfilename),Fo.a,Lb));
   
   Fo.winLen = winLen;


elseif flags.do_naive
   Fo = frameaccel(F,Lb);
end

Fo.blockalg = flags.blockalg;

function winLen = framefirlen(F)
%FRAMEFIRLEN Frame window/filter length
%
%   Function returns length of the longest FIR window/filter. The function
%   returns -1 if the frame does not have FIR windows.

winLen = -1;
info = [];
switch(F.type)
      case {'dgt','dgtreal'}
        [~, info] =  gabwin(F.g,F.a,F.M,[],F.kv.lt);
      case {'dwilt','wmdct'}
        [~, info] = wilwin(F.g,F.M,[],upper(mfilename));
      case {'filterbank','ufilterbank'}
        [~, info]  = filterbankwin(F.g,F.a);
      case {'filterbankreal','ufilterbankreal'}
        [~, info]  = filterbankwin(F.g,F.a,'real');
      case 'fwt' 
        winLen = (F.g.a(1)^F.J-1)/(F.g.a(1)-1)*(numel(F.g.g{1}.h)-1)+1; 
end;

  
if ~isempty(info) && isfield(info,'isfir') && info.isfir
   if isfield(info,'longestfilter')
      winLen = info.longestfilter;
   else
      winLen = max(info.gl);
   end
end




