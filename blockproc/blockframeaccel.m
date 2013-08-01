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
%                            of length *2L*. It is used in the slicing
%                            window approach.



definput.flags.blockalg = {'naive','sliced','segola'};
definput.keyvals.sliwin = [];
[flags,kv]=ltfatarghelper({},definput,varargin);

assert(~(~flags.do_sliced && ~isempty(kv.sliwin)),...
   '%s: Definig slicing window without setting the ''silced'' flag.',mfilename);

if flags.do_sliced 
   if isempty(kv.sliwin)
      kv.sliwin = 'hann';
   end

   if ~isnumeric(kv.sliwin)
      kv.sliwin = fftshift(firwin(kv.sliwin,2*Lb));
   else
      if numel(kv.sliwin)~=2*Lb
         error('%s: The slicing window length has to be 2*Lb=%i.',upper(mfilename),2*Lb);
      end
   end

   Fo = frameaccel(F,2*Lb);
   Fo.sliwin = kv.sliwin;
elseif flags.do_segola
   Fo = frameaccel(F,Lb);

   if ~isfield(Fo,'winLen')
      error(['%s: Segment overlap cannot be used with this analysis frame.,'...
             ' It does not have FIR windows.'],upper(mfilename));
   end
    
   switch(Fo.type) 
      case {'dgt','dgtreal'}
         Fo = frameaccel(F,Lb+Fo.winLen-1+Fo.a);
         assert(Fo.a <= Lb ,['%s: Time step is bigger than the',...
      ' block length.'],mfilename);

      case {'dwilt'}
         Fo = frameaccel(F,Lb+Fao.winLen-1+2*Fo.M);
         Fo.a = 2*Fo.M;
      case {'wmdct'}
         Fo = frameaccel(F,Lb+Fao.winLen-1+Fo.M);
         Fo.a = Fo.M;
   end
   


elseif flags.do_naive
   Fo = frameaccel(F,Lb);
end

Fo.blokalg = flags.blockalg;