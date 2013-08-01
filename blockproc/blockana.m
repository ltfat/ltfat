function [c, fola] = blockana(F, f, fola)
%BLOCKANA Blockwise analysis interface
%   Usage: c=blockana(F, f)
%
%   Input parameters:
%      Fa   : Analysis frame object.    
%      f    : Block of signal.
%      fola : Explicitly defined overlap
%   Output parameters:
%      c    : Block coefficients.
%      fola : Stored overlap
%
%   `c=blockana(Fa,f)` calculates the coefficients *c* of the signal block *f* using 
%   the frame defined by *F*. The block overlaps are handled according to the 
%   `F.blokalg`. Assuming |blockana| is called in the loop only once, *fola*
%   can be omitted and the overlaps are handled in the background
%   automatically.    
%
%   See also: block, blocksyn, blockplay
    
    if nargin<2
        error('%s: Too few input parameters.',upper(mfilename));
    end;
    
    if ~isstruct(F)
        error('%s: First agument must be a frame definition structure.',upper(mfilename));
    end;
    
    if nargin<3
       fola = [];
    end
    
    % Block length
    Lb = size(f,1);
    % Next block index start (from a global point of view, starting with zero)
    nextSb = block_interface('getPos');
    % Block index start (from a global point of view, starting with zero)
    Sb = nextSb-Lb;
    
    if strcmp(F.blokalg,'sliced')
        % General processing
        % Equal block length assumtion
        % Append the previous block
        fext = [loadOverlap(Lb,size(f,2),fola);f];
        % Save the current block
        if nargout>1
           fola = storeOverlap(fext,Lb);
        else
           storeOverlap(fext,Lb);
        end
        % Multiply by the slicing window (all channels)
        fwin = bsxfun(@times,F.sliwin,fext);
        fwin = [fwin; zeros(F.L-size(fwin,1),size(fwin,2))];
        % Apply transform
        c = F.frana(fwin);
    elseif strcmp(F.blokalg,'segola')
      if ~isfield(F,'winLen') 
         error('%s: Frame does not have FIR windows.',upper(mfilename));
      end
      Lw = F.winLen;
      
      switch(F.type)
         case 'fwt'
           J = F.J;
           w = F.g;
           m = numel(w.h{1}.h);
           a = w.a(1);
           if Lb<a^J
               error('%s: Minimum block length is %i.',upper(mfilename),a^J);
           end
           rred = (a^J-1)/(a-1)*(m-a);
           Sbolen = rred + mod(Sb,a^J);
           nextSbolen = rred + mod(nextSb,a^J);
           
           % Extend from the left side
           fext = [loadOverlap(Sbolen,size(f,2),fola);f];
           % Save the overlap for the next block
           if nargout>1
              fola = storeOverlap(fext,nextSbolen);
           else
              storeOverlap(fext,nextSbolen);
           end
           
           c = block_fwt(fext,w,J);
         % Uniform window lengths  
         case {'dgtreal','dgt'}
%            if Lb<F.a
%                error('%s: Minimum block length is %i.',upper(mfilename),F.a);
%            end 
           a = F.a; 
           Lwl = floor(Lw/2);

           % Overlap length 
           Sbolen = ceil((Lw-1)/a)*a + mod(Sb,a);
           % Next block overlap length
           nextSbolen = ceil((Lw-1)/a)*a + mod(nextSb,a);
           Lext = Sbolen + Lb - nextSbolen + Lwl;
           startc = ceil(Lwl/a)+1;
           endc = ceil((Lext)/a);
           
           
           % Extend from the left side
           fext = [loadOverlap(Sbolen,size(f,2),fola);f];
           % Save the overlap for the next block
           if nargout>1
              fola = storeOverlap(fext,nextSbolen);
           else
              storeOverlap(fext,nextSbolen);
           end

           % Pad with zeros to comply with the frame requirements
           fext = [fext; zeros(F.L-size(fext,1),size(fext,2))];
           c = F.frana(fext(1:F.L,:));
           % Pick just valid coefficients
           cc = F.coef2native(c,size(c));
           cc = cc(:,startc:endc,:);
           c = F.native2coef(cc);
         case {'filterbank','filterbankreal','ufilterbank','ufilterbankreal'}
            
         otherwise
           error('%s: Unsupported frame.',upper(mfilename));
      end
    elseif strcmp(F.blokalg,'naive')
       f = [f; zeros(F.L-size(f,1),size(f,2))];
       c = F.frana(f);
    else
       error('%s: Frame was not created with blockaccel.',upper(mfilename));
    end

end % BLOCKANA

function overlap = loadOverlap(L,chan,overlap)
%LOADOVERLAP Loads overlap
%
%
    if isempty(overlap)
       overlap = block_interface('getAnaOverlap');
    end

    % Supply zeros if it is empty
    if isempty(overlap)
        overlap = zeros(L,chan,block_interface('getClassId'));
    end
    Lo = size(overlap,1);
    if nargin<1
        L = Lo;
    end
    % If required more than stored, do zero padding
    if L>Lo
        oTmp = zeros(L,size(overlap,2));
        oTmp(end-Lo+1:end) = oTmp(end-Lo+1:end)+overlap;
        overlap = oTmp;
    else
        overlap = overlap(end-L+1:end,:);
    end
    
end % LOADOVERLAP

function overlap = storeOverlap(fext,L)
%STOREOVERLAP Stores overlap
%
%
    if L>size(fext,1)
        error('%s: Storing more samples than passed.',upper(mfilename));
    end
    overlap = fext(end-L+1:end,:);
    
    if nargout<1
       block_interface('setAnaOverlap',overlap); 
    end
end % STOREOVERLAP

