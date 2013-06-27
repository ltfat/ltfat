function c = blockana(F, f)
%BLOCKANA Blockwise analysis interface
%   Usage: c=blockana(F, f)
%
%   Input parameters:
%      Fa : Analysis frame object.    
%      f  : Block of signal.
%   Output parameters:
%      c  : Block coefficients.
%
%   `c=blockana(Fa,f)` calculates the coefficients *c* of the signal block *f* using 
%   the frame defined by *Fa*.
%
%   See also: block, blocksyn, blockplay
    
    if nargin<2
        error('%s: Too few input parameters.',upper(mfilename));
    end;
    
    if ~isstruct(F)
        error('%s: First agument must be a frame definition structure.',upper(mfilename));
    end;
    
    % Block length
    Lb = size(f,1);
    % Next block index start (from a global point of view, starting with zero)
    nextSb = block_interface('getPos');
    % Block index start (from a global point of view, starting with zero)
    Sb = nextSb-Lb;
    
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
        fext = [loadOverlap(Sbolen,size(f,2));f];
        c = block_fwt(fext,w,J);
        storeOverlap(fext,nextSbolen);
      otherwise
        % General processing
        % Equal block length assumtion
        % Slicing window
        g = fftshift(firwin('hann',2*Lb));
        % Append the previous block
        fext = [loadOverlap(Lb,size(f,2));f];
        % Save the current block
        storeOverlap(fext,Lb);
        % Multiply by the slicing window (all channels)
        fwin = bsxfun(@times,g,fext);
        % Apply transform
        c = F.frana(fwin);
    end

end % BLOCKANA

function overlap = loadOverlap(L,chan)
%LOADOVERLAP Loads overlap
%
%
    overlap = block_interface('getAnaOverlap');
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

function storeOverlap(fext,L)
%STOREOVERLAP Stores overlap
%
%
    if L>size(fext,1)
        error('%s: Storing more samples than passed.',upper(mfilename));
    end
    block_interface('setAnaOverlap',fext(end-L+1:end,:)); 
end % STOREOVERLAP

