function out = up_conv_td(in,outLen,filts,up,skip,ext,filtUps)
%UP_CONV_TD Upsampling followed by time-domain convolution
%   Usage:  out = up_conv_td(in,outLen,filts,up,skip,ext,filtUps)
%
%   Input parameters:
%         in       : length(filts) cell-array containing column vectors.
%         outLen   : Required output vectors lengths.
%         filts    : Cell-array containing impulse responses.
%         sub      : Subsampling factor.
%         skip     : Number of initial input samples skipped.
%         ext      : Boundary extension type.
%         filtUps  : Filter upsampling factor.
%
%   Output parameters:
%         out      : output collumn vector
%
%   Calculates length(filts)-channel synthesis filterbank response preceed by
%   upsampling with the factor *up*. Filter impulse
%   responses in *filts* are assumed to be of a equal length, initial
%   filter position in the upsampled input is given by *skip*. Optionally *filtUps* introduces
%   filter impulse response upsampling.


noOfFilts = length(filts);
fLen = length(ups(filts{1},filtUps,1));
inLen = length(in{1});
out = zeros(outLen,1);

inExt = zeros(inLen + 2*(fLen-1),1);

for ff=1:noOfFilts
    if(ext)
      inExt = extendBoundary(in{ff},fLen-1,'per');  
    else  
      inExt(fLen:end-fLen+1) = in{ff};
    end

    outTemp = conv(ups(inExt,up,1),ups(filts{ff}(:),filtUps,1));
    outTemp = outTemp(1+up*(fLen-1)+skip:end);
    toWrite = min([length(outTemp),outLen]);
    out(1:toWrite) = out(1:toWrite) + outTemp(1:toWrite);
end
    
    
