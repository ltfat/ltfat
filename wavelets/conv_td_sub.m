function out = conv_td_sub(in,outLen,filts,sub,skip,ext,filtUps)
%CONV_TD_SUB Time-domain convolution followed by subsampling
%   Usage:  out = conv_td_sub(in,outLen,filts,sub,skip,ext,filtUps)
%
%   Input parameters:
%         in       : Input column vector.
%         outLen   : Required output vectors lengths.
%         filts    : Cell-array containing impulse responses.
%         sub      : Subsampling factor.
%         skip     : Number of initial input samples skipped.
%         ext      : Boundary extension type.
%         filtUps  : Filter upsampling factor.
%
%   Output parameters:
%         out      : length(filts) cell-array containing outputs
%                    if just 1 filter is provided, the output is vector.
%
%   Calculates length(filts)-channel analysis filterbank response followed by
%   subsampling with the factor *sub*. Filter impulse
%   responses in *filts* are assumed to be of a equal length, initial
%   filter position is given by *skip*. Optionally *filtUps* introduces
%   filter impulse response upsampling.


noOfFilts = length(filts);
fLen = length(ups(filts{1},filtUps,1));
inLen = length(in);

if(noOfFilts>1)
  out = cell(noOfFilts,1);
  for ff=1:noOfFilts
     out{ff} = zeros(outLen,1);
  end
else
   out = zeros(outLen,1);
end

inExt = extendBoundary(in,fLen-1,ext);

for ff=1:noOfFilts
    outTemp = conv(inExt,ups(filts{ff}(:),filtUps,1));
    outTemp = downs(outTemp(fLen+skip:end-fLen+1),sub,1);
    toWrite = min([length(outTemp),outLen]);
    if(noOfFilts>1)
       out{ff}(1:toWrite) = outTemp(1:toWrite);
    else
       out(1:toWrite) =  outTemp(1:toWrite);
    end
end

