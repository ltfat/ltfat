function [y,state] = partconv_execute( x, state )
%PARTCONV_EXECUTE Execute single FDL partitioned convolution
%   Usage: [y,state] = partconv_init( x, state )
%
%   Input parameters:
%       x       : Input signal
%       state   : Convolution state
%
%   Output parameters:
%       y       : Output signal
%       state   : Output convolution state
%
%   [y, state] = PARTCONV_EXECUTE( x, state) executes partitioned convolution. The size of *x* must be 
%   *state.B x state.W*. The dimesnsions of the output *y* are identical. If multiple impulse responses 
%   were passed to the *_init function, the respective outputs are stacked along the 3rd dimension.
%   

%   AUTHOR : Zdenek Prusa (2023)

[B,W] = size(x);
[hw] = size(state.H_parts,3);
y = zeros([B,W,hw]);

state.FDL = circshift(state.FDL,[0,1]);

for w=1:W
    state.FDL(:,1,w) = fftreal([x(:,w);state.prevblock(:,w)]);
    y(:,w,:) = postpad( ifftreal( sum( bsxfun(@times, state.H_parts, state.FDL(:,:,w)), 2), 2*B), B );
end

state.prevblock = x;

end % function
