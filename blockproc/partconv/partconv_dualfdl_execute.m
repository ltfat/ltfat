function [y,state] = partconv_dualfdl_execute( x, state )
%PARTCONV_DUALFDL_EXECUTE Compute dual FDL partitioned convolution
%   Usage: [y, state] = partconv_dualfdl_execute( x, state )
%
%   Input parameters:
%       x       : Input signal
%       state   : Input convolution state
%
%   Output parameters:
%       y       : Output signal
%       state   : Output convolution state
%
%   [y,state] = PARTCONV_DUALFDL_EXECUTE( x, state ) computes dual-length FDL partitioned covolution. The size of
%   *x* must be *state.B_short x state.W*. The long FDL is computed at a reduced block rate and is independent of the short FDL.
%   It is therefore bound to produce a comp. spike every state.B_long/state.B_short block.
%
%   If multiple impulse responses were passed to the *_init function, the respective outputs are stacked along the 3rd dimension.
%

%   AUTHOR : Zdenek Prusa (2023)

[B,W] = size(x);
[hw] = size(state.H_parts_long,3);

if B ~= state.B_short
    error('Wrong buffer length');
end

y = zeros([B,W,hw]);

% Short FDL
state.FDL_short = circshift(state.FDL_short,[0,1]);

for w =1:W
    state.FDL_short(:,1,w) = fftreal([ x(:,w); state.prevblock_short(:,w)]);

    y(:,w,:) = postpad( ifftreal( sum( bsxfun(@times,state.H_parts_short,state.FDL_short(:,:,w)),2),2*state.B_short), state.B_short ) + ...
             state.FDL_buf_out( state.FDL_counter*state.B_short+1:(state.FDL_counter+1)*state.B_short,w,:,state.FDL_buf_out_ridx);
end

state.prevblock_short = x;

state.FDL_buf_in( state.FDL_counter*state.B_short+1:(state.FDL_counter+1)*state.B_short,:,state.FDL_buf_in_widx) = x;

state.FDL_counter = state.FDL_counter + 1;

if state.FDL_counter == state.B_long/state.B_short
    [state.FDL_buf_in_ridx, state.FDL_buf_in_widx ]   = swap(state.FDL_buf_in_ridx, state.FDL_buf_in_widx );

    % Long FDL
    state.FDL_long = circshift(state.FDL_long,[0,1]);
    
    for w = 1:W
        state.FDL_long(:,1,w) = ...
            fftreal([state.FDL_buf_in(:,w,state.FDL_buf_in_ridx);state.prevblock_long(:,w)]);

        state.FDL_buf_out(:,w,:,state.FDL_buf_out_widx) = ...
            postpad( ifftreal( sum( bsxfun(@times,state.H_parts_long,state.FDL_long(:,:,w)),2),2*state.B_long), state.B_long );
    end

    state.prevblock_long = state.FDL_buf_in(:,:,state.FDL_buf_in_ridx);
    state.FDL_counter = 0;

    [state.FDL_buf_out_ridx, state.FDL_buf_out_widx ] = swap(state.FDL_buf_out_ridx, state.FDL_buf_out_widx );
        
end
end %function

function [b, a] = swap(a, b)
%Empty
end %function
