function failed = test_partconv_dualfdl

failed = 1;

[f,fs] = gspi;
ftmp = f(1:fs,:);
f = ftmp;
hl = fs;
B_long  = 8192;
B_short = B_long/8;

Warr  = [1,1,2,2];
HWarr = [1,2,1,2];

for ii = 1:numel(Warr)
    f = repmat(ftmp,1,Warr(ii));
    f = bsxfun(@times,f,1:Warr(ii));
    hw = HWarr(ii);

    [L,W] = size(f);
    h = randn(hl,hw);
    Lpad = ceil((L + hl)/B_long)*B_long;

    state = partconv_dualfdl_init( B_short, B_long, W, h);

    f = postpad(f,Lpad);

    f_blocks = zeros( B_short, W, Lpad/B_short);

    for w = 1:W
        f_blocks(:,w,:) = reshape( f(:,w), B_short, Lpad/B_short );
    end

    fout_blocks = zeros( B_short, W, hw, Lpad/B_short );

    block_count = Lpad/B_short;

    for b = 1:block_count
        [ fout_blocks(:,:,:,b), state ] =  partconv_dualfdl_execute( f_blocks(:,:,b), state );
    end

    fout = reshape(permute(fout_blocks,[1,4,2,3]),Lpad,W,hw);

    fout_ref = zeros(Lpad,W,hw);
    for w = 1:W
        for n = 1:hw
            fout_ref(:,w,n) = postpad(conv(f(:,w),h(:,n),'full'),Lpad);
        end
    end

    % plot([fout,fout_ref]);

    failed = norm(fout(:)-fout_ref(:), 'inf') > 10^-10;

    if failed, break; end
end



