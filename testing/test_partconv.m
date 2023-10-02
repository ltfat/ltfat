function failed = test_partconv

failed = 1;

[f,fs] = gspi;
ftmp = f(1:fs,:);
hl = fs;
B = 2048;

Warr  = [1,1,2,2];
HWarr = [1,2,1,2];

for ii = 1:numel(Warr)
    f = repmat(ftmp,1,Warr(ii));
    f = bsxfun(@times,f,1:Warr(ii));
    hw = HWarr(ii);

    [L,W] = size(f);
    h = randn(hl,hw);
    Lpad = ceil((L + hl)/B)*B;

    state = partconv_init( B, W, h);
    f = postpad(f, Lpad);

    f_blocks = zeros( B, W, Lpad/B);

    for w = 1:W
        f_blocks(:,w,:) = reshape( f(:,w), B, Lpad/B );
    end

    fout_blocks = zeros( B, W, hw, Lpad/B );

    block_count = Lpad/B;

    for b = 1:block_count
        [ fout_blocks(:,:,:,b), state ] =  partconv_execute( f_blocks(:,:,b), state );
    end

    fout = reshape(permute(fout_blocks,[1,4,2,3]),Lpad,W,hw);

    fout_ref = zeros(Lpad,W,hw);

    for w = 1:W
        for n = 1:hw
            fout_ref(:,w,n) = postpad(conv(f(:,w),h(:,n),'full'),Lpad);
        end
    end

    % plot(fout-fout_ref);

    failed = norm(fout(:)-fout_ref(:), 'inf') > 10^-10;

    if failed, break; end
end

