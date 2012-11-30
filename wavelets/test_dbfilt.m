function test_dbfilt

for ord=1:25

    [lo,hi,lo_r,hi_r] = wfilters(sprintf('db%d',ord));
    [h,g] = dbfilt(ord);
    err=cmp_filt(lo,hi,lo_r,hi_r,h,g);
    if(err>1e-6)
        figure(1);
        subplot(2,2,1); stem([lo(:),h{1}(:)]); title('Analysis lowpass'); legend('wfilters','dbfilt');
        subplot(2,2,2); stem([hi(:),h{2}(:)]); title('Analysis highpass');
        subplot(2,2,3); stem([lo_r(:),g{1}(:)]); title('Synthesis lowpass');
        subplot(2,2,4); stem([hi_r(:),g{2}(:)]); title('Synthesis highpass');
        error('Filters are not equal! Error is %g',err);
       pause;
    end



end


function err = cmp_filt(lo,hi,lo_r,hi_r,h,g)

err = 0;
err = err + norm(lo-h{1});
err = err + norm(hi-h{2});
err = err + norm(lo_r-g{1});
err = err + norm(hi_r-g{2});

