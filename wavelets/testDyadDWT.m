function c = testDyadDWT
%f = traindoppler;
f = expchirp(50000,0.1,1);
undec = 1;

% call to Wavelet Toolbox to get wavelet filters
wavelet = 'sym2';
[lo,hi,lo_r,hi_r] = wfilters(wavelet);
J = 10;


gan = makeMultirateIdentity(lo,hi,J);
grec = makeMultirateIdentity(lo_r,hi_r,J);
if undec == 1
 a = ones(J+1,1);
else
 a = zeros(J+1,1);
 for j=1:J ,  a(j) = 2^j; end
 a(J+1) = 2^J; 
end

c=filterbank(f,gan,a);


plotfilterbank(c,a,'surf');
%fhat=ifilterbank(c,grec,a);

%stem(fhat);
