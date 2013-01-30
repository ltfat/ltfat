function c = wtfft(f,H,a,varargin)
%WTFFT Wavelet Transform via FFT



%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end

definput.import = {'wtfft'};
definput.keyvals.L = [];
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);




[Hr,Hc] = size(H);
if(isempty(a))
    a = ones(Hc,1);
end

if(isempty(L))
    L = filterbanklength(Ls,a);
end

if(strcmpi(H{1},'fir'))
    H = wtfftfreqz(H,L);
end



N=L./a;

c=cell(Hc,1);
for m=1:Hc
  c{m}=zeros(N(m),W);
end;


%% ----- step 2 : Run computation

for w=1:W
   F = fft(f(:,w),L);
   for hh=1:Hc
      c{hh}(:,w) = real( ifft(sum(reshape(F.*H(:,hh),N(hh),a(hh)),2))/a(hh) );
   end
end


