function f = iwtfft(c,G,a,varargin)
%IWTFFT  Inverse Wavelet Transform in the frequency-domain
%   Usage: f=iwtfft(c,G,a);
%          f=iwtfft(c,G,a,...);
%
%   `iwtfft(c,G,a)` computes XXX.
%
%   See also: wtfft

    if nargin<2
      error('%s: Too few input parameters.',upper(mfilename));
    end;


%% PARSE INPUT
definput.keyvals.L=[];    
definput.import = {'wtfft'};

[flags,kv,Ls]=ltfatarghelper({'L'},definput,varargin);

[Gr,Gc] = size(G);

if(isempty(a))
    a = ones(Gc,1);
end

L=filterbanklengthcoef(c,a);
[sigHalfLen,W] = size(c{end});
f = zeros(L,W);


N = zeros(Gc,1);
for gg=1:Gc
    N(gg)=size(c{gg},1);
end


for w=1:W
   for gg=1:Gc
      f(:,w)=f(:,w)+ifft(repmat(fft(c{gg}(:,w)),a(gg),1).*G(:,gg));
   end
end

if ~isempty(Ls)
  f=postpad(f,Ls);
else
  Ls=L;
end;

