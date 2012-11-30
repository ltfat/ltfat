function c = wtfft(f,H,a,varargin)
%WTFFT Wavelet Transform via FFT

% if nargin<3
%   error('%s: Too few input parameters.',upper(mfilename));
% end;
% 
% if ~isnumeric(J) || ~isscalar(J)
%   error('%s: "J" must be a scalar.',upper(mfilename));
% end;
% 
% if(J<1 && rem(a,1)~=0)
%    error('%s: J must be a positive integer.',upper(mfilename)); 
% end
% 
% do_definedfb = 0;
% if(iscell(h))
%     if(length(h)<2 || length(h) > 2)
%        error('%s: h is expected to be a cell array containing two wavelet filters.',upper(mfilename)); 
%     end
% 
%     if(length(h{1})< 2 && length(h{2})< 2)
%         error('%s: Wavelet filters should have at least two coefficients.',upper(mfilename)); 
%     end
% 
%     if(length(h{1})~=length(h{2}))
%         error('%s: Wavelet filters have to have equal length.',upper(mfilename));
%     end
% elseif(isstruct(h))
%     do_definedfb = 1;
% elseif(ischar(h))
%     h = waveletfb(h);
%     do_definedfb = 1;
% else
%    error('%s: Unrecognized Wavelet filters definition.',upper(mfilename)); 
% end

if isa(H, 'function_handle')
    disp('it is function handle!');
    return;
end


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


