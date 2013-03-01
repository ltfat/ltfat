function fdowns = comp_downs(f,varargin)
%COMP_DOWNS Downsampling
%   Usage: fdowns = comp_downs(f,a,type) 
%          fdowns = comp_downs(f,a,type,'dim',dim) 
%   
%   Input parameters:
%         f     : Input vector/matrix.
%         a     : Downsampling factor.
%         type  : Type of the downsampling.
%         dim   : Direction of downsampling.
%   Output parameters:
%         fdowns  : Downsampled vector/matrix.
%
%   Downsamples input *f* by a factor *a* (leaves every *a*th sample) along
%   dimension *dim*. If *dim* is not specified, first non-singleton
%   dimension is used. Parameter *type* (integer from [0:a-1]) specifies how
%   many samples to skip from the beginning.
%
%   Examples:
%   ---------
%
%   The outcome of the default upsampling type is equal to the subsampling performed
%   directly in the frequency domain using reshape and sum:::
%
%      f = 1:9;
%      a = 3;
%      fupsTD = comp_downs(f,a)
%      fupsFD = real(ifft(sum(reshape(fft(f),length(f)/a,a),2).'))/a;
%


definput.keyvals.dim = [];
definput.keyvals.a = 2;
definput.keyvals.type = 0;
[flags,kv,a,type]=ltfatarghelper({'a','type'},definput,varargin);

% a have to be positive integer
if(a<1)
    a = 1;
end
if(a<0 || rem(a,1)~=0)
    error('%s: Parameter *a* have to be a positive integer.',upper(mfilename));
end

% supported type are [0--a-1]
if(type<0||type>(a-1))
    error('%s: Unsupported downsampling type.',upper(mfilename));
end

if(ndims(f)>2)
    error('%s: Multidimensional signals (d>2) are not supported.',upper(mfilename));
end

% ----- Verify f and determine its length -------
[f,L,Ls,~,dim,permutedsize,order]=assert_sigreshape_pre(f,[],kv.dim,upper(mfilename));

% Actual computation
fdowns = f(1+type:a:end,:);

permutedSizeAlt = size(fdowns);
fdowns=assert_sigreshape_post(fdowns,dim,permutedSizeAlt,order);