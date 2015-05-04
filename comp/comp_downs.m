function fdowns = comp_downs(f,a,varargin)
%COMP_DOWNS Downsampling
%   Usage: fdowns = comp_downs(f,a) 
%          fdowns = comp_downs(f,a,skip,L,'dim',dim) 
%   
%   Input parameters:
%         f     : Input vector/matrix.
%         a     : Downsampling factor.
%         skip  : Skipped initial samples.
%         L     : Length of the portion of the input to be used.
%         dim   : Direction of downsampling.
%   Output parameters:
%         fdowns  : Downsampled vector/matrix.
%
%   Downsamples input *f* by a factor *a* (leaves every *a*th sample) along
%   dimension *dim*. If *dim* is not specified, first non-singleton
%   dimension is used. Parameter *skip* (integer) specifies how
%   many samples to skip from the beginning and *L* defines how many
%   elements of the input data are to be used starting at index 1+skip. 
%
%   Examples:
%   ---------
%
%   The default behavior is equal to the subsampling performed
%   in the frequency domain using reshape and sum:::
%
%      f = 1:9;
%      a = 3;
%      fupsTD = comp_downs(f,a)
%      fupsFD = real(ifft(sum(reshape(fft(f),length(f)/a,a),2).'))/a
%


if(nargin<2)
    error('%s: Too few input parameters.',upper(mfilename));
end

definput.keyvals.dim = [];
definput.keyvals.skip = 0;
definput.keyvals.L = [];
[flags,kv,skip,L]=ltfatarghelper({'skip','L','dim'},definput,varargin);

% a have to be a positive integer
if(a<1)
    a = 1;
end
if(a<0 || rem(a,1)~=0)
    error('%s: Parameter *a* have to be a positive integer.',upper(mfilename));
end

% supported type are [0--a-1]
% if(type<0||type>(a-1))
%     error('%s: Unsupported downsampling type.',upper(mfilename));
% end

if(ndims(f)>2)
    error('%s: Multidimensional signals (d>2) are not supported.',upper(mfilename));
end

% ----- Verify f and determine its length -------
[f,Lreq,Ls,~,dim,permutedsize,order]=assert_sigreshape_pre(f,L,kv.dim,upper(mfilename));

if(skip>=Ls)
    error('%s: Parameter *skip* have to be less than the input length.',upper(mfilename));
end

if(~isempty(L))
   if(Lreq+skip>Ls)
       error('%s: Input length is less than required samples count: L+skip>Ls.',upper(mfilename)); 
   end 
   Ls = Lreq+skip;
end


% Actual computation
fdowns = f(1+skip:a:Ls,:);   


permutedSizeAlt = size(fdowns);
fdowns=assert_sigreshape_post(fdowns,dim,permutedSizeAlt,order);
