function c = fwt(f,h,J,varargin)
%FWT   Fast Wavelet Transform 
%   Usage:  c = fwt(f,h,J);
%           c = fwt(f,h,J,...);
%           c = fwt(f,h,J,a,...);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%         a     : Subsampling factor(s).
%
%   Output parameters:
%         c      : Coefficients stored in a cell-array.
%
%   `c=fwt(f,h,J)` computes wavelet coefficients *c* of the input signal *f*
%   using a basis (or frame) constructed from the filters specified in *h* and the depth *J* using
%   the MRA principle.  For computing the coefficients the fast wavelet
%   transform algorithm (or Mallat's algorithm) is emplyed. If *f* is a
%   matrix, the transformation is applied to each of *W* columns.
%
%   The wavelet filters $h$ should be a structure produced by the waveletfb
%   function. The .h field of the structure passed is a cell-array as in
%   the following case.
%
%   The wavelet filters can be also passed directly as a collumn cell-array in which each entry
%   contain single filter impulse response. Functions with wfilt_ prefix generate
%   such cell arrays. The number of filters will be refered to as $length(h)$.
%
%   The coefficients in the cell array $c\{j\}$ for $j=1,\ldots,J*(length(h)-1)+1$ are ordered
%   with inceasing central frequency of the equivalent filter frequency
%   response or equivalently with decreasing wavelet scale.
%
%   If the input *f* is matrix with *W* columns, each element of the cell
%   array $c\{j\}$ is a matrix with *W* columns with coefficients belonging
%   to the appropriate input channel.
%
%   The following flag groups are supported:
%
%         'dec','undec'
%                Type of the wavelet transform.
%
%         'per','zpd','sym','symw','asym','asymw','ppd','sp0'
%                Type of the boundary handling.
%
%   Time-invariant wavelet tranform:
%   --------------------------------
%
%   `c=fwt(f,h,J,'undec')` computes redundant time (or shift) invariant
%   wavelet representation of the input signal *f* using the "a-trous"
%   algorithm. `length(c{j})=length(f)` for all *j*.
%
%   Other names for this version of the wavelet transform are: the
%   undecimated wavelet transform, the stationary wavelet transform or even
%   the "continuous" (as the time step is one sample) wavelet transform.  The
%   redundancy is exactly *J(length(h)-1)+1*.
%
%   Boundary handling:
%   ------------------
%
%   `fwt(f,h,J,ext)` computes slightly redundant wavelet representation of
%   the input signal *f* with the chosen boundary extension *ext*.
%
%   The default periodic extension at the signal boundaries can result in
%   "false" high wavelet coefficients near the boundaries due to the
%   possible discontinuity introduced by the periodic extension. Using
%   different kind of boundary extensions comes with a price of a slight
%   redundancy of the wavelet representation.
%
%   See also: ifwt, waveletfb, wfilt_db
%
%   Demos: demo_fwt
%
%   References: ma98  



if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1 && rem(a,1)~=0)
   error('%s: J must be a positive integer.',upper(mfilename)); 
end

do_definedfb = 0;
if(iscell(h))
    if(length(h)<2)
       error('%s: h is expected to be a cell array containing two or more filters.',upper(mfilename)); 
    end

    for ii=2:numel(h)
     if(length(h{1})~=length(h{ii}))
        error('%s: Wavelet filters have to have equal length.',upper(mfilename));
     end
    end
    
    if(length(h{1})< 2)
        error('%s: Wavelet filters should have at least two coefficients.',upper(mfilename)); 
    end
elseif(isstruct(h))
    do_definedfb = 1;
elseif(ischar(h))
    h = waveletfb(h);
    do_definedfb = 1;
elseif(isnumeric(h))
    % TO DO: Does it make sense to accept filter coefficients as matrix?
else
   error('%s: Unrecognized Wavelet filters definition.',upper(mfilename)); 
end

%% ----- step 0 : Check inputs -------
definput.keyvals.a = [];
definput.import = {'fwt'};
[flags,kv,a]=ltfatarghelper({'a'},definput,varargin);
anotequal = 0;
if(length(a)>1)
   if(do_definedfb)
       if(length(a)~=length(h.h))
          anotequal=1;
       end
   else
       if(length(a)~=length(h))
          anotequal=1;
       end
   end
end

if(anotequal)
   error('%s: Number of the subsampling factors is not equal to the number of filters.',upper(mfilename)); 
end

if(length(a)==1)
   if(do_definedfb)
       a = a*ones(length(h.h),1); 
   else
       a = a*ones(length(h),1); 
   end
end

if(do_definedfb)
    if(flags.do_type_null)
       flags.type = h.type; 
    end

    if(flags.do_ext_null)
       flags.ext = h.ext; 
    end
    
    if(isempty(a))
      a = h.a;
    end
    h = h.h;
else
    % setting defaults
    if(flags.do_type_null)
       flags.type = 'dec'; 
    end

    if(flags.do_ext_null)
       flags.ext = 'per'; 
    end

    if(isempty(a))
      a = length(h)*ones(length(h),1);
    end
end




%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
% TO DO: if elements of a are not equal and the flags are 'dec' and 'per',
% do zero padding to the next multiple of 2^J
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);
if(Ls<2)
   error('%s: Input signal seems not to be a vector of length > 1.',upper(mfilename));  
end
 


%% ----- step 2 : Check whether the input signal is long enough
% input signal length is not restricted for expansive wavelet transform (extension type other than the default 'per')
flen = length(h{1});
if(strcmp(flags.ext,'per'))
   if(strcmp(flags.type,'dec'))
     minLs = (a(1)^J-1)*(flen-1); % length of the longest equivalent filter -1
   else
     minLs = (a(1)^(J-1))*(flen-1); % length of the longest upsampled filter - 1
   end
   if Ls<minLs
     error('%s: Input signal length is %d. Minimum signal length is %d or use %s flag instead. \n',upper(mfilename),Ls,minLs,'''ppd''');
   end;
end



%% ----- step 3 : Run computation
 c = comp_fwt_all(f,h,J,a,flags.type,flags.ext);


