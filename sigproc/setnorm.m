function [f,fnorm]=setnorm(f,varargin)
%SETNORM  Set a norm of the input signal to a specified value
%   Usage:  h=setnorm(f,...);
%
% 
%   `setnorm(f,...)` will set the input signal *f* to a specified norm.
%
%   `[f,fnorm]=setnorm(f,...)` does the same thing, but in addition
%   returns norm *fnorm* of a signal *f*.
%
%   The norm is specified as a string and may be one of:
%
%     '1'       Normalize the $l^1$ norm to be *1*.
%
%     'area'    Normalize the area of the signal to be *1*. This is exactly the same as `'1'`.
%
%     '2'       Normalize the $l^2$ norm to be *1*.
%
%     'energy'  Normalize the energy of the signal to be *1*. This is exactly
%               the same as `'2'`.
%
%     'inf'     Normalize the $l^{\inf}$ norm to be *1*.
%
%     'peak'    Normalize the peak value of the signal to be *1*. This is exactly
%               the same as `'inf'`.
%
%     'rms'     Normalize the Root Mean Square (RMS) norm of the
%               signal to be *1*.
%
%     's0'      Normalize the S0-norm to be *1*.
%
%     'wav'     Normalize to the $l^{\inf}$ norm to be *0.99* to avoid 
%               possible clipping introduced by the quantization procedure 
%               when saving as a wav file. This only works with floating
%               point data types.
%
%     'null'    Do NOT normalize, output is identical to input.
%
%
%   It is possible to specify the value to which the norm is set and the dimension:
%
%     'val',v    Value to set the norm to
%
%      'dim',d  
%                Work along specified dimension. The default value of `[]`
%                means to work along the first non-singleton one.
%
%   See also: rms, s0norm
  
if ~isnumeric(f) 
  error('%s: Input must be numerical.',upper(mfilename));
end;

if nargin<1
  error('%s: Too few input parameters.',upper(mfilename));
end;


definput.import={'setnorm'};
definput.keyvals.dim=[];
definput.keyvals.val=[];

[flags,kv]=ltfatarghelper({},definput,varargin,'setnorm');

if flags.do_null || flags.do_norm_notset || isempty(f);
  return
end;

if isa(f,'integer') && ~flags.do_wav
   error('%s: Integer data types are unsupported.',upper(mfilename)); 
end

if isempty(kv.val)
    kv.val=1;
end
%% ------ Computation --------------------------
 
[f,L,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(f,[],kv.dim, ...
                                                  upper(mfilename));
fnorm = zeros(W,1);

for ii=1:W  
  
  if flags.do_1 || flags.do_area
    fnorm(ii) =  norm(f(:,ii),1);
  end;

  if flags.do_2 || flags.do_energy 
    fnorm(ii) = norm(f(:,ii),2);
  end;

  if flags.do_inf || flags.do_peak
    fnorm(ii) = norm(f(:,ii),Inf);
  end;

  if flags.do_rms
    fnorm(ii) = rms(f(:,ii));
  end;
  
  if flags.do_s0 
    fnorm(ii) = s0norm(f(:,ii));
  end;
  
  if flags.do_wav
    if isa(f,'float')
       fnorm(ii) = (1/0.99)*norm(f(:,ii),Inf); 
    else
       error(['%s: TO DO: Normalizing integer data types not supported ',...
              'yet.'],upper(mfilename));
    end
  end;
  
  fnorm = 1/kv.val*fnorm;
  
  if fnorm(ii) > 0
     f(:,ii)=f(:,ii)/fnorm(ii);
  end
end;


f=assert_sigreshape_post(f,kv.dim,permutedsize,order);
