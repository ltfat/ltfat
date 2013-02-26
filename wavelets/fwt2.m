function c = fwt2(f,h,J,varargin)
%FWT2   Fast Wavelet Transform 2D
%   Usage:  c = fwt2(f,h,J);
%           c = fwt2(f,h,J,...);
%
%   Input parameters:
%         f     : Input data.
%         h     : Analysis Wavelet Filterbank. 
%         J     : Number of filterbank iterations.
%
%   Output parameters:
%         c      : Coefficients stored in a cell-array.
%
%   `c=fwt2(f,h,J)` returns wavelet coefficients *c* of the input matrix *f*
%   using *J* iterations of the basic wavelet filterbank defined by *h*.
%   
%   
%
%   Examples:
%   ---------
%   
%   A simple example of calling the |fwt2|_ function:::
% 
%     figure(1);
%     image(cameraman);axis equal;axis off;colormap Gray;
%     c = fwt2(cameraman,{'db',8},4);
%     figure(2);
%     imagesc(20*log10(abs(c)));axis equal;axis off;colormap jet;
%
%   See also: ifwt2, fwtinit
%
%   References: ma98  


if nargin<3
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isnumeric(J) || ~isscalar(J)
  error('%s: "J" must be a scalar.',upper(mfilename));
end;

if(J<1 || rem(J,1)~=0)
   error('%s: J must be a positive integer.',upper(mfilename)); 
end

[M,N]=size(f);
if(M==1||N==1)
   error('%s: The input data is vector.',upper(mfilename)); 
end

% Initialize the wavelet filters structure
h = fwtinit(h,'ana');

if(~all(h.a==length(h.filts)))
   error('%s: Non-critically subsampled filterbanks not supported.',upper(mfilename));  
end


%% ----- step 0 : Check inputs -------
definput.import = {'fwt','fwt2'};
[flags,kv]=ltfatarghelper({},definput,varargin);
nFilts = numel(h.filts);

Lcrows = fwtclength(size(f,1),h,J);
Lccols = fwtclength(size(f,2),h,J);

if(flags.do_standard)
   Jstep = 1;
   c = fwt(f,h,Jstep,1,'per');
   c = fwt(c,h,Jstep,2,'per');
   for jj=1:J-1
      colRange = 1:Lcrows(end-jj*(nFilts-1)+1);
      rowRange = 1:Lccols(end-jj*(nFilts-1)+1);
      c(colRange,rowRange) = fwt(c(colRange,rowRange),h,Jstep,1,'per');
      c(colRange,rowRange) = fwt(c(colRange,rowRange),h,Jstep,2,'per');
   end
elseif(flags.do_tensor)
   c = fwt(f,h,J,'dim',1,'per');
   c = fwt(c,h,J,'dim',2,'per');
else
    error('%s: Should not get here. Bug somewhere else.',upper(mfilename))
end


%% ----- step 1 : Run calc -----------
% filtNo = length(h.filts);
% subbNo = filtNo^2-1;
% c = cell(J*subbNo+1,1);
% cJidx = J*subbNo+1;
% cTmp = f;
% for jj=1:J
%      cCols = comp_fwt_all(cTmp,h.filts,1,h.a,'dec',flags.ext);
%      [cColsPack,rows] = wavcell2pack(cCols); 
%      cRows = comp_fwt_all(cColsPack.',h.filts,1,h.a,'dec',flags.ext);
%      [cRowsPack,cols] = wavcell2pack(cRows); 
%      cTmp = cRowsPack(1:cols(1),1:rows(1)).';
%      
%      cJidxTmp = cJidx - jj*subbNo+1;
% 
%      for cc= 1:filtNo
%         for rr = 1:filtNo
%            if(cc==1&&rr==1), continue; end; 
%            c{cJidxTmp} = cRowsPack(sum(cols(1:cc-1))+1:sum(cols(1:cc)),sum(rows(1:rr-1))+1:sum(rows(1:rr))).';
%            cJidxTmp = cJidxTmp + 1;
%         end
%      end
% end
% c{1} = cTmp;


