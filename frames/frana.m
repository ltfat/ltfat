function outsig=frana(F,insig);
%FRANA  Frame analysis operator
%   Usage: c=frana(F,f);
%
%   `c=frana(F,f)` computes the frame coefficients *c* of the input
%   signal *f* using the frame *F*. The frame object *F* must have been
%   created using |frame| or |framepair|.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   The output coefficients are stored as columns. This is usually
%   **not** the same format as the 'native' format of the frame. As an
%   examples, the output from |frana| for a gabor frame cannot be
%   passed to |idgt| without a reshape.
%
%   See also: frame, framepair, frsyn, plotframe

if nargin<2
  error('%s: Too few input parameters.',upper(mfilename));
end;

if ~isstruct(F)
  error('%s: First agument must be a frame definition structure.',upper(mfilename));
end;


%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[insig,~,Ls,W,dim,permutedsize,order]=assert_sigreshape_pre(insig,[],[],upper(mfilename));
 
F=frameaccel(F,Ls);

insig=postpad(insig,F.L);

%% ----- do the computation ----

if isfield(F,'frana')

    outsig=F.frana(insig);

else

    switch(F.type)
      case 'dgt'
        outsig=framenative2coef(F,comp_dgt(insig,F.g,F.a,F.M,F.kv.lt,F.flags.do_timeinv,0,0));
      case 'dgtreal'
        outsig=framenative2coef(F,comp_dgtreal(insig,F.g,F.a,F.M,F.kv.lt,F.flags.do_timeinv));
      case 'dwilt'
        outsig=framenative2coef(F,comp_dwilt(insig,F.g,F.M,F.L));
      case 'wmdct'
        outsig=framenative2coef(F,comp_dwiltiii(insig,F.g,F.M,F.L));
        
      case 'filterbank'
        outsig=framenative2coef(F,filterbank(insig,F.g,F.a));
      case 'filterbankreal'
        outsig=framenative2coef(F,filterbank(insig,F.g,F.a));
      case 'ufilterbank'
        outsig=framenative2coef(F,ufilterbank(insig,F.g,F.a));
      case 'ufilterbankreal'
        outsig=framenative2coef(F,ufilterbank(insig,F.g,F.a));
        
      case 'nsdgt'
        outsig=framenative2coef(F,nsdgt(insig,F.g,F.a,F.M));
      case 'unsdgt'
        outsig=framenative2coef(F,unsdgt(insig,F.g,F.a,F.M));
      case 'nsdgtreal'
        outsig=framenative2coef(F,nsdgtreal(insig,F.g,F.a,F.M));
      case 'unsdgtreal'
        outsig=framenative2coef(F,unsdgtreal(insig,F.g,F.a,F.M));
        
      case 'fwt'
        outsig=fwt(insig,F.g,F.J);
      case 'fusion'
        % All frames must use the same length signal.
        L=framelength(F,size(insig,1));
        insig=postpad(insig,L);
        
        coefs = cell(F.Nframes,1);
        for ii=1:F.Nframes
            coefs(ii)={F.w(ii)*frana(F.frames{ii},insig)};
        end;
        outsig=cell2mat(coefs);
      case 'tensor'
        outsig=frana(F.frames{1},insig);
        perm=[circshift((1:F.Nframes).',-1);
              F.Nframes+1:ndims(insig)];
        for ii=2:F.Nframes
            outsig=permute(outsig,perm);
            outsig=frana(F.frames{ii},outsig);
        end;
        outsig=permute(outsig,perm);
    end;

end;

%% --- cleanup -----

permutedsize=[size(outsig,1),permutedsize(2:end)];

outsig=assert_sigreshape_post(outsig,dim,permutedsize,order);

  
