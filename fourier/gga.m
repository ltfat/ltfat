function c = gga(f,fvec,fs,dim)
%GGA Generalized Goertzel algorithm
%   Usage:  c = gga(x,indvec)
%
%   Input parameters:
%         x      : Input data.
%         fvec   : Indices to calculate. 
%         fs     : Sampling frequency.
%
%   Output parameters:
%         c      : Coefficient vector.
%
%   `c=gga(f,fvec)` computes the discrete-time fourier transform DTFT of
%   *f* at frequencies in `fvec` as $c(k)=F(2*pi*fvec(k))$ where
%   $F=DTFT(f)$, $k=1,\dots\K$ and `K=length(fvec)` using the generalized
%   second-order Goertzel algorithm. Thanks to the generalization, values
%   in `fvec` can be arbitrary numbers in range $0-1$ and not restricted to
%   $l/Ls$, $l=0,\dots\Ls-1$ (usual DFT samples) as the original Goertzel 
%   algorithm is. *Ls* is the length of the first non-singleton dimension
%   of *f*. If `indvec` is empty or ommited, `indvec` is assumed to be
%   `(0:Ls-1)/L` and results in the same output as `fft`.
%
%   `c=gga(f,fvec,fs)` computes the same with frequencies relative to *fs*.
%
%   The input *f* is processed along the first non-singleton dimension or
%   along dimension *dim* if specified.
%
%   **Remark:**
%   Besides the generalization the algorithm is also shortened by one
%   iteration compared to the conventional Goertzel.
%
%   Examples:
%   ---------
%   
%   Calculating DTFT samples of interest:::
% 
%     % Generate input signal
%     fs = 8000;
%     L = 2^10;
%     k = (0:L-1).';
%     freq = [400,510,620,680,825];
%     phase = [pi/4,-pi/4,-pi/8,pi/4,-pi/3];
%     amp = [5,3,4,1,2];
%     f = arrayfun(@(a,f,p) a*sin(2*pi*k*f/fs+p),...
%                  amp,freq,phase,'UniformOutput',0);
%     f = sum(cell2mat(f),2);
% 
%     % This is equal to fft(f)
%     ck = gga(f);
% 
%     %GGA to FFT error:
%     norm(ck-fft(f))
% 
%     % DTFT samples at 400,510,620,680,825 Hz
%     ckchzt = gga(f,freq,fs);
% 
%     % Plot modulus of coefficients
%     figure(1);
%     hold on;
%     stem(k/L*fs,abs(ck),'k');
%     stem(freq,abs(ckchzt),'r:');
%     set(gca,'XLim',[freq(1)-50,freq(end)+50]);
%     set(gca,'YLim',[0 3*L]);
%     xlabel('f[Hz]');
%     ylabel('|ck|');
% 
%     % Plot phase of coefficients
%     figure(2);
%     hold on;
%     stem(k/L*fs,angle(ck),'k');
%     stem(freq,angle(ckchzt),'r:');
%     set(gca,'XLim',[freq(1)-50,freq(end)+50]);
%     set(gca,'YLim',[-pi pi]);
%     xlabel('f[Hz]');
%     ylabel('angle(ck)');
%
%   See also: chirpzt
%
%   References: syra2012goertzel
       
% The original copyright goes to
% 2013 Pavel Rajmic, Brno University of Technology, Czech Rep.


%% Check the input arguments
if nargin < 1
    error('%s: Not enough input arguments.',upper(mfilename))
end

if isempty(f)
    error('%s: X must be a nonempty vector or a matrix.',upper(mfilename))
end

if nargin<4
  dim=[];  
end;

if nargin<3 || isempty(fs)
  fs=1;  
end;

[f,~,Ls,~,dim,permutedsize,order]=assert_sigreshape_pre(f,[],dim,'GGA');

if nargin > 1 && ~isempty(fvec)
   if ~isreal(fvec) || ~isvector(fvec)
      error('%s: INDVEC must be a real vector.',upper(mfilename))
   end
else
   fvec = (0:Ls-1)/Ls;
end

c = comp_gga(f,fvec/fs*Ls);

permutedsize(1)=numel(fvec);

c=assert_sigreshape_post(c,dim,permutedsize,order);

