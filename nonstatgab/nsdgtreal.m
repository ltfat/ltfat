function [c,Ls] = nsdgtreal(f,g,a,M)
%NSDGTREAL  Non-stationary Discrete Gabor transform for real valued signals
%   Usage:  c=nsdgtreal(f,g,a,M);
%           [c,Ls]=nsdgtreal(f,g,a,M);
%
%   Input parameters:
%         f     : Input signal.
%         g     : Cell array of window functions.
%         a     : Vector of time positions of windows.
%         M     : Vector of numbers of frequency channels.
%   Output parameters:
%         c     : Cell array of coefficients.
%         Ls    : Length of input signal.
%
%   `nsdgtreal(f,g,a,M)` computes the non-stationary Gabor coefficients of the
%   input signal *f*. The signal *f* can be a multichannel signal, given in
%   the form of a 2D matrix of size $Ls \times W$, with *Ls* the signal
%   length and *W* the number of signal channels.
%
%   As opposed to |nsdgt| only the coefficients of the positive frequencies
%   of the output are returned. `nsdgtreal` will refuse to work for complex
%   valued input signals.
%
%   The non-stationary Gabor theory extends standard Gabor theory by
%   enabling the evolution of the window over time. It is therefor necessary
%   to specify a set of windows instead of a single window.  This is done by
%   using a cell array for *g*. In this cell array, the n'th element `g{n}`
%   is a row vector specifying the n'th window.
%
%   The resulting coefficients also require a storage in a cell array, as
%   the number of frequency channels is not constant over time. More
%   precisely, the n'th cell of *c*, `c{n}`, is a 2D matrix of size
%   $M(n)/2+1 \times W$ and containing the complex local spectra of the
%   signal channels windowed by the n'th window `g{n}` shifted in time at
%   position $a(n)$.  `c{n}(m,l)` is thus the value of the coefficient for
%   time index *n*, frequency index *m* and signal channel *l*.
%
%   The variable *a* contains the distance in samples between two
%   consequtive blocks of coefficients. The variable *M* contains the
%   number of channels for each block of coefficients. Both *a* and *M* are
%   vectors of integers.
%
%   The variables *g*, *a* and *M* must have the same length, and the result *c*
%   will also have the same length.
%   
%   The time positions of the coefficients blocks can be obtained by the
%   following code. A value of 0 correspond to the first sample of the
%   signal::
%
%     timepos = cumsum(a)-a(1);
%
%   `[c,Ls]=nsdgtreal(f,g,a,M)` additionally returns the length *Ls* of the input 
%   signal *f*. This is handy for reconstruction::
%
%     [c,Ls]=nsdgtreal(f,g,a,M);
%     fr=insdgtreal(c,gd,a,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is, 
%   provided that *gd* are dual windows of *g*.
%
%   Notes:
%   ------
%
%   `nsdgtreal` uses circular border conditions, that is to say that the signal is
%   considered as periodic for windows overlapping the beginning or the 
%   end of the signal.
%
%   The phaselocking convention used in `nsdgtreal` is different from the
%   convention used in the |dgt| function. `nsdgtreal` results are phaselocked (a
%   phase reference moving with the window is used), whereas |dgt| results are
%   not phaselocked (a fixed phase reference corresponding to time 0 of the
%   signal is used). See the help on |phaselock| for more details on
%   phaselocking conventions.
%
%   See also:  nsdgt, insdgtreal, nsgabdual, nsgabtight, phaselock
%
%   Demos:  demo_nsdgt
%
%   References: ltfatnote018
  
%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGTREAL
%   REFERENCE: 

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(mfilename));
end;

if ~isnumeric(M)
  error('%s: M must be numeric.',upper(mfilename));
end;


%% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,upper(mfilename),0);

L=nsdgtlength(Ls,a);
f=postpad(f,L);

[g,info]=nsgabwin(g,a,M);

if L<info.gl
  error('%s: Window is too long.',upper(mfilename));
end;

timepos=cumsum(a)-a(1);

N=length(a); % Number of time positions

c=cell(N,1); % Initialisation of the result
   
for ii = 1:N
    Lg = length(g{ii});
    gt = g{ii}; gt = gt([end-floor(Lg/2)+1:end,1:ceil(Lg/2)]);
    win_range = mod(timepos(ii)+(-floor(Lg/2):ceil(Lg/2)-1),L)+1;
    if M(ii) < Lg 
        % if the number of frequency channels is too small, aliasing is introduced
        col = ceil(Lg/M(ii));
        temp = zeros(col*M(ii),W,assert_classname(f,g{1}));
        temp([col*M(ii)-floor(Lg/2)+1:end,1:ceil(Lg/2)],:) = bsxfun(@times,f(win_range,:),gt);
        temp = reshape(temp,M(ii),col,W);
        
        c{ii}=squeeze(fftreal(sum(temp,2)));
    else
        temp = zeros(M(ii),W,assert_classname(f,g{1}));
        temp([end-floor(Lg/2)+1:end,1:ceil(Lg/2)],:) = bsxfun(@times,f(win_range,:),gt);
        c{ii} = fftreal(temp);
    end       
end

if 0
for ii=1:N
  shift=floor(length(g{ii})/2);
  temp=zeros(M(ii),W,assert_classname(f,g{1}));
  
  % Windowing of the signal.
  % Possible improvements: The following could be computed faster by 
  % explicitely computing the indexes instead of using modulo and the 
  % repmat is not needed if the number of signal channels W=1 (but the time 
  % difference when removing it whould be really small)
  temp(1:length(g{ii}))=f(mod((1:length(g{ii}))+timepos(ii)-shift-1,L)+1,:).*...
    repmat(conj(circshift(g{ii},shift)),1,W);
  
  temp=circshift(temp,-shift);
  if M(ii)<length(g{ii}) 
    % Fft size is smaller than window length, some aliasing is needed
    x=floor(length(g{ii})/M(ii));
    y=length(g{ii})-x*M(ii);
    % Possible improvements: the following could probably be computed 
    % faster using matrix manipulation (reshape, sum...)
    temp1=temp;
    temp=zeros(M(ii),size(temp,2),assert_classname(f,g{1}));
    for jj=0:x-1
      temp=temp+temp1(jj*M(ii)+(1:M(ii)),:);
    end
    temp(1:y,:)=temp(1:y,:)+temp1(x*M(ii)+(1:y),:);
  end
  
  c{ii}=fftreal(temp); % FFT of the windowed signal
end
end;
