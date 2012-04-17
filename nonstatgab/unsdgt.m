function [c,Ls] = unsdgt(f,g,a,M)
%UNSDGT  Uniform Nonstationary Discrete Gabor transform
%   Usage:  c=unsdgt(f,g,a,M);
%           [c,Ls]=unsdgt(f,g,a,M);
%
%   Input parameters:
%         f     : Input signal.
%         g     : Cell array of window functions.
%         a     : Vector of time positions of windows.
%         M     : Numbers of frequency channels.
%   Output parameters:
%         c     : Cell array of coefficients.
%         Ls    : Length of input signal.
%
%   `unsdgt(f,g,a,M)` computes the uniform nonstationary Gabor coefficients
%   of the input signal *f*. The signal f can be a multichannel signal,
%   given in the form of a 2D matrix of size $Ls \times W$, with *Ls* being
%   the signal length and *W* the number of signal channels.
%
%   The non-stationary Gabor theory extends standard Gabor theory by
%   enabling the evolution of the window over time. It is therefor necessary
%   to specify a set of windows instead of a single window.  This is done by
%   using a cell array for *g*. In this cell array, the n'th element `g{n}`
%   is a row vector specifying the n'th window. However, the uniformity
%   means that the number of channels is fixed.
%
%   The resulting coefficients is stored as a $M \times N \times W$
%   array. `c(m,n,l)` is thus the value of the coefficient for time index *n*,
%   frequency index *m* and signal channel *l*.
%
%   The variable *a* contains the distance in samples between two consequtive
%   blocks of coefficients. *a* is a vectors of integers. The variables *g* and
%   *a* must have the same length.
%   
%   The time positions of the coefficients blocks can be obtained by the
%   following code. A value of 0 correspond to the first sample of the
%   signal::
%
%     timepos = cumsum(a)-a(1);
%
%   `[c,Ls]=nsdgt(f,g,a,M)` additionally returns the length *Ls* of the input 
%   signal *f*. This is handy for reconstruction::
%
%     [c,Ls]=unsdgt(f,g,a,M);
%     fr=iunsdgt(c,gd,a,Ls);
%
%   will reconstruct the signal *f* no matter what the length of *f* is, 
%   provided that *gd* are dual windows of *g*.
%
%   Notes:
%   ------
%
%   `unsdgt` uses circular border conditions, that is to say that the signal is
%   considered as periodic for windows overlapping the beginning or the 
%   end of the signal.
%
%   The phaselocking convention used in `unsdgt` is different from the
%   convention used in the |dgt|_ function. `unsdgt` results are phaselocked
%   (a phase reference moving with the window is used), whereas |dgt|_ results
%   are not phaselocked (a fixed phase reference corresponding to time 0 of
%   the signal is used). See the help on |phaselock|_ for more details on
%   phaselocking conventions.
%
%   See also:  insdgt, nsgabdual, nsgabtight, phaselock
%
%   Demos:  demo_nsdgt
%
%   References: ltfatnote010
  
%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGT
%   REFERENCE: 

if ~isnumeric(a)
  error('%s: a must be numeric.',upper(callfun));
end;

if ~isnumeric(M)
  error('%s: M must be numeric.',upper(callfun));
end;

L=sum(a);

[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,'UNSDGT',0);
f=postpad(f,L);

[g,info]=nsgabwin(g,a,M);

timepos=cumsum(a)-a(1);
  
N=length(a); % Number of time positions

c=zeros(M,N,W); % Initialisation of the result

for ii=1:N
  shift=floor(length(g{ii})/2);
  temp=zeros(M,W);
  
  % Windowing of the signal.
  % Possible improvements: The following could be computed faster by 
  % explicitely computing the indexes instead of using modulo and the 
  % repmat is not needed if the number of signal channels W=1 (but the time 
  % difference when removing it whould be really small)
  temp(1:length(g{ii}))=f(mod((1:length(g{ii}))+timepos(ii)-shift-1,L)+1,:).*...
    repmat(conj(circshift(g{ii},shift)),1,W);
  
  temp=circshift(temp,-shift);
  if M<length(g{ii}) 
    % Fft size is smaller than window length, some aliasing is needed
    x=floor(length(g{ii})/M);
    y=length(g{ii})-x*M;
    % Possible improvements: the following could probably be computed 
    % faster using matrix manipulation (reshape, sum...)
    temp1=temp;
    temp=zeros(M,size(temp,2));
    for jj=0:x-1
      temp=temp+temp1(jj*M+(1:M),:);
    end
    temp(1:y,:)=temp(1:y,:)+temp1(x*M+(1:y),:);
  end
  
  % FFT of the windowed signal
  c(:,ii,:) = reshape(fft(temp),M,1,W); 
end
