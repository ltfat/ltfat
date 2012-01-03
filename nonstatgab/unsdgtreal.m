function [c,Ls] = nsdgtreal(f,g,a,M)
%NSDGTREAL  Nonsationary Discrete Gabor transform.
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
%   NSDGTREAL(f,g,a,M) computes the nonstationary Gabor coefficients of the 
%   input signal f. The signal f can be a multichannel signal, given in the
%   form of a 2D matrix of size Ls x W, with Ls the signal length and W the
%   number of signal channels.
%
%   The nonstationnary Gabor theory extends standard Gabor theory by 
%   enabling the evolution of the window over time. It is therefor
%   necessary to specify a set of windows instead of a single window. 
%   This is done by using a cell array for g. In this cell array, the nth 
%   element g{n} is a row vetor specifying the nth window.
%
%   The resulting coefficients also require a storage in a cell array, as 
%   the number of frequency channels is not constant over time. More 
%   precisely, the nth cell of c, c{n}, is a 2D matrix of size M(n) x W 
%   and containing the complex local spectra of the signal channels 
%   windowed by the nth window g{n} shifted in time at position timepos(n). 
%   c{n}(m,l) is thus the value of the coefficient for time index n, 
%   frequency index m and signal channel l.
%
%   The variable _a contains the distance in samples between two
%   consequtive blocks of coefficients. The variable M contains the
%   number of channels for each block of coefficients. Both _a and M are
%   vectors of integers.
%
%   The variables g, a and M must have the same length, and the result c 
%   will also have the same length.
%   
%   The time positions of the coefficients blocks can be obtained by the
%   following code. A value of 0 correspond to the first sample of the
%   signal:
%
%C     timepos = cumsum(a)-a(1);
%
%   [c,Ls]=nsdgtreal(f,g,a,M) additionally returns the length Ls of the input 
%   signal f. This is handy for reconstruction:
%
%      [c,Ls]=nsdgtreal(f,g,a,M);
%      fr=insdgtreal(c,gd,a,Ls);
%
%   will reconstruct the signal f no matter what the length of f is, 
%   provided that gd are dual windows of g.
%
%   Notes: 
%   nsdgtreal uses circular border conditions, that is to say that the signal is
%   considered as periodic for windows overlapping the beginning or the 
%   end of the signal.
%
%   The phaselocking convention used in NSDGTREAL is different from the
%   convention used in the DGT function. NSDGTREAL results are phaselocked (a
%   phase reference moving with the window is used), whereas DGT results are
%   not phaselocked (a fixed phase reference corresponding to time 0 of the
%   signal is used). See the help on PHASELOCK for more details on
%   phaselocking conventions.
%
%   See also:  nsdgt, insdgtreal, nsgabdual, nsgabtight, phaselock
%
%   Demos:  demo_nsdgtreal
%
%R  ltfatnote010
  
%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGTREAL
%   REFERENCE: 

timepos=cumsum(a)-a(1);
  
Ls=length(f);

N=length(a); % Number of time positions

W=size(f,2); % Number of signal channels

c=cell(N,1); % Initialisation of the result

for ii=1:N
  shift=floor(length(g{ii})/2);
  temp=zeros(M(ii),W);
  
  % Windowing of the signal.
  % Possible improvements: The following could be computed faster by 
  % explicitely computing the indexes instead of using modulo and the 
  % repmat is not needed if the number of signal channels W=1 (but the time 
  % difference when removing it whould be really small)
  temp(1:length(g{ii}))=f(mod((1:length(g{ii}))+timepos(ii)-shift-1,Ls)+1,:).*...
    repmat(conj(circshift(g{ii},shift)),1,W);
  
  temp=circshift(temp,-shift);
  if M(ii)<length(g{ii}) 
    % Fft size is smaller than window length, some aliasing is needed
    x=floor(length(g{ii})/M(ii));
    y=length(g{ii})-x*M(ii);
    % Possible improvements: the following could probably be computed 
    % faster using matrix manipulation (reshape, sum...)
    temp1=temp;
    temp=zeros(M(ii),size(temp,2));
    for jj=0:x-1
      temp=temp+temp1(jj*M(ii)+(1:M(ii)),:);
    end
    temp(1:y,:)=temp(1:y,:)+temp1(x*M(ii)+(1:y),:);
  end
  
  c{ii}=fftreal(temp); % FFT of the windowed signal
end


%OLDFORMAT
