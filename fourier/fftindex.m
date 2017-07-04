function n=fftindex(N,nyquistzero)
%FFTINDEX  Frequency index of FFT modulations
%   Usage: n=fftindex(N);
%
%   `fftindex(N)` returns the index of the frequencies of the standard FFT of
%   length *N* as they are ordered in the output from the `fft` routine. The
%   numbers returned are in the range `-ceil(N/2)+1:floor(N/2)`
%
%   `fftindex(N,0)` does as above, but sets the Nyquist frequency to zero.
%
%   See also: dft

%   AUTHOR : Peter L. Søndergaard
%   TESTING: OK
%   REFERENCE: OK

complainif_argnonotinrange(nargin,1,2,mfilename);

if nargin ==1
    if rem(N,2)==0
        n=[0:N/2,-N/2+1:-1].';
  else
      n=[0:(N-1)/2,-(N-1)/2:-1].';
  end;
else
    if rem(N,2)==0
        n=[0:N/2-1,0,-N/2+1:-1].';
  else
      n=[0:(N-1)/2,-(N-1)/2:-1].';
  end;
end;
