function h=pheaviside(L)
%PHEAVISIDE  Periodic Heaviside function
%   Usage: h=pheaviside(L);
%
%   `pheaviside(L)` returns a periodic Heaviside function. The periodic
%   Heaviside function takes on the value 1 for indices corresponding to
%   positive frequencies, 0 corresponding to negative frequencies and the
%   value .5 for the zero and Nyquist frequencies.
%
%   To get a function that weights the negative frequencies by 1 and the
%   positive by 0, use `involute(pheaviside(L))`
%
%   As an example, the `pheaviside` function can be use to calculate the
%   Hilbert transform for a column vector *f*::
%
%     h=2*ifft(fft(f).*pheaviside(length(f)));
%
%   See also: middlepad, involute, fftindex

%   AUTHOR : Peter L. Søndergaard.
%   REFERENCE: OK
%   TESTING: OK

complainif_argnonotinrange(nargin,1,1,mfilename);

h=zeros(L,1);
if L>0
    % First term is .5
    h(1)=.5;

    % Set positive frequencies to 1.
    h(2:ceil(L/2))=1;

    % Last term (Nyquist frequency) is also .5, if it exists.
    if rem(L,2)==0
        h(L/2+1)=.5;
  end;
end;

