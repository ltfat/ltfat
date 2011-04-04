function nfft=nextfastfft(n,varargin)
%NEXTNICEFFT  Next higher number with a fast FFT
%   Usage: nfft=nextfastfft(n);
%
%   NEXTFASTFFT(n) will return the next number greater than or equal to
%   n, for which the computation of an FFT is fast. Such a number is
%   solely comprised of small prime-factors.
%
%   NEXTFASTFFT is intended as a replacement of NEXTPOW2, which if often
%   used for the same purpose. However, a modern FFT implementation (like FFTW)
%   usually performs well for size which are powers or 2,3,5 and 7.
%  
  
  persistent table;
  
  maxval=1e6;

  if isempty(table)
    for i2=0:floor(log(maxval)/log(2))
      for i3=0:floor(log(maxval)/log(3))
        for i5=0:floor(log(maxval)/log(5))
          for i7=0:floor(log(maxval)/log(7))
            candidate=2^i2*3^i3*5^i5*7^i7;
            if candidate<=maxval
              table=[table;candidate];
            end;
          end;
        end;
      end;
    end;
    table=sort(table);
  end;

  nfft=zeros(size(n));
  % Handle input of any shape by Fortran indexing.
  for ii=1:numel(n)    
    ind=find(table>=n(ii));
    nfft(ii)=table(ind(1));  
  end;