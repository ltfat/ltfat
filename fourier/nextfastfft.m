function [nfft,tableout]=nextfastfft(n)
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
%   The algorithm will look up the best size in a table, which is
%   computed the first time the function is run. If the input size is
%   larger than the largest value in the table, the input size will be
%   reduced by factors of 2, until it is in range.
%
%   [n,nfft]=NEXTFASTFFT(n) additionally returns the table used for lookup.
%
%   Demos: demo_nextfastfft
  
%   AUTHOR: Peter L. Soendergaard and Johan Sebastian Rosenkilde Nielsen
  
  
  persistent table;
  
  maxval=2^20;

  if isempty(table)
    % Compute the table for the first time, it is empty.
    l2=log(2);
    l3=log(3);
    l5=log(5);
    l7=log(7);
    tables=zeros(1286,1);
    ii=1;
    for i2=0:floor(log(maxval)/l2)
      for i3=0:floor((log(maxval)-i2*l2)/l3)
        for i5=0:floor((log(maxval)-i2*l2-i3*l3)/l5)
          for i7=0:floor((log(maxval)-i2*l2-i3*l3-i5*l5)/l7)
            table(ii)=2^i2*3^i3*5^i5*7^i7;
            ii=ii+1;
          end;
        end;
      end;
    end;
    table=sort(table);
  end;

  % Copy input to output. This allows us to efficiently work in-place.
  nfft=n;
    
  % Handle input of any shape by Fortran indexing.
  for ii=1:numel(n)
    n2reduce=0;
    
    if n(ii)>maxval
      % Reduce by factors of 2 to get below maxval
      n2reduce=ceil(log2(nfft(ii)/maxval));
      nfft(ii)=nfft(ii)/2^n2reduce;
    end;
    
    % Use a simple bisection method to find the answer in the table.
    from=1;
    to=numel(table);
    while from<=to
      mid = round((from + to)/2);    
      diff = table(mid)-nfft(ii);
      if diff<0
        from=mid+1;
      else
        to=mid-1;                       
      end
    end
    nfft(ii)=table(from);

    % Add back the missing factors of 2 (if any)
    nfft(ii)=nfft(ii)*2^n2reduce;

  end;
  
  tableout=table;