function [a,M,L,N,Ngood]=gabimagepars(Ls,x,y)
%GABIMAGEPARS  Find Gabor parameters to generate image
%   Usage: [a,M,L,N,Ngood]=gabimagepars(Ls,x,y);
%
%   `[a,M,L,N,Ngood]=gabimagepars(Ls,x,y)` will compute a reasonable set of
%   parameters *a*, *M* and *L* to produce a nice Gabor 'image' of a signal
%   of length *Ls*. The approximate number of pixels in the time direction is
%   given as *x* and the number of pixels in the frequency direction is given
%   as *y*.
%
%   The output parameter *Ngood* contains the number of time steps (columns
%   in the coefficients matrix) that contains relevant information. The
%   columns from *Ngood* until *N* only contains information from a
%   zero-extension of the signal.
%
%   If you use this function to calculate a grid size for analysis of a
%   real-valued signal (using |dgtreal|), please input twice of the desired
%   size *y*. This is because |dgtreal| only returns half as many
%   coefficients in the frequency direction as |dgt|.
%
%   An example: We wish to compute a Gabor image of a real valued signal *f*
%   of length $7500$. The image should have an approximate resolution of
%   $600 \times 800$ pixels:::
%
%     [f,fs]=linus; f=f(4001:4000+7500);
%     [a,M,L,N,Ngood] = gabimagepars(7500,800,2*600);
%     c = dgtreal(f,'gauss',a,M);
%     plotdgtreal(c,a,M,fs,90);
%
%   The size of c is $(M/2)+1 \times N$ equal to $601 \times 700$ pixels. 
%
%   For this function to work properly, the specified numbers for *x* and
%   *y* must not be large prime numbers.
%  
%   See also: dgt, dgtreal, sgram

if min(x,y)>Ls
  % Small values case, just do an STFT
  M=Ls;
  N=Ls;
  a=1;
  Ngood=N;
  L=Ls;
else

  % Set M and N to be what the user specified
  M=y;
  N=x;

  % Determine the minimum transform size.
  K=lcm(M,N);
    
  % This L is good, but is it not the same as DGT will choose.
  Llong=ceil(Ls/K)*K;
  
  % Fix a from the long L
  a=Llong/N;
  
  % Now we have fixed a and M, so we can use the standard method of choosing L
  Lsmallest=lcm(a,M);
  L=ceil(Ls/Lsmallest)*Lsmallest;
  
  % We did not get N as desired.
  N=L/a;
  
  % Number of columns to display
  Ngood=ceil(Ls/a);
  
  if M<=a
    error('LTFAT:noframe',['Cannot generate a frame, the signal is too long as compared ' ...
           'to the size of the image. Increase x and y.']);
  end;
  
end;
