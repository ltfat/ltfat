function outsig=framet(insig,F);
%FRAMET  Frame transform
%   Usage: c=framet(f,F);
%
%   `c=framet(f,F)` computes the frame coefficients *c* of the input
%   signal *f* using the frame *F*. The frame object *F* must have been
%   created using |newframe|_.
%
%   If *f* is a matrix, the transform will be applied along the columns
%   of *f*. If *f* is an N-D array, the transform will be applied along
%   the first non-singleton dimension.
%
%   The output coefficients are stored as columns. This is usually
%   **not** the same format as the 'native' format of the frame. As an
%   examples, the output from |framet|_ for a gabor frame cannot be
%   passed to |idgt|_ without a reshape.
%
%   See also: newframe, iframet, plotframe
  
switch(F.type)
 case 'dgt'
  outsig=dgt(insig,F.g,F.a,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 case 'dgtreal'
  outsig=dgtreal(insig,F.g,F.a,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 case 'dwilt'
  outsig=dwilt(insig,F.g,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
 case 'wmdct'
  outsig=wmdct(insig,F.g,F.M);
  [M,N,W]=size(outsig);
  outsig=reshape(outsig,[M*N,W]);
end;

  