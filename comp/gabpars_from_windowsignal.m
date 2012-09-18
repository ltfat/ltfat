function [f,g,L,Ls,W,info] = gabpars_from_windowsignal(f,g,a,M,L,lt,callfun)
%GABPARS_FROM_WINDOWSIGNAL  Compute g and L from window and signal
%   Usage: [g,g.info,L] = gabpars_from_windowsignal(f,g,a,M,L);
%
%   Use this function if you know an input signal, a window and a lattice
%   for the DGT. The function will calculate a transform length L and
%   evaluate the window g into numerical form. The signal will be padded and
%   returned as a column vector.
%
%   If the transform length is unknown (as it usually is unless explicitly
%   specified by the user), set L to be [] in the input to this function.
  
if nargin<6
  stacknames=dbstack;  
  callfun=stacknames(2).name;
end;

% ----- step 1 : Verify f and determine its length -------
% Change f to correct shape.
[f,Ls,W,wasrow,remembershape]=comp_sigreshape_pre(f,callfun,0);


if isempty(L)

    % ----- step 2b : Verify a, M and get L from the signal length f----------
    L=dgtlength(Ls,a,M);

else

    % ----- step 2a : Verify a, M and get L
    Luser=dgtlength(L,a,M);
    if Luser~=L
        error(['%s: Incorrect transform length L=%i specified. Next valid length ' ...
               'is L=%i.'],callfun,L,Luser)
    end;

end;

% ----- step 3 : Determine the window 

[g,info]=gabwin(g,a,M,L,'callfun',callfun);

if L<info.gl
  error('%s: Window is too long.',callfun);
end;

% ----- final cleanup ---------------

f=postpad(f,L);

% If the signal is single precision, make the window single precision as
% well to avoid mismatches.
if isa(f,'single')
  g=single(g);
end;




