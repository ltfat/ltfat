function [gamma]=wildual(g,M,L)
%WILDUAL  Wilson dual window
%   Usage:  gamma=wildual(g,M);
%           gamma=wildual(g,M,L);
%
%   Input parameters:
%         g     : Gabor window.
%         M     : Number of modulations.
%         L     : Length of window. (optional)
%   Output parameters:
%         gamma : Canonical dual window.
%
%   `wildual(g,M)` returns the dual window of the Wilson or WMDCT basis with
%   window *g*, parameter *M* and length equal to the length of the window *g*.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |wilwin| for more details.
%
%   If the length of *g* is equal to $2\cdot M$ then the input window is
%   assumed to be an FIR window. In this case, the dual window also has
%   length of $2\cdot M$. Otherwise the smallest possible transform length is
%   chosen as the window length.
%
%   `wildual(g,M,L)` does the same, but now *L* is used as the length of the
%   Wilson basis.
%
%   The input window *g* must be real and whole-point even. If *g* is not
%   whole-point even, then reconstruction using the dual window will not be
%   perfect. For a random window *g*, the window closest to *g* that satisfies
%   these restrictions can be found by ::
%
%     g_wpe = real(peven(g));
%
%   All windows in the toolbox satisfies these restrictions unless
%   clearly stated otherwise.
%
%   See also:  dwilt, wilwin, wmdct, wilorth, isevenfunction

%   AUTHOR : Peter L. Søndergaard.
%   TESTING: TEST_DWILT
%   REFERENCE: OK

complainif_argnonotinrange(nargin,2,3,mfilename);

if nargin==2
    L=[];
end;

%% ------ step 2: Verify a, M and L
if isempty(L)
    if isnumeric(g)
        % Use the window length
        Ls=length(g);
    else
        % Use the smallest possible length
        Ls=1;
    end;

    % ----- step 2b : Verify M and get L from the window length ----------
    L=dwiltlength(Ls,M);

else

    % ----- step 2a : Verify M and get L

    Luser=dwiltlength(L,M);
    if Luser~=L
        error(['%s: Incorrect transform length L=%i specified. Next valid length ' ...
               'is L=%i. See the help of DWILTLENGTH for the requirements.'],...
              upper(mfilename),L,Luser);
    end;

end;


%% ----- step 3 : Determine the window 

[g,info]=wilwin(g,M,L,upper(mfilename));

if L<info.gl
  error('%s: Window is too long.',upper(mfilename));
end;

%% ----- call gabdual ----------------
a=M;

g=fir2long(g,L);
gamma=2*comp_gabdual_long(g,a,2*M);
