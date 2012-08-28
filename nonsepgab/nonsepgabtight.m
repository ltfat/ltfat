function gt=nonsepgabtight(p1,p2,p3,p4,varargin)
%NONSEPGABTIGHT  Canonical tight window for non-separable lattices
%   Usage:  gt=nonsepgabtight(a,M,lt,L);
%           gt=nonsepgabtight(g,a,M,lt);
%           gt=nonsepgabtight(g,a,M,lt,L);
%
%   Input parameters:
%         g     : Gabor window.
%         a     : Length of time shift.
%         M     : Number of modulations.
%         lt    : Lattice type
%         L     : Length of window. (optional)
%   Output parameters:
%         gt    : Canonical tight window, column vector.
%
%   `nonsepgabtight(a,M,L,lt)` computes a nice tight window of length *L* for
%   a lattice with parameters *a*, *M* and *lt*. The window is not an FIR
%   window, meaning that it will only generate a tight system if the system
%   length is equal to *L*.
%
%   `nonsepgabtight(g,a,M,lt)` computes the canonical tight window of the Gabor
%   frame with window *g* and parameters *a*, *M* and *lt*.
%
%   The window *g* may be a vector of numerical values, a text string or a
%   cell array. See the help of |gabwin|_ for more details.
%  
%   If the length of *g* is equal to *M*, then the input window is assumed to
%   be a FIR window. In this case, the canonical dual window also has
%   length of *M*. Otherwise the smallest possible transform length is
%   chosen as the window length.
%
%   `nonsepgabtight(g,a,M,lt,L)` returns a window that is tight for a system of
%   length *L*. Unless the input window *g* is a FIR window, the returned
%   tight window will have length *L*.
%
%   If $a>M$ then an orthonormal window of the Gabor Riesz sequence with
%   window *g* and parameters *a* and *M* will be calculated.
%
%   See also:  gabdual, gabwin, fir2long, dgt

%   AUTHOR : Peter Soendergaard.
%   TESTING: TEST_DGT
%   REFERENCE: OK

%% ------------ decode input parameters ------------
  
if numel(p1)==1
  % First argument is a scalar.

  error(nargchk(4,4,nargin));

  a=p1;
  M=p2;
  lt=p3;
  varargin=[L,varargin];

  g='gauss';
  
else    
  % First argument assumed to be a vector.
      
  g=p1;
  a=p2;
  M=p3;
  lt=p4;
    
end;

definput.keyvals.L=[];
definput.flags.nsalg={'multiwin','smith','shear'};
[flags,kv,L]=ltfatarghelper({'L'},definput,varargin);

[g,L,info] = nonsepgabpars_from_window(g,a,M,lt,L);

% -------- Are we in the Riesz sequence of in the frame case

scale=1;
if a>M
  % Handle the Riesz basis (dual lattice) case.
  % Swap a and M, and scale differently.
  scale=sqrt(a/M);
  tmp=a;
  a=M;
  M=tmp;
end;

%% -------- Compute ------------- 

% Just in case, otherwise the call is harmless. 
g=fir2long(g,L);

if flags.do_multiwin

    mwin=comp_nonsepwin2multi(g,a,M,lt);
    
    gtfull=comp_gabtight_long(mwin,a*lt(2),M)*scale;
    
    % We need just the first vector
    gt=gtfull(:,1);
    
end;

if flags.do_smith

    % Convert the lattice type representation to matrix type representation
    latm = latticetype2matrix(L,a,M,lt);
    
    % and convert this representation to Smith normal form
    [U,S,V] = smithnf(latm);
    
    % Convert the S matrix back to lattice type, for use with 'gabdual'
    [a0,M0,~]=matrix2latticetype(L,S);
    
    g0 = metaplecop(g,U,'inv');
    gt0 = gabtight(g0,a0,M0);
    gt  = metaplecop(gt0,U);    
    
end;

if flags.do_shear
    b=L/M;
    s=b*lt(1)/lt(2);
    
    [s0,s1,X] = shearfind(a,b,s,L);
    
    if s0 ~= 0
        g = ifft(pchirp(L,-s0).*fft(g));
    end
    
    if s1 ~= 0
        g = pchirp(L,s1).*g;
    end
    
    Mr = L/X;
    ar = a*b/X;
    
    gt=gabtight(g,ar,Mr,L);
    
    if s1 ~= 0
        gt = pchirp(L,-s1).*gt;
    end

    if s0 ~= 0
        gt = ifft(pchirp(L,s0).*fft(gt));
    end
    
end;


%% --------- post process result -------
      
if info.wasrow
  gt=gt.';
end;

if isreal(g) && (lt(2)==1 || lt(2)==2)
  % If g is real and the lattice is either rectangular or quinqux, then
  % the output is known to be real.
  gt=real(gt);
end;