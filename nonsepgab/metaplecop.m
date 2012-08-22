function h = metaplecop(f,A,varargin)
%METAPLECOP  Metaplectic operator
%   Usage: h = metaplecop(f,A);
%
%   Input parameters:  
%          f     : Column vector of length *L*.
%          A     : Symplectic $2\times 2$ matrix over $\mathbb{Z}_L$.
%  
%   Output parameters: 
%          h     : The metaplectic operator `U`, associated to the
%                  symplectic matrix `A`, applied to the vector `f`.
%
%   `h=metaplecop(f,A)` returns the result `h` of applying the
%   metaplectic operator *U* associated to `A` to the vector `f`. `A` must
%   be a $2\times 2$ matrix `A` with integer entries.
%
%   `metaplec(f,A,'inv')` instead applies the inverse operator of *U* to `f`.
%
%   To accomplish this, the first step computes the Weyl decomposition of `A`,
%   representing `A` as the product of 6 matrices of the forms
%   
%   ..    S_c = [1 0; c 1]
%         F   = [0 1; -1 0]
%         D_a = [a 0, 0 a^{-1}]
%
%   .. math:: S_c = \left(\begin{array}{cc}
%     1&0\\c&1\end{array}\right),
%     \quad F = \left(\begin{array}{cc} 0&1\\ -1&0 \end{array}\right), \quad
%     D_a = \left(\begin{array}{cc}
%     a &  0 \\ 
%     0 &  a^{-1}
%     \end{array}\right),
%
%   i.e. shear matrices, standard symplectic form/flip matrix and
%   dilation/permutation or their inverses. Note the $a$ in $D_a$ must be
%   invertible over $\mathbb{Z}_L$, where $L$ is the length of the input
%   signal *f*.
%  
%   Each of these matrices corresponds to an elementary metaplectic operator
%   and consequently the operator *U* associated to `A` can be written as the
%   composition of 6 elementary metaplectic operators. The correspondences
%   are
%  
%     * $S_c$: Multiplication by a chirp with steepness $c$.
%           
%     * $F$:   Fourier transform.
%      
%     * $D_a$: Dilation by $a$, with $a$ invertible over $\mathbb{Z}_L$.
%  
%   Note that every $2\times 2$ matrix with determinant `1` is a symplectic
%   matrix (not true in general for larger matrices). Therefore, if
%   `M=A*D*V`, with `D` being the smith normal form of `M`, then `A` is
%   symplectic.
%  
%   For more details see Feichtinger et al. 2008. This code was originally
%   part of the Metaplectic Toolbox compiled by E. Matusiak (NuHAG)
%   
%   See also: smithnf
% 
%   References: feichtinger2008metaplectic
    
% Authors : Ewa Matusiak and Nicki Holighaus 21.08.12
    
%Check input
  
  if nargin<2
      error('%s: Too few input parameters.',upper(mfilename));
  end;
    
  definput.flags.oppower={'forward','inv'};
  [flags,kv]=ltfatarghelper({},definput,varargin);
  
  if size(A,1) ~= size(A,2) || size(A,1) ~= 2 || sum(sum(A ~= round(A)))>0
      error('A must be an integer valued 2x2 matrix');
  elseif round(det(A)) ~= 1 % rounding errors might occur in det(A)
      error('A is not symplectic, det(A) != 1');
  end
    
  if ~isvector(f)
      error('f must be a vector.');
  end
  
  L = length(f);
  
  % Determine the Weyl decomposition of A. Note that the matrices themselves
  % are unimportant, only the shear and dilation parameters must be
  % determined.
  
  a = A(1,1);
  b = A(1,2);
  c = A(2,1);
  d = A(2,2);
  t = theta(a,L);
  
  a0 = mod(a + t*b,L);    % Invertible element equivalent to a
  c0 = mod(c + t*d,L);    % c has to be changed accordingly
  
  % Compute the inverse of a0 
  [~,c] = gcd(a0,L);
  a0_inv = mod(c,L);
  
  % Determine shear parameters 
  ca = mod(c0*a0_inv,L);  
  ab = mod(a0_inv*b,L);
  
  % Determine dilation parameters
  sub = mod((0:a0_inv:a0_inv*(L-1)).',L) + 1;
    
  if flags.do_forward
      % Compute chirps (corresponding to the shear parameters)
      c_t = pchirp(L,-t);
      c_ab = pchirp(L,-ab);
      c_ca = pchirp(L,ca);

      % Apply metaplectic operator to v
      f = c_t .* f;
      f = fft(f);
      f = c_ab .* f;
      f = ifft(f);
      f = f(sub);
      h = c_ca .* f;
      
  else
      
      % Compute chirps (corresponding to the shear parameters)
      c_t = pchirp(L,t);
      c_ab = pchirp(L,ab);
      c_ca = pchirp(L,-ca);
      
      % Apply inverse metaplectic operator to v
      f = f .* c_ca;
      f(sub) = f;
      f = fft(f);
      f = f .* c_ab;
      f = ifft(f);
      h = f .* c_t;
  end
  
end

function t = theta(a,L)

% THETA.M - Ewa Matusiak, edited by Nicki Holighaus 21.08.12
%
% t = theta(a,L)
%
% Find theta as in Lemma 3 in Feichtinger et al. 2008
%
% This theta is necessary to determine from a=A(1,1) in the matrix A an 
% equivalent invertible element a0 = mod(a + theta*b,L)
%
    p = unique(factor(L));
    sub = find(mod(a,p));
    t = prod(p(sub));
    t = mod(t,L);
    
end