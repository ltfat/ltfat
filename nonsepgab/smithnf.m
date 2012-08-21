function [U,S,V] = smithnormalform(A);
%SMITHNF Smith normal form of an integer matrix
% 
%   `[U,S,V] = smithnf(A)` returns integer matrices `U`, `S`, and
%   `V` such that::
%
%     A = U*S*V',
%
%   `S` is diagonal and nonnegative, `S(i,i)` divides `S(i+1,i+1)` for
%   all `i`, $\det U = \pm 1$, and $det V =\pm 1$.
%
%   `S = smithnf(A)` just returns `diag(S)`.
% 
%   This function is in some ways analogous to SVD for integer matrices.
%
%   An example:::
%
%     [U,S,V] = smithnf([10 5; 0 10])
%
%   For more information, see `<http://en.wikipedia.org/wiki/Smith_normal_form>`_.
%
%   See also: hermitenf, metaplecop

%   Authors: John Gilbert, Ewa Matusiak

%  This function was originally written by John Gilbert, Xerox Corporation,
%  December 1993.

if round(A) ~= A
    error('Requires integer input.');
end

% This looks much like an SVD algorithm
% that first bidiagonalizes A by Givens rotations and then chases zeros,
% except for the construction of the 2 by 2 elementary transformation.

[m,n] = size(A);
S = A;
U = eye(m);
V = eye(n);

% Bidiagonalize S with elementary Hermite transforms.

for j = 1:min(m,n)
    % Zero column j below the diagonal.
    for i = j+1:m
        if S(i,j)
            % Construct an elementary Hermite transformation E
            % to zero S(i,j) by combining rows i and j.
            E = ehermite(S(j,j),S(i,j));
            % Apply the transform to S and U.
            S([j i],:) = E * S([j i],:);
            U(:,[j i]) = U(:,[j i]) / E;
        end;
    end;
    % Zero row j after the superdiagonal.
    for i = j+2:n
        if S(j,i)
            % Construct an elementary Hermite transformation E
            % to zero S(j,i) by combining columns j+1 and i.
            E = ehermite(S(j,j+1),S(j,i));
            % Apply the transform to S and V.
            S(:,[j+1 i]) = S(:,[j+1 i]) * E';
            V(:,[j+1 i]) = V(:,[j+1 i]) / E;
        end;
    end;
end;

% Now S is upper bidiagonal.
% Chase the superdiagonal nonzeros away.

D = diag(S,1);
while any(D)
    b = min(find(D));
    % Start chasing a bulge at the first nonzero superdiagonal element.
    
    % To guarantee reduction in S(b,b), first make S(b,b) positive
    % and make S(b,b+1) nonnegative and less than S(b,b).
    if S(b,b) < 0
        S(b,:) = -S(b,:);
        U(:,b) = -U(:,b);
    end;
    q = floor(S(b,b+1)/S(b,b));
    E = [1 0 ; -q 1];
    S(:,[b b+1]) = S(:,[b b+1]) * E';
    V(:,[b b+1]) = V(:,[b b+1]) / E;
    
    if S(b,b+1)
        
        % Zero the first nonzero superdiagonal element
        % using columns b and b+1, to start the bulge at S(b+1,b).
        E = ehermite(S(b,b),S(b,b+1));
        S(:,[b b+1]) = S(:,[b b+1]) * E';
        V(:,[b b+1]) = V(:,[b b+1]) / E;
        for j = 1:min(m,n)
            if j+1 <= m
                % Zero S(j+1,j) using rows j and j+1.
                E = ehermite(S(j,j),S(j+1,j));
                S([j j+1],:) = E * S([j j+1],:);
                U(:,[j j+1]) = U(:,[j j+1]) / E;
            end
            if j+2 <= n
                % Zero S(j,j+2) using columns j+1 and j+2.
                E = ehermite(S(j,j+1),S(j,j+2));
                S(:,[j+1 j+2]) = S(:,[j+1 j+2]) * E';
                V(:,[j+1 j+2]) = V(:,[j+1 j+2]) / E;
            end;
        end;
    end;
    D = diag(S,1);
end;

% Now S is diagonal. Make it nonnegative.

for j = 1:min(m,n)
    if S(j,j) < 0
        S(j,:) = -S(j,:);
        U(:,j) = -U(:,j);
    end;
end;

% Squeeze factors to the lower right to enforce the divisibility condition.

for i = 1 : min(m,n)
    for j = i+1 : min(m,n)
        % Replace S(i,i) and S(j,j) by their gcd and lcm respectively.
        a = S(i,i);
        b = S(j,j);
        [g,c,d] = gcd(a,b);
        E = [ 1 d ; -b/g a*c/g];
        F = [ c 1 ; -b*d/g a/g];
        S([i j],[i j]) = E * S([i j],[i j]) * F';
        U(:,[i j]) = U(:,[i j]) / E;
        V(:,[i j]) = V(:,[i j]) / F;
    end;
end;

U = round(U);
V = round(V);
if nargout <= 1
    U = diag(S);
end;


function E = ehermite(a,b);
% EHERMITE : Elementary Hermite tranformation.
%
% For integers a and b, E = ehermite(a,b) returns
% an integer matrix with determinant 1 such that E * [a;b] = [g;0],
% where g is the gcd of a and b.
%
% This function is in some ways analogous to GIVENS.

[g,c,d] = gcd(a,b);
if g
    E = [c d ; -b/g a/g];
else
    E = [1 0 ; 0 1];
end

