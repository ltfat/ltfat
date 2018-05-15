function x = ref_pcg(A, b, x, tol)
if nargin < 3
    x = b;
    tol = 1e-10;
end
if nargin < 4
    tol = 1e-10;
end

if isa(A,'function_handle')
r = b - A(x);
else
r = b - A * x;
end
p = r;
rsold = norm(r)^2;
if rsold < tol^2, return; end

disp('Starting iterations')
for i = 1:length(b)
    if isa(A,'function_handle')
        Ap = A(p);
    else
        Ap = A * p;
    end
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = norm(r)^2;
    if rsnew < tol^2, break; end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
    fprintf('Iter %d, err=%.6f\n',i,rsnew);
end