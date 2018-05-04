function x = ref_pcg(A, b, x, tol)
if nargin < 3
    x = zeros(size(b));
    tol = 1e-10;
end
if nargin < 4
    tol = 1e-10;
end

r = b - A * x;
p = r;
rsold = norm(r)^2;

for i = 1:length(b)
    Ap = A * p;
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = norm(r)^2;
    if rsnew < tol^2, break; end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
end