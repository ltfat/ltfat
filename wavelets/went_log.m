function E = went_log(x)
%u = x(:)./norm(x(:));
%u = u(u~=0);
x = x(x~=0);
E = sum(log(x(:).^2));