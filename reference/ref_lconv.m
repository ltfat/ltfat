function h = ref_lconv(f,g,ctype)
%REF_LCONV  Reference linear convolution
%   Usage:  h=ref_lconv(f,g)
%
%   PCONV(f,g) computes the linear convolution of f and g.

% AUTHOR: Jordy van Velthoven

Lf = length(f);
Lg = length(g);

Lh = Lf+Lg-1;

f = [f; zeros(Lh - Lf, 1)];
g = [g; zeros(Lh - Lg, 1)];

h = zeros(Lf+Lg-1, 1);

switch(lower(ctype))
	case {'default'}
    for ii = 0 : Lh-1
      for jj = 0 : Lh-1
        h(ii+1)=h(ii+1)+f(jj+1)*g(mod(ii-jj,Lh)+1);
      end
    end
  case {'r'}
    for ii=0:Lh-1
      for jj=0:Lh-1
	      h(ii+1)=h(ii+1)+f(jj+1)*conj(g(mod(jj-ii, Lh)+1));
      end;
   	end;
  case {'rr'}
    for ii=0:Lh-1
      for jj=0:Lh-1
	      h(ii+1)=h(ii+1)+conj(f(mod(-jj, Lh)+1))*conj(g(mod(jj-ii,Lh)+1));
      end;
    end;
end
      
