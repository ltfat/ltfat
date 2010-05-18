function h=ref_pconv(f,g,ctype)
%REF_PCONV  Reference PCONV
%   Usage:  h=ref_pconv(f,g)
%
%   PCONV(f,g) computes the periodic convolution of f and g.

% AUTHOR: Peter Soendergaard

L=length(f);
h=zeros(L,1);

if nargin==2
  ctype='';
end;

% The following is a HACK to work around broken support for switch
% statements of empty strings in some Octave versions.
if strcmp(ctype,'')
  ctype=' ';
end;

switch(lower(ctype))
  case {' '}    
    for ii=0:L-1
      for jj=0:L-1
	h(ii+1)=h(ii+1)+f(jj+1)*g(mod(ii-jj,L)+1);
      end;
    end;
  case {'r'}
    for ii=0:L-1
      for jj=0:L-1
	h(ii+1)=h(ii+1)+f(jj+1)*conj(g(mod(jj-ii,L)+1));
      end;
    end;
  case {'rr'}
    for ii=0:L-1
      for jj=0:L-1
	h(ii+1)=h(ii+1)+conj(f(mod(-jj,L)+1))*conj(g(mod(jj-ii,L)+1));
      end;
    end;
end;

% Clean output if input was real-valued
if isreal(f) && isreal(g)
  h=real(h);
end;