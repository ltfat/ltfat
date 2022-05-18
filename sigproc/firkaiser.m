function g=firkaiser(L,beta,varargin)
%FIRKAISER  Kaiser-Bessel window
%   Usage:  g=firkaiser(L,beta);
%           g=firkaiser(L,beta,...);
%
%   `firkaiser(L,beta)` computes the Kaiser-Bessel window of length *L* with
%   parameter *beta*. The smallest element of the window is set to zero when
%   the window has an even length. This gives the window perfect whole-point
%   even symmetry, and makes it possible to use the window for a Wilson
%   basis.
%
%   `firkaiser` takes the following flags at the end of the input arguments:
%
%     'normal'   Normal Kaiser-Bessel window. This is the default.
%
%     'derived'  Derived Kaiser-Bessel window.
%
%     'wp'       Generate a whole point even window. This is the default.
%
%     'hp'       Generate half point even window.
%  
%   Additionally, `firkaiser` accepts flags to normalize the output. Please
%   see the help of |setnorm|. Default is to use `'null'` normalization.
%
%   Note that odd-length Derived Kaiser-Bessel windows are not
%   mathematically defined, yet they are supported by this code.
%
%   See also: firwin, setnorm
%
%   References: opsc89

%   AUTHOR: unknown. Additions by Clara Hollomey

if nargin<2
  error('Too few input arguments.');
end;

if numel(beta)>1
  error('beta must be a scalar.');
end;

% Define initial value for flags and key/value pairs.
definput.import={'setnorm'};
definput.importdefaults={'null'};
definput.flags.centering={'wp','hp'};
definput.flags.stype={'normal','derived'};

[flags,keyvals]=ltfatarghelper({},definput,varargin);

cent=0;
if flags.do_hp
  cent=.5;
end;

if flags.do_normal
 
  if (L == 1)
    g = 1;
  else
    m = L - 1;
    k = (0:m)'+rem(L,2)/2-.5+cent;
    k = 2*beta/(m)*sqrt(k.*(m-k));
    g = besseli(0,k)/besseli(0,beta);
  end;

  g=ifftshift(g);
 
  if ((flags.do_wp && rem(L,2)==0) || ...
      (flags.do_hp && rem(L,2)==1))
   
    % Explicitly zero last element. This is done to get the right
    % symmetry, and because that element sometimes turns negative.
    g(floor(L/2)+1)=0;
  end;
 
else
 
  %if rem(L,2)==1    
  %  error('The length of the choosen window must be even.');
  %end;
 
  if flags.do_wp
    %if rem(L,4)==0
    %  L2=L/2+2;
    %else
      L2=floor(L/2+1);
    %end;
  else
    L2=floor((L+1)/2);
  end;
 
  % Compute a normal Kaiser window
  g_normal=fftshift(firkaiser(L2,beta,flags.centering));
 
  g1=sqrt(cumsum(g_normal(1:L2))./sum(g_normal(1:L2)));
 
  if flags.do_wp
    if rem(L,2)==0
      g=[flipud(g1);...
         g1(2:L/2)];
    else
      g=[flipud(g1);...
         g1(1:floor(L/2))];
    end;    
  else
      if rem(L,2)==0
         g=[flipud(g1);0;...
         g1(1:end-1)];
      else
         g=[flipud(g1);...
         g1(1:end-1)];
      end
  end;

  if ((flags.do_wp && rem(L,2)==0)) %|| ...
      %(flags.do_hp && rem(L,2)==1))
   
    % Explicitly zero last element. This is done to get the right
    % symmetry, and because that element sometimes turn negative.
    g(floor(L/2)+1)=0;
  end;
 
end;

% The besseli computation sometimes generates a zero imaginary component.
g=real(g);

g=setnorm(g,flags.norm);
