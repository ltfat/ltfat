function gt=nsgabtight(g,a,Ls)
%NSGABTIGHT  Canonical tight window for nsionnary Gabor frames
%   Usage:  gt=nsgabtight(g,a,Ls)
%
%   Input parameters:
%         g     : Cell array of windows
%         a     : Vector of time positions of windows.
%         Ls    : Length of analyzed signal.
%   Output parameters:
%         gt : Cell array of canonical tight windows
%
%   NSGABTIGHT(g,a,Ls) computes the canonical tight windows of the 
%   nonstationary discrete Gabor frame defined by windows given in g an  
%   time-shifts given by _a.
%   
%   NSGABTIGHT is designed to be used with functions NSDGT and INSDGT.  Read
%   the help on NSDGT for more details about the variables structure.
%
%   The computed tight windows are only valid for the 'painless case', that
%   is to say that they ensure perfect reconstruction only if for each 
%   window the number of frequency channels used for computation of NSDGT is
%   greater than or equal to the window length. This correspond to cases
%   for which the frame operator is diagonal.
%
%   See also:  nsgabtight, nsdgt, insdgt
%
%R  ltfatnote010

%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGT
%   REFERENCE:

timepos=cumsum(a)-a(1);

N=length(a); % Number of time positions
f=zeros(Ls,1); % Diagonal of the frame operator

% Compute the diagonal of the frame operator:
% sum up in time (overlap-add) all the contributions of the windows as if 
% we where using windows in g as analysis and synthesis windows
for ii=1:N
  shift=floor(length(g{ii})/2);
  temp=abs(circshift(g{ii},shift)).^2*length(g{ii});
  tempind=mod((1:length(g{ii}))+timepos(ii)-shift-1,Ls)+1;
  f(tempind)=f(tempind)+temp;
end

% As we want tight frame, we will use the sqrt of the operator
f=sqrt(f);

% Initialize the result with g
gt=g;

% Correct each window to ensure perfect reconstrution
for ii=1:N
  shift=floor(length(g{ii})/2);
  tempind=mod((1:length(g{ii}))+timepos(ii)-shift-1,Ls)+1;
  gt{ii}(:)=circshift(circshift(g{ii},shift)./f(tempind),-shift);
end
