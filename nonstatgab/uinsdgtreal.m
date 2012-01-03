function f=insdgtreal(c,g,a,M,Ls)
%INSDGTREAL  Inverse nonstationary discrete Gabor transform
%   Usage:  f=insdgtreal(c,g,a,M,Ls);
%
%   Input parameters:
%         c     : Cell array of coefficients.
%         g     : Cell array of window functions.
%         a     : Vector of time positions of windows.
%         M     : Vector of numbers of frequency channels.
%         Ls    : Length of input signal.
%   Output parameters:
%         f     : Signal.
%
%   INSDGTREAL(c,g,a,M,Ls) computes the nonstationary Gabor expansion of the 
%   input coefficients c.
%
%   INSDGTREAL is used to invert the function NSDGTREAL. Read the help of NSDGTREAL
%   for details of variables format and usage.
%
%   For perfect reconstruction, the windows used must be dual windows of 
%   the ones used to generate the coefficients. The windows can be
%   generated unsing NSGABDUAL.
%
%   See also:  nsdgt, nsgabdual, nsgabtight
%
%   Demos:  demo_nsdgt
%
%   References: ltfatnote010

%   AUTHOR : Florent Jaillet
%   TESTING: TEST_NSDGT
%   REFERENCE: 
%   Last changed 2009-05


timepos=cumsum(a)-a(1);
  
N=length(c); % Number of time positions

W=size(c{1},2); % Number of signal channels

f=zeros(Ls,W); % Initialisation of the result

for ii=1:N
  shift=floor(length(g{ii})/2);
  % Note: the *size(c{ii},1) in the following is here to ensure the 
  % coherence of the ifft normalisation convention with the function dgt 
  % of ltfat
  temp=ifftreal(c{ii}*M(ii),M(ii),1); 
    
  if size(c{ii},1)<length(g{ii}) 
    % The number of frequency channels is smaller than window length, 
    % we need to periodize.
    % We have to periodize symetrically around time 0 which requires
    % some heavy indexing because time zero is index 1 and negative times
    % are at the end of the vector.
    temp=circshift(temp,mod(floor(length(g{ii})/2),length(g{ii})));
    x=floor(length(g{ii})/size(c{ii},1));
    y=length(g{ii})-x*size(c{ii},1);
    temp=[repmat(temp,x,1);temp(1:y,:)];
  else
    temp=circshift(temp,shift);
  end
  
  % Windowing with the synthesis window
  % Possible improvement: the following repmat is not needed if W=1 
  % (but the time difference when removing it should be really small)
  temp=temp(1:length(g{ii}),:).*repmat(circshift(g{ii},shift),1,W);
  
  % Possible improvement: the following could be computed faster by 
  % explicitely computing the indexes instead of using modulo
  tempind=mod((1:length(g{ii}))+timepos(ii)-shift-1,Ls)+1;
  f(tempind,:)=f(tempind,:)+temp;
end

%OLDFORMAT
