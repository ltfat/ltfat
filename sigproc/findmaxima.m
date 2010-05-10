function [maximaout,yout] = findmaxima(x,F)
%FINDMAXIMA  Find location of local maxima
%  Usage: [maximaout,yout] = findmaxima(x,F);
%
%  F is number of maxima to find.
%  

%  Google for local maxima matlab

% Unwrap to vector
x = x(:);
% Identify whether signal is rising or falling
upordown = sign(diff(x));
% Find points where signal is rising before, falling after
maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
maxima   = find(maxflags);

% Sort them with the largest one first.
values=x(maxima);
[Y,ii]=sort(values);
ii=flipud(ii);

% tlen handles the case where fewer than F local maxima are found.
tlen=min(length(maxima),F);

% Keep the largest tlen local maxima, sort them by index.
ikeep=sort(maxima(ii(1:tlen)));

maximaout=zeros(F,1);
yout=zeros(F,1);

% Handle the case where fewer than F local maxima are found.
maximaout(1:tlen)=ikeep;
yout(1:tlen)=x(ikeep);
