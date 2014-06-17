function tgrad = pgrpdelay(g,L,varargin)


[g,info] = comp_fourierwindow(g,L,upper(mfilename));

H = (comp_transferfunction(g,L));
Harg = angle(H);
Hmod = abs(H);

% Forward approximation
tgrad_1 = Harg-circshift(Harg,-1);
tgrad_1 = tgrad_1 - 2*pi*round(tgrad_1/(2*pi));
% Backward approximation
tgrad_2 = circshift(Harg,1)-Harg;
tgrad_2 = tgrad_2 - 2*pi*round(tgrad_2/(2*pi));
% Average
tgrad = (tgrad_1+tgrad_2)/2;
 
tgrad = tgrad/(2*pi)*L;


%plot(tgrad);
