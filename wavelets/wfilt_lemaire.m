function [h,g,a]=wfilt_lemaire(num_coefs)

%LEMARIE    Generates the quadrature filters given by Battle and Lemarie.
%
%	    [H,G,RH,RG] = LEMARIE (NUM_COEFS) returns the coeficients of
%	    orthonormal  Battle-Lemarie wavelets. NUM_COEFS specifies the
%           number of coefficients. 
%
%	    H is the analysis lowpass filter, RH the synthesis lowpass 
%	    filter, G the analysis higthpass filter and RG the synthesis
%	    highpass filter.
%
%	    References: S. Mallat, "A Theory for Multiresolution Signal
%	    		Decomposition: The Wavelet Representation", IEEE Trans.
%			on Patt. An. and Mach. Intell., July 1989

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es



L = 1024;
H = wfreq_lemaire(L);
hh=real(ifft(H,L));
hh=[ hh(L-floor(num_coefs/2)+1:L) hh(1:ceil(num_coefs/2))];
hh=sqrt(2)/sum(hh)*hh;

[g{1},g{2},h{1},h{2}]=rh2rg(fliplr(hh));
a= [2;2];




