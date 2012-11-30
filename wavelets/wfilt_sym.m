function [h,g,a]=wfilt_sym(num_coefs)

%SYMLETS    Generates the "least asymmetric" Daubechies' wavelets or "symlets".
%
%	    [H,G,RH,RG] = SYMLETS (NUM_COEFS) generates  the quadrature 
%	    filters, given by Daubechies, with zeros of the trigonometrical 
%	    polynomial selected alternatively inside and outside the
%	    unit circle. NUM_COEFS specifies the number of coefficients.
%
%	    In this way "least asymmetric" compactly supported wavelets are
%	    obtained. H is the analysis lowpass filter, RH the synthesis lowpass 
%	    filter, G the analysis higthpass filter and RG the synthesis
%	    highpass filter.
%
%	    References: I. Daubechies
%			"Orthonormal bases of compactly supported wavelets.
%		        II.Variations on a theme.", SIAM
%			March 1993

%--------------------------------------------------------
% Copyright (C) 1994, 1995, 1996, by Universidad de Vigo 
% Author: Jose Martin Garcia
% e-mail: Uvi_Wave@tsc.uvigo.es



if rem(num_coefs,2)
   error( 'Error: NUM_COEFS must be even!!!')
end

if num_coefs==2	    % Haar filters
	rh=0.5*[sqrt(2) sqrt(2)];
	[rh,rg,h,g]=rh2rg(rh);
	return
end

N=num_coefs/2;

poly=trigpol(N);    %Calculate trigonometric polynomial 

ceros=roots(poly);  %Calculate roots

realzeros=[];
imagzeros=[];
numrealzeros=0;
numimagzeros=0;

for i=1:2*(N-1)
	if (imag(ceros(i))==0)
		numrealzeros=numrealzeros+1;
		realzeros(numrealzeros)=ceros(i);
	else
		numimagzeros=numimagzeros+1;
		imagzeros(numimagzeros)=ceros(i);
	end
end


%% complex zeros are grouped together
i=0;
for cont=1:numimagzeros/4
	modulos(cont)=abs(imagzeros(cont+i));
	alfa(cont)=angle(imagzeros(cont+i));
	i=i+1;
end

%% Calculate phase contribution of complex and real zeros for all the
%% combination of these zeros. Each group of zeros is identified with a binary
%% number.

indice=2^(numimagzeros/4+numrealzeros/2);
fase=zeros(indice,1001);
for cont=0:indice-1,
	bin=dec2bina(cont,log2(indice));
   	for i=1:length(bin)-numrealzeros/2
		if bin(i)
			R=1/modulos(i);
		else
			R=modulos(i);
		end
		alf=alfa(i);
		fase(cont+1,:)=fase(cont+1,:)+atang1(R,alf);
	end
	ind=1;
	for i=length(bin)-numrealzeros/2+1:length(bin)
		if bin(i)
			R=realzeros(ind+1);		
			R=realzeros(ind+1);
		else
			R=realzeros(ind);
		end
		ind=ind+2;
	 	fase(cont+1,:)=fase(cont+1,:)+atang2(R);

	end	
end

%% To retain only the non linear part of the phase.

fas=linealiz(fase);

imagzeros=[];
zerosreales=[];


%% To see which phase is closer to zero we select the one with minimun variance

[maximo,pos]=min(sum(fas'.^2));  

bin=dec2bina(pos-1,log2(indice));

for i=1:length(bin)-numrealzeros/2
	if bin(i)
		z1=1/modulos(i)*exp(j*alfa(i));
	else
		z1=modulos(i)*exp(j*alfa(i));	
	end
	imagzeros=[imagzeros z1 conj(z1)];
end

ind=1;
for i=length(bin)-numrealzeros/2+1:length(bin)
	if bin(i)
		zerosreales=[zerosreales realzeros(ind+1)];
	else
		zerosreales=[zerosreales realzeros(ind)];
	end
	ind=ind+2;
end

% Construction of rh from its zeros

numrealzeros=numrealzeros/2;
numimagzeros=numimagzeros/2;

rh=[1 1];

for i=2:N
	rh=conv(rh,[1 1]);
end

for i=1:numrealzeros
	rh=conv(rh,[1 -zerosreales(i)]);
end

for i=1:2:numimagzeros
	rh=conv(rh,[1 -2*real(imagzeros(i)) abs(imagzeros(i))^2]);
end

% Once ho is factorized in its zeros, it must be normalized multiplying by "cte".

cte=sqrt(2)/sum(rh);
rh=cte*rh;
[g{1},g{2},h{1},h{2}]=rh2rg(rh);
a= [2;2];




