clf;
hold all;

L=30;
dr=110;
magresp(firwin('hanning',L,'1'),'fir','dynrange',dr);
magresp(firwin('hamming',L,'1'),'fir','dynrange',dr);
magresp(firwin('blackman',L,'1'),'fir','dynrange',dr);
magresp(firwin('nuttall',L,'1'),'fir','dynrange',dr);
magresp(firwin('itersine',L,'1'),'fir','dynrange',dr);

legend('Hann','Hamming','Blackman','Nuttall','Itersine');

