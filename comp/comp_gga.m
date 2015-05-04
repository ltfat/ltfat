function c = comp_gga(f,indvec)
%COMP_GGA Generalized Goertzel Algorithm
%


%% Initialization
L = size(f,1);
W = size(f,2);
no_freq = length(indvec); %number of frequencies to compute
classname = assert_classname(f);
c = zeros(no_freq,W,classname); %memory allocation for the output coefficients

%% Computation via second-order system
% loop over the particular frequencies

for cnt_freq = 1:no_freq
    
    %for a single frequency:
    %a/ precompute the constants
    pik_term = 2*pi*(indvec(cnt_freq))/(L);
    cos_pik_term2 = cos(pik_term) * 2;
    cc = exp(-1i*pik_term); % complex constant
    for w=1:W
    %b/ state variables
       s0 = 0; 
       s1 = 0;
       s2 = 0;
       %c/ 'main' loop
       for ind = 1:L-1 %number of iterations is (by one) less than the length of signal
          %new state
          s0 = f(ind,w) + cos_pik_term2 * s1 - s2;  % (*)
          %shifting the state variables
          s2 = s1;
          s1 = s0;
       end
       %d/ final computations
       s0 = f(L,w) + cos_pik_term2 * s1 - s2; %correspond to one extra performing of (*)
       c(cnt_freq,w) = s0 - s1*cc; %resultant complex coefficient
    
       %complex multiplication substituting the last iteration
       %and correcting the phase for (potentially) non-integer valued
       %frequencies at the same time
       c(cnt_freq,w) = c(cnt_freq,w) * exp(-1i*pik_term*(L-1));
    end
end
