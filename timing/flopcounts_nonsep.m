function out = flopcounts_nonsep(L,a,M,lt)
%
%   Compute flopcounts for the nonseparable algorithms.

% Constants for the rectangular case, as if lt=[0 1]
N=L/a;
b=L/M;
c=gcd(a,M);
p=a/c;
q=M/c;
d=N/q;

% Constants for the multiwin factorization, M does not change. For q and
% p it holds that q_mw/p_mw=q/(p*lt(2)) 
a_mw=a*lt(2);
N_mw=L/a_mw;
c_mw=gcd(a_mw,M);
p_mw=a_mw/c_mw;
q_mw=M/c_mw;
d_mw=N_mw/q_mw;

% Constants for the Smith factorization, p and q does not change
% The product c*d is unchanged because of this.
latm = latticetype2matrix(L,a,M,lt);
[U,S,V] = smithnf(latm);
[a_sm,M_sm,~]=matrix2latticetype(L,S);

N_sm=L/a_sm;
c_sm=M_sm/q;
d_sm=N_sm/q;


% Constants for the Shear factorization, p and q does not change
% The product c*d is unchanged because of this.
s=b*lt(1)/lt(2);
[s0,s1,X] = shearfind(a,b,s,L);
M_sh = L/X;
a_sh = a*b/X;    
N_sh = L/a_sh;
c_sh=M_sh/q;
d_sh=N_sh/q;

k_freq=(s0>0);
k_time=(s1>0);

fc_mw_dual = 16*L*p_mw*lt(2)+4/3*c_mw*d_mw*p_mw^3+8*L*lt(2)*log2(d_mw);
fc_sm_dual = 16*L*p+4/3*c*d*p^3+8*L*log2(d_sm)+36*L+16*L*log2(L);
fc_sh_dual = 16*L*p+4/3*c*d*p^3+8*L*log2(d_sh)+(k_time+k_freq)*12*L+k_freq*16*L*log2(L);


fc_mw_dgt  = L*8*q_mw*lt(2)+4*L*(1+q/p)*log2(d_mw)*lt(2)+4*M*N*log2(M)+6*M*N;

fc_sm_dgt  = L*8*q+4*L*(1+q/p)*log2(d_sm)+4*M*N*log2(M)+18*L+8*L*log2(L)+6*M*N;

if k_freq
    fc_sh_dgt  = L*8*q+4*L*(1+q/p)*log2(c_sh)+4*M*N*log2(N)+...
        (k_time+1)*6*L+4*L*log2(L)+6*M*N;
else
    fc_sh_dgt  = L*8*q+4*L*(1+q/p)*log2(d_sh)+4*M*N*log2(M)+...
        k_time*6*L+6*M*N;    
end;

[c,d]

% pack it up
out = [fc_mw_dual,fc_sm_dual,fc_sh_dual,...
       fc_mw_dgt, fc_sm_dgt, fc_sh_dgt];