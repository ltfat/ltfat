function gf=filterbankresponse(g,a,L,varargin)
%FILTERBANKRESPONSE  Response of filterbank as function of frequency
%   Usage:  gf=filterbankresponse(g,a,L);
%      
%   `filterbankresponse(g,a,L)` computes the total response in frequency of
%   a filterbank specified by *g* and *a* for a signal length of
%   *L*. This corresponds to summing up all channels. The output is a
%   usefull tool to investigate the behaviour of the windows, as peaks
%   indicate that a frequency is overrepresented in the filterbank, while
%   a dip indicates that it is not well represented.
%
%   In mathematical terms, this function computes the diagonal of the
%   Fourier transform of the frame operator when the filterbank is painless.
%
%   `filterbankresponse(g,a,L,'real')` does the same for a filterbank
%   intended for positive-only filterbank.
%
%   See also: filterbank, filterbankbounds
  
definput.flags.ctype={'complex','real'};
definput.flags.plottype={'noplot','plot'};
[flags,kv]=ltfatarghelper({},definput,varargin);


% G1 is done this way just so that we can determine the data type.
G1=comp_transferfunction(g{1},L);
gf=abs(G1).^2*(L/a(1));

M=numel(g);
  
for m=2:M
  gf=gf+abs(comp_transferfunction(g{m},L)).^2*(L/a(m));
end;
  
if flags.do_real
  gf=gf+involute(gf);   
end;

gf=gf/L;

if flags.do_plot
    if flags.do_real
        plotfftreal(gf(1:floor(L/2)+1),'lin');
    else
        plotfft(gf,'lin');
    end;
end;
