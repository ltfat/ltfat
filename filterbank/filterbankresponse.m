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
%   Fourier transform of the frame operator.
%
%   `filterbankresponse(g,a,L,'real')` does the same for a filterbank
%   intended for positive-only filterbank.
%
%   See also: filterbank, filterbankbounds
  
definput.flags.ctype={'complex','real'};
definput.flags.plottype={'plot','noplot'};
[flags,kv]=ltfatarghelper({},definput,varargin);

gf=zeros(L,1);
M=numel(g);
  
for m=1:M
  gf=gf+abs(fft(middlepad(g{m},L))).^2;
end;
  
if flags.do_real
  gf=gf+involute(gf);   
end;

if flags.do_plot
  plotfft(gf,'lin');
end;