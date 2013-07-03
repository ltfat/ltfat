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
%   `filterbankresponse(g,a,L,fs)` specifies the sampling rate *fs*. This
%   is only used for plotting purposes.
%
%   `filterbankresponse` takes the following optional parameters:
%
%      'fs',fs    Sampling rate, used only for plotting.
%
%      'complex'  Assume that the filters cover the entire frequency
%                 range. This is the default.
%
%      'real'     Assume that the filters only cover the positive
%                 frequencies (and is intended to work with real-valued
%                 signals only).
%
%      'noplot'   Don't plot the response, just return it.
%
%      'plot'     Plot the response using |plotfftreal|.
%
%   See also: filterbank, filterbankbounds
  
definput.flags.ctype={'complex','real'};
definput.flags.plottype={'noplot','plot'};
definput.keyvals.fs=[];
[flags,kv,fs]=ltfatarghelper({'fs'},definput,varargin);

[g,info]=filterbankwin(g,a,L,'normal');
M=info.M;

gf=comp_filterbankresponse(g,info.a,L,flags.do_real);

if flags.do_plot
    if flags.do_real
        plotfftreal(gf(1:floor(L/2)+1),fs,'lin');
    else
        plotfft(gf,fs,'lin');
    end;
end;
