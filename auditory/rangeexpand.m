function outsig = rangeexpand(insig,varargin);
%RANGEEXPAND  Expand the dynamic range of a signal
%   Usage:  sig = rangeexpand(insig,mu,sigweight);
%
%   `rangeexpand(insig,mu,sigweight)` inverts a previously
%   applied $\mu$-law companding to the signal *insig*. The parameters
%   *mu* and *sigweight* must match those from the call to |rangecompress|
%
%   `rangeexpand` takes the following optional arguments:
%
%     'mulaw'   Do mu-law compression, this is the default.
%
%     'alaw'    Do A-law compression.
%
%     'mu',mu   $\mu$-law parameter. Default value is 255.
%  
%   References: jano90

% AUTHOR: Bruno Torresani and Peter L. Søndergaard

definput.flags.method={'mulaw','alaw'};
definput.keyvals.mu=255;
definput.keyvals.A=87.7;
[flags,kv]=ltfatarghelper({},definput,varargin);

if flags.do_mulaw

  cst = (1+kv.mu);
  outsig = cst.^(abs(insig));
  outsig = sign(insig) .* (outsig-1);
  outsig = outsig/kv.mu;

end;

if flags.do_alaw
  absx=abs(insig);
  tmp=1+log(kv.A);
  mask=absx<1/tmp;

  outsig = sign(insig).*(mask.*(absx*tmp/kv.A)+(1-mask).*exp(absx*tmp-1)/kv.A);
end;

