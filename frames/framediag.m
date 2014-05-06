function d=framediag(F,L);
%FRAMEDIAG  Compute the diagonal of the frame operator
%   Usage: d=framediag(F,L);
%
%   `framediag(F,L)` computes the diagonal of the frame operator for a
%   frame of type *F* of length *L*.
%
%   The diagonal of the frame operator can for instance be used as a
%   preconditioner.
%
%   See also: franaiter, frsyniter

callfun = upper(mfilename);
complainif_notenoughargs(nargin,2,callfun);
complainif_notposint(L,'L',callfun);
complainif_notvalidframeobj(F,callfun);

% Standard response, works for all tight and orthonormal systems
d=ones(L,1);

switch(F.type)
  case 'gen'
    d=diag(F.g*F.g');
    
  case {'dgt','dgtreal'}
    d=gabframediag(F.g,F.a,F.M,L,F.vargs{:});  
  
  case {'dwilt','wmdct'}
    d=wilframediag(F.g,F.M,L);
    
  case {'filterbank','ufilterbank','filterbankreal','ufilterbankreal'}
    error('%s: TO DO: Not supported yet',upper(mfilename));

  case {'nsdgt','unsdgt','nsdgtreal','unsdgtreal'}
    d=nsgabframediag(F.g,F.a,F.M);
    
  case 'fusion'
    error('Not implemented yet.');
    
  case {'fwt','ufwt','wfbt','uwfbt','wpfbt','uwpfbt'}
    error('%s: TO DO: Not supported yet',upper(mfilename));  
  
end;


