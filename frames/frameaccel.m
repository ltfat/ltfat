function F=frameaccel(F,Ls);  
%FRAMEACCEL  Precompute structures
%   Usage: F=frameaccel(F,L);
%
%   `F=frameaccel(F,Ls)` precomputes certain structures that makes the basic
%   frame operations |frana| and |frsyn| faster (like instantiating the
%   window from a textual description). If you only need to call the
%   routines once, calling `frameaccel` first will not provide any total
%   gain, but if you are repeatedly calling these routines, for instance in
%   an iterative algorithm, is will be a benefit.
%
%   Notice that you need to input the signal length *Ls*, so this routines
%   is only a benefit if *Ls* stays fixed.
%
%   If `frameaccel` is called twice for the same transform length, no
%   additional computations will be done.
%
%   See also: frame, frana, framelength, framelengthcoef
  
L=framelength(F,Ls);

if (isfield(F,'L') && (L==F.L))
  % Quick return, we have already accelerated
  return
end;

F.L=L;

if strcmp(F.type,'fusion')
    for ii=1:F.Nframes
        accel_frames{ii}=frameaccel(F.frames{ii},Ls);
    end;
    F=frame('fusion',F.w,accel_frames{:});
    F.L=L;
    return;
end;

if strcmp(F.type,'tensor')
    for ii=1:F.Nframes
        accel_frames{ii}=frameaccel(F.frames{ii},Ls);
    end;
    F=frame('tensor',accel_frames{:});
    F.L=L;
    return;
end;


if ~isfield(F,'g')
  % Quick exit, the frame does not use a window. In this case, the frame
  % always has a factorization
  
  % Default values for a lot of transforms
  F.isfac=1;

  return;
end;
  
% From this point and on, we are sure that F.g exists

if ~isempty(F.g)
    
    switch(F.type)
      case 'gen'
        F.isfac=~issparse(F.g);  
      case {'dgt','dgtreal'}
        F = frame(F.type,gabwin(F.g,F.a,F.M,L,F.kv.lt),F.origargs{2:end});
        F.isfac=1;
      case {'dwilt','wmdct'}
        F = frame(F.type,wilwin(F.g,F.M,L,upper(mfilename)),F.origargs{2:end});
        F.isfac=1;
      case {'filterbank','ufilterbank'}
        [F.g,F.g_info]  = filterbankwin(F.g,F.a,L);
        F.isfac=F.g_info.isfac;
      case {'filterbankreal','ufilterbankreal'}
        [F.g,F.g_info]  = filterbankwin(F.g,F.a,L,'real');
        F.isfac=F.g_info.isfac;
      case {'nsdgt','unsdgt','nsdgtreal','unsdgtreal'}
        [F.g,F.g_info]  = nsgabwin(F.g,F.a,F.M);
        F.isfac=F.g_info.isfac;
      case 'fwt'
        F.isfac=1;
  end;
  
end;

F.L=L;
