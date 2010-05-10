function [flags,keyvals,varargout]  = ltfatarghelper(posdepnames,defnopos,arglist,callfun)
%LTFATARGHELPER  Parse arguments for LTFAT
%   Usage: [flags,varargout]  = ltfatarghelper(posdepnames,posdepdef,arglist);
%
%   Input parameters:
%      posdepnames : Names of the position dependant parameters.
%      posdepdef   : Default values for pos. dep. pars. (cell array)
%      defnopos    : Definitions of the pos. independent pars (cell array)
%      arglist     : Commandline of the calling function (varargin)
%      callfun     : Name of calling function
%
%   Output parameters:
%      nopos     : Cell array containing parsed pos. independent pars.
%      varargout : The position dependant pars. properly initialized
%
%   [nopos,varargout]=LTFATARGHELPER(posdepnames,posdepdef,arglist) assist in
%   parsing input parameters for a function in LTFAT. Parameters come in
%   two categories:
%  
%      * Position dependant parameters. These must not be strings. These are
%      the first parameters passed to a function. If they are not present a
%      default value is given.
%
%      * Position independant parameters. These come in key,value pairs
%      at the end of the parameter list of the calling function.
%
%   The typical way of calling LTFATARGHELPER is as follows:
%  
%C    FIXME
%
%   This will pass any (or a number) of the input arguments from the calling
%   function onto LTFATARGHELPER. In this case, there are 2 position
%   dependant parameters (n and betamul), which will have default values 4
%   and [] if they are not present.

global TF_CONF;

if isempty(TF_CONF)

  basepath=which('ltfatarghelper');
  % Kill the function name and comp from the path.
  basepath=basepath(1:end-22);
  % add the base path
  addpath(basepath);
  ltfatstart;

end;

if nargin<4
  f=dbstack;  
  callfun=f(2).name;
end;

nposdep=numel(posdepnames);

if isfield(defnopos,'flags')
  defflags=defnopos.flags;
else
  defflags=struct;
end;

if isfield(defnopos,'keyvals')
  defkeyvals=defnopos.keyvals;
else
  defkeyvals=struct;
end;

if isfield(defnopos,'groups')
  groups=defnopos.groups;
else
  groups=struct;
end;

total_args = numel(arglist);

  % Determine the position of the first optional argument.
  % If no optional argument is given, return nposdep+1
  first_str_pos = 1;
  while first_str_pos<=total_args && ~ischar(arglist{first_str_pos}) 
    first_str_pos = first_str_pos +1;    
  end;
    
  % If more than nposdep arguments are given, the first additional one must
  % be a string
  if (first_str_pos>nposdep+1)
    error('%s: Too many input arguments',upper(callfun));
  end;

  n_first_args=min(nposdep,first_str_pos-1);

  keyvals=defkeyvals;      
  
  % Copy the given first arguments
  for ii=1:n_first_args
    keyvals.(posdepnames{ii})=arglist{ii};
  end;

  % Initialize the position independent parameters.
  % and create reverse mapping of flag -> group
  flagnames=fieldnames(defflags);
  flags=struct;
  flagsreverse=struct;
  for ii=1:numel(flagnames)
    name=flagnames{ii};
    flaggroup=defflags.(name);
    flags.(name)=flaggroup{1};
    for jj=1:numel(flaggroup)
      flagsreverse.(flaggroup{jj})=name;
      flags.(['do_',flaggroup{jj}])=0;
    end;
    flags.(['do_',flaggroup{1}])=1;
  end;
  
  %Get the rest of the arguments
  restlist = arglist(first_str_pos:end);

  %Check for default arguments
  if isfield(TF_CONF,callfun)
    s=TF_CONF.(callfun);
    restlist={s{:},restlist{:}};
  end;

  while ~isempty(restlist)
    argname=restlist{1};
    restlist=restlist(2:end);  % pop
    found=0;
    % Is this name a flag? If so, set it
    if isfield(flagsreverse,argname)
      % Unset all other flags in this group
      flaggroup=defflags.(flagsreverse.(argname));
      for jj=1:numel(flaggroup)
        flags.(['do_',flaggroup{jj}])=0;
      end;
      
      flags.(flagsreverse.(argname))=argname;
      flags.(['do_',argname])=1;
      found=1;
    end;
      
    % Is this name the key of a key/value pair? If so, set the value.
    if isfield(defkeyvals,argname)      
      keyvals.(argname)=restlist{1};
      restlist=restlist(2:end);
      found=1;
    end;

    % Is this name a group definition? If so, put the group in front of the parameters
    if isfield(groups,argname)
      s=groups.(argname);
      restlist={s{:},restlist{:}};
      found=1;
    end;      
    
    if found==0
      error('%s: Unknown parameter: %s',upper(callfun),argname);
    end;

    ii=ii+1;
  end;

% Fill varargout

varargout={};
for ii=1:nposdep
    varargout(ii)={keyvals.(posdepnames{ii})};
end;
