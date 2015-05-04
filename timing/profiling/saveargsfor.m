function saveargsfor(matfile,args,varargin)
%SAVEARGSFOR Save input arguments for mexExecuter

if nargin<2
   error('%s: Too few input arguments.',upper(mfilename) );
end

if isempty(args)
   error('%s: No variables to be saved.',upper(mfilename) );
end

definput.keyvals.res=[];
[flags,kv]=ltfatarghelper({},definput,varargin);


varNo = numel(args);
varNames = cell(varNo,1);


for ii=1:numel(args)
   varNames{ii} = sprintf('arg%i',ii-1);
   eval([varNames{ii},' = args{ii};']);   
end

if ~isempty(kv.res)
   varNames{end+1} = 'res';
   eval([varNames{end},' = kv.res;']);
end

save(matfile,varNames{:});
