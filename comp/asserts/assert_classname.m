function classname = assert_classname(varargin)
% ASSERT_CLASSNAME 
%
% Returns name of the least "simplest" common data type.

% Array of data types to be checked. Ordered from the "simplest" to the
% most "complex".
typesToTest = {'single','double'}; 

if nargin==0 || isempty(varargin)
   classname = 'double';
   return;
end

if ~all(cellfun(@(vEl) isnumeric(vEl),varargin))
   error('%s: Parameters are not numeric types. ',upper(mfilename));
end

% Shortcut to double
if all(cellfun(@(vEl) isa(vEl,'double'),varargin))
   classname = 'double';
   return;
end

% Go trough all the types, halt if any of the inputs is of the specified
% type.
for ii=1:numel(typesToTest)
   if any(cellfun(@(vEl) isa(vEl,typesToTest{ii}),varargin))
      classname = typesToTest{ii};
      return;
   end
end
