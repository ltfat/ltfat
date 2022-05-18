function definput=arg_normalize(definput)
  
% Both 'null' and 'empty' do no scaling when normalize is called
% directly.
% When used in different functions,
% 'empty' can be set as default by definput.importdefaults={'empty'};
% to detect whether any of the other flags were set.
%
warning(['LTFAT: NORMALIZE has been deprecated and will be removed',...
         ' in the future releases, please use SETNORM instead.']);
         
definput.flags.norm={'2','1','inf','area','energy','peak',...
                     's0','rms','null','wav','norm_notset'};


