function definput = arg_fwt(definput)

% default NULL incorporated for the absence of the flag have to be known too
definput.flags.type = {'type_null','dec', 'undec', 'dtdwt','hddwt','full','cust'};
definput.flags.ext=  {'ext_null','per','zpd','sym','symw','asym','asymw','ppd','sp0'};
