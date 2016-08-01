function sstring=ltfatstatusstring(status)

[~,~,enuminfo]=libltfatprotofile;

map = structfun(@(a) a==status ,enuminfo.ltfaterr_status);
names = fieldnames(enuminfo.ltfaterr_status);
sstring = names{map};
  


