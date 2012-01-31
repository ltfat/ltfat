function gamma=ref_gabdualns_3(g,V);
%REF_GABDUALNS_3  GABDUALNS by multiwindow method.
%   Usage:  gamma=ref_gabdualns_3(g,V);
%

[gm,a,M]=ref_nonsep2multiwin(g,V);

gmd=gabdual(gm,a,M);

gamma=gmd(:,1);


