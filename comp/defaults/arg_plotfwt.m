function definput=arg_plotfwt(definput)

  definput.flags.frqbands = {'uni', 'dyad'};
  definput.flags.wavplottype = {'image', 'stem', 'surf','waterfall'};
  definput.keyvals.fc=[];
  definput.keyvals.ntickpos=10;
  definput.keyvals.tick=[];
  
  definput.groups.audtick={'tick',[0,100,250,500,1000,2000,4000,8000,16000,32000]};


