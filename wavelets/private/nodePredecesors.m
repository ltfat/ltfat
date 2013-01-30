function pred = nodePredecesors(nodeNo,treeStruct)
pred = [];
tmpNodeNo = nodeNo;
while treeStruct.parents(tmpNodeNo)~=0
   tmpNodeNo = treeStruct.parents(tmpNodeNo);
   pred(end+1) = tmpNodeNo;
end