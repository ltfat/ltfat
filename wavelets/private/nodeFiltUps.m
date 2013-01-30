function upsNo = nodeFiltUps(nodeNo,treeStruct)
tmpNodeNo =nodeNo;
upsNo = 1;
while(treeStruct.parents(tmpNodeNo))
    parentNo = treeStruct.parents(tmpNodeNo);
    upsNo=upsNo*treeStruct.nodes{parentNo}.a(find(treeStruct.children{parentNo}==tmpNodeNo));
    tmpNodeNo=parentNo;
end