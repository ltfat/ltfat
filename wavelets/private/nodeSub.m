function subNo = nodeSub(nodeNo,treeStruct)
subNo = nodeFiltUps(nodeNo,treeStruct).*treeStruct.nodes{nodeNo}.a;